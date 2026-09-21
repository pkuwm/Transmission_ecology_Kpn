# Figure 5A: Hospital strain-sharing network
rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggnewscale)

# Settings
message("Select the hospital isolate metadata Excel file.")
metadata_file <- file.choose()

message("Select the Excel file containing isolate pairs at <=21 SNPs.")
pairs_file <- file.choose()

metadata_sheet <- "医院传播菌448株"
pairs_sheet <- 1
cross_group_only <- TRUE

output_dir <- file.path(dirname(metadata_file), "Figure_5A")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read metadata
raw <- read_excel(metadata_file, sheet = metadata_sheet)

source_column <- select.list(
    names(raw),
    title = "Select the sample source column",
    graphics = TRUE
)
if (source_column == "") stop("No sample source column selected.")

meta <- raw %>%
    transmute(
        Isolate = trimws(as.character(`服务器ID`)),
        Location = trimws(as.character(`院内传播位置大类`)),
        Label = trimws(as.character(`院内传播小类`)),
        Source = trimws(as.character(.data[[source_column]]))
    ) %>%
    mutate(across(everything(), ~ na_if(.x, "")))

if (anyNA(meta) || anyDuplicated(meta$Isolate)) {
    stop("Metadata must have complete annotations and unique isolate IDs.")
}

# Edit aliases here if the metadata use different source names
source_aliases <- c(
    "Clinical isolates" = "Clinical isolates",
    "Clinical isolate" = "Clinical isolates",
    "HCW Stool" = "HCW Stool",
    "HCW stool" = "HCW Stool",
    "Environmental Surface" = "Environmental Surface",
    "Hospital Surface" = "Environmental Surface",
    "Sewage" = "Sewage",
    "Hospital Sewage" = "Sewage",
    "Air" = "Air",
    "Hospital Air" = "Air"
)

source_colors <- c(
    "Clinical isolates" = "#C39A9E",
    "HCW Stool" = "#EC8A76",
    "Environmental Surface" = "#92789C",
    "Sewage" = "#446AA2",
    "Air" = "#C7D7ED"
)

unknown <- setdiff(unique(meta$Source), names(source_aliases))
if (length(unknown)) {
    stop("Add these values to source_aliases: ",
         paste(unknown, collapse = ", "))
}

meta$Source <- unname(source_aliases[meta$Source])

# Generate nodes
nodes <- meta %>%
    count(Location, Label, Source, name = "Isolates") %>%
    arrange(Location, Label, Source) %>%
    mutate(Node = row_number())

meta <- meta %>%
    left_join(nodes, by = c("Location", "Label", "Source"))

# Read and deduplicate unordered isolate pairs
pairs <- read_excel(pairs_file, sheet = pairs_sheet)

if ("SNP_distance" %in% names(pairs)) {
    pairs <- pairs %>% filter(SNP_distance <= 21)
}

pairs <- pairs %>%
    transmute(
        A = trimws(as.character(Isolate1)),
        B = trimws(as.character(Isolate2))
    ) %>%
    filter(!is.na(A), !is.na(B), A != "", B != "", A != B) %>%
    transmute(Isolate1 = pmin(A, B), Isolate2 = pmax(A, B)) %>%
    distinct()

# Retain pairs with both isolates in the hospital metadata
links <- pairs %>%
    inner_join(
        meta %>% select(Isolate1 = Isolate, From = Node, Group1 = Location),
        by = "Isolate1"
    ) %>%
    inner_join(
        meta %>% select(Isolate2 = Isolate, To = Node, Group2 = Location),
        by = "Isolate2"
    )

if (cross_group_only) {
    links <- links %>% filter(Group1 != Group2)
}

# Within-node pairs are not drawn as loops
links <- links %>%
    filter(From != To) %>%
    transmute(From = pmin(From, To), To = pmax(From, To)) %>%
    count(From, To, name = "Pairs")

if (!nrow(links)) stop("No links remain under the selected display rule.")

# Sector layout
nodes <- nodes %>%
    mutate(Sector = paste(Location, Label, sep = " / "))

sector_order <- unique(nodes$Sector)
nodes$Sector <- factor(nodes$Sector, levels = sector_order)

gap <- 0.045
step <- (2 * pi - length(sector_order) * gap) / nrow(nodes)
if (step <= 0) stop("Reduce the sector gap.")

sectors <- nodes %>%
    group_by(Sector) %>%
    summarise(
        Location = first(Location),
        Label = first(Label),
        n = n(),
        .groups = "drop"
    ) %>%
    mutate(
        start = lag(cumsum(n * step + gap), default = 0),
        end = start + n * step,
        middle = (start + end) / 2,
        x = 1.12 * sin(middle),
        y = 1.12 * cos(middle),
        angle = -middle * 180 / pi,
        angle = if_else(
            middle > pi / 2 & middle < 3 * pi / 2,
            angle + 180, angle
        )
    )

nodes <- nodes %>%
    left_join(sectors %>% select(Sector, start), by = "Sector") %>%
    group_by(Sector) %>%
    mutate(theta = start + (row_number() - 0.5) * step) %>%
    ungroup() %>%
    mutate(x = sin(theta), y = cos(theta))

# Curved links
edges <- links %>%
    left_join(
        nodes %>% select(From = Node, x1 = x, y1 = y, Location),
        by = "From"
    ) %>%
    left_join(
        nodes %>% select(To = Node, x2 = x, y2 = y),
        by = "To"
    )

curves <- bind_rows(lapply(seq_len(nrow(edges)), function(i) {
    e <- edges[i, ]
    data.frame(
        Edge = i,
        x = c(e$x1, 0, e$x2),
        y = c(e$y1, 0, e$y2),
        Pairs = e$Pairs,
        Location = e$Location
    )
}))

location_levels <- unique(nodes$Location)
location_colors <- setNames(
    grDevices::colorRampPalette(
        c("#A8DAD2", "#A4B7DA", "#9995B7", "#EFC3C4")
    )(length(location_levels)),
    location_levels
)

# Draw and export
p <- ggplot() +
    geom_bezier(
        data = curves,
        aes(x, y, group = Edge, colour = Location, linewidth = Pairs),
        alpha = 0.65, lineend = "round"
    ) +
    scale_colour_manual(values = location_colors, guide = "none") +
    scale_linewidth_continuous(
        range = c(0.3, 3), breaks = c(10, 100, 500),
        name = "Strain-sharing pairs"
    ) +
    geom_arc_bar(
        data = sectors,
        aes(
            x0 = 0, y0 = 0, r0 = 1.20, r = 1.26,
            start = start, end = end, fill = Location
        ),
        colour = NA
    ) +
    scale_fill_manual(values = location_colors, name = "Sampling location") +
    ggnewscale::new_scale_colour() +
    geom_point(
        data = nodes,
        aes(x, y, colour = Source, size = Isolates)
    ) +
    scale_colour_manual(values = source_colors, name = "Source") +
    scale_size_continuous(
        range = c(3, 10), breaks = c(5, 20, 100),
        name = "Number of isolates"
    ) +
    geom_text(
        data = sectors,
        aes(x, y, label = Label, angle = angle),
        size = 3.5, colour = "black"
    ) +
    guides(
        linewidth = guide_legend(
            order = 1, nrow = 1,
            override.aes = list(colour = "grey40", alpha = 1)
        ),
        size = guide_legend(
            order = 2, nrow = 1,
            override.aes = list(shape = 1, colour = "grey40")
        ),
        fill = guide_legend(order = 3),
        colour = guide_legend(order = 4, override.aes = list(size = 4))
    ) +
    coord_fixed(
        xlim = c(-1.36, 1.36), ylim = c(-1.36, 1.36),
        clip = "off"
    ) +
    theme_void(base_size = 11) +
    theme(
        legend.position = "left",
        legend.box = "vertical",
        legend.title = element_text(face = "bold"),
        plot.margin = margin(15, 15, 15, 15)
    )

print(p)

ggsave(
    file.path(output_dir, "Figure_5A_hospital_transmission_network.pdf"),
    p, width = 10, height = 8,
    bg = "white", useDingbats = FALSE
)

message("Figure saved in: ", normalizePath(output_dir))
