# SNP threshold sensitivity analysis
rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(writexl)

# Settings
message("Select the annotated SNP pair file (including pairs up to 100 SNPs).")
pair_file <- file.choose()
pair_sheet <- 1
# Use "pairs_within_100_SNPs" for the previous combined output workbook.

message("Select the isolate metadata file.")
info_file <- file.choose()
info_sheet <- 1

source_column <- "Source"
bubble_thresholds <- c(1, 10, 15, 21, 25, 35)
curve_thresholds <- 0:100
key_thresholds <- c(15, 21, 25, 35)

output_dir <- file.path(
    dirname(pair_file),
    paste0("SNP_sensitivity_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read unique unordered pairs
pairs <- read_excel(pair_file, sheet = pair_sheet) %>%
    transmute(
        Isolate1 = trimws(as.character(Isolate1)),
        Isolate2 = trimws(as.character(Isolate2)),
        SNP_distance = as.numeric(SNP_distance)
    ) %>%
    mutate(
        A = pmin(Isolate1, Isolate2),
        B = pmax(Isolate1, Isolate2)
    ) %>%
    filter(A != B)

if (
    anyNA(pairs) ||
    any(!is.finite(pairs$SNP_distance)) ||
    any(pairs$SNP_distance < 0)
) {
    stop("Missing or invalid isolate IDs or SNP distances.")
}

conflicts <- pairs %>%
    group_by(A, B) %>%
    summarise(n_distances = n_distinct(SNP_distance), .groups = "drop") %>%
    filter(n_distances > 1)

if (nrow(conflicts) > 0) {
    stop("Some isolate pairs have conflicting SNP distances.")
}

pairs <- pairs %>%
    distinct(A, B, .keep_all = TRUE) %>%
    transmute(
        Isolate1 = A,
        Isolate2 = B,
        SNP_distance
    )

# Read metadata
info <- read_excel(info_file, sheet = info_sheet) %>%
    transmute(
        StrainID = trimws(as.character(StrainID)),
        Source = tolower(trimws(as.character(.data[[source_column]]))),
        Location = trimws(as.character(Location)),
        ST = sub("^ST\\s*", "", toupper(trimws(as.character(ST)))),
        KL = sub("^KL\\s*", "", toupper(trimws(as.character(KL))))
    )

if (anyDuplicated(info$StrainID)) {
    stop("Duplicate StrainID values in metadata.")
}

info1 <- info
info2 <- info
names(info1) <- c("Isolate1", "Source1", "Location1", "ST1", "KL1")
names(info2) <- c("Isolate2", "Source2", "Location2", "ST2", "KL2")

pairs <- pairs %>%
    left_join(info1, by = "Isolate1") %>%
    left_join(info2, by = "Isolate2")

required_metadata <- c(
    "Source1", "Source2", "Location1", "Location2",
    "ST1", "ST2", "KL1", "KL2"
)

if (any(vapply(
    pairs[required_metadata],
    function(x) any(is.na(x) | x == ""),
    logical(1)
))) {
    stop("Missing pair metadata. Check isolate IDs and metadata fields.")
}

if (!all(c(pairs$Source1, pairs$Source2) %in%
         c("hospital", "community", "farm"))) {
    stop("The source column must contain Hospital, Community and Farm.")
}

if (any(pairs$ST1 != pairs$ST2)) {
    stop("Pairs with different ST assignments found. Check the input data.")
}

# Define source and clone categories
source_order <- c(
    "Hospital-Hospital",
    "Community-Community",
    "Farm-Farm",
    "Community-Farm",
    "Hospital-Community",
    "Hospital-Farm"
)

clone_order <- c(
    "ST11-KL25", "ST11-KL64", "ST11-KL47", "Other",
    "ST15", "ST17", "ST37", "ST307", "ST23"
)

source_lookup <- c(
    "hospital|hospital" = "Hospital-Hospital",
    "community|community" = "Community-Community",
    "farm|farm" = "Farm-Farm",
    "community|farm" = "Community-Farm",
    "community|hospital" = "Hospital-Community",
    "farm|hospital" = "Hospital-Farm"
)

pairs <- pairs %>%
    mutate(
        source_key = paste(
            pmin(Source1, Source2),
            pmax(Source1, Source2),
            sep = "|"
        ),
        Source_category = unname(source_lookup[source_key]),
        Clone_category = case_when(
            ST1 == "11" & KL1 == "25" & KL2 == "25" ~ "ST11-KL25",
            ST1 == "11" & KL1 == "64" & KL2 == "64" ~ "ST11-KL64",
            ST1 == "11" & KL1 == "47" & KL2 == "47" ~ "ST11-KL47",
            ST1 %in% c("15", "17", "37", "307", "23") ~ paste0("ST", ST1),
            TRUE ~ "Other"
        )
    )

# Cluster and event sensitivity
curve_data <- bind_rows(lapply(curve_thresholds, function(threshold) {
    sub <- pairs %>%
        filter(SNP_distance <= threshold)
    
    if (nrow(sub) == 0) {
        return(data.frame(
            threshold = threshold,
            n_pairs = 0,
            n_clusters = 0,
            n_events = 0
        ))
    }
    
    g <- graph_from_data_frame(
        sub %>% select(Isolate1, Isolate2),
        directed = FALSE
    )
    
    comp <- components(g)
    
    sub$Cluster <- unname(comp$membership[sub$Isolate1])
    
    events <- sub %>%
        mutate(
            Institution_A = pmin(Location1, Location2),
            Institution_B = pmax(Location1, Location2)
        ) %>%
        distinct(Cluster, Institution_A, Institution_B)
    
    data.frame(
        threshold = threshold,
        n_pairs = nrow(sub),
        n_clusters = comp$no,
        n_events = nrow(events)
    )
}))

if (max(curve_data$n_clusters) == 0) {
    stop("No linked pairs found within the selected threshold range.")
}

# Exploratory cluster-event curve
cluster_color <- "#A67C7C"
event_color <- "#8699AD"

event_scale <- max(curve_data$n_events) /
    max(curve_data$n_clusters)

key_data <- curve_data %>%
    filter(threshold %in% key_thresholds) %>%
    mutate(event_y = n_events / event_scale)

curve_plot <- ggplot(curve_data, aes(x = threshold)) +
    geom_vline(
        xintercept = key_thresholds,
        linetype = "dashed",
        colour = "grey60",
        linewidth = 0.4
    ) +
    geom_line(
        aes(y = n_clusters, colour = "Clusters"),
        linewidth = 0.9
    ) +
    geom_line(
        aes(y = n_events / event_scale, colour = "Events"),
        linewidth = 0.9
    ) +
    geom_point(
        data = key_data,
        aes(y = n_clusters),
        shape = 21, fill = cluster_color, colour = "white",
        size = 3.5, stroke = 1
    ) +
    geom_point(
        data = key_data,
        aes(y = event_y),
        shape = 21, fill = event_color, colour = "white",
        size = 3.5, stroke = 1
    ) +
    geom_label_repel(
        data = key_data,
        aes(y = n_clusters, label = paste0(threshold, ": C=", n_clusters)),
        colour = "#6B4F4F", fill = "#E8DCDC",
        size = 3, direction = "y",
        nudge_y = max(curve_data$n_clusters) * 0.05,
        seed = 2026
    ) +
    geom_label_repel(
        data = key_data,
        aes(y = event_y, label = paste0(threshold, ": E=", n_events)),
        colour = "#54667A", fill = "#DCE3EC",
        size = 3, direction = "y",
        nudge_y = -max(curve_data$n_clusters) * 0.06,
        seed = 2026
    ) +
    scale_y_continuous(
        name = "Number of clusters",
        sec.axis = sec_axis(
            ~ . * event_scale,
            name = "Number of putative transmission events"
        ),
        expand = expansion(mult = c(0.05, 0.12))
    ) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_colour_manual(
        values = c(Clusters = cluster_color, Events = event_color)
    ) +
    labs(x = "SNP threshold", colour = NULL) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "bottom",
        plot.margin = margin(10, 15, 10, 10)
    )
curve_plot
ggsave(
    file.path(output_dir, "Exploratory_cluster_event_sensitivity.pdf"),
    curve_plot,
    width = 7.5,
    height = 5.5
)

# Strain-sharing pairs composition at selected thresholds
summarise_categories <- function(column, categories) {
    bind_rows(lapply(bubble_thresholds, function(threshold) {
        sub <- pairs %>% filter(SNP_distance <= threshold)
        
        counts <- sub %>%
            count(category = .data[[column]], name = "n")
        
        tibble(category = categories) %>%
            left_join(counts, by = "category") %>%
            mutate(
                threshold = threshold,
                n = replace_na(n, 0L),
                N = nrow(sub),
                percent = if_else(N > 0, 100 * n / N, NA_real_),
                label = if_else(
                    n > 0,
                    sprintf("%.2f%%\n(%d/%d)", percent, n, N),
                    ""
                )
            )
    }))
}

source_data <- summarise_categories("Source_category", source_order)
clone_data <- summarise_categories("Clone_category", clone_order)

# Figure S5C
plot_bubbles <- function(data, categories, show_x = TRUE) {
    data <- data %>%
        mutate(
            threshold = factor(threshold, levels = bubble_thresholds),
            category = factor(category, levels = rev(categories))
        )
    
    positive_data <- data %>% filter(n > 0)
    
    p <- ggplot(data, aes(x = threshold, y = category)) +
        geom_tile(
            width = 1, height = 1,
            fill = "white", colour = "#F2F4F5",
            linewidth = 0.4
        ) +
        geom_point(
            data = positive_data,
            aes(size = percent, fill = percent),
            shape = 21, colour = "black", stroke = 0.4
        ) +
        geom_text(
            data = positive_data,
            aes(label = label),
            nudge_y = -0.3,
            size = 2,
            lineheight = 1,
            colour = "black"
        ) +
        scale_size_continuous(
            range = c(1.5, 7),
            trans = "sqrt",
            limits = c(0, 100),
            guide = "none"
        ) +
        scale_fill_gradient(
            low = "#EAF1F6",
            high = "#5B7E96",
            trans = "sqrt",
            limits = c(0, 100),
            guide = "none"
        ) +
        scale_x_discrete(
            drop = FALSE,
            expand = expansion(add = 0.5)
        ) +
        scale_y_discrete(
            drop = FALSE,
            expand = expansion(add = 0.5)
        ) +
        labs(x = if (show_x) "SNP threshold" else NULL, y = NULL) +
        theme_minimal(base_family = "sans", base_size = 10) +
        theme(
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_text(colour = "black"),
            axis.title.x = element_text(margin = margin(t = 6)),
            plot.margin = margin(4, 8, 4, 6)
        )
    
    if (!show_x) {
        p <- p + theme(axis.text.x = element_blank())
    }
    
    p
}

source_plot <- plot_bubbles(source_data, source_order, show_x = FALSE)
clone_plot <- plot_bubbles(clone_data, clone_order, show_x = TRUE)

bubble_plot <- (source_plot / clone_plot) +
    plot_layout(heights = c(6, 9)) +
    plot_annotation(
        title = "C",
        theme = theme(
            plot.title = element_text(size = 15, face = "plain")
        )
    )
bubble_plot

ggsave(
    file.path(output_dir, "Figure_S5C_pair_sensitivity.pdf"),
    bubble_plot,
    width = 7,
    height = 10,
    bg = "white"
)

# Export calculated values
write_xlsx(
    list(
        cluster_event_curve = curve_data,
        source_pairs = source_data,
        clone_pairs = clone_data
    ),
    file.path(output_dir, "SNP_sensitivity_results.xlsx")
)

print(curve_data %>% filter(threshold %in% bubble_thresholds))

message("Results saved in: ", output_dir)
