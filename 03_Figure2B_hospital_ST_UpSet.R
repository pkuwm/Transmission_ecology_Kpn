# Figure 2B: ST sharing among hospital sources
rm(list = ls())
library(UpSetR)
library(readxl)
library(dplyr)

# Settings
input_file <- file.choose()

sheet <- "kpn902元数据"
st_column <- "ST绘图"

output_file <- file.path(
    dirname(input_file),
    "Figure_2B_hospital_ST_UpSet.pdf"
)

group_order <- c(
    "Clinical isolates",
    "Hospital Environment",
    "HCW Stool"
)

group_colors <- c(
    "Clinical isolates" = "#BD9296",
    "Hospital Environment" = "#84A0AF",
    "HCW Stool" = "#D8B6B1"
)

# Read data
df <- read_excel(input_file, sheet = sheet)
names(df) <- trimws(names(df))

required_columns <- c("Site", "DSampleType", st_column)
missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

# Define hospital sources
hospital <- df %>%
    mutate(
        Site = trimws(as.character(Site)),
        DSampleType = trimws(as.character(DSampleType)),
        ST_id = trimws(as.character(.data[[st_column]]))
    ) %>%
    filter(tolower(Site) == "hospital") %>%
    mutate(
        Source = case_when(
            tolower(DSampleType) == "clinical isolates" ~
                "Clinical isolates",
            tolower(DSampleType) %in% c(
                "hospital surface",
                "hospital sewage",
                "hospital air"
            ) ~ "Hospital Environment",
            tolower(DSampleType) == "hcw stool" ~
                "HCW Stool",
            TRUE ~ NA_character_
        )
    )

if (nrow(hospital) == 0) {
    stop("No hospital isolates found.")
}

if (anyNA(hospital$Source)) {
    stop(
        "Unmatched hospital sample types: ",
        paste(
            unique(hospital$DSampleType[is.na(hospital$Source)]),
            collapse = ", "
        )
    )
}

if (any(
    is.na(hospital$ST_id) |
    toupper(hospital$ST_id) %in% c("", "NA", "N/A", "NULL")
)) {
    stop("Missing ST identifiers among hospital isolates.")
}

# Build ST sets
st_lists <- setNames(
    lapply(
        group_order,
        function(group) {
            unique(hospital$ST_id[hospital$Source == group])
        }
    ),
    group_order
)

if (any(lengths(st_lists) == 0)) {
    stop(
        "Groups without isolates: ",
        paste(names(st_lists)[lengths(st_lists) == 0], collapse = ", ")
    )
}

upset_data <- UpSetR::fromList(st_lists)
plot_sets <- rev(group_order)

# Draw plot
draw_upset <- function() {
    p <- UpSetR::upset(
        upset_data,
        nsets = 3,
        sets = plot_sets,
        keep.order = TRUE,
        nintersects = NA,
        order.by = "freq",
        decreasing = TRUE,
        mainbar.y.label = "Number of STs in intersection",
        sets.x.label = "Number of STs",
        main.bar.color = "#8EA6B2",
        sets.bar.color = unname(group_colors[plot_sets]),
        matrix.color = "#414141",
        show.numbers = "yes",
        set_size.show = TRUE,
        set_size.numbers_size = 5,
        point.size = 3.2,
        line.size = 1,
        shade.color = "#F0F2F4",
        shade.alpha = 0.6,
        mb.ratio = c(0.65, 0.35),
        text.scale = c(1.2, 1.1, 1.1, 1, 1.15, 1)
    )
    
    print(p)
}

# Export PDF
pdf(
    output_file,
    width = 8,
    height = 6,
    useDingbats = FALSE
)

tryCatch(
    draw_upset(),
    finally = dev.off()
)

if (interactive()) {
    draw_upset()
}
