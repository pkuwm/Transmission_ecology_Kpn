# Figure S9A: Global distribution of K. pneumoniae KL25 isolates
rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# Settings
input_file <- file.choose()
# input_file <- "data/KL25_metadata.xlsx"

sheet <- "NCBI-KL25"

output_file <- file.path(
    dirname(input_file),
    "KL25_global_distribution_map.pdf"
)

# Read data
df <- read_excel(input_file, sheet = sheet)

if (!"Country" %in% names(df)) {
    stop("Missing column: Country")
}

country_counts <- df %>%
    mutate(Country = trimws(as.character(Country))) %>%
    filter(!is.na(Country), Country != "") %>%
    mutate(
        Country = case_when(
            Country %in% c("USA", "US", "United States") ~
                "United States of America",
            Country %in% c("UK", "Great Britain") ~
                "United Kingdom",
            Country == "Czech Republic" ~ "Czechia",
            TRUE ~ Country
        )
    ) %>%
    count(Country, name = "Isolate_count")

# Join country counts to the map
world <- ne_countries(scale = "medium", returnclass = "sf")

unmatched <- setdiff(country_counts$Country, world$name)

if (length(unmatched) > 0) {
    stop(
        "Country names not matched to the map: ",
        paste(unmatched, collapse = ", ")
    )
}

world_data <- world %>%
    left_join(country_counts, by = c("name" = "Country")) %>%
    mutate(
        Isolate_count = coalesce(Isolate_count, 0L),
        Count_group = cut(
            Isolate_count,
            breaks = c(0, 1, 10, 50, 100, 250, 500, Inf),
            labels = c(
                "0", "1–9", "10–49", "50–99",
                "100–249", "250–499", "500+"
            ),
            right = FALSE,
            include.lowest = TRUE
        )
    )

# Label countries with at least 10 isolates
world_labels <- world_data %>%
    filter(Isolate_count >= 10) %>%
    mutate(
        Country_label = if_else(
            name == "United States of America",
            "USA",
            name
        )
    )

map_colors <- c(
    "0" = "#F2F2F2",
    "1–9" = "#9BAEB9",
    "10–49" = "#7E95A4",
    "50–99" = "#6A8392",
    "100–249" = "#597180",
    "250–499" = "#4C5E6C",
    "500+" = "#354652"
)

# Draw map
p <- ggplot(world_data) +
    geom_sf(
        aes(fill = Count_group),
        colour = "grey65",
        linewidth = 0.15
    ) +
    geom_sf_text(
        data = world_labels,
        aes(label = Country_label),
        size = 3,
        colour = "black",
        check_overlap = TRUE
    ) +
    scale_fill_manual(
        values = map_colors,
        drop = FALSE,
        name = "Number of KL25 isolates"
    ) +
    coord_sf(expand = FALSE) +
    theme_void(base_family = "sans") +
    theme(
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(
        fill = guide_legend(
            nrow = 1,
            title.position = "top"
        )
    )

print(p)

# Export PDF
ggsave(
    filename = output_file,
    plot = p,
    width = 14,
    height = 7.5,
    units = "in",
    bg = "white",
    useDingbats = FALSE
)

message("Figure saved to: ", output_file)