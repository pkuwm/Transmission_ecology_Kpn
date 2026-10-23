# Figure 1B: Distribution of CRKP and CSKP across provinces and sample types
rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scatterpie)
library(ggnewscale)

# Select the input Excel file
input_file <- file.choose()

# Save figures in the same folder as the input file
output_dir <- dirname(input_file)
output_pdf <- file.path(
    output_dir,
    "Figure_2B_bubble_pie.pdf"
)

# Read data
df <- read_excel(input_file, sheet = 3)

stopifnot(all(c("Province", "DSampleType", "Site", "CRKP") %in% names(df)))
stopifnot(all(df$CRKP %in% c(0, 1)) & !anyNA(df$CRKP))

# Display order
province_levels <- c(
    "Guangdong", "Zhejiang", "Sichuan", "Henan", "Hebei"
)

sample_levels <- c(
    "Clinical isolates",
    "HCW Stool",
    "Hospital Surface",
    "Hospital Sewage",
    "Hospital Air",
    "Municipal Sewage",
    "River",
    "Community Surface",
    "Community Air",
    "Animal Stool",
    "Farm Sewage",
    "Farm Air"
)

sample_site <- tibble(
    DSampleType = sample_levels,
    Site = c(
        rep("Hospital", 5),
        rep("Community", 4),
        rep("Farm", 3)
    )
)

# Sampling layout shown in the final figure
plot_layout <- crossing(
    Province = province_levels,
    DSampleType = sample_levels
) %>%
    left_join(sample_site, by = "DSampleType") %>%
    filter(
        Site != "Farm" |
            Province %in% c("Hebei", "Henan", "Sichuan") |
            (Province == "Zhejiang" &
                 DSampleType %in% c("Farm Sewage", "Farm Air"))
    )

# Counts for each bubble
bubble_data <- df %>%
    count(
        Province,
        DSampleType,
        Site,
        wt = CRKP,
        name = "CRKP_Count"
    ) %>%
    left_join(
        df %>%
            count(Province, DSampleType, Site, name = "Total"),
        by = c("Province", "DSampleType", "Site")
    )

plot_data <- plot_layout %>%
    left_join(
        bubble_data,
        by = c("Province", "DSampleType", "Site")
    ) %>%
    mutate(
        across(c(CRKP_Count, Total), ~replace_na(.x, 0)),
        CSKP_Count = Total - CRKP_Count,
        x = match(DSampleType, sample_levels),
        y = match(Province, province_levels),
        radius = sqrt(Total / max(df %>% count(Province, DSampleType) %>% pull(n))) * 0.29,
        label = paste0("(", CRKP_Count, "/", Total, ")"),
        
        Hospital_CRKP = if_else(Site == "Hospital", CRKP_Count, 0),
        Hospital_CSKP = if_else(Site == "Hospital", CSKP_Count, 0),
        Community_CRKP = if_else(Site == "Community", CRKP_Count, 0),
        Community_CSKP = if_else(Site == "Community", CSKP_Count, 0),
        Farm_CRKP = if_else(Site == "Farm", CRKP_Count, 0),
        Farm_CSKP = if_else(Site == "Farm", CSKP_Count, 0)
    )

pie_columns <- c(
    "Hospital_CRKP", "Hospital_CSKP",
    "Community_CRKP", "Community_CSKP",
    "Farm_CRKP", "Farm_CSKP"
)

# Province totals
province_totals <- df %>%
    group_by(Province) %>%
    summarise(
        CRKP_Count = sum(CRKP),
        Total = n(),
        .groups = "drop"
    ) %>%
    mutate(
        Province = factor(Province, levels = province_levels),
        axis_label = paste0(Province, "\n(", CRKP_Count, "/", Total, ")")
    )

province_labels <- setNames(
    province_totals$axis_label,
    province_totals$Province
)

# Setting totals
site_totals <- df %>%
    group_by(Site) %>%
    summarise(
        CRKP_Count = sum(CRKP),
        Total = n(),
        .groups = "drop"
    ) %>%
    mutate(
        x = c(Hospital = 3, Community = 7.5, Farm = 11)[Site],
        label = paste0(Site, "\n(", CRKP_Count, "/", Total, ")")
    )

# Dummy layers used only to generate legends
size_legend <- tibble(x = 1, y = 1, Total = c(0, 40, 80, 120, 160))

phenotype_legend <- tibble(
    x = 1,
    y = 1,
    Phenotype = factor(c("CRKP", "CSKP"), levels = c("CRKP", "CSKP"))
)

# Plot
p <- ggplot() +
    geom_hline(
        yintercept = 1:5,
        color = "#DEDEDE",
        linewidth = 0.45
    ) +
    geom_vline(
        xintercept = 1:12,
        color = "#DEDEDE",
        linewidth = 0.45
    ) +
    geom_vline(
        xintercept = c(5.5, 9.5),
        color = "black",
        linewidth = 0.55,
        linetype = "longdash"
    ) +
    geom_scatterpie(
        data = plot_data,
        aes(x = x, y = y, r = radius),
        cols = pie_columns,
        color = "#484848",
        linewidth = 0.35
    ) +
    scale_fill_manual(
        values = c(
            Hospital_CRKP = "#65578F",
            Hospital_CSKP = "#AAA3C1",
            Community_CRKP = "#4F689E",
            Community_CSKP = "#899AC1",
            Farm_CRKP = "#718962",
            Farm_CSKP = "#AFC1A5"
        ),
        guide = "none"
    ) +
    geom_text(
        data = plot_data,
        aes(x = x, y = y - 0.32, label = label),
        size = 3.1,
        color = "black"
    ) +
    geom_text(
        data = site_totals,
        aes(x = x, y = 5.66, label = label),
        size = 4.8,
        lineheight = 0.9
    ) +
    geom_point(
        data = size_legend,
        aes(x = x, y = y, size = Total),
        shape = 1,
        alpha = 0
    ) +
    scale_size_continuous(
        name = "Number",
        range = c(1.5, 12),
        breaks = c(0, 40, 80, 120, 160),
        limits = c(0, 160)
    ) +
    guides(
        size = guide_legend(
            order = 1,
            override.aes = list(alpha = 1, shape = 1, color = "black")
        )
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(
        data = phenotype_legend,
        aes(x = x, y = y, fill = Phenotype),
        shape = 22,
        size = 7,
        alpha = 0
    ) +
    scale_fill_manual(
        name = NULL,
        values = c(CRKP = "#333132", CSKP = "#A6A4A5"),
        guide = guide_legend(
            order = 2,
            direction = "horizontal",
            override.aes = list(alpha = 1, shape = 22, size = 7)
        )
    ) +
    scale_x_continuous(
        breaks = 1:12,
        labels = sample_levels,
        limits = c(0.4, 12.6),
        expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
        breaks = 1:5,
        labels = province_labels[province_levels],
        limits = c(0.45, 6.0),
        expand = expansion(mult = 0)
    ) +
    coord_fixed(clip = "off") +
    labs(
        x = "Sample Type",
        y = "Province"
    ) +
    theme_bw(base_size = 13) +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(
            color = "#3F3F3F",
            linewidth = 0.8
        ),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            color = "black",
            size = 10.5
        ),
        axis.text.y = element_text(
            color = "black",
            size = 11,
            lineheight = 0.9
        ),
        axis.title = element_text(
            color = "black",
            size = 14
        ),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.margin = margin(10, 15, 10, 10)
    )

print(p)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

ggsave(
    filename = output_pdf,
    plot = p,
    width = 13,
    height = 8,
    device = cairo_pdf
)
