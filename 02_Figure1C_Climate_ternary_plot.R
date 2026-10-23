# Figure 1D: Climate ternary plot of sample positivity
rm(list = ls())
library(dplyr)
library(ggtern)

# Settings
input_file <- file.choose()
# input_file <- "data/sample_metadata.xlsx"

sheet <- "样本库"
precipitation_col <- "Rainfall"

out_dir <- file.path(
    dirname(input_file),
    paste0("Climate_ternary_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

province_levels <- c(
    "Hebei", "Henan", "Sichuan", "Zhejiang", "Guangdong"
)

season_levels <- c("spring", "summer", "autumn", "winter")

# Read data
dat <- as.data.frame(
    readxl::read_excel(input_file, sheet = sheet)
)

names(dat) <- trimws(names(dat))

required_columns <- c(
    "Province", "Season", "Positive",
    precipitation_col, "Temperature", "NO2"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

clean_text <- function(x) {
    x <- trimws(as.character(x))
    x[toupper(x) %in% c("", "NA", "N/A", "NULL")] <- NA_character_
    x
}

for (variable in required_columns) {
    dat[[variable]] <- clean_text(dat[[variable]])
}

dat$Precipitation <- dat[[precipitation_col]]

# Validate outcome and grouping variables
outcome <- tolower(dat$Positive)

invalid_outcomes <- setdiff(
    unique(na.omit(outcome)),
    c("yes", "no", "1", "0", "positive", "negative")
)

if (length(invalid_outcomes) > 0) {
    stop(
        "Unexpected Positive values: ",
        paste(invalid_outcomes, collapse = ", ")
    )
}

dat$positive <- NA_integer_
dat$positive[outcome %in% c("yes", "1", "positive")] <- 1L
dat$positive[outcome %in% c("no", "0", "negative")] <- 0L

dat$Season <- tolower(dat$Season)
dat$Season[dat$Season %in% "fall"] <- "autumn"

invalid_provinces <- setdiff(
    unique(na.omit(dat$Province)),
    province_levels
)

invalid_seasons <- setdiff(
    unique(na.omit(dat$Season)),
    season_levels
)

if (length(invalid_provinces) > 0) {
    stop(
        "Unexpected Province values: ",
        paste(invalid_provinces, collapse = ", ")
    )
}

if (length(invalid_seasons) > 0) {
    stop(
        "Unexpected Season values: ",
        paste(invalid_seasons, collapse = ", ")
    )
}

if (anyNA(dat[, c("Province", "Season", "positive")])) {
    stop("Missing Province, Season, or Positive values. Please check the input.")
}

dat$Province <- factor(dat$Province, levels = province_levels)
dat$Season <- factor(dat$Season, levels = season_levels)

climate_variables <- c("Precipitation", "Temperature", "NO2")

for (variable in climate_variables) {
    original <- dat[[variable]]
    numeric_value <- suppressWarnings(as.numeric(original))
    
    invalid <- !is.na(original) &
        (is.na(numeric_value) | !is.finite(numeric_value))
    
    if (any(invalid)) {
        stop(
            "Invalid numeric values in ", variable, ": ",
            paste(unique(original[invalid]), collapse = ", ")
        )
    }
    
    dat[[variable]] <- numeric_value
}

# Province-season summaries
mean_available <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

summary_data <- dat %>%
    group_by(Province, Season) %>%
    summarise(
        n_total = n(),
        n_positive = sum(positive),
        PositiveRate_pct = 100 * n_positive / n_total,
        across(
            all_of(climate_variables),
            mean_available
        ),
        .groups = "drop"
    )

ternary_data <- summary_data %>%
    filter(
        if_all(
            all_of(climate_variables),
            ~ is.finite(.x)
        )
    )

if (nrow(ternary_data) == 0) {
    stop("No province-season groups have complete climate means.")
}

if (nrow(ternary_data) < nrow(summary_data)) {
    warning(
        nrow(summary_data) - nrow(ternary_data),
        " groups excluded because climate means were missing."
    )
}

if (any(as.matrix(ternary_data[climate_variables]) < 0)) {
    stop(
        "Negative climate means are incompatible with this ",
        "maximum-scaling ternary transformation."
    )
}

# Maximum scaling followed by closure to a sum of one
climate_maxima <- vapply(
    ternary_data[climate_variables],
    max,
    numeric(1)
)

if (any(climate_maxima <= 0)) {
    stop("Each climate variable must have a positive maximum.")
}

ternary_data <- ternary_data %>%
    mutate(
        precipitation_rel = Precipitation /
            climate_maxima[["Precipitation"]],
        temperature_rel = Temperature /
            climate_maxima[["Temperature"]],
        no2_rel = NO2 / climate_maxima[["NO2"]],
        total_rel = precipitation_rel + temperature_rel + no2_rel
    )

if (any(ternary_data$total_rel <= 0)) {
    stop("A group has zero total relative exposure.")
}

ternary_data <- ternary_data %>%
    mutate(
        precipitation_tern = precipitation_rel / total_rel,
        temperature_tern = temperature_rel / total_rel,
        no2_tern = no2_rel / total_rel
    )

scaling_parameters <- data.frame(
    variable = names(climate_maxima),
    maximum = unname(climate_maxima)
)

writexl::write_xlsx(
    list(
        province_season_summary = as.data.frame(summary_data),
        ternary_coordinates = as.data.frame(ternary_data),
        scaling_parameters = scaling_parameters
    ),
    file.path(out_dir, "Figure_1D_ternary_data.xlsx")
)

# Plot
axis_breaks <- seq(0, 1, by = 0.2)

p <- ggtern(
    data = ternary_data,
    aes(
        x = no2_tern,
        y = temperature_tern,
        z = precipitation_tern
    )
) +
    geom_point(
        aes(fill = PositiveRate_pct),
        shape = 21,
        size = 4,
        stroke = 0.8,
        colour = "black"
    ) +
    theme_rgbw() +
    theme_showarrows() +
    labs(
        x = expression(NO[2]),
        y = "Temperature",
        z = "Precipitation",
        fill = "Positivity (%)",
        caption = paste(
            "Coordinates represent closed proportions of",
            "maximum-scaled climate variables."
        )
    ) +
    scale_T_continuous(
        limits = c(0, 1),
        breaks = axis_breaks,
        labels = sprintf("%.1f", axis_breaks)
    ) +
    scale_L_continuous(
        limits = c(0, 1),
        breaks = axis_breaks,
        labels = sprintf("%.1f", axis_breaks)
    ) +
    scale_R_continuous(
        limits = c(0, 1),
        breaks = axis_breaks,
        labels = sprintf("%.1f", axis_breaks)
    ) +
    scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1
    ) +
    theme(
        text = element_text(family = "sans"),
        tern.panel.grid.major = element_line(
            linewidth = 0.7,
            linetype = "dashed"
        ),
        tern.panel.grid.minor = element_line(
            colour = "grey90",
            linewidth = 0.3
        ),
        tern.axis.line = element_line(linewidth = 0.7),
        tern.axis.arrow = element_line(linewidth = 0.8),
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.caption = element_text(
            size = 9,
            colour = "grey40",
            hjust = 0
        ),
        plot.margin = margin(15, 15, 15, 15)
    )

print(p)

ggsave(
    filename = file.path(out_dir, "Figure_1D_climate_ternary.pdf"),
    plot = p,
    width = 8,
    height = 6,
    units = "in",
    bg = "white"
)
