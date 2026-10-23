# Figure S1: Sampling site maps across five provinces
rm(list = ls())
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
library(ggspatial)

# Settings
message("Select the sampling metadata Excel file.")
input_file <- file.choose()

message("Select the county boundary GeoJSON file.")
map_file <- file.choose()

sheet <- "采样点位"
provinces <- c("Hebei", "Henan", "Sichuan", "Guangdong", "Zhejiang")

output_dir <- file.path(dirname(input_file), "Figure_S1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("svglite", quietly = TRUE)) {
    stop('Please install svglite: install.packages("svglite")')
}

# Convert decimal or degree-minute-second coordinates
to_decimal <- function(x) {
    x <- trimws(as.character(x))
    decimal <- suppressWarnings(as.numeric(x))
    
    pattern <- paste0(
        "^\\s*([+-]?[0-9]+(?:\\.[0-9]+)?)\\s*[°º]\\s*",
        "([0-9]+(?:\\.[0-9]+)?)\\s*['′’]\\s*",
        "([0-9]+(?:\\.[0-9]+)?)"
    )
    
    parts <- str_match(x, pattern)
    use_dms <- is.na(decimal) & !is.na(parts[, 1])
    
    degrees <- as.numeric(parts[use_dms, 2])
    minutes <- as.numeric(parts[use_dms, 3])
    seconds <- as.numeric(parts[use_dms, 4])
    
    values <- abs(degrees) + minutes / 60 + seconds / 3600
    
    negative <- degrees < 0 |
        str_detect(toupper(x[use_dms]), "[SW]\\s*$")
    
    values[negative] <- -values[negative]
    values[minutes >= 60 | seconds >= 60] <- NA_real_
    
    decimal[use_dms] <- values
    decimal
}

# Read sampling sites
site_data <- read_excel(input_file, sheet = sheet)
names(site_data) <- trimws(names(site_data))

required_columns <- c(
    "Province", "Longitude", "Latitude", "Source", "Location"
)

missing_columns <- setdiff(required_columns, names(site_data))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

site_data <- site_data %>%
    select(all_of(required_columns)) %>%
    mutate(
        across(
            everything(),
            ~ na_if(trimws(as.character(.x)), "")
        )
    ) %>%
    filter(Province %in% provinces) %>%
    filter(!(
        is.na(Location) &
            is.na(Longitude) &
            is.na(Latitude)
    )) %>%
    mutate(
        Longitude = to_decimal(Longitude),
        Latitude = to_decimal(Latitude)
    )

if (nrow(site_data) == 0) {
    stop("No sampling sites found.")
}

if (any(
    is.na(site_data$Longitude) |
    is.na(site_data$Latitude) |
    abs(site_data$Longitude) > 180 |
    abs(site_data$Latitude) > 90
)) {
    stop("Missing or invalid coordinates for sampling sites.")
}

if (anyNA(site_data$Source)) {
    stop("Missing sampling site types.")
}

# Keep source colors consistent across maps
source_levels <- sort(unique(site_data$Source))
site_data$Source <- factor(site_data$Source, levels = source_levels)

# Read county boundaries
china_county <- st_read(map_file, quiet = TRUE)

if (is.na(st_crs(china_county))) {
    stop("The boundary file has no coordinate reference system.")
}

china_county <- china_county %>%
    st_transform(4326) %>%
    st_make_valid()

# Draw provincial sampling areas
plot_province_sites <- function(province) {
    sites <- site_data %>%
        filter(Province == province)
    
    if (nrow(sites) == 0) {
        warning("No sampling sites found for ", province)
        return(invisible(NULL))
    }
    
    sites_sf <- st_as_sf(
        sites,
        coords = c("Longitude", "Latitude"),
        crs = 4326
    )
    
    bbox <- st_bbox(sites_sf)
    bbox[c("xmin", "ymin")] <- bbox[c("xmin", "ymin")] - 0.2
    bbox[c("xmax", "ymax")] <- bbox[c("xmax", "ymax")] + 0.2
    
    keep <- lengths(
        st_intersects(china_county, st_as_sfc(bbox))
    ) > 0
    
    area_map <- china_county[keep, ]
    
    p <- ggplot() +
        geom_sf(
            data = area_map,
            fill = "white",
            colour = "#B2B2B2",
            linewidth = 0.5
        ) +
        geom_point(
            data = sites,
            aes(x = Longitude, y = Latitude, colour = Source),
            size = 1.5
        ) +
        geom_text(
            data = sites,
            aes(x = Longitude, y = Latitude, label = Location),
            vjust = -0.8,
            size = 1,
            check_overlap = TRUE,
            na.rm = TRUE
        ) +
        scale_colour_brewer(
            palette = "Set1",
            limits = source_levels,
            drop = FALSE,
            name = "Sampling site type"
        ) +
        coord_sf(
            crs = st_crs(4326),
            default_crs = st_crs(4326),
            xlim = unname(bbox[c("xmin", "xmax")]),
            ylim = unname(bbox[c("ymin", "ymax")]),
            expand = FALSE
        ) +
        annotation_scale(
            location = "bl",
            width_hint = 0.3
        ) +
        annotation_north_arrow(
            location = "tr",
            style = north_arrow_fancy_orienteering()
        ) +
        labs(
            title = province,
            x = NULL,
            y = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 14),
            legend.position = "bottom"
        )
    
    ggsave(
        filename = file.path(
            output_dir, paste0("Figure_S1_", province, ".pdf")
        ),
        plot = p,
        width = 8,
        height = 6,
        bg = "white",
        useDingbats = FALSE
    )
    
    ggsave(
        filename = file.path(
            output_dir, paste0("Figure_S1_", province, ".svg")
        ),
        plot = p,
        width = 8,
        height = 6,
        bg = "white",
        device = svglite::svglite
    )
    
    invisible(p)
}

# Generate five maps
plots <- setNames(
    lapply(provinces, plot_province_sites),
    provinces
)

message("Figures saved in: ", normalizePath(output_dir))
