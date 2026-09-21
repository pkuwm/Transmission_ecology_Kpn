# Figure S4: ARG and VFG distributions across five source groups
rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# Settings
input_file <- file.choose()
# input_file <- "data/analysis_results.xlsx"

sheet <- "kpn902元数据"
output_dir <- file.path(dirname(input_file), "Figure_S4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_order <- c(
    "Hospital human",
    "Hospital environment",
    "Community environment",
    "Farm environment",
    "Farm animal"
)

group_labels <- c(
    "Hospital\nhuman",
    "Hospital\nenvironment",
    "Community\nenvironment",
    "Farm\nenvironment",
    "Farm\nanimal"
)

group_colors <- setNames(
    c("#92849E", "#92849E", "#84A0AF", "#A7C6B0", "#A7C6B0"),
    group_order
)

# Read and group isolates
raw <- read_excel(input_file, sheet = sheet)

required <- c("Site", "Niche", "ARGnum", "VFGnum")
missing_columns <- setdiff(required, names(raw))

if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

df <- raw %>%
    mutate(
        Site_key = tolower(trimws(as.character(Site))),
        Niche_key = tolower(trimws(as.character(Niche))),
        Group = case_when(
            Site_key == "hospital" & Niche_key == "human" ~
                "Hospital human",
            Site_key == "hospital" & Niche_key == "environment" ~
                "Hospital environment",
            Site_key == "community" & Niche_key == "environment" ~
                "Community environment",
            Site_key == "farm" & Niche_key == "environment" ~
                "Farm environment",
            Site_key == "farm" & Niche_key == "animal" ~
                "Farm animal",
            TRUE ~ NA_character_
        ),
        Group = factor(Group, levels = group_order)
    )

if (anyNA(df$Group)) {
    print(df %>% filter(is.na(Group)) %>% count(Site, Niche))
    stop("Unmatched Site/Niche combinations.")
}

for (variable in c("ARGnum", "VFGnum")) {
    text <- trimws(as.character(df[[variable]]))
    missing <- is.na(text) | text == ""
    value <- suppressWarnings(as.numeric(text))
    
    if (any(!missing & (!is.finite(value) | value < 0))) {
        stop("Invalid values in ", variable)
    }
    
    value[missing] <- NA_real_
    df[[variable]] <- value
}

# Statistical tests and plotting
make_plot <- function(variable, axis_label) {
    d <- df %>%
        transmute(Group, Value = .data[[variable]]) %>%
        filter(!is.na(Value)) %>%
        mutate(id = as.integer(Group))
    
    counts <- tabulate(d$id, nbins = length(group_order))
    
    if (any(counts < 2)) {
        stop("Each group requires at least two observations: ", variable)
    }
    
    # Two-sided Wilcoxon tests with Holm adjustment
    comparisons <- combn(seq_along(group_order), 2, simplify = FALSE)
    
    tests <- bind_rows(lapply(comparisons, function(ids) {
        a <- d$Value[d$id == ids[1]]
        b <- d$Value[d$id == ids[2]]
        
        p <- if (length(unique(c(a, b))) == 1) {
            NA_real_
        } else {
            wilcox.test(
                a, b,
                alternative = "two.sided",
                paired = FALSE,
                exact = FALSE,
                correct = TRUE
            )$p.value
        }
        
        data.frame(
            i = ids[1],
            j = ids[2],
            Group1 = group_order[ids[1]],
            Group2 = group_order[ids[2]],
            P_raw = p
        )
    })) %>%
        mutate(
            P_adjusted = p.adjust(P_raw, method = "holm", n = 10),
            Significance = case_when(
                is.na(P_adjusted) ~ "NE",
                P_adjusted < 0.001 ~ "***",
                P_adjusted < 0.01 ~ "**",
                P_adjusted < 0.05 ~ "*",
                TRUE ~ "ns"
            )
        )
    
    message(variable, ": ", nrow(d), " nonmissing observations")
    print(tests %>% select(-i, -j))
    
    # Right-sided density shapes, trimmed to the observed range
    violin_data <- bind_rows(lapply(seq_along(group_order), function(i) {
        values <- d$Value[d$id == i]
        
        if (length(unique(values)) < 2) {
            return(NULL)
        }
        
        den <- density(
            values,
            from = min(values),
            to = max(values),
            n = 256
        )
        
        baseline <- i + 0.10
        
        data.frame(
            Group = group_order[i],
            x = c(
                baseline,
                baseline + 0.43 * den$y / max(den$y),
                baseline
            ),
            y = c(min(values), den$x, max(values))
        )
    }))
    
    # Significance brackets
    value_span <- max(diff(range(d$Value)), 1)
    
    annotation <- tests %>%
        arrange(j - i, i) %>%
        mutate(
            y = max(d$Value) +
                value_span * (0.10 + (row_number() - 1) * 0.085)
        )
    
    tip <- value_span * 0.018
    
    p <- ggplot(d, aes(x = id, y = Value))
    
    if (nrow(violin_data) > 0) {
        p <- p +
            geom_polygon(
                data = violin_data,
                aes(x = x, y = y, group = Group, fill = Group),
                inherit.aes = FALSE,
                alpha = 0.75,
                colour = NA
            )
    }
    
    p +
        geom_boxplot(
            aes(group = Group, fill = Group),
            width = 0.12,
            colour = "#414141",
            linewidth = 0.5,
            alpha = 0.95,
            outlier.shape = NA,
            orientation = "x"
        ) +
        geom_segment(
            data = annotation,
            aes(x = i, xend = j, y = y, yend = y),
            inherit.aes = FALSE,
            linewidth = 0.35,
            colour = "#414141"
        ) +
        geom_segment(
            data = annotation,
            aes(x = i, xend = i, y = y, yend = y - tip),
            inherit.aes = FALSE,
            linewidth = 0.35,
            colour = "#414141"
        ) +
        geom_segment(
            data = annotation,
            aes(x = j, xend = j, y = y, yend = y - tip),
            inherit.aes = FALSE,
            linewidth = 0.35,
            colour = "#414141"
        ) +
        geom_text(
            data = annotation,
            aes(
                x = (i + j) / 2,
                y = y + value_span * 0.025,
                label = Significance
            ),
            inherit.aes = FALSE,
            size = 3.5,
            vjust = 0,
            colour = "#303030"
        ) +
        scale_fill_manual(values = group_colors, guide = "none") +
        scale_x_continuous(
            breaks = seq_along(group_order),
            labels = paste0(group_labels, "\n(n = ", counts, ")"),
            limits = c(0.65, 5.65),
            expand = expansion(mult = 0)
        ) +
        scale_y_continuous(
            expand = expansion(mult = c(0.03, 0.07))
        ) +
        labs(x = NULL, y = axis_label) +
        theme_classic(base_size = 12, base_family = "sans") +
        theme(
            axis.line = element_line(
                linewidth = 0.45, colour = "#414141"
            ),
            axis.ticks = element_line(
                linewidth = 0.4, colour = "#414141"
            ),
            axis.text = element_text(colour = "#333333"),
            axis.text.x = element_text(
                size = 10.5,
                lineheight = 0.95,
                margin = margin(t = 8)
            ),
            axis.title.y = element_text(
                face = "bold",
                margin = margin(r = 10)
            ),
            plot.margin = margin(12, 15, 10, 10)
        )
}

# Generate panels
p_arg <- make_plot("ARGnum", "Number of ARGs per isolate")
p_vfg <- make_plot("VFGnum", "Number of VFGs per isolate")

combined <- (p_arg | p_vfg) +
    plot_annotation(
        tag_levels = "A",
        theme = theme(
            plot.tag = element_text(size = 20)
        )
    )

print(combined)

# Export figures
ggsave(
    file.path(output_dir, "Figure_S4_ARG_VFG.pdf"),
    combined,
    width = 14,
    height = 7,
    bg = "white",
    useDingbats = FALSE
)

message("Figures saved in: ", output_dir)