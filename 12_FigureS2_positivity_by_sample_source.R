# Figure S2: Positivity by sample source
rm(list = ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

# Settings
input_file <- file.choose()
# input_file <- "data/sample_metadata.xlsx"

sheet <- "样本库"
output_dir <- file.path(dirname(input_file), "Figure_S3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
df <- read_excel(input_file, sheet = sheet)

required_columns <- c("SampleSource", "DSampleSource", "Positive")
missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

df <- df %>%
    mutate(
        SampleSource = trimws(as.character(SampleSource)),
        DSampleSource = trimws(as.character(DSampleSource)),
        Positive = trimws(as.character(Positive))
    ) %>%
    filter(
        SampleSource %in% c("Water", "Air", "Environmental Surface"),
        !is.na(Positive), Positive != "",
        !is.na(DSampleSource), DSampleSource != ""
    )

if (!all(df$Positive %in% c("0", "1"))) {
    stop("Positive must be coded as 0 or 1.")
}

df$Positive <- as.integer(df$Positive)

# Select the test using expected cell counts
choose_test <- function(ct) {
    if (any(rowSums(ct) == 0) || any(colSums(ct) == 0)) {
        return(list(method = "Not applicable", p = NA_real_))
    }
    
    chi <- suppressWarnings(chisq.test(ct, correct = FALSE))
    use_fisher <- any(chi$expected < 1) ||
        mean(chi$expected < 5) > 0.20
    
    if (!use_fisher) {
        return(list(method = "Pearson chi-square test", p = chi$p.value))
    }
    
    if (nrow(ct) == 2 && ncol(ct) == 2) {
        test <- fisher.test(ct, alternative = "two.sided")
        method <- "Fisher exact test"
    } else {
        set.seed(12345)
        test <- fisher.test(ct, simulate.p.value = TRUE, B = 100000)
        method <- "Fisher test (Monte Carlo P value)"
    }
    
    list(method = method, p = test$p.value)
}

significance_label <- function(p) {
    case_when(
        is.na(p) ~ "NA",
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
    )
}

# Calculate positivity and draw each panel
plot_positivity <- function(data, source) {
    summary_data <- data %>%
        filter(SampleSource == source) %>%
        group_by(DSampleSource) %>%
        summarise(
            n_positive = sum(Positive),
            n_total = n(),
            .groups = "drop"
        ) %>%
        mutate(
            n_negative = n_total - n_positive,
            rate = n_positive / n_total
        ) %>%
        arrange(desc(rate), DSampleSource) %>%
        mutate(
            x = row_number(),
            label = sprintf(
                "%.1f%% (%d/%d)",
                rate * 100, n_positive, n_total
            )
        )
    
    if (nrow(summary_data) < 2) {
        stop("At least two subgroups are required for: ", source)
    }
    
    counts <- as.matrix(
        summary_data[, c("n_positive", "n_negative")]
    )
    
    overall <- choose_test(counts)
    
    comparisons <- combn(seq_len(nrow(summary_data)), 2)
    
    pairwise <- bind_rows(lapply(seq_len(ncol(comparisons)), function(k) {
        i <- comparisons[1, k]
        j <- comparisons[2, k]
        test <- choose_test(counts[c(i, j), , drop = FALSE])
        
        data.frame(
            group1 = summary_data$DSampleSource[i],
            group2 = summary_data$DSampleSource[j],
            x1 = i,
            x2 = j,
            method = test$method,
            p_value = test$p,
            symbol = significance_label(test$p)
        )
    }))
    
    message(
        source, ": ", overall$method,
        "; overall P = ", signif(overall$p, 4)
    )
    print(pairwise %>% select(group1, group2, method, p_value, symbol))
    
    max_rate <- max(max(summary_data$rate), 0.05)
    
    pairwise <- pairwise %>%
        mutate(
            y = max_rate * 1.15 +
                (row_number() - 1) * max_rate * 0.10
        )
    
    y_max <- max(pairwise$y) + max_rate * 0.12
    
    palette <- rep_len(
        c("#BBD7DF", "#E3B6A1", "#CDBED6", "#BD9495"),
        nrow(summary_data)
    )
    
    ggplot(summary_data, aes(x = x, y = rate)) +
        geom_col(
            aes(fill = factor(x)),
            width = 0.68,
            colour = "white",
            linewidth = 0.6
        ) +
        geom_text(
            aes(y = rate + max_rate * 0.025, label = label),
            size = 3.2,
            vjust = 0,
            colour = "grey25"
        ) +
        geom_segment(
            data = pairwise,
            aes(x = x1, xend = x2, y = y, yend = y),
            inherit.aes = FALSE,
            linewidth = 0.45,
            colour = "grey45"
        ) +
        geom_segment(
            data = pairwise,
            aes(
                x = x1, xend = x1,
                y = y, yend = y - max_rate * 0.015
            ),
            inherit.aes = FALSE,
            linewidth = 0.45,
            colour = "grey45"
        ) +
        geom_segment(
            data = pairwise,
            aes(
                x = x2, xend = x2,
                y = y, yend = y - max_rate * 0.015
            ),
            inherit.aes = FALSE,
            linewidth = 0.45,
            colour = "grey45"
        ) +
        geom_text(
            data = pairwise,
            aes(
                x = (x1 + x2) / 2,
                y = y + max_rate * 0.014,
                label = symbol
            ),
            inherit.aes = FALSE,
            size = 3.7,
            vjust = 0,
            colour = "grey25"
        ) +
        scale_fill_manual(values = palette, guide = "none") +
        scale_x_continuous(
            breaks = summary_data$x,
            labels = summary_data$DSampleSource,
            expand = expansion(add = 0.6)
        ) +
        scale_y_continuous(
            labels = label_percent(accuracy = 0.1),
            breaks = breaks_pretty(n = 4),
            expand = expansion(mult = c(0, 0))
        ) +
        coord_cartesian(ylim = c(0, y_max)) +
        labs(x = NULL, y = "Positive rate (%)") +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x = element_text(
                angle = 45, hjust = 1,
                colour = "grey30", size = 10.5
            ),
            axis.text.y = element_text(colour = "grey30"),
            axis.title.y = element_text(
                colour = "grey25",
                margin = margin(r = 8)
            ),
            axis.line = element_line(colour = "grey60", linewidth = 0.4),
            axis.ticks = element_line(colour = "grey60", linewidth = 0.3),
            panel.grid.major.y = element_line(
                colour = "grey90",
                linewidth = 0.3,
                linetype = "dashed"
            ),
            plot.margin = margin(10, 12, 8, 8)
        )
}

# Panels A–C
p_water <- plot_positivity(df, "Water")
p_air <- plot_positivity(df, "Air")
p_surface <- plot_positivity(df, "Environmental Surface")

combined <- (p_water | p_air | p_surface) +
    plot_layout(widths = c(4, 3, 2)) +
    plot_annotation(
        tag_levels = "A",
        theme = theme(
            plot.tag = element_text(size = 18, colour = "black")
        )
    )

print(combined)

# Export figures
ggsave(
    file.path(output_dir, "Figure_S3_positivity.pdf"),
    combined,
    width = 14,
    height = 6,
    useDingbats = FALSE
)