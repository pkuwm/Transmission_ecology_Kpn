# Figure 4D: Geographic distance and strain-sharing pair counts
rm(list = ls())
library(dplyr)
library(ggplot2)

# Data
df <- tibble(
    from = c(
        "Hospital A", "Hospital A", "Farm A", "Farm B", "Hospital A",
        "Hospital A", "Hospital B", "Hospital C", "Hospital C",
        "Community A", "Community A", "Community B",
        "Farm B", "Farm A", "Hospital A", "Hospital A",
        "Hospital B", "Hospital A", "Hospital B", "Hospital A", "Hospital B"
    ),
    to = c(
        "Hospital B", "WWTP A", "Farm C", "WWTP C", "Hospital B",
        "Hospital C", "Hospital C", "Community A", "Community B",
        "Community B", "Community C", "Community C",
        "Farm E", "Farm E", "Hospital B", "River A",
        "River A", "Hospital B", "Community B", "WWTP", "WWTP"
    ),
    freq = c(
        63, 1, 3, 2, 133, 476, 413, 6, 3, 3, 1,
        1, 1, 1, 23, 6, 1, 60, 6, 1, 10
    ),
    dist = c(
        4, 3, 3, 6, 10, 6, 5, 8, 7, 13, 7.5,
        9, 9.5, 3, 18, 2, 20, 2.5, 2, 7.5, 4
    )
)

output_dir <- choose.dir(caption = "Select the output folder")
if (is.na(output_dir)) stop("No output folder selected.")

# Classify institution pairs
df <- df %>%
    mutate(
        type = case_when(
            grepl("Hospital", from) & grepl("Hospital", to) ~
                "Hospital-Hospital",
            grepl("Hospital", from) | grepl("Hospital", to) ~
                "Hospital-Other",
            TRUE ~ "Non-hospital"
        ),
        type = factor(
            type,
            levels = c(
                "Hospital-Hospital",
                "Hospital-Other",
                "Non-hospital"
            )
        )
    )

# Negative binomial regression
model <- MASS::glm.nb(freq ~ dist + type, data = df)

print(summary(model))

coefficients <- summary(model)$coefficients
estimate <- coefficients["dist", "Estimate"]
se <- coefficients["dist", "Std. Error"]

distance_result <- data.frame(
    Term = "Distance (per km)",
    IRR = exp(estimate),
    CI_low = exp(estimate - qnorm(0.975) * se),
    CI_high = exp(estimate + qnorm(0.975) * se),
    P_value = coefficients["dist", "Pr(>|z|)"]
)

print(distance_result)

# Descriptive log-log plot
df <- df %>%
    mutate(
        log_dist = log10(dist),
        log_freq = log10(freq)
    )

p_value <- distance_result$P_value

p_label <- if (p_value < 0.001) {
    "P < 0.001"
} else {
    paste0("P = ", format.pval(p_value, digits = 2))
}

annotation <- paste0(
    "Adjusted IRR per km = ",
    sprintf("%.2f", distance_result$IRR),
    "\n", p_label
)

p <- ggplot(df, aes(x = log_dist, y = log_freq)) +
    geom_point(
        size = 4.5,
        alpha = 0.8,
        colour = "#5C8694",
        shape = 16
    ) +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        colour = "#333333",
        linewidth = 1.1
    ) +
    annotate(
        "text",
        x = min(df$log_dist) + 0.05,
        y = max(df$log_freq) - 0.1,
        label = annotation,
        hjust = 0,
        vjust = 1,
        size = 4,
        fontface = "bold"
    ) +
    labs(
        x = expression(log[10] ~ "Geographic distance (km)"),
        y = expression(log[10] ~ "Strain-sharing pair count"),
        caption = paste(
            "Line: descriptive log-log linear fit.",
            "IRR and P: negative binomial model adjusted for institution-pair type."
        )
    ) +
    theme_classic(base_size = 12, base_family = "sans") +
    theme(
        axis.title = element_text(
            size = 13, face = "bold", colour = "black"
        ),
        axis.text = element_text(size = 11, colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        plot.caption = element_text(
            size = 8, hjust = 0, colour = "grey40"
        ),
        plot.margin = margin(10, 15, 10, 10)
    )

print(p)

# Export PDF
ggsave(
    filename = file.path(
        output_dir,
        "Geographic_distance_pair_counts.pdf"
    ),
    plot = p,
    width = 8,
    height = 4.5,
    bg = "white",
    useDingbats = FALSE
)

message("Figure saved in: ", output_dir)