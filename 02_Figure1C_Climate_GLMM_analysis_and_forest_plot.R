# # Figure 1C: Climate GLMM analysis and forest plot
rm(list = ls)
library(ggplot2)
library(patchwork)

# Settings
input_file <- file.choose()
# input_file <- "data/metadata.xlsx"

sheet <- 1
sample_type_col <- "Detail SampleType"

exposures <- c(
    "Precipitation", "Temperature", "Velocity", "AQI",
    "PM2.5", "PM10", "CO", "SO2", "NO2", "O3"
)

joint_exposures <- c(
    "Precipitation", "Temperature", "PM2.5", "PM10", "NO2"
)

out_dir <- file.path(
    dirname(input_file),
    paste0("Climate_GLMM_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

control <- lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
)

# Helper functions
logs <- list()

save_log <- function() {
    if (length(logs) > 0) {
        write.csv(
            do.call(rbind, logs),
            file.path(out_dir, "00_run_log.csv"),
            row.names = FALSE
        )
    }
    
    capture.output(
        sessionInfo(),
        file = file.path(out_dir, "00_session_info.txt")
    )
}

run_step <- function(label, fun) {
    warnings <- character()
    messages <- character()
    error <- NULL
    
    result <- tryCatch(
        withCallingHandlers(
            fun(),
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            },
            message = function(m) {
                messages <<- c(messages, conditionMessage(m))
                invokeRestart("muffleMessage")
            }
        ),
        error = function(e) {
            error <<- conditionMessage(e)
            NULL
        }
    )
    
    status <- if (is.null(error)) "completed" else "failed"
    
    logs[[length(logs) + 1L]] <<- data.frame(
        step = label,
        status = status,
        warnings = paste(unique(warnings), collapse = " | "),
        messages = paste(unique(messages), collapse = " | "),
        error = if (is.null(error)) "" else error
    )
    
    save_log()
    message(label, ": ", status)
    
    if (length(warnings) > 0) {
        message("Warnings: ", paste(unique(warnings), collapse = " | "))
    }
    
    if (!is.null(error)) {
        message("Error: ", error)
    }
    
    result
}

clean_text <- function(x) {
    x <- trimws(as.character(x))
    x[toupper(x) %in% c("", "NA", "N/A", "NULL")] <- NA_character_
    x
}

fixed_table <- function(model, ratio_name = "OR") {
    cf <- summary(model)$coefficients
    critical_value <- qnorm(0.975)
    
    tab <- data.frame(
        term = rownames(cf),
        estimate = cf[, 1],
        SE = cf[, 2],
        OR = exp(cf[, 1]),
        CI_low = exp(cf[, 1] - critical_value * cf[, 2]),
        CI_high = exp(cf[, 1] + critical_value * cf[, 2]),
        P_Wald = cf[, 4],
        row.names = NULL
    )
    
    names(tab)[names(tab) == "OR"] <- ratio_name
    tab
}

model_flags <- function(model) {
    flags <- model@optinfo$conv$lme4$messages
    optimizer_code <- model@optinfo$conv$opt
    
    if (length(optimizer_code) > 0 && any(optimizer_code != 0)) {
        flags <- c(
            flags,
            paste0("optimizer_code=", paste(optimizer_code, collapse = ","))
        )
    }
    
    if (lme4::isSingular(model, tol = 1e-4)) {
        flags <- c(flags, "singular_fit")
    }
    
    paste(flags, collapse = " | ")
}

validate_model_data <- function(d, variables, groups) {
    if (nrow(d) < 2L || length(unique(d$positive)) != 2L) {
        stop("Model data must contain both outcome classes.")
    }
    
    for (variable in variables) {
        variable_sd <- sd(d[[variable]])
        
        if (!is.finite(variable_sd) || variable_sd <= 0) {
            stop("Invalid or zero variance: ", variable)
        }
    }
    
    for (group in groups) {
        n_groups <- length(unique(d[[group]]))
        
        if (n_groups < 2L || n_groups >= nrow(d)) {
            stop("Invalid number of random-effect levels: ", group)
        }
    }
}

# Read and validate data
dat <- as.data.frame(
    readxl::read_excel(input_file, sheet = sheet)
)

names(dat) <- trimws(names(dat))

if (anyDuplicated(names(dat))) {
    stop("Duplicate column names after trimming whitespace.")
}

required <- unique(c(
    "Positive", "Province", "Season",
    sample_type_col, exposures
))

missing_columns <- setdiff(required, names(dat))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

for (variable in required) {
    dat[[variable]] <- clean_text(dat[[variable]])
}

outcome <- tolower(dat$Positive)

unknown_outcomes <- setdiff(
    unique(na.omit(outcome)),
    c("yes", "no", "1", "0", "positive", "negative")
)

if (length(unknown_outcomes) > 0) {
    stop(
        "Unexpected Positive values: ",
        paste(unknown_outcomes, collapse = ", ")
    )
}

dat$positive <- NA_integer_
dat$positive[outcome %in% c("yes", "1", "positive")] <- 1L
dat$positive[outcome %in% c("no", "0", "negative")] <- 0L

for (variable in exposures) {
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

season_value <- tolower(dat$Season)
season_value[season_value %in% "fall"] <- "autumn"

season_levels <- c("spring", "summer", "autumn", "winter")
unknown_seasons <- setdiff(
    unique(na.omit(season_value)),
    season_levels
)

if (length(unknown_seasons) > 0) {
    stop(
        "Unexpected Season values: ",
        paste(unknown_seasons, collapse = ", ")
    )
}

dat$province <- factor(dat$Province)
dat$season <- factor(season_value, levels = season_levels)
dat$sample_type <- factor(dat[[sample_type_col]])

dat$province_season <- interaction(
    dat$province,
    dat$season,
    drop = TRUE,
    sep = "_"
)

data_summary <- data.frame(
    total_samples = nrow(dat),
    positive_samples = sum(dat$positive == 1, na.rm = TRUE),
    negative_samples = sum(dat$positive == 0, na.rm = TRUE),
    missing_outcomes = sum(is.na(dat$positive)),
    sample_types = nlevels(dat$sample_type),
    province_season_groups = nlevels(dat$province_season)
)

missing_summary <- data.frame(
    variable = c(
        "positive", "sample_type", "province_season", exposures
    ),
    missing_n = vapply(
        dat[c("positive", "sample_type", "province_season", exposures)],
        function(x) sum(is.na(x)),
        integer(1)
    ),
    row.names = NULL
)

writexl::write_xlsx(
    list(
        summary = data_summary,
        missing_values = missing_summary
    ),
    file.path(out_dir, "01_data_summary.xlsx")
)

# Univariable screening
screen_results <- list()
screen_scaling <- list()

for (variable in exposures) {
    result <- run_step(
        paste0("screen_", variable),
        function() {
            columns <- c("positive", "province_season", variable)
            
            ds <- droplevels(
                dat[complete.cases(dat[, columns]), , drop = FALSE]
            )
            
            validate_model_data(
                ds,
                variables = variable,
                groups = "province_season"
            )
            
            exposure_mean <- mean(ds[[variable]])
            exposure_sd <- sd(ds[[variable]])
            
            ds$exposure_z <- (
                ds[[variable]] - exposure_mean
            ) / exposure_sd
            
            model <- lme4::glmer(
                positive ~ exposure_z + (1 | province_season),
                data = ds,
                family = binomial(),
                control = control
            )
            
            tab <- fixed_table(model)
            tab <- tab[tab$term == "exposure_z", , drop = FALSE]
            
            tab <- data.frame(
                exposure = variable,
                n_samples = nrow(ds),
                positive_n = sum(ds$positive),
                tab,
                singular = lme4::isSingular(model, tol = 1e-4),
                flags = model_flags(model),
                row.names = NULL
            )
            
            list(
                coefficients = tab,
                scaling = data.frame(
                    exposure = variable,
                    mean = exposure_mean,
                    SD = exposure_sd
                )
            )
        }
    )
    
    if (!is.null(result)) {
        screen_results[[variable]] <- result$coefficients
        screen_scaling[[variable]] <- result$scaling
    }
}

if (length(screen_results) > 0) {
    screen_table <- do.call(rbind, screen_results)
    screen_scaling_table <- do.call(rbind, screen_scaling)
    rownames(screen_table) <- NULL
    rownames(screen_scaling_table) <- NULL
    
    writexl::write_xlsx(
        list(
            screening = screen_table,
            scaling = screen_scaling_table
        ),
        file.path(out_dir, "02_univariable_screening.xlsx")
    )
} else {
    warning("All univariable models failed. See the run log.")
}

# Multivariable GLMM
model_variables <- c(
    "positive", "sample_type", "province_season",
    joint_exposures
)

d <- droplevels(
    dat[complete.cases(dat[, model_variables]), , drop = FALSE]
)

validate_model_data(
    d,
    variables = joint_exposures,
    groups = c("sample_type", "province_season")
)

scaling_table <- data.frame(
    exposure = joint_exposures,
    mean = vapply(d[joint_exposures], mean, numeric(1)),
    SD = vapply(d[joint_exposures], sd, numeric(1)),
    row.names = NULL
)

for (i in seq_along(joint_exposures)) {
    variable <- joint_exposures[i]
    
    d[[paste0(variable, "_z")]] <- (
        d[[variable]] - scaling_table$mean[i]
    ) / scaling_table$SD[i]
}

z_names <- paste0(joint_exposures, "_z")

joint_formula <- as.formula(
    paste(
        "positive ~",
        paste(z_names, collapse = " + "),
        "+ (1 | sample_type) + (1 | province_season)"
    )
)

message(
    "Multivariable model: ",
    paste(deparse(joint_formula), collapse = " ")
)

fit <- run_step(
    "multivariable_GLMM",
    function() {
        lme4::glmer(
            joint_formula,
            data = d,
            family = binomial(),
            control = control
        )
    }
)

if (is.null(fit)) {
    stop("The multivariable model failed. See the run log.")
}

saveRDS(
    fit,
    file.path(out_dir, "03_final_model.rds")
)

coef_table <- fixed_table(fit, ratio_name = "aOR")
random_effects <- as.data.frame(lme4::VarCorr(fit))

model_information <- data.frame(
    formula = paste(deparse(joint_formula), collapse = " "),
    n_samples = nrow(d),
    excluded_samples = nrow(dat) - nrow(d),
    positive_n = sum(d$positive),
    sample_type_groups = nlevels(d$sample_type),
    province_season_groups = nlevels(d$province_season),
    singular = lme4::isSingular(fit, tol = 1e-4),
    flags = model_flags(fit),
    confidence_interval = "95% Wald CI",
    exposure_unit = "One standard deviation",
    lme4_version = as.character(packageVersion("lme4"))
)

# Collinearity
vif_results <- run_step(
    "collinearity",
    function() {
        as.data.frame(performance::check_collinearity(fit))
    }
)

result_sheets <- list(
    coefficients = coef_table,
    random_effects = random_effects,
    scaling = scaling_table,
    model_information = model_information
)

if (!is.null(vif_results)) {
    result_sheets$VIF <- vif_results
}

writexl::write_xlsx(
    result_sheets,
    file.path(out_dir, "04_final_results.xlsx")
)

# DHARMa diagnostics
sim <- run_step(
    "DHARMa_simulation",
    function() {
        DHARMa::simulateResiduals(
            fittedModel = fit,
            n = 1000,
            refit = FALSE,
            seed = 2026
        )
    }
)

if (!is.null(sim)) {
    saveRDS(
        sim,
        file.path(out_dir, "05_DHARMa_residuals.rds")
    )
    
    run_step(
        "DHARMa_plots",
        function() {
            pdf(
                file.path(out_dir, "05_DHARMa_diagnostics.pdf"),
                width = 8,
                height = 6
            )
            on.exit(dev.off(), add = TRUE)
            
            plot(sim)
            invisible(TRUE)
        }
    )
    
    diagnostic_tests <- list(
        uniformity = run_step(
            "DHARMa_uniformity",
            function() DHARMa::testUniformity(sim, plot = FALSE)
        ),
        dispersion = run_step(
            "DHARMa_dispersion",
            function() DHARMa::testDispersion(sim, plot = FALSE)
        ),
        sample_type = run_step(
            "DHARMa_sample_type",
            function() {
                DHARMa::testCategorical(
                    sim,
                    catPred = d$sample_type,
                    plot = FALSE
                )
            }
        )
    )
    
    saveRDS(
        diagnostic_tests,
        file.path(out_dir, "05_DHARMa_tests.rds")
    )
    
    capture.output(
        diagnostic_tests,
        file = file.path(out_dir, "05_DHARMa_tests.txt")
    )
}

# Forest plot
plot_data <- coef_table[
    match(z_names, coef_table$term),
    ,
    drop = FALSE
]

plot_data$exposure <- joint_exposures
plot_data$y <- rev(seq_along(joint_exposures))

if (
    any(!is.finite(as.matrix(
        plot_data[, c("aOR", "CI_low", "CI_high")]
    ))) ||
    any(plot_data$CI_low <= 0)
) {
    stop("Invalid estimates or confidence intervals for the forest plot.")
}

blue <- "#39758C"
red <- "#B66560"
ink <- "#30383C"

plot_data$direction <- ifelse(
    plot_data$aOR > 1, "Higher", "Lower"
)

plot_data$excludes_one <- (
    plot_data$CI_low > 1 | plot_data$CI_high < 1
)

plot_data$point_fill <- ifelse(
    plot_data$excludes_one,
    ifelse(plot_data$direction == "Higher", red, blue),
    "white"
)

plot_data$value_label <- sprintf(
    "%.2f (%.2f–%.2f)",
    plot_data$aOR,
    plot_data$CI_low,
    plot_data$CI_high
)

n_terms <- nrow(plot_data)
data_top <- n_terms + 0.35
header_y <- n_terms + 0.75
rule_y <- n_terms + 0.38
y_limits <- c(0.45, n_terms + 1.05)

label_map <- c(
    Precipitation = "Precipitation",
    Temperature = "Temperature",
    "PM2.5" = "PM[2.5]",
    PM10 = "PM[10]",
    NO2 = "NO[2]"
)

term_labels <- parse(
    text = unname(label_map[rev(joint_exposures)])
)

# Expand the axis if any confidence interval exceeds 0.5–2
log_span <- max(
    log(2),
    abs(log(c(plot_data$CI_low, plot_data$CI_high))) * 1.05
)

x_limits <- exp(c(-log_span, log_span))

if (log_span <= log(2) + 1e-8) {
    x_breaks <- c(0.5, 0.75, 1, 1.5, 2)
} else {
    x_breaks <- sort(unique(c(
        exp(pretty(log(x_limits), n = 5)),
        1
    )))
    x_breaks <- x_breaks[
        x_breaks >= x_limits[1] & x_breaks <= x_limits[2]
    ]
}

p_forest <- ggplot(plot_data, aes(x = aOR, y = y)) +
    annotate(
        "rect",
        xmin = x_limits[1], xmax = 1,
        ymin = 0.55, ymax = data_top,
        fill = "#EEF4F6", colour = NA
    ) +
    annotate(
        "rect",
        xmin = 1, xmax = x_limits[2],
        ymin = 0.55, ymax = data_top,
        fill = "#FAF0EE", colour = NA
    ) +
    geom_hline(
        yintercept = seq_len(n_terms - 1) + 0.5,
        colour = "white",
        linewidth = 0.5
    ) +
    annotate(
        "segment",
        x = 1, xend = 1,
        y = 0.55, yend = data_top,
        colour = "#899397",
        linetype = "22",
        linewidth = 0.45
    ) +
    geom_segment(
        aes(
            x = CI_low, xend = CI_high,
            yend = y, colour = direction
        ),
        linewidth = 0.65,
        lineend = "round"
    ) +
    geom_segment(
        aes(
            x = CI_low, xend = CI_low,
            y = y - 0.055, yend = y + 0.055,
            colour = direction
        ),
        linewidth = 0.55
    ) +
    geom_segment(
        aes(
            x = CI_high, xend = CI_high,
            y = y - 0.055, yend = y + 0.055,
            colour = direction
        ),
        linewidth = 0.55
    ) +
    geom_point(
        aes(colour = direction, fill = point_fill),
        shape = 21,
        size = 2.9,
        stroke = 0.8
    ) +
    annotate(
        "text",
        x = sqrt(x_limits[1]), y = header_y,
        label = "Lower odds",
        colour = blue,
        size = 3.3
    ) +
    annotate(
        "text",
        x = sqrt(x_limits[2]), y = header_y,
        label = "Higher odds",
        colour = red,
        size = 3.3
    ) +
    scale_colour_manual(
        values = c(Higher = red, Lower = blue),
        guide = "none"
    ) +
    scale_fill_identity() +
    scale_x_log10(
        limits = x_limits,
        breaks = x_breaks,
        labels = function(x) format(signif(x, 3), trim = TRUE),
        expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
        limits = y_limits,
        breaks = seq_len(n_terms),
        labels = term_labels,
        expand = expansion(mult = 0)
    ) +
    labs(x = "Adjusted odds ratio", y = NULL) +
    theme_classic(base_size = 11, base_family = "sans") +
    theme(
        axis.line.y = element_blank(),
        axis.line.x = element_line(
            colour = "#92999C", linewidth = 0.35
        ),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(
            colour = "#92999C", linewidth = 0.35
        ),
        axis.ticks.length = grid::unit(2, "mm"),
        axis.text.y = element_text(
            colour = ink, size = 11,
            margin = margin(r = 12)
        ),
        axis.text.x = element_text(
            colour = ink, size = 9.5,
            margin = margin(t = 5)
        ),
        axis.title.x = element_text(
            colour = ink, size = 11,
            margin = margin(t = 9)
        ),
        plot.margin = margin(8, 8, 6, 8)
    )

p_values <- ggplot(plot_data, aes(y = y)) +
    annotate(
        "text",
        x = 0, y = header_y,
        label = "aOR (95% CI)",
        hjust = 0,
        size = 3.3,
        fontface = "bold",
        colour = ink
    ) +
    annotate(
        "segment",
        x = 0, xend = 1,
        y = rule_y, yend = rule_y,
        colour = "#D5DADD",
        linewidth = 0.4
    ) +
    geom_text(
        aes(x = 0, label = value_label),
        hjust = 0,
        size = 3.3,
        colour = ink
    ) +
    scale_x_continuous(
        limits = c(0, 1.08),
        expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
        limits = y_limits,
        expand = expansion(mult = 0)
    ) +
    theme_void(base_family = "sans") +
    theme(
        plot.margin = margin(8, 8, 6, 8)
    )

forest_plot <- (p_forest | p_values) +
    plot_layout(widths = c(3.6, 1.45)) +
    plot_annotation(
        caption = paste(
            "Effects per 1-SD increase.",
            "Filled markers indicate 95% CIs excluding 1."
        ),
        theme = theme(
            plot.caption = element_text(
                size = 8,
                colour = "#697378",
                hjust = 0,
                margin = margin(t = 5)
            ),
            plot.background = element_rect(
                fill = "white", colour = NA
            )
        )
    )

print(forest_plot)

ggsave(
    file.path(out_dir, "Figure_1C_climate_forest_plot.pdf"),
    plot = forest_plot,
    width = 5.5,
    height = 2.8,
    units = "in",
    bg = "white"
)

save_log()
print(coef_table)
message("Results saved in: ", out_dir)
