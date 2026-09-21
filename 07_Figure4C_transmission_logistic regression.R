# Figure 4C: Screening, mixed-effects logistic regression and forest plot
rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(patchwork)
library(writexl)

# Settings
input_file <- file.choose()
sheet <- "大-回归分析原始数据表"

output_dir <- file.path(
    dirname(input_file),
    paste0("Figure_4C_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.1
candidate_vars <- c("Region", "StrainType", "Niche", "ST")

niche_levels <- c(
    "Farm animal",
    "Hospital human",
    "Hospital environment",
    "Community environment",
    "Farm environment"
)

# Read and prepare data
raw <- read_excel(input_file, sheet = sheet)

dat <- raw %>%
    transmute(
        outcome = as.character(`传播等级二分类`),
        Region = as.character(Region),
        StrainType = as.character(`菌株类型`),
        Site = tolower(trimws(as.character(Site))),
        OriginalNiche = tolower(trimws(as.character(Niche))),
        ST = as.character(ST)
    ) %>%
    mutate(across(everything(), ~ na_if(trimws(.x), "")))

if (
    any(!na.omit(dat$outcome) %in% c("0", "1")) ||
    any(!na.omit(dat$Region) %in% c("South", "North")) ||
    any(!na.omit(dat$StrainType) %in% c("1", "2", "3", "4"))
) {
    stop("Unexpected outcome, Region or StrainType values.")
}

dat <- dat %>%
    mutate(
        Niche = case_when(
            Site == "hospital" & OriginalNiche == "human" ~
                "Hospital human",
            Site == "hospital" & OriginalNiche == "environment" ~
                "Hospital environment",
            Site == "community" & OriginalNiche == "environment" ~
                "Community environment",
            Site == "farm" & OriginalNiche == "environment" ~
                "Farm environment",
            Site == "farm" & OriginalNiche == "animal" ~
                "Farm animal",
            TRUE ~ NA_character_
        )
    )

if (any(
    !is.na(dat$Site) & !is.na(dat$OriginalNiche) & is.na(dat$Niche)
)) {
    stop("Unmatched Site/Niche combinations.")
}

df <- dat %>%
    dplyr::select(-OriginalNiche) %>%
    drop_na() %>%
    mutate(
        outcome = as.integer(outcome),
        Region = factor(Region, levels = c("South", "North")),
        StrainType = factor(
            StrainType,
            levels = c("1", "2", "3", "4"),
            labels = c("cKP", "CRKP", "hvKP", "hv-CRKP")
        ),
        Niche = factor(Niche, levels = niche_levels),
        Site = factor(Site, levels = c("community", "farm", "hospital")),
        ST = factor(ST)
    ) %>%
    droplevels()

if (n_distinct(df$outcome) != 2) {
    stop("Both outcome classes are required.")
}

references <- c(
    Region = "South",
    StrainType = "cKP",
    Niche = "Farm animal",
    ST = if ("Others" %in% levels(df$ST)) "Others" else levels(df$ST)[1]
)

for (v in candidate_vars) {
    if (nlevels(df[[v]]) < 2 || !references[[v]] %in% levels(df[[v]])) {
        stop("Insufficient levels or missing reference: ", v)
    }
    
    df[[v]] <- relevel(df[[v]], ref = references[[v]])
    contrasts(df[[v]]) <- contr.treatment(levels(df[[v]]), base = 1)
}

reference_table <- data.frame(
    Variable = names(references),
    Reference = unname(references)
)

# Map coefficients to factor levels
term_map <- bind_rows(lapply(candidate_vars, function(v) {
    mm <- model.matrix(reformulate(v), data = df)
    
    data.frame(
        Term = colnames(mm)[-1],
        Variable = v,
        Level = levels(df[[v]])[-1],
        Reference = levels(df[[v]])[1]
    )
}))

# Overall association tests
global_results <- bind_rows(lapply(candidate_vars, function(v) {
    tab <- table(df[[v]], df$outcome)
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
    
    use_fisher <- any(chi$expected < 1) ||
        mean(chi$expected < 5) > 0.20
    
    test <- if (use_fisher) {
        fisher.test(tab, workspace = 2e7)
    } else {
        chi
    }
    
    data.frame(
        Variable = v,
        Method = if (use_fisher) "Fisher exact test" else "Pearson chi-square",
        P_value = test$p.value
    )
}))

screened_vars <- global_results$Variable[
    global_results$P_value < alpha
]

# Univariable logistic regression
univariable <- bind_rows(lapply(screened_vars, function(v) {
    fit <- glm(
        reformulate(v, response = "outcome"),
        data = df,
        family = binomial()
    )
    
    cf <- summary(fit)$coefficients
    
    ci <- tryCatch(
        suppressMessages(confint(fit)),
        error = function(e) {
            warning("Profile confidence intervals failed for ", v)
            matrix(
                NA_real_, nrow(cf), 2,
                dimnames = list(rownames(cf), NULL)
            )
        }
    )
    
    data.frame(
        Term = rownames(cf),
        Estimate = cf[, 1],
        SE = cf[, 2],
        OR = exp(cf[, 1]),
        Lower = exp(ci[rownames(cf), 1]),
        Upper = exp(ci[rownames(cf), 2]),
        P_value = cf[, 4],
        row.names = NULL
    ) %>%
        filter(Term != "(Intercept)") %>%
        left_join(term_map, by = "Term")
}))

write_xlsx(
    list(
        references = reference_table,
        global_tests = global_results,
        univariable = univariable
    ),
    file.path(output_dir, "01_screening_results.xlsx")
)

if (nrow(univariable) == 0) {
    stop("No variables passed the overall association tests.")
}

selected_vars <- univariable %>%
    filter(!is.na(P_value), P_value < alpha) %>%
    pull(Variable) %>%
    unique()

if (length(selected_vars) == 0) {
    stop("No variables passed univariable screening.")
}

# Retain the whole factor if any nonreference level passes screening
model_formula <- as.formula(paste(
    "outcome ~",
    paste(selected_vars, collapse = " + "),
    "+ (1 | Region/Site)"
))

fixed_matrix <- model.matrix(reformulate(selected_vars), data = df)

if (qr(fixed_matrix)$rank < ncol(fixed_matrix)) {
    stop("The fixed-effects design matrix is rank deficient.")
}

message("Final formula: ", paste(deparse(model_formula), collapse = " "))

# Frequentist mixed-effects model
fit_warnings <- character()

fit_freq <- tryCatch(
    withCallingHandlers(
        glmer(
            model_formula,
            data = df,
            family = binomial(),
            control = glmerControl(
                optimizer = "bobyqa",
                optCtrl = list(maxfun = 5e5)
            )
        ),
        warning = function(w) {
            fit_warnings <<- c(fit_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    ),
    error = function(e) e
)

freq_ok <- FALSE

if (inherits(fit_freq, "error")) {
    fit_notes <- conditionMessage(fit_freq)
} else {
    fit_notes <- c(
        fit_warnings,
        fit_freq@optinfo$conv$lme4$messages,
        if (isSingular(fit_freq)) "Singular fit"
    )
    
    optimizer_ok <- all(fit_freq@optinfo$conv$opt == 0)
    freq_ok <- length(fit_notes) == 0 && optimizer_ok
    
    if (!optimizer_ok) {
        fit_notes <- c(fit_notes, "Non-zero optimizer code")
    }
    
    saveRDS(fit_freq, file.path(output_dir, "02_frequentist_model.rds"))
}

writeLines(
    c(paste("Accepted:", freq_ok), fit_notes),
    file.path(output_dir, "02_frequentist_diagnostics.txt")
)

# Bayesian fallback
if (freq_ok) {
    cf <- summary(fit_freq)$coefficients
    
    result <- data.frame(
        Term = rownames(cf),
        Estimate = cf[, 1],
        SE = cf[, 2],
        aOR = exp(cf[, 1]),
        Lower = exp(cf[, 1] - 1.96 * cf[, 2]),
        Upper = exp(cf[, 1] + 1.96 * cf[, 2]),
        P_value = cf[, 4],
        row.names = NULL
    ) %>%
        filter(Term != "(Intercept)")
    
    interval <- "95% CI"
    model_type <- "Frequentist mixed-effects logistic regression"
    diagnostics_passed <- TRUE
    
} else {
    fit_bayes <- brms::brm(
        model_formula,
        data = df,
        family = brms::bernoulli(link = "logit"),
        prior = brms::set_prior("normal(0, 1)", class = "b"),
        chains = 4,
        cores = 4,
        iter = 2000,
        warmup = 1000,
        seed = 123,
        control = list(adapt_delta = 0.95, max_treedepth = 10),
        refresh = 100
    )
    
    saveRDS(fit_bayes, file.path(output_dir, "03_bayesian_model.rds"))
    
    diagnostics <- as.data.frame(
        posterior::summarise_draws(
            posterior::as_draws_array(fit_bayes)
        )
    )
    
    write.csv(
        diagnostics,
        file.path(output_dir, "03_bayesian_diagnostics.csv"),
        row.names = FALSE
    )
    
    pars <- diagnostics %>%
        filter(grepl("^(b_|sd_|r_)", variable))
    
    nuts <- brms::nuts_params(fit_bayes)
    
    diagnostics_passed <- nrow(pars) > 0 &&
        all(
            is.finite(pars$rhat) & pars$rhat < 1.05 &
                is.finite(pars$ess_bulk) & pars$ess_bulk > 200 &
                is.finite(pars$ess_tail) & pars$ess_tail > 200
        ) &&
        !any(nuts$Parameter == "divergent__" & nuts$Value == 1) &&
        !any(nuts$Parameter == "treedepth__" & nuts$Value >= 10)
    
    cf <- brms::fixef(fit_bayes, probs = c(0.025, 0.975))
    
    result <- data.frame(
        Term = rownames(cf),
        Estimate = cf[, "Estimate"],
        SE = cf[, "Est.Error"],
        aOR = exp(cf[, "Estimate"]),
        Lower = exp(cf[, "Q2.5"]),
        Upper = exp(cf[, "Q97.5"]),
        row.names = NULL
    ) %>%
        filter(Term != "Intercept")
    
    interval <- "95% CrI"
    model_type <- "Bayesian mixed-effects logistic regression"
}

# Match model coefficients to labels
write.csv(
    result,
    file.path(output_dir, "04_raw_coefficients.csv"),
    row.names = FALSE
)

selected_map <- term_map %>%
    filter(Variable %in% selected_vars)

normalize_term <- function(x) gsub("[^[:alnum:]_]", "", x)

index <- match(result$Term, selected_map$Term)

for (i in which(is.na(index))) {
    matches <- which(
        normalize_term(selected_map$Term) == normalize_term(result$Term[i])
    )
    if (length(matches) == 1) index[i] <- matches
}

if (anyNA(index) || anyDuplicated(index)) {
    stop("Coefficient labels could not be matched uniquely. Check raw coefficients.")
}

result <- bind_cols(
    result,
    selected_map[index, c("Variable", "Level", "Reference")]
) %>%
    arrange(
        match(Variable, candidate_vars),
        match(Term, result$Term)
    )

model_info <- data.frame(
    Model = model_type,
    Formula = paste(deparse(model_formula), collapse = " "),
    Original_n = nrow(raw),
    Complete_n = nrow(df),
    Outcome_positive_n = sum(df$outcome),
    Interval = interval,
    Diagnostics_passed = diagnostics_passed
)

write_xlsx(
    list(
        model_info = model_info,
        references = reference_table,
        global_tests = global_results,
        univariable = univariable,
        multivariable = result
    ),
    file.path(output_dir, "04_regression_results.xlsx")
)

capture.output(
    sessionInfo(),
    file = file.path(output_dir, "session_info.txt")
)

if (!diagnostics_passed) {
    stop("Bayesian diagnostics failed. Results saved; review the fit before plotting.")
}

# Forest plot
category_names <- c(
    Region = "Region",
    StrainType = "Pathotype",
    Niche = "Source",
    ST = "Clone type"
)

colors <- c(
    Region = "#8499AD",
    Pathotype = "#A87987",
    Source = "#7E9F95",
    "Clone type" = "#B78D67"
)

plot_data <- result %>%
    mutate(
        Category = factor(
            unname(category_names[Variable]),
            levels = names(colors)
        ),
        y = rev(seq_len(n())),
        Excludes_one = Lower > 1 | Upper < 1,
        Label = paste0(Level, " (vs ", Reference, ")"),
        Result_text = sprintf("%.2f (%.2f–%.2f)", aOR, Lower, Upper)
    )

if (any(
    !is.finite(as.matrix(plot_data[, c("aOR", "Lower", "Upper")])) |
    plot_data$Lower <= 0
)) {
    stop("Invalid estimates or intervals for the forest plot.")
}

n_terms <- nrow(plot_data)
y_limits <- c(0.3, n_terms + 1.3)
x_limits <- c(
    min(0.1, min(plot_data$Lower) / 1.2),
    max(25, max(plot_data$Upper) * 1.2)
)

background <- plot_data[seq(1, n_terms, by = 2), ]

p_forest <- ggplot(plot_data, aes(y = y)) +
    geom_rect(
        data = background,
        aes(
            xmin = x_limits[1], xmax = x_limits[2],
            ymin = y - 0.43, ymax = y + 0.43
        ),
        inherit.aes = FALSE,
        fill = "#F5F6F7",
        colour = NA
    ) +
    annotate(
        "segment",
        x = 1, xend = 1, y = 0.5, yend = n_terms + 0.45,
        linetype = "dashed", colour = "#8B8B8B", linewidth = 0.5
    ) +
    geom_segment(
        aes(x = Lower, xend = Upper, yend = y, colour = Category),
        linewidth = 0.85
    ) +
    geom_segment(
        aes(
            x = Lower, xend = Lower,
            y = y - 0.09, yend = y + 0.09, colour = Category
        ),
        linewidth = 0.6
    ) +
    geom_segment(
        aes(
            x = Upper, xend = Upper,
            y = y - 0.09, yend = y + 0.09, colour = Category
        ),
        linewidth = 0.6
    ) +
    geom_point(
        aes(x = aOR, colour = Category, shape = Excludes_one),
        fill = "white", size = 3.4, stroke = 1
    ) +
    scale_colour_manual(values = colors, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16), guide = "none") +
    scale_x_log10(limits = x_limits) +
    scale_y_continuous(
        limits = y_limits,
        breaks = plot_data$y,
        labels = plot_data$Label,
        expand = expansion(mult = 0)
    ) +
    labs(x = "Adjusted odds ratio", y = NULL) +
    theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "#333333", size = 10.5)
    )

p_values <- ggplot(plot_data, aes(y = y)) +
    geom_rect(
        data = background,
        aes(xmin = 0, xmax = 1, ymin = y - 0.43, ymax = y + 0.43),
        inherit.aes = FALSE,
        fill = "#F5F6F7",
        colour = NA
    ) +
    geom_text(
        aes(x = 0.02, label = Result_text),
        hjust = 0, size = 3.5
    ) +
    annotate(
        "text",
        x = 0.02, y = n_terms + 1,
        label = paste0("aOR (", interval, ")"),
        hjust = 0, fontface = "bold", size = 3.7
    ) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = 0)) +
    theme_void()

p <- (p_forest | p_values) +
    plot_layout(widths = c(3.8, 1.5)) +
    plot_annotation(
        caption = paste0(
            "Bars indicate ", interval,
            ". Filled circles: interval excludes 1; open circles: interval includes 1."
        ),
        theme = theme(
            plot.caption = element_text(size = 9, hjust = 0, colour = "grey40")
        )
    )

print(p)

ggsave(
    file.path(output_dir, "Figure_4C_transmission_forest.pdf"),
    p,
    width = 10,
    height = max(4, 2 + 0.45 * n_terms),
    bg = "white",
    useDingbats = FALSE
)

message("Results saved in: ", output_dir)
