# Figure 5B: IPC factors associated with putative transmission
rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(broom)
library(brms)
library(ggplot2)
library(writexl)

# Settings
message("Select the IPC metadata Excel file.")
input_file <- file.choose()
output_dir <- file.path(dirname(input_file), "Figure_5B")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

continuous_vars <- c(
    "Number of Beds",
    "Annual Average Number of Admitted Patients",
    "Doctor-to-Nurse Ratio",
    "Bed-to-Nurse Ratio",
    "Average Day of Hospital Stay",
    "Department Area",
    "Bed Spacing",
    "Frequency of Bedding Change",
    "Family Visiting Number Limit"
)

surfaces <- c(
    "Bed Rails, Adjusters and Nightstands",
    "Fixed Medical Equipment",
    "Mobile Medical Equipment",
    "Ward Toilets",
    "Ward Furniture",
    "Door Handles and Light Switches",
    "Public Sinks",
    "Public Door Handles, Handrails and Elevator Buttons",
    "Utility Carts",
    "Computer Mice and Keyboards"
)

references <- c(
    setNames(rep("5", 10), paste("Disinfection Method for", surfaces)),
    setNames(rep("0", 10), paste("Disinfection Frequency for", surfaces)),
    "Air Ventilation Method for Wards" = "1",
    "Air Ventilation Frequency for Wards" = "0",
    "Type of Hand Hygiene Facilities in Public Area" = "1",
    "Whether Patient Beds Are Fixed" = "0",
    "Family Visiting Management Method" = "0",
    "Family Caregiver Presence" = "0",
    "Patient Activity Range" = "4",
    "Department Type" = "1"
)

variables <- c(continuous_vars, names(references))

# Read data
raw <- read_excel(input_file, sheet = "重编码-英文表头")
names(raw) <- trimws(names(raw))

required <- c("Transmission Status", variables)
missing_columns <- setdiff(required, names(raw))

if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

to_numeric <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "N/A")] <- NA_character_
    value <- suppressWarnings(as.numeric(x))
    
    if (any(!is.na(x) & is.na(value))) {
        stop("Non-numeric values found in a coded numeric column.")
    }
    
    value[!is.finite(value)] <- NA_real_
    value
}

dat <- raw %>%
    select(all_of(required)) %>%
    mutate(across(everything(), to_numeric)) %>%
    rename(outcome = `Transmission Status`) %>%
    filter(!is.na(outcome))

if (!setequal(unique(dat$outcome), c(0, 1))) {
    stop("Transmission Status must contain both 0 and 1 only.")
}

# Univariable logistic regression
run_univariable <- function(v) {
    d <- data.frame(outcome = dat$outcome, x = dat[[v]]) %>%
        drop_na()
    
    if (!nrow(d)) {
        return(tibble(Variable = v, Status = "No complete observations"))
    }
    
    categorical <- v %in% names(references)
    reference <- NA_character_
    
    if (categorical) {
        d$x <- factor(d$x)
        reference <- references[[v]]
        
        if (!reference %in% levels(d$x)) {
            if (v == "Department Type") {
                stop("Department Type reference level 1 is missing.")
            }
            reference <- levels(d$x)[1]
            warning(v, ": using reference level ", reference)
        }
        
        d$x <- relevel(d$x, ref = reference)
        if (nlevels(d$x) > 1) {
            contrasts(d$x) <- contr.treatment(levels(d$x))
        }
    }
    
    counts <- d %>%
        mutate(Level = if (categorical) as.character(x) else "Continuous") %>%
        group_by(Level) %>%
        summarise(
            Total_n = n(),
            Transmission_n = sum(outcome),
            Transmission_percent = 100 * mean(outcome),
            .groups = "drop"
        )
    
    counts$Is_reference <- categorical & counts$Level == reference
    
    fit_result <- tryCatch({
        if (length(unique(d$x)) < 2 ||
            length(unique(d$outcome)) < 2) {
            stop("Predictor or outcome has fewer than two observed values.")
        }
        
        fit <- glm(outcome ~ x, data = d, family = binomial())
        
        tidy(fit, conf.int = TRUE) %>%
            filter(term != "(Intercept)") %>%
            transmute(
                Level = if (categorical) sub("^x", "", term) else "Continuous",
                Estimate = estimate,
                SE = std.error,
                OR = exp(estimate),
                CI_low = exp(conf.low),
                CI_high = exp(conf.high),
                P_value = p.value,
                Status = if (fit$converged) "Fitted" else "Not converged"
            )
    }, error = function(e) {
        tibble(
            Level = counts$Level,
            Estimate = NA_real_, SE = NA_real_,
            OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
            P_value = NA_real_,
            Status = conditionMessage(e)
        )
    })
    
    counts %>%
        left_join(fit_result, by = "Level") %>%
        mutate(
            Variable = v,
            Reference = reference,
            OR = if_else(Is_reference, 1, OR),
            Status = if_else(Is_reference, "Reference", Status)
        ) %>%
        select(Variable, Level, Reference, everything())
}

univariable <- bind_rows(lapply(variables, run_univariable))

write_xlsx(
    univariable,
    file.path(output_dir, "01_univariable_results.xlsx")
)

# Fixed six-variable multivariable model
model_data <- dat %>%
    transmute(
        outcome,
        Disinfection = `Disinfection Frequency for Mobile Medical Equipment`,
        Activity = `Patient Activity Range`,
        Department = `Department Type`,
        BedNurseRatio = `Bed-to-Nurse Ratio`,
        BedSpacing = `Bed Spacing`,
        VisitingLimit = `Family Visiting Number Limit`
    ) %>%
    drop_na()

factor_levels <- list(
    Disinfection = c("0", "1", "2", "3"),
    Activity = c("4", "0", "1", "2", "3"),
    Department = c("1", "2")
)

for (v in names(factor_levels)) {
    observed <- as.character(model_data[[v]])
    allowed <- factor_levels[[v]]
    
    if (any(!observed %in% allowed) || !allowed[1] %in% observed) {
        stop("Invalid codes or missing reference level in ", v)
    }
    
    model_data[[v]] <- droplevels(factor(observed, levels = allowed))
    
    if (nlevels(model_data[[v]]) < 2) {
        stop("Fewer than two observed levels in ", v)
    }
    
    contrasts(model_data[[v]]) <- contr.treatment(levels(model_data[[v]]))
}

if (!setequal(unique(model_data$outcome), c(0, 1))) {
    stop("The complete-case outcome must contain both 0 and 1.")
}

model_formula <- outcome ~
    Disinfection + Activity + Department +
    BedNurseRatio + BedSpacing + VisitingLimit

mm <- model.matrix(model_formula, model_data)
if (qr(mm)$rank < ncol(mm)) {
    stop("The multivariable design matrix is rank deficient.")
}

fit <- brm(
    formula = model_formula,
    data = model_data,
    family = bernoulli(link = "logit"),
    prior = c(
        prior(normal(0, 1), class = "b"),
        prior(normal(0, 2.5), class = "Intercept")
    ),
    chains = 4,
    cores = 4,
    iter = 2000,
    warmup = 1000,
    seed = 123,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE),
    refresh = 100
)

saveRDS(fit, file.path(output_dir, "02_bayesian_model.rds"))

# Posterior estimates and diagnostics
posterior_table <- as.data.frame(fixef(fit, probs = c(0.025, 0.975))) %>%
    tibble::rownames_to_column("Term") %>%
    mutate(
        aOR = exp(Estimate),
        Lower = exp(Q2.5),
        Upper = exp(Q97.5)
    )

parameter_diagnostics <- as.data.frame(
    posterior::summarise_draws(posterior::as_draws_array(fit))
)

nuts <- nuts_params(fit)
fixed_diagnostics <- parameter_diagnostics %>%
    filter(grepl("^b_", variable))

diagnostics <- tibble(
    N = nrow(model_data),
    Transmission_n = sum(model_data$outcome),
    Maximum_Rhat = max(fixed_diagnostics$rhat),
    Minimum_ESS_bulk = min(fixed_diagnostics$ess_bulk),
    Minimum_ESS_tail = min(fixed_diagnostics$ess_tail),
    Divergences = sum(nuts$Parameter == "divergent__" & nuts$Value == 1),
    Treedepth_hits = sum(nuts$Parameter == "treedepth__" & nuts$Value >= 15)
)

# Match model terms to table labels
variable_labels <- c(
    Disinfection = "Mobile equipment disinfection frequency",
    Activity = "Patient activity range",
    Department = "Department type",
    BedNurseRatio = "Bed-to-nurse ratio",
    BedSpacing = "Bed spacing",
    VisitingLimit = "Family visiting number limit"
)

level_labels <- list(
    Disinfection = c(
        "0" = "Irregular", "1" = "Once daily",
        "2" = "Twice daily", "3" = "Three times daily"
    ),
    Activity = setNames(paste("Code", 0:4), as.character(0:4)),
    Department = c("1" = "Type 1", "2" = "Type 2")
)

layout <- bind_rows(lapply(names(variable_labels), function(v) {
    x <- model_data[[v]]
    categorical <- is.factor(x)
    levels_v <- if (categorical) levels(x) else "Continuous"
    reference <- if (categorical) levels_v[1] else NA_character_
    
    tibble(
        Key = v,
        Variable = variable_labels[[v]],
        Raw_level = levels_v,
        Level = if (categorical) {
            unname(level_labels[[v]][levels_v])
        } else if (v == "VisitingLimit") {
            "Per 1-person increase"
        } else {
            "Per 1-unit increase"
        },
        Is_reference = if (categorical) levels_v == reference else FALSE,
        Term = if (categorical) {
            c(NA_character_, paste0(v, levels_v[-1]))
        } else {
            v
        },
        Total_n = vapply(levels_v, function(z) {
            if (categorical) sum(x == z) else length(x)
        }, integer(1)),
        Transmission_n = vapply(levels_v, function(z) {
            sum(model_data$outcome[if (categorical) x == z else rep(TRUE, length(x))])
        }, numeric(1))
    )
}))

results <- layout %>%
    left_join(posterior_table, by = "Term") %>%
    mutate(
        aOR = if_else(Is_reference, 1, aOR),
        Transmission = sprintf(
            "%d (%.1f%%)", Transmission_n, 100 * Transmission_n / Total_n
        ),
        Interval = if_else(
            Is_reference, "Reference", sprintf("%.2f–%.2f", Lower, Upper)
        )
    )

if (any(!results$Is_reference & is.na(results$aOR))) {
    stop("Some model terms could not be matched to table labels.")
}

publication_table <- results %>%
    transmute(
        Variable,
        Level = if_else(Is_reference, paste0(Level, " (ref)"), Level),
        `Total (n)` = Total_n,
        `Transmission n (%)` = Transmission,
        aOR = round(aOR, 2),
        `95% CrI` = Interval
    )

write_xlsx(
    list(
        Publication_table = publication_table,
        Coefficients = results,
        Model_diagnostics = diagnostics,
        Parameter_diagnostics = parameter_diagnostics
    ),
    file.path(output_dir, "03_multivariable_results.xlsx")
)

capture.output(
    sessionInfo(),
    file = file.path(output_dir, "sessionInfo.txt")
)

print(publication_table, n = Inf)
print(diagnostics)

if (
    any(!is.finite(fixed_diagnostics$rhat)) ||
    any(fixed_diagnostics$rhat >= 1.05) ||
    any(!is.finite(fixed_diagnostics$ess_bulk)) ||
    any(fixed_diagnostics$ess_bulk <= 200) ||
    any(!is.finite(fixed_diagnostics$ess_tail)) ||
    any(fixed_diagnostics$ess_tail <= 200) ||
    diagnostics$Divergences > 0 ||
    diagnostics$Treedepth_hits > 0
) {
    warning("Sampling diagnostics require review; inspect the saved results.")
}

# Figure 5B
forest_data <- results %>%
    filter(!Is_reference) %>%
    mutate(
        Label = paste0(Variable, ": ", Level),
        Label = factor(Label, levels = rev(Label))
    )

p <- ggplot(forest_data, aes(x = Label, y = aOR)) +
    geom_errorbar(
        aes(ymin = Lower, ymax = Upper),
        width = 0.18, linewidth = 0.55, colour = "#3F3F3F"
    ) +
    geom_point(
        shape = 21, size = 4,
        fill = "#71869A", colour = "#3F3F3F", stroke = 0.5
    ) +
    geom_hline(
        yintercept = 1, linetype = "dashed",
        linewidth = 0.5, colour = "#8A8A8A"
    ) +
    scale_y_log10() +
    coord_flip() +
    labs(x = NULL, y = "Adjusted odds ratio (95% CrI)") +
    theme_classic(base_size = 11, base_family = "sans") +
    theme(
        axis.text = element_text(colour = "black", size = 9.5),
        axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(10, 15, 10, 10)
    )

print(p)

ggsave(
    file.path(output_dir, "Figure_5B_IPC_forest.pdf"),
    p, width = 8, height = 6.5,
    bg = "white", useDingbats = FALSE
)

message("Results saved in: ", normalizePath(output_dir))

