# Figure 4D: geographic distance and strain-sharing pair counts
# Run from the repository root:
# Rscript Figure4D_distance_decay.R
# Or supply input and output paths:
# Rscript Figure4D_distance_decay.R pairs.csv metadata.xlsx results/distance_decay
# See README_distance_decay.md for input schemas and interpretation.

main <- function() {
    required <- c("readr", "readxl", "dplyr", "tidyr", "ggplot2",
                  "glmmTMB", "DHARMa", "MASS", "scales")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop("Install missing packages: install.packages(c(",
             paste(sprintf('"%s"', missing), collapse = ", "), "))")
    }
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(ggplot2)
    })
    
    # Configuration: edit these paths when running interactively.
    pair_file <- "data/institution_pairs_within_province.csv"
    metadata_file <- "data/isolate_metadata.xlsx"
    out_dir <- "results/distance_decay"
    metadata_sheet <- NULL  # NULL detects the original kpn902 metadata sheet.
    expected_n_pairs <- 137L
    expected_n_zero <- 116L
    seed <- 2026L
    
    args <- if (interactive()) character() else commandArgs(trailingOnly = TRUE)
    if (length(args) && length(args) != 3L) {
        stop("Usage: Rscript Figure4D_distance_decay.R [pairs.csv metadata.xlsx output_dir]")
    }
    if (length(args) == 3L) {
        pair_file <- args[1]
        metadata_file <- args[2]
        out_dir <- args[3]
    }
    if (!all(file.exists(c(pair_file, metadata_file)))) {
        stop("Input files not found. Set the paths above or supply all three arguments.")
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
    set.seed(seed)
    
    require_columns <- function(d, columns, label) {
        absent <- setdiff(columns, names(d))
        if (length(absent)) stop(label, " lacks columns: ", paste(absent, collapse = ", "))
    }
    numeric_column <- function(x, label) {
        z <- suppressWarnings(as.numeric(as.character(x)))
        if (any(!is.finite(z))) stop(label, " contains missing or non-finite values.")
        z
    }
    
    # 1. Read all within-province pairs, including observed zero counts.
    pair_raw <- readr::read_csv(
        pair_file, col_types = readr::cols(.default = readr::col_character()),
        name_repair = "minimal", show_col_types = FALSE
    )
    require_columns(pair_raw, c("province", "instA", "instB",
                                "Strain-sharing pairs", "Distance"), "Pair table")
    pair <- pair_raw %>%
        transmute(
            province = trimws(.data[["province"]]),
            instA = trimws(.data[["instA"]]),
            instB = trimws(.data[["instB"]]),
            shared_pairs = numeric_column(.data[["Strain-sharing pairs"]], "Pair count"),
            distance_km = numeric_column(.data[["Distance"]], "Distance")
        ) %>%
        mutate(pair_id = paste(province, pmin(instA, instB), pmax(instA, instB), sep = "___"))
    if (anyNA(pair) || any(pair$province == "") || any(pair$instA == "") ||
        any(pair$instB == "")) stop("Missing institution or province identifiers.")
    if (any(pair$shared_pairs < 0 | pair$shared_pairs != floor(pair$shared_pairs)) ||
        any(pair$distance_km < 0) || any(pair$instA == pair$instB)) {
        stop("Counts must be non-negative integers; distances must be non-negative; no self-pairs.")
    }
    if (anyDuplicated(pair$pair_id)) stop("Duplicated unordered institution pairs.")
    if (nrow(pair) != expected_n_pairs || sum(pair$shared_pairs == 0) != expected_n_zero) {
        stop("Published dataset requires 137 pairs, including 116 zero-count pairs.")
    }
    
    sheets <- readxl::excel_sheets(metadata_file)
    if (is.null(metadata_sheet)) {
        candidates <- sheets[grepl("kpn902\\s*元数据", sheets, ignore.case = TRUE)]
        if (length(candidates) == 1L) {
            metadata_sheet <- candidates
        } else if (length(sheets) == 1L) {
            metadata_sheet <- sheets
        } else {
            stop("Set metadata_sheet explicitly; a unique metadata sheet could not be detected.")
        }
    }
    meta_raw <- readxl::read_excel(metadata_file, sheet = metadata_sheet)
    site_column <- if ("unit_code" %in% names(meta_raw)) "unit_code" else "单位代码"
    require_columns(meta_raw, site_column, "Metadata")
    # Preserve the original analysis: one metadata row is one sequenced isolate.
    # Do not replace these counts with counts from the IPC analysis subset.
    site_counts <- meta_raw %>%
        transmute(unit_code = trimws(as.character(.data[[site_column]]))) %>%
        filter(!is.na(unit_code), unit_code != "") %>%
        count(unit_code, name = "n_iso")
    unmatched <- setdiff(unique(c(pair$instA, pair$instB)), site_counts$unit_code)
    if (length(unmatched)) stop("Unmatched institution codes: ", paste(unmatched, collapse = ", "))
    pair <- pair %>%
        left_join(transmute(site_counts, instA = unit_code, nA = n_iso), by = "instA") %>%
        left_join(transmute(site_counts, instB = unit_code, nB = n_iso), by = "instB") %>%
        mutate(possible_pairs = as.double(nA) * nB)
    if (anyNA(pair$possible_pairs) || any(pair$possible_pairs <= 0) ||
        any(pair$shared_pairs > pair$possible_pairs)) stop("Invalid isolate-count denominators.")
    
    institutions <- bind_rows(
        transmute(pair, province, institution = instA),
        transmute(pair, province, institution = instB)
    ) %>% distinct()
    if (anyDuplicated(institutions$institution)) {
        stop("Institution codes must be unique across provinces.")
    }
    expected_pairs <- institutions %>% count(province) %>%
        summarise(total = sum(n * (n - 1) / 2)) %>% pull(total)
    if (expected_pairs != nrow(pair)) stop("Missing within-province institution combinations.")
    # This check uses the institutions represented in the input table. The input
    # must also include institutions with no sharing links at all.
    
    # 2. Construct model terms. Preserve endpoint ordering from the supplied CSV.
    is_hospital <- function(x) grepl("Hospital", x, ignore.case = TRUE)
    pair <- pair %>% mutate(
        pair_type = case_when(
            is_hospital(instA) & is_hospital(instB) ~ "Hospital-Hospital",
            xor(is_hospital(instA), is_hospital(instB)) ~ "Hospital-Other",
            TRUE ~ "Non-hospital"
        ),
        pair_type = factor(pair_type, levels = c("Hospital-Hospital", "Hospital-Other", "Non-hospital")),
        province = factor(province),
        distance_10km = distance_km / 10
    )
    unit_levels <- sort(unique(c(pair$instA, pair$instB)))
    pair <- pair %>% mutate(
        instA = factor(instA, levels = unit_levels),
        instB = factor(instB, levels = unit_levels)
    )
    readr::write_csv(pair, file.path(out_dir, "analysis_pairs.csv"))
    readr::write_csv(site_counts, file.path(out_dir, "institution_isolate_counts.csv"))
    readr::write_csv(tibble(
        n_pairs = nrow(pair), n_zero = sum(pair$shared_pairs == 0),
        n_positive = sum(pair$shared_pairs > 0),
        min_distance_km = min(pair$distance_km), max_distance_km = max(pair$distance_km)
    ), file.path(out_dir, "data_summary.csv"))
    
    # 3. Main model: NB1, log link, offset for all possible cross-site isolate pairs.
    # The two endpoint random-intercept terms are independent, as in the revised
    # analysis. They do not impose a shared institution effect across endpoint roles.
    fit_nb1 <- glmmTMB::glmmTMB(
        shared_pairs ~ distance_10km + pair_type + province +
            offset(log(possible_pairs)) + (1 | instA) + (1 | instB),
        family = glmmTMB::nbinom1(link = "log"), data = pair
    )
    check_fit <- function(fit, label) {
        if (!isTRUE(fit$fit$convergence == 0L) || !isTRUE(fit$sdr$pdHess)) {
            stop(label, " did not pass optimizer/Hessian checks; inspect the saved model.")
        }
    }
    saveRDS(fit_nb1, file.path(out_dir, "main_model.rds"))
    check_fit(fit_nb1, "Main model")
    writeLines(capture.output(summary(fit_nb1)), file.path(out_dir, "main_model_summary.txt"))
    
    distance_result <- function(fit, label) {
        co <- summary(fit)$coefficients$cond
        beta <- co["distance_10km", "Estimate"]
        se <- co["distance_10km", "Std. Error"]
        z <- qnorm(0.975)
        data.frame(
            model = label, IRR_per_km = exp(beta / 10),
            CI_lower = exp((beta - z * se) / 10),
            CI_upper = exp((beta + z * se) / 10),
            P_value = co["distance_10km", "Pr(>|z|)"]
        )
    }
    main_result <- distance_result(fit_nb1, "Main: all institution pairs")
    
    # 4. Simulation-based diagnostics. NB1 is retained; no zero-inflation term
    # is automatically added solely because the dataset contains many zeros.
    save_pdf <- function(filename, draw, width = 8, height = 6) {
        grDevices::pdf(file.path(out_dir, filename), width = width, height = height,
                       useDingbats = FALSE)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw()
    }
    sim_nb1 <- DHARMa::simulateResiduals(fittedModel = fit_nb1, n = 1000)
    save_pdf("main_DHARMa_residuals.pdf", function() plot(sim_nb1))
    save_pdf("main_residuals_by_distance.pdf", function() {
        DHARMa::plotResiduals(sim_nb1, form = pair$distance_10km)
    })
    # These DHARMa test functions can draw their own diagnostic figures.
    tests <- NULL
    save_pdf("main_DHARMa_tests.pdf", function() {
        tests <<- list(
            uniformity = DHARMa::testUniformity(sim_nb1),
            dispersion = DHARMa::testDispersion(sim_nb1),
            zero_inflation = DHARMa::testZeroInflation(sim_nb1),
            distance_quantiles = DHARMa::testQuantiles(sim_nb1, predictor = pair$distance_10km),
            pair_type = DHARMa::testCategorical(sim_nb1, catPred = pair$pair_type)
        )
    })
    saveRDS(tests, file.path(out_dir, "main_diagnostic_tests.rds"))
    writeLines(capture.output(print(tests)), file.path(out_dir, "main_diagnostic_tests.txt"))
    
    # 5. Sensitivity analysis: omit the three within-Hebei Hospital A/B/C pairs.
    is_abc_pair <- as.character(pair$province) == "Hebei" &
        grepl("Hospital [ABC]$", as.character(pair$instA)) &
        grepl("Hospital [ABC]$", as.character(pair$instB))
    if (sum(is_abc_pair) != 3L) stop("Expected three Hebei Hospital A/B/C pairs; check institution codes.")
    pair_no_abc <- pair[!is_abc_pair, ]
    readr::write_csv(pair[is_abc_pair, ], file.path(out_dir, "sensitivity_excluded_pairs.csv"))
    fit_no_abc <- update(fit_nb1, data = pair_no_abc)
    saveRDS(fit_no_abc, file.path(out_dir, "sensitivity_model.rds"))
    check_fit(fit_no_abc, "Sensitivity model")
    writeLines(capture.output(summary(fit_no_abc)), file.path(out_dir, "sensitivity_model_summary.txt"))
    sensitivity_result <- distance_result(fit_no_abc, "Sensitivity: excluding three Hebei pairs")
    results <- bind_rows(main_result, sensitivity_result)
    readr::write_csv(results, file.path(out_dir, "Table_S6_distance_effects.csv"))
    print(results)
    sim_no_abc <- DHARMa::simulateResiduals(fittedModel = fit_no_abc, n = 1000)
    sensitivity_tests <- NULL
    save_pdf("sensitivity_DHARMa_tests.pdf", function() {
        sensitivity_tests <<- list(
            dispersion = DHARMa::testDispersion(sim_no_abc),
            zero_inflation = DHARMa::testZeroInflation(sim_no_abc)
        )
    })
    writeLines(capture.output(print(sensitivity_tests)),
               file.path(out_dir, "sensitivity_diagnostic_tests.txt"))
    saveRDS(sensitivity_tests, file.path(out_dir, "sensitivity_diagnostic_tests.rds"))
    
    # 6. Figure 4D: observed positive rates and zero-count rug marks.
    # The curve averages fixed-effect predictions over the observed distribution
    # of province and pair type. Random effects are set to zero, not integrated out.
    # Bands are approximate 95% confidence intervals from fixed-effect coefficient
    # simulations; they are not prediction intervals for individual pairs.
    pair_plot <- pair %>% mutate(rate_per_10000 = shared_pairs / possible_pairs * 10000)
    distance_grid <- seq(min(pair$distance_km), max(pair$distance_km), length.out = 100)
    profiles <- pair %>% mutate(profile_id = row_number()) %>%
        select(profile_id, province, pair_type)
    pred_data <- tidyr::crossing(distance_km = distance_grid, profile_id = profiles$profile_id) %>%
        left_join(profiles, by = "profile_id") %>%
        mutate(
            distance_10km = distance_km / 10,
            province = factor(province, levels = levels(pair$province)),
            pair_type = factor(pair_type, levels = levels(pair$pair_type))
        )
    X <- model.matrix(~ distance_10km + pair_type + province, data = pred_data)
    beta <- glmmTMB::fixef(fit_nb1)$cond
    V <- as.matrix(vcov(fit_nb1)$cond)
    if (!all(names(beta) %in% colnames(X))) stop("Prediction matrix does not match fixed effects.")
    X <- X[, names(beta), drop = FALSE]
    V <- V[names(beta), names(beta), drop = FALSE]
    set.seed(seed)
    beta_draws <- MASS::mvrnorm(n = 3000, mu = beta, Sigma = V)
    pred_summary <- bind_rows(lapply(distance_grid, function(d) {
        X_d <- X[pred_data$distance_km == d, , drop = FALSE]
        point <- mean(exp(X_d %*% beta + log(10000)))
        draws <- colMeans(exp(X_d %*% t(beta_draws) + log(10000)))
        data.frame(distance_km = d, predicted_rate = point,
                   lower = unname(quantile(draws, 0.025)),
                   upper = unname(quantile(draws, 0.975)))
    }))
    readr::write_csv(pred_summary, file.path(out_dir, "Figure4D_curve.csv"))
    p_text <- if (main_result$P_value < 0.001) "< 0.001" else
        paste0("= ", sprintf("%.3f", main_result$P_value))
    effect_label <- sprintf(
        "Adjusted IRR per km = %.2f (95%% CI %.2f-%.2f)\nP %s",
        main_result$IRR_per_km, main_result$CI_lower, main_result$CI_upper, p_text
    )
    p <- ggplot() +
        geom_ribbon(data = pred_summary, aes(x = distance_km, ymin = lower, ymax = upper),
                    fill = "grey60", alpha = 0.30) +
        geom_line(data = pred_summary, aes(distance_km, predicted_rate),
                  colour = "#222222", linewidth = 1) +
        geom_point(data = filter(pair_plot, shared_pairs > 0),
                   aes(distance_km, rate_per_10000, colour = pair_type), size = 2.4, alpha = 0.85) +
        geom_rug(data = filter(pair_plot, shared_pairs == 0), aes(x = distance_km),
                 sides = "b", colour = "grey45", alpha = 0.45,
                 length = grid::unit(0.035, "npc")) +
        annotate("label", x = Inf, y = Inf, label = effect_label,
                 hjust = 1.05, vjust = 1.1, size = 3.3, fill = "white") +
        scale_colour_manual(values = c("Hospital-Hospital" = "#B24C4C",
                                       "Hospital-Other" = "#D69A37", "Non-hospital" = "#4C78A8")) +
        scale_x_continuous(limits = range(pair$distance_km), breaks = scales::pretty_breaks(n = 6)) +
        scale_y_log10() +
        labs(x = "Geographic distance between institutions (km)",
             y = "Sharing pairs per 10,000 possible pairs", colour = "Institution-pair type",
             subtitle = paste0("All within-province pairs: ", nrow(pair),
                               "; zero-sharing pairs: ", sum(pair$shared_pairs == 0)),
             caption = paste("Curve: standardized fixed-effect predictions with approximate 95% CI.",
                             "\nLog10 y-axis; bottom rug marks indicate zero-count pairs.")) +
        theme_classic(base_size = 12) +
        theme(legend.position = "top", axis.title = element_text(colour = "black"),
              axis.text = element_text(colour = "black"), plot.subtitle = element_text(colour = "grey30"),
              plot.caption = element_text(colour = "grey35", hjust = 0))
    ggsave(file.path(out_dir, "Figure4D_distance_decay_NB1.pdf"), p,
           width = 8, height = 5.5, bg = "white", useDingbats = FALSE)
    ggsave(file.path(out_dir, "Figure4D_distance_decay_NB1.png"), p,
           width = 8, height = 5.5, dpi = 300, bg = "white")
    if (interactive()) print(p)
    writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
    message("Analysis complete. Outputs: ", normalizePath(out_dir))
    invisible(list(main = fit_nb1, sensitivity = fit_no_abc, results = results, plot = p))
}

main()
