calc_aicc <- function(model_obj, n_obs) {
  k_params <- length(coef(model_obj)) + 1
  aic_val <- AIC(model_obj)
  if ((n_obs - k_params - 1) <= 0) {
    return(NA_real_)
  }
  aic_val + (2 * k_params * (k_params + 1)) / (n_obs - k_params - 1)
}

calc_loocv_stats <- function(model_formula, model_df, min_n = 6) {
  n <- nrow(model_df)
  if (n < min_n) {
    return(list(rmse = NA_real_, mae = NA_real_, r2 = NA_real_))
  }

  response <- all.vars(model_formula)[1]
  obs <- model_df[[response]]
  pred <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    fit_i <- tryCatch(lm(model_formula, data = model_df[-i, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(fit_i)) {
      pred[i] <- tryCatch(
        as.numeric(predict(fit_i, newdata = model_df[i, , drop = FALSE])),
        error = function(e) NA_real_
      )
    }
  }

  valid <- is.finite(obs) & is.finite(pred)
  if (sum(valid) < 3) {
    return(list(rmse = NA_real_, mae = NA_real_, r2 = NA_real_))
  }

  errs <- obs[valid] - pred[valid]
  sst <- sum((obs[valid] - mean(obs[valid], na.rm = TRUE))^2, na.rm = TRUE)
  sse <- sum(errs^2, na.rm = TRUE)

  list(
    rmse = sqrt(mean(errs^2, na.rm = TRUE)),
    mae = mean(abs(errs), na.rm = TRUE),
    r2 = ifelse(sst > 0, 1 - sse / sst, NA_real_)
  )
}

compute_residual_diagnostics <- function(model_obj) {
  pull_scalar <- function(obj, key) {
    val <- tryCatch(obj[[key]], error = function(e) NULL)
    num <- suppressWarnings(as.numeric(val))
    if (length(num) >= 1 && is.finite(num[1])) {
      return(num[1])
    }
    NA_real_
  }

  resid_vals <- residuals(model_obj)
  resid_vals <- resid_vals[is.finite(resid_vals)]

  shapiro_w <- NA_real_
  shapiro_p <- NA_real_
  if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
    sh <- tryCatch(shapiro.test(resid_vals), error = function(e) NULL)
    if (!is.null(sh)) {
      shapiro_w <- suppressWarnings(as.numeric(unname(sh$statistic)))
      shapiro_p <- suppressWarnings(as.numeric(sh$p.value))
    }
  }

  ncv <- tryCatch(car::ncvTest(model_obj), error = function(e) NULL)
  ncv_chisq <- if (!is.null(ncv)) pull_scalar(ncv, "Chisquare") else NA_real_
  if (!is.finite(ncv_chisq) && !is.null(ncv)) {
    ncv_chisq <- pull_scalar(ncv, "ChiSquare")
  }
  ncv_p <- if (!is.null(ncv)) pull_scalar(ncv, "p") else NA_real_

  tibble::tibble(
    n_residuals = length(resid_vals),
    shapiro_W = shapiro_w,
    shapiro_p = shapiro_p,
    ncv_chisq = ncv_chisq,
    ncv_p = ncv_p,
    normality_pass_p05 = ifelse(is.finite(shapiro_p), shapiro_p > 0.05, NA),
    homoscedasticity_pass_p05 = ifelse(is.finite(ncv_p), ncv_p > 0.05, NA)
  )
}

round_export_cols <- function(df, cols, digits = 3) {
  keep <- intersect(cols, names(df))
  if (length(keep) == 0) {
    return(df)
  }
  dplyr::mutate(df, dplyr::across(dplyr::all_of(keep), ~ signif(.x, digits)))
}

format_export_outcome <- function(x, strip_mean = TRUE, drop_underscores = TRUE) {
  out <- as.character(x)
  if (isTRUE(strip_mean)) {
    out <- gsub("_mean$", "", out)
  }
  if (isTRUE(drop_underscores)) {
    out <- gsub("_", "", out, fixed = TRUE)
  }
  out
}

apply_vif_filter <- function(model_obj, outcome, model_df, threshold, protected = character()) {
  fit <- model_obj
  protected <- as.character(protected)

  repeat {
    retained <- setdiff(names(coef(fit)), "(Intercept)")
    drop_pool <- setdiff(retained, protected)
    if (length(drop_pool) == 0 || length(retained) <= 1) {
      break
    }

    vif_vals <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (is.null(vif_vals) || max(vif_vals, na.rm = TRUE) <= threshold) {
      break
    }

    ordered <- names(sort(vif_vals, decreasing = TRUE))
    drop_var <- ordered[ordered %in% drop_pool][1]
    if (is.na(drop_var) || !nzchar(drop_var)) {
      break
    }

    keep <- setdiff(retained, drop_var)
    if (length(keep) == 0) {
      return(NULL)
    }

    fit <- tryCatch(
      lm(as.formula(paste(outcome, "~", paste(keep, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- tryCatch(MASS::stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)
  }

  fit
}

apply_correlated_predictor_rules_catchment <- function(model_obj, outcome, model_df) {
  refit_drop_var <- function(drop_var) {
    keep <- setdiff(setdiff(names(coef(model_obj_curr)), "(Intercept)"), drop_var)
    if (length(keep) == 0) {
      return(NULL)
    }
    fit <- tryCatch(
      lm(as.formula(paste(outcome, "~", paste(keep, collapse = " + "))), data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    fit <- tryCatch(MASS::stepAIC(fit, direction = "backward", trace = 0), error = function(e) fit)
    keep_fit <- setdiff(names(coef(fit)), "(Intercept)")
    if (length(keep_fit) == 0) {
      return(NULL)
    }
    list(model = fit, aicc = calc_aicc(fit, nrow(model_df)))
  }

  model_obj_curr <- model_obj
  repeat {
    retained <- setdiff(names(coef(model_obj_curr)), "(Intercept)")
    has_ash <- "Ash_Per" %in% retained
    has_lava <- any(c("Lava1_per", "Lava2_per") %in% retained)
    has_both_landslide <- all(c("Landslide_Total", "Landslide_Young") %in% retained)

    if (!(has_ash && has_lava) && !has_both_landslide) {
      break
    }

    options <- list()
    if (has_ash && has_lava) {
      drop_candidates <- c("Ash_Per", intersect(c("Lava1_per", "Lava2_per"), retained))
      options <- lapply(drop_candidates, refit_drop_var)
    } else if (has_both_landslide) {
      options <- lapply(c("Landslide_Total", "Landslide_Young"), refit_drop_var)
    }

    valid <- which(vapply(options, function(x) !is.null(x), logical(1)))
    if (length(valid) == 0) {
      break
    }

    best <- valid[which.min(vapply(options[valid], function(x) x$aicc, numeric(1)))]
    model_obj_curr <- options[[best]]$model
  }

  model_obj_curr
}
