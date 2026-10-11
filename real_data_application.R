#!/usr/bin/env Rscript

# CEASE-SQR real-data application
# Dataset: 2024 ACS 1-year California person-level PUMS
#
# This standalone file downloads the official public-use archive when needed,
# constructs the analytic cohort, produces descriptive statistics, fits the
# nonprivate and private estimators, and creates the manuscript figures.
#
# Requirements:
#   R >= 4.5.3
#   install.packages(c("ggplot2", "patchwork"))
#
# Usage:
#   Rscript --vanilla real_data_application.R help
#   Rscript --vanilla real_data_application.R prepare
#   Rscript --vanilla real_data_application.R analyze 20
#   Rscript --vanilla real_data_application.R all 20
#   Rscript --vanilla real_data_application.R figures
#
# The number following "analyze" or "all" is the number of independent
# privacy-noise draws per epsilon value; the manuscript setting is 20.
# When sourced in RStudio, the analysis does not start automatically. Use
# real_data_main(c("all", "20")) after sourcing.
#
# Official data:
# https://www2.census.gov/programs-surveys/acs/data/pums/2024/1-Year/csv_pca.zip

options(stringsAsFactors = FALSE, warn = 1)

real_data_script_directory <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_line, value = TRUE)
  if (length(file_argument)) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_argument[[1L]]),
      winslash = "/", mustWork = FALSE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

app_dir <- real_data_script_directory()
data_dir <- file.path(app_dir, "data")
result_dir <- file.path(app_dir, "results")
dp_repetitions <- 20L

zip_file <- file.path(data_dir, "csv_pca.zip")
cohort_file <- file.path(data_dir, "acs2024_ca_wage_cohort.rds")
data_url <- paste0(
  "https://www2.census.gov/programs-surveys/acs/data/pums/2024/1-Year/",
  "csv_pca.zip"
)

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  cw <- cumsum(w[ord]) / sum(w)
  vapply(probs, function(pr) x[which(cw >= pr)[1L]], numeric(1L))
}

prepare_cohort <- function() {
  if (!file.exists(zip_file)) {
    message("Downloading the official 2024 California ACS PUMS archive ...")
    download.file(data_url, zip_file, mode = "wb", quiet = FALSE)
  }

  csv_name <- "psam_p06.csv"
  zfiles <- unzip(zip_file, list = TRUE)$Name
  if (!csv_name %in% zfiles) stop("Expected file psam_p06.csv was not found.")

  wanted <- c(
    "SERIALNO", "SPORDER", "ADJINC", "PWGTP", "AGEP", "CIT", "COW",
    "MAR", "SCHL", "SEX", "WAGP", "WKHP", "WKWN", "DIS", "ESR",
    "HISP", "RAC1P"
  )
  con_header <- unz(zip_file, csv_name, open = "rt")
  header <- strsplit(readLines(con_header, n = 1L), ",", fixed = TRUE)[[1L]]
  close(con_header)
  if (!all(wanted %in% header)) {
    stop("Missing required PUMS variables: ", paste(setdiff(wanted, header), collapse = ", "))
  }

  classes <- rep("NULL", length(header))
  names(classes) <- header
  classes[wanted] <- "numeric"
  classes["SERIALNO"] <- "character"

  message("Reading 17 selected fields from the compressed PUMS file ...")
  con <- unz(zip_file, csv_name, open = "rt")
  raw <- read.csv(
    con, header = TRUE, colClasses = unname(classes),
    na.strings = c("", "b", "bb", "bbb", "bbbb", "bbbbb", "bbbbbb"),
    check.names = FALSE
  )
  close(con)

  initial_n <- nrow(raw)
  duplicate_id <- duplicated(raw[c("SERIALNO", "SPORDER")])
  complete_fields <- c(
    "ADJINC", "PWGTP", "AGEP", "CIT", "COW", "MAR", "SCHL", "SEX",
    "WAGP", "WKHP", "WKWN", "DIS", "ESR", "HISP", "RAC1P"
  )
  eligible <- !duplicate_id & complete.cases(raw[complete_fields]) &
    raw$AGEP >= 25 & raw$AGEP <= 64 &
    raw$ESR %in% c(1, 2) & raw$WAGP > 0 &
    raw$WKHP >= 1 & raw$WKWN >= 1 &
    raw$SCHL >= 1 & raw$SCHL <= 24 &
    raw$COW >= 1 & raw$COW <= 8
  d <- raw[eligible, , drop = FALSE]

  d$wage_2024 <- d$WAGP * d$ADJINC / 1e6
  d$log_wage <- log(d$wage_2024)
  d$education <- cut(
    d$SCHL, breaks = c(0, 15, 17, 20, 21, 24),
    labels = c("Less than high school", "High school/GED", "Some college/associate",
               "Bachelor's", "Graduate/professional")
  )
  d$race_group <- c(
    "White", "Black", "Other", "Other", "Other", "Asian", "Other", "Other", "Other"
  )[d$RAC1P]
  d$sex_group <- ifelse(d$SEX == 2, "Female", "Male")

  attrition <- data.frame(
    stage = c(
      "All California person records",
      "Unique person identifiers",
      "Age 25--64",
      "Civilian employed",
      "Positive wage, hours, and weeks",
      "Complete model covariates (analytic cohort)"
    ),
    n = c(
      initial_n,
      sum(!duplicate_id),
      sum(!duplicate_id & raw$AGEP >= 25 & raw$AGEP <= 64, na.rm = TRUE),
      sum(!duplicate_id & raw$AGEP >= 25 & raw$AGEP <= 64 & raw$ESR %in% c(1, 2), na.rm = TRUE),
      sum(!duplicate_id & raw$AGEP >= 25 & raw$AGEP <= 64 & raw$ESR %in% c(1, 2) &
            raw$WAGP > 0 & raw$WKHP >= 1 & raw$WKWN >= 1, na.rm = TRUE),
      nrow(d)
    )
  )

  keep <- c(
    "SERIALNO", "SPORDER", "PWGTP", "AGEP", "CIT", "COW", "MAR", "SCHL",
    "SEX", "WKHP", "WKWN", "DIS", "HISP", "RAC1P", "wage_2024",
    "log_wage", "education", "race_group", "sex_group"
  )
  saveRDS(list(data = d[keep], attrition = attrition), cohort_file, compress = "xz")
  write.csv(attrition, file.path(result_dir, "cohort_attrition.csv"), row.names = FALSE)
  message("Prepared analytic cohort: ", format(nrow(d), big.mark = ","), " records.")
  invisible(d)
}

project_ball <- function(beta, radius) {
  nr <- sqrt(sum(beta^2))
  if (is.finite(nr) && nr > radius) beta * radius / nr else beta
}

rho_from_eps_delta <- function(epsilon, delta) {
  (sqrt(log(1 / delta) + epsilon) - sqrt(log(1 / delta)))^2
}

smooth_loss <- function(beta, X, y, tau, h) {
  r <- y - drop(X %*% beta)
  u <- r / h
  mean(r * (tau - 1 + pnorm(u)) + h * dnorm(u))
}

smooth_gradient <- function(beta, X, y, tau, h) {
  r <- y - drop(X %*% beta)
  -drop(crossprod(X, tau - 1 + pnorm(r / h))) / nrow(X)
}

fit_sqr <- function(X, y, tau, h, beta_init, radius, linear_term = NULL,
                    ridge = 0, max_iter = 12L, grad_tol = 1e-6) {
  p <- ncol(X)
  if (is.null(linear_term)) linear_term <- numeric(p)
  fn <- function(beta) {
    smooth_loss(beta, X, y, tau, h) + sum(linear_term * beta) +
      0.5 * ridge * sum(beta^2)
  }
  beta <- project_ball(beta_init, radius)
  for (iter in seq_len(max_iter)) {
    r <- y - drop(X %*% beta)
    u <- r / h
    cdf_u <- pnorm(u)
    density_u <- dnorm(u)
    g <- -drop(crossprod(X, tau - 1 + cdf_u)) / nrow(X) +
      linear_term + ridge * beta
    if (max(abs(g)) < grad_tol) break
    H <- crossprod(X, X * (density_u / h)) / nrow(X) +
      diag(ridge + 1e-8, p)
    step <- tryCatch(solve(H, g), error = function(e) g)
    directional <- sum(g * step)
    old <- mean(r * (tau - 1 + cdf_u) + h * density_u) +
      sum(linear_term * beta) + 0.5 * ridge * sum(beta^2)
    step_size <- 1
    proposal <- beta
    for (ls in seq_len(10L)) {
      proposal <- project_ball(beta - step_size * step, radius)
      if (fn(proposal) <= old - 1e-4 * step_size * directional) break
      step_size <- step_size / 2
    }
    if (sqrt(sum((proposal - beta)^2)) < 1e-7) {
      beta <- proposal
      break
    }
    beta <- proposal
  }
  project_ball(beta, radius)
}

split_workers <- function(X, y, K, seed = 20260814L) {
  set.seed(seed)
  ord <- sample.int(nrow(X))
  sizes <- rep.int(nrow(X) %/% K, K)
  sizes[seq_len(nrow(X) %% K)] <- sizes[seq_len(nrow(X) %% K)] + 1L
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  ids <- Map(function(a, b) ord[seq.int(a, b)], starts, ends)
  workers <- lapply(ids, function(ii) list(X = X[ii, , drop = FALSE], y = y[ii]))
  attr(workers, "X_all") <- X
  attr(workers, "y_all") <- y
  attr(workers, "weights") <- sizes / nrow(X)
  workers
}

worker_average <- function(estimates, workers) {
  drop(do.call(cbind, estimates) %*% attr(workers, "weights"))
}

global_gradient <- function(workers, beta, tau, h) {
  smooth_gradient(beta, attr(workers, "X_all"), attr(workers, "y_all"), tau, h)
}

fit_local_average <- function(workers, tau, b, radius) {
  p <- ncol(workers[[1L]]$X)
  fits <- lapply(workers, function(w) {
    init <- c(median(w$y), numeric(p - 1L))
    fit_sqr(w$X, w$y, tau, b, init, radius, max_iter = 12L)
  })
  worker_average(fits, workers)
}

fit_cease_path <- function(workers, tau, h, alpha, rounds, radius, beta_init,
                           local_iter = 5L) {
  p <- length(beta_init)
  path <- matrix(NA_real_, p, rounds)
  beta <- project_ball(beta_init, radius)
  for (tt in seq_len(rounds)) {
    gbar <- global_gradient(workers, beta, tau, h)
    fits <- lapply(workers, function(w) {
      gk <- smooth_gradient(beta, w$X, w$y, tau, h)
      fit_sqr(
        w$X, w$y, tau, h, beta, radius,
        linear_term = gbar - gk - alpha * beta,
        ridge = alpha, max_iter = local_iter, grad_tol = 5e-6
      )
    })
    beta <- project_ball(worker_average(fits, workers), radius)
    path[, tt] <- beta
  }
  path
}

fit_cease_two_bandwidth_path <- function(
    workers, tau, h, b, alpha, rounds, radius, beta_init, local_iter = 5L
) {
  p <- length(beta_init)
  path <- matrix(NA_real_, p, rounds)
  beta <- project_ball(beta_init, radius)
  for (tt in seq_len(rounds)) {
    gbar_h <- global_gradient(workers, beta, tau, h)
    fits <- lapply(workers, function(w) {
      gk_b <- smooth_gradient(beta, w$X, w$y, tau, b)
      fit_sqr(
        w$X, w$y, tau, b, beta, radius,
        linear_term = gbar_h - gk_b - alpha * beta,
        ridge = alpha, max_iter = local_iter, grad_tol = 5e-6
      )
    })
    beta <- project_ball(worker_average(fits, workers), radius)
    path[, tt] <- beta
  }
  path
}

fit_csl_path <- function(workers, tau, h, b, rounds, radius, beta_init) {
  p <- length(beta_init)
  path <- matrix(NA_real_, p, rounds)
  beta <- project_ball(beta_init, radius)
  master <- workers[[1L]]
  for (tt in seq_len(rounds)) {
    gbar <- global_gradient(workers, beta, tau, h)
    g1 <- smooth_gradient(beta, master$X, master$y, tau, b)
    beta <- fit_sqr(
      master$X, master$y, tau, b, beta, radius,
      linear_term = gbar - g1, max_iter = 12L
    )
    path[, tt] <- beta
  }
  path
}

fit_dp_cease_path <- function(workers, tau, h, alpha, rounds, radius,
                              epsilon, delta, Bx, beta_init = NULL,
                              local_iter = 5L) {
  p <- ncol(workers[[1L]]$X)
  N <- nrow(attr(workers, "X_all"))
  path <- matrix(NA_real_, p, rounds)
  if (is.null(beta_init)) beta_init <- numeric(p)
  beta <- project_ball(beta_init, radius) # data-independent initialization
  rho_release <- rho_from_eps_delta(epsilon, delta) / (2 * rounds)
  sigma_g <- (2 * Bx / N) / sqrt(2 * rho_release)
  sigma_beta <- (4 * Bx / (N * alpha)) / sqrt(2 * rho_release)
  for (tt in seq_len(rounds)) {
    gbar <- global_gradient(workers, beta, tau, h) + rnorm(p, sd = sigma_g)
    fits <- lapply(workers, function(w) {
      gk <- smooth_gradient(beta, w$X, w$y, tau, h)
      fit_sqr(
        w$X, w$y, tau, h, beta, radius,
        linear_term = gbar - gk - alpha * beta,
        ridge = alpha, max_iter = local_iter, grad_tol = 5e-6
      )
    })
    beta <- project_ball(
      worker_average(fits, workers) + rnorm(p, sd = sigma_beta), radius
    )
    path[, tt] <- beta
  }
  attr(path, "sigma_gradient") <- sigma_g
  attr(path, "sigma_parameter") <- sigma_beta
  path
}

fit_dp_cease_two_bandwidth_path <- function(
    workers, tau, h, b, alpha, rounds, radius, epsilon, delta, Bx,
    beta_init = NULL, local_iter = 5L
) {
  p <- ncol(workers[[1L]]$X)
  N <- nrow(attr(workers, "X_all"))
  path <- matrix(NA_real_, p, rounds)
  if (is.null(beta_init)) beta_init <- numeric(p)
  beta <- project_ball(beta_init, radius)
  rho_release <- rho_from_eps_delta(epsilon, delta) / (2 * rounds)
  sigma_g <- (2 * Bx / N) / sqrt(2 * rho_release)
  sigma_beta <- (4 * Bx / (N * alpha)) / sqrt(2 * rho_release)
  for (tt in seq_len(rounds)) {
    gbar_h <- global_gradient(workers, beta, tau, h) +
      rnorm(p, sd = sigma_g)
    fits <- lapply(workers, function(w) {
      gk_b <- smooth_gradient(beta, w$X, w$y, tau, b)
      fit_sqr(
        w$X, w$y, tau, b, beta, radius,
        linear_term = gbar_h - gk_b - alpha * beta,
        ridge = alpha, max_iter = local_iter, grad_tol = 5e-6
      )
    })
    beta <- project_ball(
      worker_average(fits, workers) + rnorm(p, sd = sigma_beta), radius
    )
    path[, tt] <- beta
  }
  attr(path, "sigma_gradient") <- sigma_g
  attr(path, "sigma_parameter") <- sigma_beta
  path
}

fit_dp_cease_pair <- function(
    workers, tau, h, b, alpha, rounds, radius, epsilon, delta, Bx,
    beta_init = NULL, local_iter = 5L
) {
  rng_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  hh <- fit_dp_cease_path(
    workers, tau, h, alpha, rounds, radius, epsilon, delta, Bx,
    beta_init = beta_init, local_iter = local_iter
  )
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  hb <- fit_dp_cease_two_bandwidth_path(
    workers, tau, h, b, alpha, rounds, radius, epsilon, delta, Bx,
    beta_init = beta_init, local_iter = local_iter
  )
  list(hh = hh, hb = hb)
}

pinball <- function(y, pred, tau) {
  r <- y - pred
  mean(r * (tau - (r < 0)))
}

evaluate_fit <- function(beta, X, y, tau, beta_global, method, round,
                         epsilon = NA_real_, repetition = NA_integer_) {
  # The model is fitted to the fixed public transformation
  # y_model = {log(wage) - 10}/2. Report predictive losses on the original
  # log-wage scale for direct interpretation.
  y_log <- 10 + 2 * y
  pred <- 10 + 2 * drop(X %*% beta)
  data.frame(
    method = method, round = round, epsilon = epsilon, repetition = repetition,
    pinball = pinball(y_log, pred, tau),
    mae = mean(abs(y_log - pred)),
    rmse = sqrt(mean((y_log - pred)^2)),
    calibration = mean(y_log <= pred),
    coefficient_distance = mean((beta - beta_global)^2)
  )
}

make_design <- function(d) {
  # Every transformation below is fixed in advance and bounded. This avoids
  # spending privacy budget on data-dependent centering or scaling.
  X <- cbind(
    `(Intercept)` = 1,
    age = (d$AGEP - 44.5) / 19.5,
    age_sq = ((d$AGEP - 44.5) / 19.5)^2,
    hours = (d$WKHP - 50) / 49,
    weeks = (d$WKWN - 26.5) / 25.5,
    female = as.numeric(d$SEX == 2),
    married = as.numeric(d$MAR == 1),
    some_college = as.numeric(d$SCHL >= 18 & d$SCHL <= 20),
    bachelors = as.numeric(d$SCHL == 21),
    graduate = as.numeric(d$SCHL >= 22),
    less_high_school = as.numeric(d$SCHL <= 15),
    government = as.numeric(d$COW %in% 3:5),
    self_employed = as.numeric(d$COW %in% 6:7),
    black = as.numeric(d$RAC1P == 2),
    asian = as.numeric(d$RAC1P == 6),
    other_race = as.numeric(!d$RAC1P %in% c(1, 2, 6)),
    hispanic = as.numeric(d$HISP != 1),
    noncitizen = as.numeric(d$CIT == 5),
    disability = as.numeric(d$DIS == 1)
  )
  storage.mode(X) <- "double"
  X
}

descriptive_outputs <- function(d, train_n, test_n) {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  wage_q <- quantile(d$wage_2024, probs = probs, names = FALSE)
  wage_wq <- weighted_quantile(d$wage_2024, d$PWGTP, probs)
  desc <- data.frame(
    variable = c(
      "Annual wage (2024 dollars)", "Log annual wage", "Age (years)",
      "Usual weekly hours", "Weeks worked", "Female", "Married",
      "Hispanic/Latino", "Non-U.S. citizen", "With a disability"
    ),
    mean_or_percent = c(
      mean(d$wage_2024), mean(d$log_wage), mean(d$AGEP), mean(d$WKHP),
      mean(d$WKWN), 100 * mean(d$SEX == 2), 100 * mean(d$MAR == 1),
      100 * mean(d$HISP != 1), 100 * mean(d$CIT == 5), 100 * mean(d$DIS == 1)
    ),
    sd = c(
      sd(d$wage_2024), sd(d$log_wage), sd(d$AGEP), sd(d$WKHP), sd(d$WKWN),
      rep(NA_real_, 5)
    ),
    median_or_percent = c(
      median(d$wage_2024), median(d$log_wage), median(d$AGEP), median(d$WKHP),
      median(d$WKWN), 100 * mean(d$SEX == 2), 100 * mean(d$MAR == 1),
      100 * mean(d$HISP != 1), 100 * mean(d$CIT == 5), 100 * mean(d$DIS == 1)
    )
  )
  write.csv(desc, file.path(result_dir, "descriptive_statistics.csv"), row.names = FALSE)
  write.csv(
    data.frame(probability = probs, unweighted = wage_q, survey_weighted = wage_wq),
    file.path(result_dir, "wage_quantiles.csv"), row.names = FALSE
  )

  groups <- do.call(rbind, lapply(c("sex_group", "education", "race_group"), function(v) {
    lev <- unique(as.character(d[[v]]))
    do.call(rbind, lapply(lev, function(z) {
      ii <- as.character(d[[v]]) == z
      data.frame(
        characteristic = v, group = z, n = sum(ii),
        mean_wage = mean(d$wage_2024[ii]),
        median_wage = median(d$wage_2024[ii]),
        weighted_median_wage = weighted_quantile(d$wage_2024[ii], d$PWGTP[ii], 0.5)
      )
    }))
  }))
  write.csv(groups, file.path(result_dir, "descriptive_groups.csv"), row.names = FALSE)

  writeLines(c(
    paste0("analytic_n=", nrow(d)),
    paste0("training_n=", train_n),
    paste0("test_n=", test_n),
    paste0("wage_mean=", signif(mean(d$wage_2024), 8)),
    paste0("wage_median=", signif(median(d$wage_2024), 8)),
    paste0("weighted_wage_median=", signif(weighted_quantile(d$wage_2024, d$PWGTP, 0.5), 8))
  ), file.path(result_dir, "analysis_counts.txt"))
}

make_figures <- function(d, metrics, dp_metrics) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    warning("ggplot2/patchwork unavailable; CSV results were still written.")
    return(invisible(NULL))
  }
  library(ggplot2)
  library(patchwork)
  theme_paper <- theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95")
    )

  p1 <- ggplot(d, aes(log_wage)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60,
                   fill = "#4C78A8", color = "white", linewidth = 0.15) +
    geom_density(color = "#B22222", linewidth = 0.8) +
    labs(x = "Log annual wage (2024 dollars)", y = "Density",
         title = "(a) Distribution of the response") + theme_paper
  plot_group <- subset(
    read.csv(file.path(result_dir, "descriptive_groups.csv")),
    characteristic == "education"
  )
  plot_group$group <- factor(
    plot_group$group,
    levels = c("Less than high school", "High school/GED", "Some college/associate",
               "Bachelor's", "Graduate/professional")
  )
  p2 <- ggplot(plot_group, aes(group, median_wage / 1000)) +
    geom_col(fill = "#59A14F", width = 0.72) +
    labs(x = NULL, y = "Median annual wage (thousand dollars)",
          title = "(b) Median wage by education") +
    scale_x_discrete(labels = c(
      "Less than high school" = "Less than\nhigh school",
      "High school/GED" = "High school/\nGED",
      "Some college/associate" = "Some college/\nassociate",
      "Bachelor's" = "Bachelor's",
      "Graduate/professional" = "Graduate/\nprofessional"
    )) + theme_paper +
    theme(axis.text.x = element_text(size = 9, lineheight = 0.95))
  descriptive <- p1 + p2 + plot_layout(widths = c(1.1, 1))
  ggsave(file.path(result_dir, "acs_descriptive.pdf"), descriptive,
         width = 11, height = 5.0, device = cairo_pdf)
  ggsave(file.path(result_dir, "acs_descriptive.png"), descriptive,
         width = 11, height = 5.0, dpi = 320)

  curve_data <- subset(
    metrics,
    method %in% c(
      "CEASE-SQR (h,h)", "CEASE-SQR (h,b)",
      "CSL-SQR (h,h)", "CSL-SQR (h,b)"
    )
  )
  for (benchmark in c("Global SQR", "One-shot average SQR")) {
    base_row <- subset(metrics, method == benchmark)[1L, , drop = FALSE]
    base_row <- base_row[rep(1L, 5L), , drop = FALSE]
    base_row$round <- seq_len(5L)
    curve_data <- rbind(curve_data, base_row)
  }
  p3 <- ggplot(curve_data, aes(round, pinball, color = method, shape = method)) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.2) +
    scale_x_continuous(breaks = 0:5) +
    labs(x = "Communication iterations", y = "Test check loss",
         color = NULL, shape = NULL, title = "(a) Nonprivate convergence") +
    theme_paper +
    theme(legend.text = element_text(size = 8.5)) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE),
           shape = guide_legend(ncol = 2, byrow = TRUE))

  dp_sum <- aggregate(
    cbind(pinball, mae, calibration) ~ epsilon + method,
    data = dp_metrics, FUN = mean
  )
  dp_sd <- aggregate(
    cbind(pinball, mae, calibration) ~ epsilon + method,
    data = dp_metrics, FUN = sd
  )
  names(dp_sd)[3:5] <- paste0(names(dp_sd)[3:5], "_sd")
  dp_plot <- merge(dp_sum, dp_sd, by = c("epsilon", "method"), all.x = TRUE)
  dp_plot$pinball_sd[is.na(dp_plot$pinball_sd)] <- 0
  dp_plot$mae_sd[is.na(dp_plot$mae_sd)] <- 0
  dp_plot$calibration_sd[is.na(dp_plot$calibration_sd)] <- 0
  write.csv(dp_plot, file.path(result_dir, "private_summary.csv"), row.names = FALSE)
  private_labels <- c(
    "CEASE-SQR (h,b), private tuning" = "Clean CEASE-SQR (h,b)",
    "CEASE-SQR (h,h), private tuning" = "Clean CEASE-SQR (h,h)",
    "DP-CEASE-SQR (h,b)" = "DP-CEASE-SQR (h,b)",
    "DP-CEASE-SQR (h,h)" = "DP-CEASE-SQR (h,h)"
  )
  p4 <- ggplot(dp_plot, aes(epsilon, pinball, color = method, shape = method)) +
    geom_errorbar(aes(ymin = pinball - pinball_sd, ymax = pinball + pinball_sd),
                  width = 0.08, linewidth = 0.45) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.3) +
    scale_color_discrete(labels = private_labels) +
    scale_shape_discrete(labels = private_labels) +
    scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
    labs(x = "Privacy budget epsilon", y = "Test check loss",
         color = NULL, shape = NULL, title = "(b) Privacy-utility trade-off") +
    theme_paper +
    theme(legend.text = element_text(size = 8.5)) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE),
           shape = guide_legend(ncol = 2, byrow = TRUE))
  methods <- p3 + p4
  ggsave(file.path(result_dir, "acs_method_comparison.pdf"), methods,
         width = 11, height = 6.0, device = cairo_pdf)
  ggsave(file.path(result_dir, "acs_method_comparison.png"), methods,
         width = 11, height = 6.0, dpi = 320)
  invisible(NULL)
}

run_analysis <- function() {
  if (!file.exists(cohort_file)) prepare_cohort()
  obj <- readRDS(cohort_file)
  d <- obj$data
  write.csv(obj$attrition, file.path(result_dir, "cohort_attrition.csv"), row.names = FALSE)

  set.seed(20260814L)
  train_id <- sample.int(nrow(d), floor(0.8 * nrow(d)))
  test_id <- setdiff(seq_len(nrow(d)), train_id)
  X <- make_design(d)
  # Fixed public transformation: unlike sample centering/scaling, this does
  # not create an unaccounted data-dependent release. It also puts the
  # data-independent zero initialization used by DP-CEASE-SQR on a sensible
  # numerical scale. All reported predictive losses are transformed back to
  # log-wage units in evaluate_fit().
  y <- (d$log_wage - 10) / 2
  X_train <- X[train_id, , drop = FALSE]
  y_train <- y[train_id]
  X_test <- X[test_id, , drop = FALSE]
  y_test <- y[test_id]

  tau <- 0.5
  K <- 50L
  rounds <- 5L
  radius <- 8
  p <- ncol(X_train)
  Bx <- sqrt(p) # all p entries are in [-1,1]
  N <- nrow(X_train)
  n_local <- floor(N / K)
  alpha_nonprivate <- 0.15 * p / n_local
  alpha_private <- 0.20
  h <- 1.4 * ((p + log(N)) / N)^(1 / 3)
  b <- max(h, 1.4 * ((p + log(n_local)) / n_local)^(1 / 3))
  delta <- N^(-2)

  descriptive_outputs(d, length(train_id), length(test_id))
  workers <- split_workers(X_train, y_train, K, seed = 20260815L)
  message("Fitting global and distributed nonprivate estimators ...")
  global_init <- c(median(y_train), numeric(p - 1L))
  beta_global <- fit_sqr(
    X_train, y_train, tau, h, global_init, radius,
    max_iter = 20L, grad_tol = 2e-7
  )
  beta_average <- fit_local_average(workers, tau, b, radius)
  master <- workers[[1L]]
  beta_master <- fit_sqr(
    master$X, master$y, tau, b,
    c(median(master$y), numeric(p - 1L)), radius, max_iter = 12L
  )
  cease_hh_path <- fit_cease_path(
    workers, tau, h, alpha_nonprivate, rounds, radius, beta_master,
    local_iter = 5L
  )
  cease_hb_path <- fit_cease_two_bandwidth_path(
    workers, tau, h, b, alpha_nonprivate, rounds, radius, beta_master,
    local_iter = 5L
  )
  csl_hh_path <- fit_csl_path(
    workers, tau, h, h, rounds, radius, beta_master
  )
  csl_hb_path <- fit_csl_path(
    workers, tau, h, b, rounds, radius, beta_master
  )

  metrics <- list(
    evaluate_fit(beta_global, X_test, y_test, tau, beta_global,
                 "Global SQR", 0L),
    evaluate_fit(beta_average, X_test, y_test, tau, beta_global,
                 "One-shot average SQR", 0L)
  )
  for (tt in seq_len(rounds)) {
    metrics[[length(metrics) + 1L]] <- evaluate_fit(
      cease_hh_path[, tt], X_test, y_test, tau, beta_global,
      "CEASE-SQR (h,h)", tt
    )
    metrics[[length(metrics) + 1L]] <- evaluate_fit(
      cease_hb_path[, tt], X_test, y_test, tau, beta_global,
      "CEASE-SQR (h,b)", tt
    )
    metrics[[length(metrics) + 1L]] <- evaluate_fit(
      csl_hh_path[, tt], X_test, y_test, tau, beta_global,
      "CSL-SQR (h,h)", tt
    )
    metrics[[length(metrics) + 1L]] <- evaluate_fit(
      csl_hb_path[, tt], X_test, y_test, tau, beta_global,
      "CSL-SQR (h,b)", tt
    )
  }
  metrics <- do.call(rbind, metrics)
  write.csv(metrics, file.path(result_dir, "nonprivate_metrics.csv"), row.names = FALSE)

  epsilon_grid <- c(0.5, 1, 2, 4, 8)
  dp_metrics <- vector("list", 2L * length(epsilon_grid) * dp_repetitions)
  z <- 1L
  message("Fitting both DP-CEASE-SQR bandwidth variants: ", dp_repetitions,
          " paired privacy-noise draws at each of ", length(epsilon_grid),
          " budgets ...")
  for (ee in seq_along(epsilon_grid)) {
    eps <- epsilon_grid[ee]
    for (rr in seq_len(dp_repetitions)) {
      set.seed(202608140L + 1000L * ee + rr)
      pair <- fit_dp_cease_pair(
        workers, tau, h, b, alpha_private, rounds, radius, eps, delta,
        Bx, beta_init = numeric(p), local_iter = 3L
      )
      row_hh <- evaluate_fit(
        pair$hh[, rounds], X_test, y_test, tau, beta_global,
        "DP-CEASE-SQR (h,h)", rounds, eps, rr
      )
      row_hh$sigma_gradient <- attr(pair$hh, "sigma_gradient")
      row_hh$sigma_parameter <- attr(pair$hh, "sigma_parameter")
      dp_metrics[[z]] <- row_hh
      z <- z + 1L
      row_hb <- evaluate_fit(
        pair$hb[, rounds], X_test, y_test, tau, beta_global,
        "DP-CEASE-SQR (h,b)", rounds, eps, rr
      )
      row_hb$sigma_gradient <- attr(pair$hb, "sigma_gradient")
      row_hb$sigma_parameter <- attr(pair$hb, "sigma_parameter")
      dp_metrics[[z]] <- row_hb
      z <- z + 1L
    }
    message("  completed epsilon = ", eps)
  }
  # Clean comparators use the same public alpha and zero initialization as
  # their private counterparts, isolating the effect of Gaussian releases.
  zero <- numeric(p)
  clean_hh_path <- fit_cease_path(
    workers, tau, h, alpha_private, rounds, radius, zero, local_iter = 5L
  )
  clean_hb_path <- fit_cease_two_bandwidth_path(
    workers, tau, h, b, alpha_private, rounds, radius, zero, local_iter = 5L
  )
  clean_hh <- evaluate_fit(
    clean_hh_path[, rounds], X_test, y_test, tau, beta_global,
    "CEASE-SQR (h,h), private tuning", rounds
  )
  clean_hb <- evaluate_fit(
    clean_hb_path[, rounds], X_test, y_test, tau, beta_global,
    "CEASE-SQR (h,b), private tuning", rounds
  )
  comparator <- do.call(rbind, lapply(epsilon_grid, function(eps) {
    out_hh <- clean_hh
    out_hb <- clean_hb
    out_hh$epsilon <- out_hb$epsilon <- eps
    out_hh$repetition <- out_hb$repetition <- 1L
    out_hh$sigma_gradient <- out_hb$sigma_gradient <- 0
    out_hh$sigma_parameter <- out_hb$sigma_parameter <- 0
    rbind(out_hh, out_hb)
  }))
  dp_metrics <- rbind(do.call(rbind, dp_metrics), comparator)
  write.csv(dp_metrics, file.path(result_dir, "private_metrics.csv"), row.names = FALSE)

  coefficient_table <- data.frame(
    term = colnames(X_train), global_sqr = beta_global,
    one_shot = beta_average,
    cease_sqr_hh = cease_hh_path[, rounds],
    cease_sqr_hb = cease_hb_path[, rounds],
    csl_sqr_hh = csl_hh_path[, rounds],
    csl_sqr_hb = csl_hb_path[, rounds],
    private_tuning_clean_hh = clean_hh_path[, rounds],
    private_tuning_clean_hb = clean_hb_path[, rounds]
  )
  write.csv(coefficient_table, file.path(result_dir, "coefficient_estimates.csv"), row.names = FALSE)
  settings <- data.frame(
    dataset = "2024 ACS 1-year California person PUMS",
    quantile = tau, analytic_n = nrow(d), training_n = N, test_n = length(test_id),
    predictors_including_intercept = p, machines = K,
    min_local_n = floor(N / K), max_local_n = ceiling(N / K),
    rounds = rounds, bandwidth_h = h, master_bandwidth_b = b,
    proximal_alpha_nonprivate = alpha_nonprivate,
    proximal_alpha_private = alpha_private, projection_radius = radius,
    covariate_bound_Bx = Bx, delta = delta,
    dp_repetitions_per_epsilon = dp_repetitions,
    seed_train_test = 20260814L, seed_machine_partition = 20260815L
  )
  write.csv(settings, file.path(result_dir, "analysis_settings.csv"), row.names = FALSE)
  make_figures(d, metrics, dp_metrics)
  saveRDS(
    list(settings = settings, metrics = metrics, private_metrics = dp_metrics,
         coefficients = coefficient_table),
    file.path(result_dir, "acs_real_data_results.rds"), compress = "xz"
  )
  message("Real-data analysis complete: ", normalizePath(result_dir))
  invisible(list(metrics = metrics, private_metrics = dp_metrics, settings = settings))
}

real_data_print_help <- function() {
  cat(paste(
    "",
    "CEASE-SQR real-data reproducibility code",
    "",
    "Usage:",
    "  Rscript --vanilla real_data_application.R prepare",
    "  Rscript --vanilla real_data_application.R analyze [PRIVACY_REPETITIONS]",
    "  Rscript --vanilla real_data_application.R all [PRIVACY_REPETITIONS]",
    "  Rscript --vanilla real_data_application.R figures",
    "",
    "The manuscript uses 20 independent privacy-noise draws per epsilon.",
    "The official ACS archive is downloaded only if it is not already cached.",
    sep = "\n"
  ))
  invisible(NULL)
}

real_data_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (!length(arguments)) return(real_data_print_help())
  mode <- tolower(arguments[[1L]])
  if (mode %in% c("help", "-h", "--help")) return(real_data_print_help())
  if (!mode %in% c("all", "prepare", "analyze", "figures")) {
    stop("MODE must be all, prepare, analyze, or figures.", call. = FALSE)
  }

  repetitions <- if (length(arguments) >= 2L) {
    suppressWarnings(as.integer(arguments[[2L]]))
  } else {
    20L
  }
  if (!is.finite(repetitions) || repetitions < 1L) {
    stop("PRIVACY_REPETITIONS must be a positive integer.", call. = FALSE)
  }
  dp_repetitions <<- repetitions
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  if (mode == "prepare") return(prepare_cohort())
  if (mode %in% c("all", "analyze")) return(run_analysis())

  if (!file.exists(cohort_file)) {
    stop("Prepared cohort is missing; run prepare or analyze first.",
         call. = FALSE)
  }
  metric_files <- file.path(
    result_dir, c("nonprivate_metrics.csv", "private_metrics.csv")
  )
  if (any(!file.exists(metric_files))) {
    stop("Analysis results are missing; run analyze first.", call. = FALSE)
  }
  obj <- readRDS(cohort_file)
  metrics <- read.csv(metric_files[[1L]], stringsAsFactors = FALSE)
  dp_metrics <- read.csv(metric_files[[2L]], stringsAsFactors = FALSE)
  make_figures(obj$data, metrics, dp_metrics)
  message("Figures regenerated: ", normalizePath(
    result_dir, winslash = "/", mustWork = TRUE
  ))
  invisible(list(metrics = metrics, private_metrics = dp_metrics))
}

if (sys.nframe() == 0L) {
  invisible(real_data_main())
}

