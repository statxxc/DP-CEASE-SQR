#!/usr/bin/env Rscript

# CEASE-SQR simulation code: nonprivate and private methods
#
# This standalone file contains all Monte Carlo procedures and ggplot2 code
# used for the simulation section of the manuscript. It does not source any
# other R file.
#
# Requirements:
#   R >= 4.5.3
#   install.packages(c("ggplot2", "patchwork"))
#
# Usage:
#   Rscript --vanilla simulation_methods.R help
#   Rscript --vanilla simulation_methods.R both       quick 1 5
#   Rscript --vanilla simulation_methods.R both       full  6 200
#   Rscript --vanilla simulation_methods.R nonprivate full  6 200
#   Rscript --vanilla simulation_methods.R private    full  6 200
#   Rscript --vanilla simulation_methods.R both       plot_only 1
#
# For backward-compatible convenience, omitting "both" is allowed:
#   Rscript --vanilla simulation_methods.R full 6 200
#
# When this file is sourced in RStudio, no simulation starts automatically.
# To run from the console after sourcing, use:
#   simulation_main(c("both", "full", "6", "200"))

options(stringsAsFactors = FALSE, warn = 1)

journal_script_directory <- function() {
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

journal_project_dir <- journal_script_directory()
simulation_result_dir <- file.path(journal_project_dir, "simulation_results")
journal_load_plot_packages <- function() {
  library_candidates <- c(
    file.path(journal_project_dir, ".Rlib"),
    file.path(journal_project_dir, ".rlib")
  )
  library_candidates <- library_candidates[dir.exists(library_candidates)]
  if (length(library_candidates)) {
    .libPaths(c(library_candidates, .libPaths()))
  }
  required <- c("ggplot2", "patchwork")
  missing <- required[
    !vapply(required, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Required package(s) not installed: ", paste(missing, collapse = ", "),
      ". Install them with install.packages(c(\"ggplot2\", \"patchwork\")).",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
  })
  invisible(TRUE)
}

journal_integer <- function(x, default, name) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) return(as.integer(default))
  value <- suppressWarnings(as.integer(x))
  if (!is.finite(value) || value < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  value
}

journal_with_project_directory <- function(expr) {
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(journal_project_dir)
  force(expr)
}
simulation_environment <- new.env(parent = .GlobalEnv)
simulation_environment$run_mode <- "quick"
simulation_environment$n_rep <- 5L
simulation_environment$n_cores <- 1L
simulation_environment$script_dir <- journal_project_dir

evalq({
  required_columns <- c(
    "study", "panel", "x_value", "error", "method", "mean_mse"
  )
  
  validate_summary <- function(data, name) {
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns)) {
      stop(
        name, " is missing required columns: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    if (any(!is.finite(data$x_value)) ||
        any(!is.finite(data$mean_mse)) || any(data$mean_mse <= 0)) {
      stop(name, " contains invalid x values or nonpositive MSE values.",
           call. = FALSE)
    }
    invisible(data)
  }
  
  pretty_error <- function(x) switch(
    x,
    normal = "Normal",
    student_t = "Student t",
    cauchy = "Cauchy",
    mixed_normal = "Mixed normal",
    x
  )
  
  x_axis_label <- function(panel) switch(
    panel,
    `Communication rounds` = "Communication iterations T",
    `Number of machines` = "Number of machines K",
    `Local sample size` = "Local sample size n",
    `Privacy budget` = "Privacy budget epsilon (log2 scale)",
    Bandwidth = "Bandwidth h",
    panel
  )
  
  focused_limits <- function(x) {
    value_range <- range(x, finite = TRUE)
    spread <- diff(value_range)
    padding <- if (spread > 0) max(0.06 * spread, 0.025) else 0.10
    value_range + c(-padding, padding)
  }
  
  x_scale_for <- function(panel, values) {
    values <- sort(unique(values))
    padding <- if (length(values) > 1L) 0.035 * diff(range(values)) else 0.5
    if (identical(panel, "Bandwidth")) {
      scale_x_continuous(
        breaks = values,
        labels = function(x) formatC(x, format = "f", digits = 3L),
        limits = range(values) + c(-padding, padding)
      )
    } else if (identical(panel, "Privacy budget")) {
      scale_x_continuous(
        transform = "log2",
        breaks = values,
        labels = function(x) format(x, trim = TRUE, scientific = FALSE),
        limits = range(values) * c(0.92, 1.08)
      )
    } else {
      scale_x_continuous(
        breaks = values,
        labels = function(x) format(x, trim = TRUE, scientific = FALSE,
                                    big.mark = ","),
        limits = range(values) + c(-padding, padding)
      )
    }
  }
  
  method_colors <- c(
    "CEASE-SQR" = "#009E73",
    "CSL-SQR" = "#0072B2",
    "CEASE-SQR (h,h)" = "#009E73",
    "CEASE-SQR (h,b)" = "#009E73",
    "CSL-SQR (h,h)" = "#0072B2",
    "CSL-SQR (h,b)" = "#0072B2",
    "CEASE-LR" = "#D55E00",
    "CSL-LR" = "#CC79A7",
    "DP-CEASE-SQR" = "#E69F00",
    "DP-CEASE-SQR (h,h)" = "#E69F00",
    "DP-CEASE-SQR (h,b)" = "#E69F00"
  )
  method_shapes <- c(
    "CEASE-SQR" = 21,
    "CSL-SQR" = 21,
    "CEASE-SQR (h,h)" = 21,
    "CEASE-SQR (h,b)" = 24,
    "CSL-SQR (h,h)" = 21,
    "CSL-SQR (h,b)" = 24,
    "CEASE-LR" = 22,
    "CSL-LR" = 22,
    "DP-CEASE-SQR" = 23,
    "DP-CEASE-SQR (h,h)" = 23,
    "DP-CEASE-SQR (h,b)" = 24
  )
  method_fills <- c(
    "CEASE-SQR" = method_colors[["CEASE-SQR"]],
    "CSL-SQR" = "white",
    "CEASE-SQR (h,h)" = method_colors[["CEASE-SQR (h,h)"]],
    "CEASE-SQR (h,b)" = "white",
    "CSL-SQR (h,h)" = method_colors[["CSL-SQR (h,h)"]],
    "CSL-SQR (h,b)" = "white",
    "CEASE-LR" = method_colors[["CEASE-LR"]],
    "CSL-LR" = "white",
    "DP-CEASE-SQR" = method_colors[["DP-CEASE-SQR"]],
    "DP-CEASE-SQR (h,h)" = method_colors[["DP-CEASE-SQR (h,h)"]],
    "DP-CEASE-SQR (h,b)" = "white"
  )
  method_linetypes <- c(
    "CEASE-SQR" = "solid",
    "CSL-SQR" = "dashed",
    "CEASE-SQR (h,h)" = "solid",
    "CEASE-SQR (h,b)" = "dashed",
    "CSL-SQR (h,h)" = "solid",
    "CSL-SQR (h,b)" = "dashed",
    "CEASE-LR" = "solid",
    "CSL-LR" = "dashed",
    "DP-CEASE-SQR" = "twodash",
    "DP-CEASE-SQR (h,h)" = "solid",
    "DP-CEASE-SQR (h,b)" = "dashed"
  )
  
  make_panel_plot <- function(
      data, panel_name, error_name, methods, row_tag = NULL,
      show_column_title = FALSE, show_y_title = FALSE
  ) {
    d <- data[
      data$panel == panel_name & data$error == error_name &
        data$method %in% methods,
      , drop = FALSE
    ]
    if (!nrow(d)) {
      stop("No rows for panel '", panel_name, "' and error '", error_name,
           "'.", call. = FALSE)
    }
    d$method <- factor(d$method, levels = methods)
    d$log10_mse <- log10(d$mean_mse)
    d <- d[order(d$method, d$x_value), , drop = FALSE]
    y_limits <- focused_limits(d$log10_mse)
  
    ggplot(
      d,
      aes(
        x = x_value,
        y = log10_mse,
        color = method,
        fill = method,
        shape = method,
        linetype = method,
        group = method
      )
    ) +
      geom_line(linewidth = 1.05) +
      geom_point(size = 3.05, stroke = 1.05) +
      x_scale_for(panel_name, d$x_value) +
      scale_y_continuous(
        breaks = pretty(y_limits, n = 5L),
        labels = function(x) formatC(x, format = "f", digits = 2L)
      ) +
      coord_cartesian(ylim = y_limits, clip = "off") +
      scale_color_manual(
        values = method_colors[methods],
        breaks = methods,
        limits = methods,
        drop = FALSE
      ) +
      scale_shape_manual(
        values = method_shapes[methods],
        breaks = methods,
        limits = methods,
        drop = FALSE
      ) +
      scale_fill_manual(
        values = method_fills[methods],
        breaks = methods,
        limits = methods,
        drop = FALSE,
        guide = "none"
      ) +
      scale_linetype_manual(
        values = method_linetypes[methods],
        breaks = methods,
        limits = methods,
        drop = FALSE
      ) +
      labs(
        title = if (show_column_title) pretty_error(error_name) else NULL,
        tag = row_tag,
        x = x_axis_label(panel_name),
        y = if (show_y_title) expression(log[10](MSE)) else NULL,
        color = NULL,
        shape = NULL,
        linetype = NULL
      ) +
      guides(
        color = guide_legend(
          nrow = 1L,
          byrow = TRUE,
          override.aes = list(
            colour = unname(method_colors[methods]),
            fill = unname(method_fills[methods]),
            shape = unname(method_shapes[methods]),
            linetype = unname(method_linetypes[methods]),
            linewidth = rep(1.25, length(methods)),
            size = rep(4.0, length(methods)),
            stroke = rep(1.05, length(methods))
          )
        ),
        shape = "none",
        linetype = "none"
      ) +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          color = "#222222",
          margin = margin(b = 8)
        ),
        plot.tag = element_text(size = 17, face = "bold", color = "#111111"),
        plot.tag.position = c(0.015, 0.985),
        axis.title.x = element_text(size = 14.5, face = "bold",
                                    margin = margin(t = 8)),
        axis.title.y = element_text(size = 14.5, face = "bold",
                                    margin = margin(r = 8)),
        axis.text = element_text(size = 12.5, color = "#222222"),
        axis.text.x = if (identical(panel_name, "Local sample size")) {
          element_text(
            angle = 30, hjust = 1, vjust = 1,
            margin = margin(t = 4)
          )
        } else {
          element_text(margin = margin(t = 3))
        },
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#E2E5E9", linewidth = 0.45),
        panel.border = element_rect(
          color = "#333333", fill = NA, linewidth = 0.65
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12.5),
        legend.key.width = grid::unit(1.15, "cm"),
        legend.margin = margin(t = 6),
        plot.margin = margin(9, 10, 8, 10)
      )
  }
  
  make_plot_matrix <- function(data, row_panels, column_errors, methods) {
    validate_summary(data, "Simulation summary")
    expected <- expand.grid(
      panel = row_panels,
      error = column_errors,
      method = methods,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    observed <- unique(data[c("panel", "error", "method")])
    observed$present <- TRUE
    coverage <- merge(
      expected, observed,
      by = c("panel", "error", "method"),
      all.x = TRUE
    )
    if (nrow(coverage) != nrow(expected) || any(is.na(coverage$present))) {
      stop("The summary does not cover the requested plot matrix.", call. = FALSE)
    }
  
    panel_list <- vector("list", length(row_panels) * length(column_errors))
    position <- 0L
    for (ri in seq_along(row_panels)) {
      for (ci in seq_along(column_errors)) {
        position <- position + 1L
        panel_list[[position]] <- make_panel_plot(
          data = data,
          panel_name = row_panels[[ri]],
          error_name = column_errors[[ci]],
          methods = methods,
          row_tag = if (ci == 1L) paste0("(", letters[[ri]], ")") else NULL,
          show_column_title = ri == 1L,
          show_y_title = ci == 1L
        )
      }
    }
  
    wrap_plots(
      panel_list,
      ncol = length(column_errors),
      guides = "collect"
    ) & theme(legend.position = "bottom")
  }
  
  save_matrix <- function(plot_object, pdf_file, width, height) {
    dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
    preview_file <- sub("\\.pdf$", "_preview.png", pdf_file,
                        ignore.case = TRUE)
    ggsave(
      filename = pdf_file,
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      device = grDevices::pdf,
      useDingbats = FALSE,
      limitsize = FALSE
    )
    ggsave(
      filename = preview_file,
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = 180,
      bg = "white",
      limitsize = FALSE
    )
    invisible(c(pdf = pdf_file, preview = preview_file))
  }
  
  plot_simulation_figures <- function(
      nonprivate_summary, private_summary, nonprivate_pdf, private_pdf
  ) {
    nonprivate_plot <- make_plot_matrix(
      nonprivate_summary,
      row_panels = c(
        "Communication rounds", "Number of machines", "Local sample size"
      ),
      column_errors = c("normal", "student_t", "cauchy", "mixed_normal"),
      methods = c(
        "CEASE-SQR (h,h)", "CEASE-SQR (h,b)",
        "CSL-SQR (h,h)", "CSL-SQR (h,b)",
        "CEASE-LR", "CSL-LR"
      )
    )
    private_plot <- make_plot_matrix(
      private_summary,
      row_panels = c(
        "Communication rounds", "Number of machines", "Local sample size",
        "Privacy budget", "Bandwidth"
      ),
      column_errors = c("normal", "student_t", "cauchy"),
      methods = c(
        "CEASE-SQR (h,h)", "CEASE-SQR (h,b)",
        "DP-CEASE-SQR (h,h)", "DP-CEASE-SQR (h,b)"
      )
    )
  
    save_matrix(nonprivate_plot, nonprivate_pdf, width = 16.5, height = 11.5)
    save_matrix(private_plot, private_pdf, width = 14.5, height = 17.5)
    invisible(list(nonprivate = nonprivate_plot, private = private_plot))
  }

  np_errors <- c("normal", "student_t", "cauchy", "mixed_normal")
  
  cfg <- list(
    p = 5L, tau = 0.5, total_N = 100000L, K = 50L, local_n = 2000L,
    rounds = 5L, allocation_rounds = 2L, epsilon = 1, delta = 1e-6,
    h_constant = 1.4, b_constant = 1.4,
    radius = 6, Bx = 4, local_solver_max_iter = 4L,
    machine_grid_nonprivate = c(25L, 50L, 75L, 100L, 125L),
    machine_grid_private = c(20L, 50L, 75L, 100L, 125L),
    local_n_grid_nonprivate = c(500L, 1000L, 1500L, 2000L, 2500L),
    local_n_grid_private = c(500L, 1000L, 1500L, 2000L, 2500L)
  )
  beta_true <- c(0.5, 1, -0.8, 0.6, -0.5)
  
  # Fan, Guo, and Wang (2023, JASA) practical CEASE specification.
  cease_alpha <- function(local_n) 0.15 * cfg$p / local_n
  
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
  
  fit_sqr <- function(X, y, tau, h, beta_init, radius,
                      linear_term = NULL, ridge = 0, max_iter = 12L,
                      grad_tol = 1e-6) {
    p <- ncol(X)
    if (is.null(linear_term)) linear_term <- numeric(p)
    fn <- function(beta) {
      smooth_loss(beta, X, y, tau, h) +
        sum(linear_term * beta) + 0.5 * ridge * sum(beta^2)
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
      weights <- density_u / h
      H <- crossprod(X, X * weights) / nrow(X) +
        diag(ridge + 1e-9, p)
      step <- tryCatch(solve(H, g), error = function(e) g)
      directional <- sum(g * step)
      old_value <- mean(
        r * (tau - 1 + cdf_u) + h * density_u
      ) + sum(linear_term * beta) + 0.5 * ridge * sum(beta^2)
      step_size <- 1
      accepted <- FALSE
      for (ls in seq_len(8L)) {
        proposal <- project_ball(beta - step_size * step, radius)
        if (fn(proposal) <= old_value - 1e-4 * step_size * directional) {
          accepted <- TRUE
          break
        }
        step_size <- step_size / 2
      }
      if (!accepted) proposal <- project_ball(beta - 1e-3 * g, radius)
      if (sqrt(sum((proposal - beta)^2)) < 1e-7) {
        beta <- proposal
        break
      }
      beta <- proposal
    }
    project_ball(beta, radius)
  }
  
  split_workers <- function(X, y, K) {
    N <- nrow(X)
    stopifnot(N >= K)
    worker_sizes <- rep.int(N %/% K, K)
    remainder <- N %% K
    if (remainder > 0L) {
      worker_sizes[seq_len(remainder)] <-
        worker_sizes[seq_len(remainder)] + 1L
    }
    end_index <- cumsum(worker_sizes)
    start_index <- c(1L, head(end_index, -1L) + 1L)
    ids <- Map(seq.int, start_index, end_index)
    out <- lapply(ids, function(ii) {
      list(X = X[ii, , drop = FALSE], y = y[ii])
    })
    attr(out, "X_all") <- X
    attr(out, "y_all") <- y
    attr(out, "worker_sizes") <- worker_sizes
    attr(out, "worker_weights") <- worker_sizes / N
    out
  }
  
  weighted_worker_average <- function(estimates, workers) {
    weights <- attr(workers, "worker_weights")
    estimate_matrix <- do.call(cbind, estimates)
    drop(estimate_matrix %*% weights)
  }
  
  global_sqr_gradient <- function(workers, beta, tau, h) {
    smooth_gradient(
      beta, attr(workers, "X_all"), attr(workers, "y_all"), tau, h
    )
  }
  
  global_lr_gradient <- function(workers, beta) {
    X <- attr(workers, "X_all")
    y <- attr(workers, "y_all")
    drop(crossprod(X, drop(X %*% beta) - y)) / nrow(X)
  }
  
  master_sqr_initialization <- function(workers, tau, h, radius) {
    p <- ncol(workers[[1L]]$X)
    master <- workers[[1L]]
    fit_sqr(
      master$X, master$y, tau, h, numeric(p), radius,
      max_iter = 8L, grad_tol = 1e-6
    )
  }
  
  master_lr_initialization <- function(workers, radius) {
    master <- workers[[1L]]
    p <- ncol(master$X)
    H <- crossprod(master$X) / nrow(master$X) + diag(1e-9, p)
    rhs <- drop(crossprod(master$X, master$y)) / nrow(master$X)
    project_ball(solve(H, rhs), radius)
  }
  
  # CEASE--SQR path.  Every worker solves its gradient-corrected proximal
  # surrogate and the server averages the local solutions, as in Algorithm 1.
  fit_cease_sqr_path <- function(workers, tau, h, alpha, max_rounds,
                                 radius, beta_init) {
    p <- length(beta_init)
    path <- matrix(NA_real_, p, max_rounds)
    beta <- project_ball(beta_init, radius)
    for (tt in seq_len(max_rounds)) {
      gbar <- global_sqr_gradient(workers, beta, tau, h)
      local_solutions <- lapply(workers, function(worker) {
        gk <- smooth_gradient(beta, worker$X, worker$y, tau, h)
        fit_sqr(
          worker$X, worker$y, tau, h, beta, radius,
          linear_term = gbar - gk - alpha * beta,
          ridge = alpha,
          max_iter = cfg$local_solver_max_iter,
          grad_tol = 5e-6
        )
      })
      beta <- project_ball(
        weighted_worker_average(local_solutions, workers), radius
      )
      path[, tt] <- beta
    }
    path
  }
  
  # Two-bandwidth CEASE--SQR.  The centralized gradient uses h, whereas every
  # worker supplies curvature through its b-bandwidth loss and local gradient.
  fit_cease_sqr_two_bandwidth_path <- function(
      workers, tau, h, b, alpha, max_rounds, radius, beta_init
  ) {
    p <- length(beta_init)
    path <- matrix(NA_real_, p, max_rounds)
    beta <- project_ball(beta_init, radius)
    for (tt in seq_len(max_rounds)) {
      gbar_h <- global_sqr_gradient(workers, beta, tau, h)
      local_solutions <- lapply(workers, function(worker) {
        gk_b <- smooth_gradient(beta, worker$X, worker$y, tau, b)
        fit_sqr(
          worker$X, worker$y, tau, b, beta, radius,
          linear_term = gbar_h - gk_b - alpha * beta,
          ridge = alpha,
          max_iter = cfg$local_solver_max_iter,
          grad_tol = 5e-6
        )
      })
      beta <- project_ball(
        weighted_worker_average(local_solutions, workers), radius
      )
      path[, tt] <- beta
    }
    path
  }
  
  # CSL--SQR: Algorithm 1 of Tan, Battey, and Zhou (2022).
  fit_csl_sqr_path <- function(workers, tau, h, b, max_rounds,
                               radius, beta_init) {
    p <- length(beta_init)
    path <- matrix(NA_real_, p, max_rounds)
    beta <- project_ball(beta_init, radius)
    master <- workers[[1L]]
    for (tt in seq_len(max_rounds)) {
      gbar <- global_sqr_gradient(workers, beta, tau, h)
      g1 <- smooth_gradient(beta, master$X, master$y, tau, b)
      beta <- fit_sqr(
        master$X, master$y, tau, b, beta, radius,
        linear_term = gbar - g1
      )
      path[, tt] <- beta
    }
    path
  }
  
  # Squared-loss CEASE path.  All workers solve their corrected local quadratic
  # problems and the server averages the solutions.
  fit_cease_lr_path <- function(workers, alpha, max_rounds, radius, beta_init) {
    p <- length(beta_init)
    path <- matrix(NA_real_, p, max_rounds)
    beta <- project_ball(beta_init, radius)
    local_stats <- lapply(workers, function(worker) {
      list(
        H = crossprod(worker$X) / nrow(worker$X),
        rhs = drop(crossprod(worker$X, worker$y)) / nrow(worker$X)
      )
    })
    for (tt in seq_len(max_rounds)) {
      gbar <- global_lr_gradient(workers, beta)
      local_solutions <- lapply(local_stats, function(stat) {
        gk <- drop(stat$H %*% beta) - stat$rhs
        linear_term <- gbar - gk - alpha * beta
        project_ball(
          solve(
            stat$H + diag(alpha + 1e-9, p),
            stat$rhs - linear_term
          ),
          radius
        )
      })
      beta <- project_ball(
        weighted_worker_average(local_solutions, workers), radius
      )
      path[, tt] <- beta
    }
    path
  }
  
  # CSL--LR of Jordan, Lee, and Yang (2019): the master solves one shifted
  # local quadratic problem at each round.
  fit_csl_lr_path <- function(workers, max_rounds, radius, beta_init) {
    p <- length(beta_init)
    path <- matrix(NA_real_, p, max_rounds)
    beta <- project_ball(beta_init, radius)
    master <- workers[[1L]]
    H1 <- crossprod(master$X) / nrow(master$X)
    rhs1 <- drop(crossprod(master$X, master$y)) / nrow(master$X)
    for (tt in seq_len(max_rounds)) {
      gbar <- global_lr_gradient(workers, beta)
      g1 <- drop(H1 %*% beta) - rhs1
      beta <- project_ball(
        solve(H1 + diag(1e-9, p), rhs1 - (gbar - g1)), radius
      )
      path[, tt] <- beta
    }
    path
  }
  
  draw_errors <- function(n, family, tau) {
    raw <- switch(
      family,
      normal = rnorm(n),
      student_t = rt(n, 3) / sqrt(3),
      cauchy = rcauchy(n),
      mixed_normal = ifelse(runif(n) < 0.9, rnorm(n), rnorm(n, sd = 8))
    )
    qtau <- switch(
      family,
      normal = qnorm(tau),
      student_t = qt(tau, 3) / sqrt(3),
      cauchy = qcauchy(tau),
      mixed_normal = uniroot(
        function(x) 0.9 * pnorm(x) + 0.1 * pnorm(x / 8) - tau,
        c(-100, 100)
      )$root
    )
    raw - qtau
  }
  
  generate_data <- function(local_n = NULL, K, family, total_N = NULL) {
    if (is.null(total_N)) {
      stopifnot(!is.null(local_n), local_n >= 1)
      N <- as.integer(local_n * K)
    } else {
      N <- as.integer(total_N)
    }
    X <- cbind(1, matrix(rnorm(N * (cfg$p - 1L)), N))
    X <- X / pmax(1, sqrt(rowSums(X^2)) / cfg$Bx)
    y <- drop(X %*% beta_true) + draw_errors(N, family, cfg$tau)
    split_workers(X, y, K)
  }
  
  bandwidths <- function(local_n, K, total_N = NULL) {
    N <- if (is.null(total_N)) local_n * K else total_N
    h <- cfg$h_constant * ((cfg$p + log(N)) / N)^(1 / 3)
    b <- cfg$b_constant * ((cfg$p + log(local_n)) / local_n)^(1 / 3)
    list(h = h, b = max(h, b))
  }
  
  mse <- function(beta) mean((beta - beta_true)^2)
  
  make_rows <- function(study, panel, x_value, family, local_n, K, rounds,
                        epsilon, bandwidth, replication, fits) {
    data.frame(
      study = study, panel = panel, x_value = x_value, error = family,
      local_n = local_n, machines = K, rounds = rounds,
      epsilon = epsilon, bandwidth = bandwidth, replication = replication,
      method = names(fits),
      mse = vapply(fits, mse, numeric(1))
    )
  }
  
  fit_nonprivate_paths <- function(workers, bw, max_rounds) {
    local_n <- mean(vapply(workers, function(worker) nrow(worker$X), numeric(1)))
    alpha <- cease_alpha(local_n)
    init_sqr <- master_sqr_initialization(
      workers, cfg$tau, bw$b, cfg$radius
    )
    init_lr <- master_lr_initialization(workers, cfg$radius)
    list(
      "CEASE-SQR (h,h)" = fit_cease_sqr_path(
        workers, cfg$tau, bw$h, alpha, max_rounds, cfg$radius, init_sqr
      ),
      "CEASE-SQR (h,b)" = fit_cease_sqr_two_bandwidth_path(
        workers, cfg$tau, bw$h, bw$b, alpha, max_rounds,
        cfg$radius, init_sqr
      ),
      "CSL-SQR (h,h)" = fit_csl_sqr_path(
        workers, cfg$tau, bw$h, bw$h, max_rounds, cfg$radius, init_sqr
      ),
      "CSL-SQR (h,b)" = fit_csl_sqr_path(
        workers, cfg$tau, bw$h, bw$b, max_rounds, cfg$radius, init_sqr
      ),
      "CEASE-LR" = fit_cease_lr_path(
        workers, alpha, max_rounds, cfg$radius, init_lr
      ),
      "CSL-LR" = fit_csl_lr_path(
        workers, max_rounds, cfg$radius, init_lr
      )
    )
  }
  
  select_path_round <- function(paths, round) {
    lapply(paths, function(path) path[, round])
  }
  
  one_bundle <- function(family, replication) {
    family_index <- match(family, np_errors)
    set.seed(20260731 + 100000L * family_index + replication)
    np <- list()
    ni <- 1L
  
    # Baseline data and complete round paths.
    n0 <- cfg$local_n
    K0 <- cfg$K
    alpha0 <- cease_alpha(n0)
    bw0 <- bandwidths(n0, K0)
    w0 <- generate_data(n0, K0, family)
    init_sqr0 <- master_sqr_initialization(
      w0, cfg$tau, bw0$b, cfg$radius
    )
    init_lr0 <- master_lr_initialization(w0, cfg$radius)
    cease_sqr_hh10 <- fit_cease_sqr_path(
      w0, cfg$tau, bw0$h, alpha0, 10L, cfg$radius, init_sqr0
    )
    cease_sqr_hb10 <- fit_cease_sqr_two_bandwidth_path(
      w0, cfg$tau, bw0$h, bw0$b, alpha0, 10L, cfg$radius, init_sqr0
    )
    csl_sqr_hh10 <- fit_csl_sqr_path(
      w0, cfg$tau, bw0$h, bw0$h, 10L, cfg$radius, init_sqr0
    )
    csl_sqr_hb10 <- fit_csl_sqr_path(
      w0, cfg$tau, bw0$h, bw0$b, 10L, cfg$radius, init_sqr0
    )
    cease_lr10 <- fit_cease_lr_path(
      w0, alpha0, 10L, cfg$radius, init_lr0
    )
    csl_lr10 <- fit_csl_lr_path(w0, 10L, cfg$radius, init_lr0)
    baseline_paths <- list(
      "CEASE-SQR (h,h)" = cease_sqr_hh10,
      "CEASE-SQR (h,b)" = cease_sqr_hb10,
      "CSL-SQR (h,h)" = csl_sqr_hh10,
      "CSL-SQR (h,b)" = csl_sqr_hb10,
      "CEASE-LR" = cease_lr10,
      "CSL-LR" = csl_lr10
    )
    for (tt in 1:10) {
      np[[ni]] <- make_rows(
        "Non-private", "Communication rounds", tt, family, n0, K0, tt,
        -1, bw0$h, replication,
        list(
          "CEASE-SQR (h,h)" = cease_sqr_hh10[, tt],
          "CEASE-SQR (h,b)" = cease_sqr_hb10[, tt],
          "CSL-SQR (h,h)" = csl_sqr_hh10[, tt],
          "CSL-SQR (h,b)" = csl_sqr_hb10[, tt],
          "CEASE-LR" = cease_lr10[, tt],
          "CSL-LR" = csl_lr10[, tt]
        )
      )
      ni <- ni + 1L
    }
  
    # Number of machines: N is fixed at 100,000; balanced partitions differ by
    # at most one observation.
    for (K in cfg$machine_grid_nonprivate) {
      n <- cfg$total_N / K
      if (K == K0) {
        w <- w0
        bw <- bw0
        paths <- baseline_paths
      } else {
        w <- generate_data(K = K, family = family, total_N = cfg$total_N)
        bw <- bandwidths(n, K, total_N = cfg$total_N)
        paths <- fit_nonprivate_paths(w, bw, cfg$allocation_rounds)
      }
      fits <- select_path_round(paths, cfg$allocation_rounds)
      np[[ni]] <- make_rows(
        "Non-private", "Number of machines", K, family, n, K,
        cfg$allocation_rounds, -1, bw$h, replication, fits
      )
      ni <- ni + 1L
    }
  
    # Local sample size: K=50 is fixed, so total N varies.
    for (n in cfg$local_n_grid_nonprivate) {
      K <- K0
      if (n == n0) {
        w <- w0
        bw <- bw0
        paths <- baseline_paths
      } else {
        w <- generate_data(n, K, family)
        bw <- bandwidths(n, K)
        paths <- fit_nonprivate_paths(w, bw, cfg$allocation_rounds)
      }
      fits <- select_path_round(paths, cfg$allocation_rounds)
      np[[ni]] <- make_rows(
        "Non-private", "Local sample size", n, family, n, K,
        cfg$allocation_rounds, -1, bw$h, replication, fits
      )
      ni <- ni + 1L
    }
  
    do.call(rbind, np)
  }
  
  run_bundles <- function() {
    tasks <- expand.grid(
      replication = seq_len(n_rep), error = np_errors,
      stringsAsFactors = FALSE
    )
    work <- function(j, tasks) {
      one_bundle(tasks$error[j], tasks$replication[j])
    }
    task_ids <- seq_len(nrow(tasks))
    if (run_mode == "full") {
      checkpoint_dir <- file.path(
        "simulation_results", "nonprivate_methods_checkpoints"
      )
      dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
      task_chunks <- split(
        task_ids,
        ceiling(seq_along(task_ids) / max(n_cores, 1L))
      )
    } else {
      checkpoint_dir <- NULL
      task_chunks <- list(task_ids)
    }
  
    if (n_cores > 1L) {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      export_env <- environment(run_bundles)
      parallel::clusterExport(
        cl, ls(envir = export_env, all.names = TRUE), envir = export_env
      )
    }
  
    out <- vector("list", length(task_ids))
    for (chunk_number in seq_along(task_chunks)) {
      ids <- task_chunks[[chunk_number]]
      checkpoint_file <- if (is.null(checkpoint_dir)) {
        NULL
      } else {
        file.path(
          checkpoint_dir,
          sprintf("chunk_%03d_of_%03d.rds", chunk_number, length(task_chunks))
        )
      }
      if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
        chunk_out <- readRDS(checkpoint_file)
        if (length(chunk_out) != length(ids)) {
          stop("Invalid checkpoint: ", checkpoint_file)
        }
      } else if (n_cores == 1L) {
        chunk_out <- lapply(ids, work, tasks = tasks)
      } else {
        chunk_out <- parallel::parLapplyLB(cl, ids, work, tasks = tasks)
      }
      if (!is.null(checkpoint_file) && !file.exists(checkpoint_file)) {
        saveRDS(chunk_out, checkpoint_file, compress = "gzip")
      }
      out[ids] <- chunk_out
      message(sprintf(
        "Completed task chunk %d/%d (%d/%d bundles).",
        chunk_number, length(task_chunks), max(ids), nrow(tasks)
      ))
    }
    do.call(rbind, out)
  }
  
  summarize_mse <- function(raw) {
    keys <- c(
      "study", "panel", "x_value", "error", "local_n", "machines",
      "rounds", "epsilon", "bandwidth", "method"
    )
    means <- aggregate(raw$mse, raw[keys], mean)
    sds <- aggregate(raw$mse, raw[keys], sd)
    names(means)[ncol(means)] <- "mean_mse"
    names(sds)[ncol(sds)] <- "sd_mse"
    merge(means, sds, by = keys, sort = FALSE)
  }

  private_errors <- c("normal", "student_t", "cauchy")
  cfg$alpha_private <- 0.20
  
  # Transcript-private one-bandwidth path.  The total zCDP budget is divided
  # equally over the gradient and parameter releases in all released rounds.
  fit_dp_cease_path <- function(workers, tau, h, alpha, max_rounds, radius,
                                epsilon, delta, Bx) {
    p <- ncol(workers[[1L]]$X)
    N <- nrow(attr(workers, "X_all"))
    path <- matrix(NA_real_, p, max_rounds)
    beta <- numeric(p)
    rho_release <- rho_from_eps_delta(epsilon, delta) / (2 * max_rounds)
    sigma_g <- (2 * Bx / N) / sqrt(2 * rho_release)
    sigma_beta <- (4 * Bx / (N * alpha)) / sqrt(2 * rho_release)
    for (tt in seq_len(max_rounds)) {
      gbar <- global_sqr_gradient(workers, beta, tau, h) +
        rnorm(p, sd = sigma_g)
      local_solutions <- lapply(workers, function(worker) {
        gk <- smooth_gradient(beta, worker$X, worker$y, tau, h)
        fit_sqr(
          worker$X, worker$y, tau, h, beta, radius,
          linear_term = gbar - gk - alpha * beta,
          ridge = alpha, max_iter = cfg$local_solver_max_iter,
          grad_tol = 5e-6
        )
      })
      clean_aggregate <- weighted_worker_average(local_solutions, workers)
      beta <- project_ball(
        clean_aggregate + rnorm(p, sd = sigma_beta), radius
      )
      path[, tt] <- beta
    }
    path
  }
  
  # Two-bandwidth private path: the released gradient uses h, while the local
  # surrogate loss and local correction gradient use b.
  fit_dp_cease_two_bandwidth_path <- function(
      workers, tau, h, b, alpha, max_rounds, radius, epsilon, delta, Bx
  ) {
    p <- ncol(workers[[1L]]$X)
    N <- nrow(attr(workers, "X_all"))
    path <- matrix(NA_real_, p, max_rounds)
    beta <- numeric(p)
    rho_release <- rho_from_eps_delta(epsilon, delta) / (2 * max_rounds)
    sigma_g <- (2 * Bx / N) / sqrt(2 * rho_release)
    sigma_beta <- (4 * Bx / (N * alpha)) / sqrt(2 * rho_release)
    for (tt in seq_len(max_rounds)) {
      gbar_h <- global_sqr_gradient(workers, beta, tau, h) +
        rnorm(p, sd = sigma_g)
      local_solutions <- lapply(workers, function(worker) {
        gk_b <- smooth_gradient(beta, worker$X, worker$y, tau, b)
        fit_sqr(
          worker$X, worker$y, tau, b, beta, radius,
          linear_term = gbar_h - gk_b - alpha * beta,
          ridge = alpha, max_iter = cfg$local_solver_max_iter,
          grad_tol = 5e-6
        )
      })
      clean_aggregate <- weighted_worker_average(local_solutions, workers)
      beta <- project_ball(
        clean_aggregate + rnorm(p, sd = sigma_beta), radius
      )
      path[, tt] <- beta
    }
    path
  }
  
  # Common random numbers make the head-to-head bandwidth comparison more
  # precise without changing either marginal DP mechanism.
  fit_dp_cease_bandwidth_pair <- function(
      workers, tau, h, b, alpha, max_rounds, radius, epsilon, delta, Bx
  ) {
    rng_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    hh <- fit_dp_cease_path(
      workers, tau, h, alpha, max_rounds, radius, epsilon, delta, Bx
    )
    assign(".Random.seed", rng_state, envir = .GlobalEnv)
    hb <- fit_dp_cease_two_bandwidth_path(
      workers, tau, h, b, alpha, max_rounds, radius, epsilon, delta, Bx
    )
    list(hh = hh, hb = hb)
  }
  
  # The private comparison uses a data-independent zero initialization for both
  # clean and private paths.  It also uses the same public alpha for all four
  # methods so that the displayed difference isolates Gaussian privacy noise.
  fit_private_clean_paths <- function(workers, bw, max_rounds) {
    zero <- numeric(ncol(workers[[1L]]$X))
    list(
      "CEASE-SQR (h,h)" = fit_cease_sqr_path(
        workers, cfg$tau, bw$h, cfg$alpha_private, max_rounds,
        cfg$radius, zero
      ),
      "CEASE-SQR (h,b)" = fit_cease_sqr_two_bandwidth_path(
        workers, cfg$tau, bw$h, bw$b, cfg$alpha_private, max_rounds,
        cfg$radius, zero
      )
    )
  }
  
  private_rows <- function(panel, x_value, family, local_n, K, rounds,
                           epsilon, bandwidth, replication,
                           clean_hh, clean_hb, private_hh, private_hb) {
    make_rows(
      "Private", panel, x_value, family, local_n, K, rounds,
      epsilon, bandwidth, replication,
      list(
        "CEASE-SQR (h,h)" = clean_hh,
        "CEASE-SQR (h,b)" = clean_hb,
        "DP-CEASE-SQR (h,h)" = private_hh,
        "DP-CEASE-SQR (h,b)" = private_hb
      )
    )
  }
  
  one_private_bundle <- function(family, replication) {
    family_index <- match(family, private_errors)
    set.seed(20260820 + 100000L * family_index + replication)
    ans <- list()
    ii <- 1L
  
    n0 <- cfg$local_n
    K0 <- cfg$K
    bw0 <- bandwidths(n0, K0)
    w0 <- generate_data(n0, K0, family)
  
    # Communication rounds.  The DP mechanism is calibrated once for the full
    # ten-round transcript; displayed intermediate iterates are post-processing.
    clean10 <- fit_private_clean_paths(w0, bw0, 10L)
    private10 <- fit_dp_cease_bandwidth_pair(
      w0, cfg$tau, bw0$h, bw0$b, cfg$alpha_private, 10L,
      cfg$radius, cfg$epsilon, cfg$delta, cfg$Bx
    )
    for (tt in 1:10) {
      ans[[ii]] <- private_rows(
        "Communication rounds", tt, family, n0, K0, tt,
        cfg$epsilon, bw0$h, replication,
        clean10[["CEASE-SQR (h,h)"]][, tt],
        clean10[["CEASE-SQR (h,b)"]][, tt],
        private10$hh[, tt], private10$hb[, tt]
      )
      ii <- ii + 1L
    }
  
    # Number of machines: N is fixed at 100,000.  For K=75, worker sizes differ
    # by at most one observation and weighted aggregation is used.
    for (K in cfg$machine_grid_private) {
      n <- cfg$total_N / K
      if (K == K0) {
        w <- w0
        bw <- bw0
        clean <- lapply(clean10, function(path) path[, cfg$rounds, drop = FALSE])
      } else {
        w <- generate_data(K = K, family = family, total_N = cfg$total_N)
        bw <- bandwidths(n, K, total_N = cfg$total_N)
        clean <- fit_private_clean_paths(w, bw, cfg$rounds)
      }
      private <- fit_dp_cease_bandwidth_pair(
        w, cfg$tau, bw$h, bw$b, cfg$alpha_private, cfg$rounds,
        cfg$radius, cfg$epsilon, cfg$delta, cfg$Bx
      )
      ans[[ii]] <- private_rows(
        "Number of machines", K, family, n, K, cfg$rounds,
        cfg$epsilon, bw$h, replication,
        drop(clean[["CEASE-SQR (h,h)"]][, ncol(clean[[1L]])]),
        drop(clean[["CEASE-SQR (h,b)"]][, ncol(clean[[2L]])]),
        private$hh[, cfg$rounds], private$hb[, cfg$rounds]
      )
      ii <- ii + 1L
    }
  
    # Local sample size: K=50 is fixed, so N=K*n varies.
    for (n in cfg$local_n_grid_private) {
      K <- K0
      if (n == n0) {
        w <- w0
        bw <- bw0
        clean <- lapply(clean10, function(path) path[, cfg$rounds, drop = FALSE])
      } else {
        w <- generate_data(n, K, family)
        bw <- bandwidths(n, K)
        clean <- fit_private_clean_paths(w, bw, cfg$rounds)
      }
      private <- fit_dp_cease_bandwidth_pair(
        w, cfg$tau, bw$h, bw$b, cfg$alpha_private, cfg$rounds,
        cfg$radius, cfg$epsilon, cfg$delta, cfg$Bx
      )
      ans[[ii]] <- private_rows(
        "Local sample size", n, family, n, K, cfg$rounds,
        cfg$epsilon, bw$h, replication,
        drop(clean[["CEASE-SQR (h,h)"]][, ncol(clean[[1L]])]),
        drop(clean[["CEASE-SQR (h,b)"]][, ncol(clean[[2L]])]),
        private$hh[, cfg$rounds], private$hb[, cfg$rounds]
      )
      ii <- ii + 1L
    }
  
    # Privacy budget: the clean reference is constant within a replication.
    clean_hh_ref <- clean10[["CEASE-SQR (h,h)"]][, cfg$rounds]
    clean_hb_ref <- clean10[["CEASE-SQR (h,b)"]][, cfg$rounds]
    for (epsilon in c(0.25, 0.5, 1, 2, 4, 8)) {
      private <- fit_dp_cease_bandwidth_pair(
        w0, cfg$tau, bw0$h, bw0$b, cfg$alpha_private, cfg$rounds,
        cfg$radius, epsilon, cfg$delta, cfg$Bx
      )
      ans[[ii]] <- private_rows(
        "Privacy budget", epsilon, family, n0, K0, cfg$rounds,
        epsilon, bw0$h, replication,
        clean_hh_ref, clean_hb_ref,
        private$hh[, cfg$rounds], private$hb[, cfg$rounds]
      )
      ii <- ii + 1L
    }
  
    # Bandwidth: h varies around the Tan--Battey--Zhou global rate, and the local
    # bandwidth is b=max(h,b0).  The zero initializer is held fixed.
    for (multiplier in c(0.5, 0.75, 1, 1.5, 2)) {
      h <- bw0$h * multiplier
      b <- max(h, bw0$b)
      bw <- list(h = h, b = b)
      clean <- fit_private_clean_paths(w0, bw, cfg$rounds)
      private <- fit_dp_cease_bandwidth_pair(
        w0, cfg$tau, h, b, cfg$alpha_private, cfg$rounds,
        cfg$radius, cfg$epsilon, cfg$delta, cfg$Bx
      )
      ans[[ii]] <- private_rows(
        "Bandwidth", h, family, n0, K0, cfg$rounds,
        cfg$epsilon, h, replication,
        clean[["CEASE-SQR (h,h)"]][, cfg$rounds],
        clean[["CEASE-SQR (h,b)"]][, cfg$rounds],
        private$hh[, cfg$rounds], private$hb[, cfg$rounds]
      )
      ii <- ii + 1L
    }
  
    do.call(rbind, ans)
  }
  
  run_private_bundles <- function() {
    tasks <- expand.grid(
      replication = seq_len(n_rep), error = private_errors,
      stringsAsFactors = FALSE
    )
    work <- function(j, tasks) {
      one_private_bundle(tasks$error[j], tasks$replication[j])
    }
    task_ids <- seq_len(nrow(tasks))
  
    if (run_mode == "full") {
      checkpoint_dir <- file.path(
        "simulation_results", "private_alpha020_checkpoints"
      )
      dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
      task_chunks <- split(
        task_ids, ceiling(seq_along(task_ids) / max(n_cores, 1L))
      )
    } else {
      checkpoint_dir <- NULL
      task_chunks <- list(task_ids)
    }
  
    if (n_cores > 1L) {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      export_env <- environment(run_private_bundles)
      parallel::clusterExport(
        cl, ls(envir = export_env, all.names = TRUE), envir = export_env
      )
    }
  
    out <- vector("list", length(task_ids))
    for (chunk_number in seq_along(task_chunks)) {
      ids <- task_chunks[[chunk_number]]
      checkpoint_file <- if (is.null(checkpoint_dir)) NULL else file.path(
        checkpoint_dir,
        sprintf("chunk_%03d_of_%03d.rds", chunk_number, length(task_chunks))
      )
      if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
        chunk_out <- readRDS(checkpoint_file)
        if (length(chunk_out) != length(ids)) {
          stop("Invalid checkpoint: ", checkpoint_file)
        }
      } else if (n_cores == 1L) {
        chunk_out <- lapply(ids, work, tasks = tasks)
      } else {
        chunk_out <- parallel::parLapplyLB(cl, ids, work, tasks = tasks)
      }
      if (!is.null(checkpoint_file) && !file.exists(checkpoint_file)) {
        saveRDS(chunk_out, checkpoint_file, compress = "gzip")
      }
      out[ids] <- chunk_out
      message(sprintf(
        "Completed private task chunk %d/%d (%d/%d bundles).",
        chunk_number, length(task_chunks), max(ids), nrow(tasks)
      ))
    }
    do.call(rbind, out)
  }
}, envir = simulation_environment)
journal_run_nonprivate <- function(mode = "full", cores = 6L, replications = NULL) {
  mode <- tolower(mode)
  if (!mode %in% c("quick", "full", "plot_only")) {
    stop("Simulation mode must be quick, full, or plot_only.", call. = FALSE)
  }
  journal_load_plot_packages()
  dir.create(simulation_result_dir, recursive = TRUE, showWarnings = FALSE)

  sim <- simulation_environment
  sim$run_mode <- mode
  sim$n_rep <- if (is.null(replications)) {
    if (mode == "full") 200L else 5L
  } else {
    journal_integer(replications, 5L, "replications")
  }
  sim$n_cores <- journal_integer(
    cores, if (mode == "full") min(6L, parallel::detectCores()) else 1L,
    "worker processes"
  )

  tag <- if (mode == "full" || mode == "plot_only") "full" else "quick"
  base_name <- paste0("bandwidth_compare_", tag, "_nonprivate")
  summary_file <- file.path(
    simulation_result_dir, paste0(base_name, "_summary.csv")
  )
  raw_file <- file.path(simulation_result_dir, paste0(base_name, "_raw.csv"))
  pdf_file <- file.path(simulation_result_dir, paste0(base_name, ".pdf"))

  if (mode == "plot_only") {
    if (!file.exists(summary_file)) {
      stop("Cannot find completed summary: ", summary_file, call. = FALSE)
    }
    np_sum <- read.csv(summary_file, stringsAsFactors = FALSE)
  } else {
    message(sprintf(
      "Running %d nonprivate replications with %d worker process(es).",
      sim$n_rep, sim$n_cores
    ))
    np_raw <- journal_with_project_directory(sim$run_bundles())
    np_sum <- sim$summarize_mse(np_raw)
    write.csv(np_raw, raw_file, row.names = FALSE)
    write.csv(np_sum, summary_file, row.names = FALSE)
  }

  np_plot <- sim$make_plot_matrix(
    np_sum,
    c("Communication rounds", "Number of machines", "Local sample size"),
    c("normal", "student_t", "cauchy", "mixed_normal"),
    c("CEASE-SQR (h,h)", "CEASE-SQR (h,b)", "CSL-SQR (h,h)",
      "CSL-SQR (h,b)", "CEASE-LR", "CSL-LR")
  )
  sim$save_matrix(np_plot, pdf_file, width = 16.5, height = 11.5)

  if (mode != "plot_only") {
    writeLines(c(
      paste("timestamp:", format(Sys.time(), tz = "Asia/Shanghai")),
      paste("run_mode:", mode),
      paste("replications:", sim$n_rep),
      paste("parallel_workers:", sim$n_cores),
      "rounds panel: T=1,...,10",
      "K and n panels: T=2",
      "baseline: N=100000, K=50, n=2000",
      "initialization: first-worker SQR at b; first-worker OLS for LR",
      "machine grid: {25,50,75,100,125}, N=100000 fixed",
      "local-n grid: {500,1000,1500,2000,2500}, K=50 fixed",
      "errors: normal, Student t, Cauchy, mixed normal",
      "kernel: Gaussian",
      "CEASE proximal parameter: alpha=0.15*p/n",
      paste(
        "methods: CEASE-SQR(h,h), CEASE-SQR(h,b), CSL-SQR(h,h),",
        "CSL-SQR(h,b), CEASE-LR, CSL-LR"
      ),
      "figure: ggplot2 log10(mean MSE), independent y-axis per panel"
    ), file.path(simulation_result_dir, paste0(base_name, "_metadata.txt")))
  }
  message("Nonprivate results written to: ", normalizePath(
    simulation_result_dir, winslash = "/", mustWork = TRUE
  ))
  invisible(list(raw = if (exists("np_raw")) np_raw else NULL,
                 summary = np_sum, plot = np_plot))
}

journal_run_private <- function(mode = "full", cores = 6L, replications = NULL) {
  mode <- tolower(mode)
  if (!mode %in% c("quick", "full", "plot_only")) {
    stop("Simulation mode must be quick, full, or plot_only.", call. = FALSE)
  }
  journal_load_plot_packages()
  dir.create(simulation_result_dir, recursive = TRUE, showWarnings = FALSE)

  sim <- simulation_environment
  sim$run_mode <- mode
  sim$n_rep <- if (is.null(replications)) {
    if (mode == "full") 200L else 5L
  } else {
    journal_integer(replications, 5L, "replications")
  }
  sim$n_cores <- journal_integer(
    cores, if (mode == "full") min(6L, parallel::detectCores()) else 1L,
    "worker processes"
  )

  tag <- if (mode == "full" || mode == "plot_only") "full" else "quick"
  base_name <- paste0("private_methods_", tag)
  summary_file <- file.path(
    simulation_result_dir, paste0(base_name, "_summary.csv")
  )
  raw_file <- file.path(simulation_result_dir, paste0(base_name, "_raw.csv"))
  pdf_file <- file.path(
    simulation_result_dir, paste0("bandwidth_compare_", tag, "_private.pdf")
  )

  if (mode == "plot_only") {
    if (!file.exists(summary_file)) {
      stop("Cannot find completed summary: ", summary_file, call. = FALSE)
    }
    private_sum <- read.csv(summary_file, stringsAsFactors = FALSE)
  } else {
    message(sprintf(
      "Running %d private replications with %d worker process(es).",
      sim$n_rep, sim$n_cores
    ))
    private_raw <- journal_with_project_directory(sim$run_private_bundles())
    private_sum <- sim$summarize_mse(private_raw)
    write.csv(private_raw, raw_file, row.names = FALSE)
    write.csv(private_sum, summary_file, row.names = FALSE)
  }

  private_plot <- sim$make_plot_matrix(
    private_sum,
    c("Communication rounds", "Number of machines", "Local sample size",
      "Privacy budget", "Bandwidth"),
    c("normal", "student_t", "cauchy"),
    c("CEASE-SQR (h,h)", "CEASE-SQR (h,b)",
      "DP-CEASE-SQR (h,h)", "DP-CEASE-SQR (h,b)")
  )
  sim$save_matrix(private_plot, pdf_file, width = 14.5, height = 17.5)

  if (mode != "plot_only") {
    writeLines(c(
      paste("timestamp:", format(Sys.time(), tz = "Asia/Shanghai")),
      paste("run_mode:", mode),
      paste("replications:", sim$n_rep),
      paste("parallel_workers:", sim$n_cores),
      "baseline: N=100000, K=50, n=2000, T=5, epsilon=1, delta=1e-6",
      "initialization: data-independent zero for all four methods",
      "machine grid: {20,50,75,100,125}, N=100000 fixed",
      "local-n grid: {500,1000,1500,2000,2500}, K=50 fixed",
      "errors: normal, Student t, Cauchy",
      "kernel: Gaussian",
      "public private proximal parameter: alpha_private=0.20",
      "gradient sensitivity: 2*Bx/N",
      "certified parameter sensitivity: 4*Bx/(N*alpha_private)",
      paste(
        "methods: CEASE-SQR(h,h), CEASE-SQR(h,b),",
        "DP-CEASE-SQR(h,h), DP-CEASE-SQR(h,b)"
      ),
      paste(
        "DP bandwidth pair uses common Gaussian random numbers within each",
        "dataset/setting; data and privacy noise are fresh across replications"
      ),
      "figure: ggplot2 log10(mean MSE), independent y-axis per panel"
    ), file.path(simulation_result_dir, paste0(base_name, "_metadata.txt")))
  }
  message("Private results written to: ", normalizePath(
    simulation_result_dir, winslash = "/", mustWork = TRUE
  ))
  invisible(list(raw = if (exists("private_raw")) private_raw else NULL,
                 summary = private_sum, plot = private_plot))
}

journal_run_simulation <- function(mode = "full", cores = 6L,
                                   replications = NULL) {
  list(
    nonprivate = journal_run_nonprivate(mode, cores, replications),
    private = journal_run_private(mode, cores, replications)
  )
}

simulation_print_help <- function() {
  cat(paste(
    "",
    "CEASE-SQR simulation reproducibility code",
    "",
    "Usage:",
    "  Rscript --vanilla simulation_methods.R COMPONENT MODE [WORKERS] [REPLICATIONS]",
    "",
    "COMPONENT:",
    "  both | nonprivate | private",
    "",
    "MODE:",
    "  quick (default 5 replications)",
    "  full (default 200 replications)",
    "  plot_only (read completed full summaries and redraw figures)",
    "",
    "Examples:",
    "  Rscript --vanilla simulation_methods.R both full 6 200",
    "  Rscript --vanilla simulation_methods.R nonprivate quick 1 5",
    "  Rscript --vanilla simulation_methods.R private plot_only 1",
    "",
    "The shorthand 'simulation_methods.R full 6 200' runs both components.",
    sep = "\n"
  ))
  invisible(NULL)
}

simulation_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (!length(arguments)) return(simulation_print_help())
  first <- tolower(arguments[[1L]])
  if (first %in% c("help", "-h", "--help")) return(simulation_print_help())

  modes <- c("quick", "full", "plot_only")
  components <- c("both", "nonprivate", "private")

  if (first %in% modes) {
    component <- "both"
    mode <- first
    remaining <- arguments[-1L]
  } else {
    component <- first
    if (!component %in% components) {
      stop("COMPONENT must be both, nonprivate, or private.", call. = FALSE)
    }
    mode <- if (length(arguments) >= 2L) {
      tolower(arguments[[2L]])
    } else {
      "quick"
    }
    remaining <- if (length(arguments) >= 3L) arguments[-c(1L, 2L)] else character()
  }

  if (!mode %in% modes) {
    stop("MODE must be quick, full, or plot_only.", call. = FALSE)
  }
  cores <- journal_integer(
    if (length(remaining) >= 1L) remaining[[1L]] else NULL,
    if (mode == "full") min(6L, parallel::detectCores()) else 1L,
    "worker processes"
  )
  replications <- if (length(remaining) >= 2L) remaining[[2L]] else NULL

  if (component == "both") {
    return(journal_run_simulation(mode, cores, replications))
  }
  if (component == "nonprivate") {
    return(journal_run_nonprivate(mode, cores, replications))
  }
  journal_run_private(mode, cores, replications)
}

if (sys.nframe() == 0L) {
  invisible(simulation_main())
}

