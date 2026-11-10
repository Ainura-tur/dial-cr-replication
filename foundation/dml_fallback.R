# ==============================================================================
# dml_fallback.R  -  partially linear DML with random forests, without DoubleML
# ==============================================================================
#
# PURPOSE
#   The paper's DML estimates come from DoubleML (DoubleMLPLR, score
#   "partialling out", ranger learners). This file gives the same estimator
#   for machines where DoubleML or mlr3 cannot be installed:
#     1. split the rows into K folds;
#     2. on each fold, predict Y and D from the controls with random forests
#        (ranger) trained on the other folds;
#     3. regress the residual of Y on the residual of D, pooled over folds
#        (DML2), which gives theta;
#     4. with several sample splits, report the median theta and the
#        standard error aggregated as the R package DoubleML (1.0.2) does.
#   The estimator is the same as DoubleML's; the random numbers are not, so
#   estimates differ from DoubleML's in the sampling noise of the splits and
#   forests. Scripts 06, 09 and 11 use DoubleML when it is installed and this
#   file otherwise, and they say which one they used. The forests use the
#   settings of mlr3's regr.ranger: ranger defaults apart from num.trees and
#   min.node.size.
#
# FUNCTIONS
#   dial_dml_plr_rf(d, Y, D, controls, folds = NULL, n_rep = 1L, n_folds = 5L,
#                   num.trees = 500L, min.node.size = 5L)
#     returns list(coef, se, coefs) where coefs holds one theta per split.
#   dial_dml_available()  TRUE when DoubleML, mlr3 and mlr3learners load.
#   dial_dml_engine_report(require = dial_require_doubleml())
#     prints the engine, the R version, the package versions and the library
#     they load from, and returns "DoubleML" or "dml_fallback.R". With
#     require = TRUE it stops when DoubleML does not load, naming the
#     packages that failed. Set the environment variable
#     DIAL_REQUIRE_DOUBLEML=1 to make 06, 09 and 11 refuse the fallback.
#   options(dial.ranger_threads = k) lets each forest use k threads (default 1,
#   which is safe inside parallel workers).
# ==============================================================================

dial_dml_available <- function() {
  all(vapply(c("DoubleML", "mlr3", "mlr3learners"), requireNamespace,
             logical(1), quietly = TRUE))
}

dial_require_doubleml <- function()
  tolower(Sys.getenv("DIAL_REQUIRE_DOUBLEML", "0")) %in% c("1", "true", "yes")

dial_dml_engine_report <- function(require = dial_require_doubleml()) {
  pk <- c("DoubleML", "mlr3", "mlr3learners", "ranger")
  ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
  ver <- vapply(pk, function(p) if (ok[[p]]) as.character(utils::packageVersion(p)) else "missing",
                character(1))
  engine <- if (all(ok[1:3])) "DoubleML" else "dml_fallback.R"
  cat(sprintf("DML engine: %s   (%s)\n", engine, R.version.string))
  cat("  ", paste(sprintf("%s %s", pk, ver), collapse = ", "), "\n", sep = "")
  if (ok[["DoubleML"]])
    cat("  DoubleML loads from", dirname(find.package("DoubleML")), "\n")
  cat("  R libraries:", paste(.libPaths(), collapse = " | "), "\n")
  if (engine != "DoubleML") {
    msg <- paste0("DoubleML does not load in this R (missing: ",
                  paste(pk[1:3][!ok[1:3]], collapse = ", "), "). Run the script with an R ",
                  "that has DoubleML, e.g. the full path to its Rscript, or install it.")
    if (require) stop(msg, " DIAL_REQUIRE_DOUBLEML is set, so the script stops.", call. = FALSE)
    cat("  NOTE:", msg, "\n  Using foundation/dml_fallback.R (same estimator, other random numbers).\n")
  }
  invisible(engine)
}

dial_dml_plr_rf <- function(d, Y, D, controls, folds = NULL, n_rep = 1L,
                            n_folds = 5L, num.trees = 500L, min.node.size = 5L) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("dial_dml_plr_rf() needs the ranger package")
  d  <- as.data.frame(d)
  X  <- d[, controls, drop = FALSE]
  y  <- d[[Y]]; dd <- d[[D]]
  n  <- nrow(d)
  fit_pred <- function(target, train, test) {
    tr <- data.frame(.t = target[train], X[train, , drop = FALSE])
    m  <- ranger::ranger(.t ~ ., data = tr, num.trees = num.trees,
                         min.node.size = min.node.size,
                         num.threads = getOption("dial.ranger_threads", 1L))
    stats::predict(m, data = X[test, , drop = FALSE])$predictions
  }
  one_split <- function(fold) {
    K <- max(fold)
    u <- v <- numeric(n)
    for (k in seq_len(K)) {
      te <- which(fold == k); tr <- which(fold != k)
      u[te] <- y[te]  - fit_pred(y,  tr, te)
      v[te] <- dd[te] - fit_pred(dd, tr, te)
    }
    theta <- sum(u * v) / sum(v * v)
    psi   <- (u - theta * v) * v
    J     <- -mean(v * v)
    c(theta = theta, se = sqrt(mean(psi^2) / J^2 / n))
  }
  res <- if (!is.null(folds)) {
    matrix(one_split(folds), nrow = 1L)
  } else {
    t(vapply(seq_len(n_rep), function(r)
      one_split(sample(rep_len(seq_len(n_folds), n))), numeric(2)))
  }
  colnames(res) <- c("theta", "se")
  theta <- stats::median(res[, "theta"])
  # Aggregation across splits as in the R package DoubleML 1.0.2
  # (agg_cross_fit): the median of n se^2 + (theta_s - theta)^2, divided by n
  se <- sqrt(stats::median(n * res[, "se"]^2 + (res[, "theta"] - theta)^2) / n)
  list(coef = theta, se = se, coefs = res[, "theta"])
}
