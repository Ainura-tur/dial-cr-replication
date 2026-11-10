#================================================================
# 03_phase2_empirical.R  —  Phase 2 empirical analyses for DIAL
#
# NOTE  This script computes the MCUB sets of tab:empirical_full. Where its
#       comments mention the bootstrap validity rule, the tau_V rule or
#       Stage 2 of Algorithm 1, they mean the MCUB-share rule of app:mc; the
#       paper's Algorithm 1 uses the validity gate of 12_rstar_gate.R. The
#       paper does not include the domain-sensitivity results
#       (tab:domain_sens).
#
# PURPOSE
#   Provides run_*_analysis_fixed() functions for the four application
#   families in tab:empirical_full and tab:domain_sens that are not
#   covered by Phase 1 (02_phase1_empirical.R). All functions follow the
#   identical residualise -> project -> check_compatibility_simple
#   template of the Phase 1 functions in Fin_Empirical5_clean.R.
#
#   Tables are referenced by LaTeX label, not by number, so a float move
#   in the paper cannot make these comments stale.
#
# REQUIRES (source before this file; assumes cwd = <repo>/scripts/):
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("../foundation/Fin_sim3_clean.R")
#   source("../foundation/Fin_Empirical5_clean.R")   # provides check_compatibility_simple
#
# Or simply run 00_setup.R once per R session before any script.
#
# FUNCTIONS
#   run_movielens_item_level()       — 3 instruments, unit = movie (used by the paper)
#   run_movielens_analysis_fixed()   — 3 instruments, unit = rating (for comparison)
#   run_papke_analysis_fixed()       — 2 instruments, canonical + within-sign wide
#   run_nhanes_analysis_fixed()      — 2 instruments (canonical = widest)
#   run_card_analysis_fixed()        — 1 instrument,  canonical + within-sign wide
#   run_all_phase2_applications()    — combines all four into one table
#
#   Every output row carries domain_label = "canonical" or
#   "wide_canonical" (same convention as 02_phase1_empirical.R).
#   The canonical and wide runs of an instrument use the same seed index
#   i, so they share the bootstrap covariance of the correlations and
#   differ only in the domain D.
#
# DATASETS
#   MovieLens 100K  — GroupLens public release (GroupLens Research, U of MN)
#                     https://files.grouplens.org/datasets/movielens/ml-100k/
#   Papke 401(k)    — wooldridge::k401k  (Papke 1995, n = 1,534)
#   NHANES          — nhanes_mercury.csv built by 01_nhanes_build.R
#                     OR CrossScreening::nhanes.fish directly (n = 1,107)
#   Card            — wooldridge::card   (Card 1995; 3,010 rows, 2,320 with
#                     fatheduc recorded)
#
# TARGET PAPER NUMBERS (tab:empirical_full, tab:domain_sens)
#   MCUB sets (CI_Bei, before clipping to [-1, 1]) and the bootstrap share
#   p_hat_B (B = 500), raw-data bootstrap, from the run of 2 November 2026.
#   Seed indices as in the functions below.
#
#   Instrument             D           CI_Bei               p_hat_B  DIAL
#   MovieLens num_genres   (0, 0.8)    [-0.400,  0.695]     1.000    B (sign)  n = 1,682  F = 163.1  R = 0.62
#   MovieLens release_year (-0.8, 0)   [-1.042, -0.244]     0.016    C (weak)  n = 1,681  F =   1.7
#   MovieLens is_sequel    (-0.8, 0)   [-1.072, -0.524]     0.000    C (weak)  n = 1,682  F =   1.3
#     (num_genres: beta_2SLS = -0.038, beta_DML = 0.061. The estimates have
#      opposite signs, so the sign gate of Algorithm 1 returns B (sign);
#      R = 0.62 < 1 / tau_R would also give B under the symmetric rule.)
#   Papke sole_plan        (0, 0.3)    [ 0.101,  0.964]     0.284    C
#   Papke sole_plan  wide  (0, 0.8)    [ 0.101,  1.082]     0.258    C
#   Papke plan_age         (0, 0.3)    [ 0.663,  1.021]     0.000    C
#   Papke plan_age   wide  (0, 0.8)    [ 0.685,  1.023]     0.000    C
#   NHANES income          (0, 0.8)    [ 0.511,  1.053]     0.000    C
#   NHANES education       (0, 0.8)    [ 0.157,  1.015]     0.062    C
#   Card fatheduc          (0, 0.4)    [-0.071,  0.584]     0.814    C         n = 2,320
#   Card fatheduc    wide  (0, 0.8)    [-0.051,  0.904]     0.798    C         n = 2,320
#
#   The paper clips endpoints beyond +/-1 to +/-1.000. Card's sets contain
#   zero, but p_hat_B < tau_V = 0.95, so DIAL returns C; Card is the one
#   instrument on which the rule and the single-sample check disagree.
#   Card's n is the 2,320 complete cases on fatheduc; 3,010 is the raw row
#   count.
#
#   A covariance bootstrap that resamples the prepared rows gives narrower
#   sets (for example Card [0.063, 0.465]), because it understates the
#   standard errors (see RAW-DATA BOOTSTRAP below).
#
# BOOTSTRAP VALIDITY RULE (prop:finite_sample)
#   This file sources ../foundation/bootstrap_rule.R, which wraps
#   check_compatibility_simple(). Every .run_cr() row therefore also
#   carries the tau_V rule of Algorithm 1: B_rule = 500 bootstrap
#   resamples, one MCUB set per resample,
#     p_hat_B    = share of resamples whose MCUB set contains zero
#     Valid_tau2 = "Valid" if p_hat_B >= tau_V = 0.95, else "Invalid"
#   with the strict verdict (zero_tol = 0). The DIAL column of the
#   MovieLens function and the summary below use Valid_tau2; the
#   single-sample check (zero in the full-sample set) stays in
#   Zero_in_CI and DIAL_single. Resamples run in parallel through
#   dial_boot_plan() and are cached in results/boot_rule/.
#   To skip the rule: options(dial.boot_rule = FALSE) before sourcing.
#
# RAW-DATA BOOTSTRAP
#   Every bootstrap resample redraws the rows of the analysis data and
#   redoes the residualise-and-project step: the rule's resamples, the
#   covariance bootstrap inside each resample's MCUB set, and the
#   covariance bootstrap behind the full-sample MCUB set (CI_Bei).
#   .residualise_and_project() attaches the analysis data to the prepared
#   frame for this purpose. options(dial.boot_unit = "rows") before
#   sourcing restores the earlier behaviour (resample the prepared rows).
#
# NOTE ON n (MovieLens)
#   The paper reports an item-level analysis (n = 1,682 movies;
#   tab:empirical_full, footnote dagger). run_movielens_item_level() does
#   that analysis and is the default in run_all_phase2_applications().
#   The rating-level function run_movielens_analysis_fixed() is kept for
#   comparison only: it uses n = 99,991 ratings for num_genres and a
#   10,000-rating subsample for the other two, and its i.i.d. bootstrap
#   treats ratings of the same movie as independent.
#   With the prepared-rows covariance, the rating-level version gives
#   [-0.109, 0.731], [-1.006, -0.566] and [-0.978, -0.405].
#================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(Matrix)
})

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# Bootstrap validity rule (see header)
.br_file <- file.path("..", "foundation", "bootstrap_rule.R")
if (!file.exists(.br_file) && exists("FOUNDATION_DIR"))
  .br_file <- file.path(FOUNDATION_DIR, "bootstrap_rule.R")
if (!file.exists(.br_file))
  stop("03_phase2_empirical.R: foundation/bootstrap_rule.R not found ",
       "(expected cwd = <repo>/scripts/).", call. = FALSE)
options(dial.zero_tol = 0)          # strict verdict in the empirical tables
source(.br_file)
dial_boot_plan()
rm(.br_file)

#================================================================
# HELPERS
#================================================================

# Shared residualise-and-project pipeline (same logic as the inline code
# in the Phase 1 functions).
.residualise_and_project <- function(data, Y, X, Z, controls) {
  vars <- c(Y, X, Z, controls)
  d    <- data[complete.cases(data[, vars, drop = FALSE]), vars, drop = FALSE]
  fml  <- function(lhs)
    stats::as.formula(paste(lhs, "~", paste(controls, collapse = "+")))
  y    <- stats::resid(stats::lm(fml(Y), data = d))
  x    <- stats::resid(stats::lm(fml(X), data = d))
  z0   <- stats::resid(stats::lm(fml(Z), data = d))
  z    <- stats::predict(stats::lm(z0 ~ x + y))
  df   <- data.frame(x = x, y = y, z = z)
  # Analysis data for the raw-data bootstrap (bootstrap_rule.R)
  if (exists("dial_attach_raw", mode = "function"))
    df <- dial_attach_raw(df, d, Y, X, Z, controls)
  list(df = df, n = nrow(d))
}

# Run the CR test on a prepared frame and tag the result.
.run_cr <- function(rp, i, alpha, rxu, instrument, domain_label) {
  out <- check_compatibility_simple(rp$df,
                                    i         = i,
                                    alpha     = alpha,
                                    rxu_range = rxu)
  if (!is.null(out)) {
    out$Instrument   <- instrument
    out$n_sample     <- rp$n
    out$rxu_range    <- sprintf("c(%.1f, %.1f)", rxu[1], rxu[2])
    out$domain_label <- domain_label
  }
  out
}

#================================================================
# APPLICATION 1 — MOVIELENS 100K
#
# Structural story
#   Y  = user rating of a movie (1–5 stars)
#   X  = log(1 + total ratings received by movie) — popularity proxy;
#        endogenous because popular movies attract both fans and
#        casual viewers, creating selection on unobserved quality
#   Z1 = num_genres — number of genre tags on the movie; more genres
#        → broader recommended audience → higher exposure (first stage)
#        but genre count does not directly determine perceived quality
#   Z2 = release_year — older films are more widely known and have
#        longer recommendation history; age of release ≠ content quality
#   Z3 = is_sequel — sequels attract fans of the original franchise,
#        creating selection in who watches; sequel status does not change
#        the quality of the content itself
#   Control = rater leniency (user_mean_rating, each user's mean rating)
#
# Sign priors (rxu_range)
#   num_genres   → positive domain (0, 0.8)
#   release_year → negative domain (-0.8, 0)
#   is_sequel    → negative domain (-0.8, 0)
#   Each canonical domain is already the widest admissible subinterval
#   of its elicited half-line, so there is no within-sign wide run.
#
# TWO VERSIONS
#   run_movielens_item_level()      — unit = movie. This is the analysis
#       tab:empirical_full describes: all three instruments are movie
#       features, so the movie is the unit of variation.
#         Y       = mean rating of the movie
#         X       = log(1 + number of ratings of the movie)
#         control = mean leniency (user_mean_rating) of the movie's raters
#       n = 1,682 movies (1,681 for release_year, which is missing for
#       one movie). It also computes the first-stage F, beta_IV (2SLS),
#       beta_DML (DoubleML PLR with ranger, 500 trees, min node size 5,
#       5 folds, 3 repetitions, as in the other applications) and
#       R = |beta_IV| / |beta_DML|, and applies Algorithm 1 in order:
#       F < 10 -> C (weak); zero outside MCUB -> C (invalid); else A or B.
#   run_movielens_analysis_fixed()  — unit = rating (n = 99,991 ratings
#       with a release year for num_genres; a 10,000-rating subsample
#       for release_year and is_sequel). Its i.i.d.
#       bootstrap treats ratings of the same movie as independent.
#
# DATA
#   Both versions read u.data and u.item from ml_dir if given (the
#   unzipped ml-100k folder); otherwise they download both files from
#   GroupLens.
#================================================================

.load_ml100k <- function(ml_dir = NULL) {
  if (!is.null(ml_dir)) {
    # Windows "Extract All" nests the files one level down
    # (ml-100k/ml-100k/u.data), so look there too.
    cand <- c(ml_dir, file.path(ml_dir, "ml-100k"))
    hit  <- cand[file.exists(file.path(cand, "u.data")) &
                 file.exists(file.path(cand, "u.item"))]
    if (!length(hit))
      stop("u.data / u.item not found in ml_dir = ", ml_dir,
           " (or its ml-100k subfolder). Current folder: ", getwd())
    ml_dir <- hit[1]
    f_data <- file.path(ml_dir, "u.data")
    f_item <- file.path(ml_dir, "u.item")
    cat(sprintf("  Reading MovieLens 100K from %s\n", ml_dir))
    ratings_src <- f_data
    items_raw   <- readLines(f_item, encoding = "latin1")
  } else {
    base_url <- "https://files.grouplens.org/datasets/movielens/ml-100k/"
    cat("  Downloading u.data and u.item from GroupLens...\n")
    ratings_src <- url(paste0(base_url, "u.data"))
    items_raw   <- tryCatch(
      readLines(url(paste0(base_url, "u.item")), encoding = "latin1"),
      error = function(e) stop("MovieLens item download failed: ", conditionMessage(e))
    )
  }

  # u.data: user_id <tab> movie_id <tab> rating <tab> timestamp, in file order
  ratings <- tryCatch(
    utils::read.table(ratings_src, sep = "\t", header = FALSE,
                      col.names  = c("user_id", "movie_id", "rating", "timestamp"),
                      colClasses = c("integer", "integer", "numeric", "numeric")),
    error = function(e) stop("MovieLens ratings read failed: ", conditionMessage(e))
  )
  ratings <- ratings[, c("user_id", "movie_id", "rating")]

  items <- .parse_ml100k_items(items_raw)
  cat(sprintf("  %d ratings, %d movies (sequels detected: %d, %.1f%%)\n",
              nrow(ratings), nrow(items), sum(items$is_sequel, na.rm = TRUE),
              100 * mean(items$is_sequel, na.rm = TRUE)))
  list(ratings = ratings, items = items)
}

# First-stage F, 2SLS and DML estimates, and R = |beta_IV| / |beta_DML|.
# DML follows the specification used for the other applications:
# partially linear model, ranger (500 trees, min node size 5), 5 folds,
# 3 repetitions.
.iv_dml_ratio <- function(d, Y, X, Z, controls, seed = 42L) {
  vars <- c(Y, X, Z, controls)
  d    <- d[stats::complete.cases(d[, vars, drop = FALSE]), vars, drop = FALSE]
  ctrl <- paste(controls, collapse = " + ")

  fml_iv  <- stats::as.formula(paste(Y, "~", X, "+", ctrl, "|", Z, "+", ctrl))
  iv_mod  <- ivreg::ivreg(fml_iv, data = d)
  beta_iv <- unname(stats::coef(iv_mod)[X])
  F_first <- tryCatch(
    summary(iv_mod, diagnostics = TRUE)$diagnostics["Weak instruments", "statistic"],
    error = function(e) {
      full <- stats::lm(stats::as.formula(paste(X, "~", Z, "+", ctrl)), data = d)
      rest <- stats::lm(stats::as.formula(paste(X, "~", ctrl)), data = d)
      stats::anova(rest, full)$F[2]
    }
  )

  if (!requireNamespace("mlr3learners", quietly = TRUE))
    stop("mlr3learners is required for the DML estimate (run 00_setup.R).")
  # Hide the per-fold INFO messages from mlr3; restore the level on exit
  if (requireNamespace("lgr", quietly = TRUE)) {
    lg  <- lgr::get_logger("mlr3")
    old <- lg$threshold
    lg$set_threshold("warn")
    on.exit(lg$set_threshold(old), add = TRUE)
  }
  set.seed(seed)
  dd <- DoubleML::DoubleMLData$new(
    data   = data.table::as.data.table(d[, c(Y, X, controls)]),
    y_col  = Y, d_cols = X, x_cols = controls
  )
  dm <- DoubleML::DoubleMLPLR$new(
    dd,
    ml_l    = mlr3::lrn("regr.ranger", num.trees = 500, min.node.size = 5),
    ml_m    = mlr3::lrn("regr.ranger", num.trees = 500, min.node.size = 5),
    n_folds = 5, n_rep = 3
  )
  dm$fit()
  beta_dml <- as.numeric(dm$coef)

  list(F = unname(F_first), beta_iv = beta_iv, beta_dml = beta_dml,
       R = abs(beta_iv) / abs(beta_dml), n = nrow(d))
}

.parse_ml100k_items <- function(item_lines) {
  # u.item is pipe-delimited; cols 6-24 are binary genre flags
  genre_names <- c("unknown","Action","Adventure","Animation",
                   "Childrens","Comedy","Crime","Documentary","Drama",
                   "Fantasy","FilmNoir","Horror","Musical","Mystery",
                   "Romance","SciFi","Thriller","War","Western")

  parsed <- lapply(item_lines, function(ln) {
    parts <- strsplit(ln, "|", fixed = TRUE)[[1]]
    if (length(parts) < 24) return(NULL)
    list(
      movie_id     = as.integer(parts[1]),
      title        = trimws(parts[2]),
      release_date = trimws(parts[3]),
      genres       = as.integer(parts[6:24])
    )
  })
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]

  items <- data.frame(
    movie_id     = vapply(parsed, `[[`, integer(1), "movie_id"),
    title        = vapply(parsed, `[[`, character(1), "title"),
    release_date = vapply(parsed, `[[`, character(1), "release_date"),
    stringsAsFactors = FALSE
  )

  # Genre count
  genre_mat <- do.call(rbind, lapply(parsed, `[[`, "genres"))
  items$num_genres <- rowSums(genre_mat, na.rm = TRUE)

  # Release year from "DD-Mon-YYYY" format
  items$release_year <- suppressWarnings(
    as.integer(sub(".*-(\\d{4})$", "\\1", items$release_date))
  )

  # is_sequel heuristic: title patterns that reliably indicate a sequel
  clean_title <- sub("\\s*\\(\\d{4}\\)\\s*$", "", items$title)
  sequel_rx   <- paste(
    "\\bII\\b", "\\bIII\\b", "\\bIV\\b", "\\bVI\\b", "\\bVII\\b",
    "\\bVIII\\b", "\\bIX\\b",
    # Space-bounded digit 2-9 (avoids matching years or addresses)
    "(?<![0-9])\\b[2-9]\\b(?![0-9])",
    "Part\\s+(2|Two|II|III|IV)",
    "Chapter\\s+(2|Two|II)",
    "Volume\\s*(2|Two|II)",
    "Episode\\s+(II|III|IV|V|VI|VII|VIII|IX|[2-9])",
    "Returns?\\b", "Revenge\\b",
    "Strikes\\s+Back", "Return\\s+of\\s+the",
    sep = "|"
  )
  items$is_sequel <- as.integer(
    grepl(sequel_rx, clean_title, ignore.case = FALSE, perl = TRUE)
  )

  items
}

run_movielens_analysis_fixed <- function(
    instruments = c("num_genres", "release_year", "is_sequel"),
    alpha       = 0.05,
    # max_n_movie_features: if not NULL, randomly sample this many ratings
    # for release_year and is_sequel. Set NULL to use the full dataset for
    # all instruments.
    max_n_movie_features = 10000L,
    seed        = 42L,
    i_start     = 1L,
    ml_dir      = NULL
) {
  cat("\n=== MovieLens 100K Analysis (rating level) ===\n")
  set.seed(seed)

  # ---- Load u.data and u.item (local folder or GroupLens) ---------------
  ml      <- .load_ml100k(ml_dir)
  ratings <- ml$ratings
  items   <- ml$items

  # ---- Build analytic frame --------------------------------------------
  # Movie-level: log popularity
  movie_pop <- ratings %>%
    group_by(movie_id) %>%
    summarise(log_pop = log(1 + n()), .groups = "drop")

  # User-level: mean rating (control for rating leniency)
  user_ctrl <- ratings %>%
    group_by(user_id) %>%
    summarise(user_mean_rating = mean(rating, na.rm = TRUE), .groups = "drop")

  dat <- ratings %>%
    left_join(items[, c("movie_id","num_genres","release_year","is_sequel")],
              by = "movie_id") %>%
    left_join(movie_pop,  by = "movie_id") %>%
    left_join(user_ctrl,  by = "user_id")  %>%
    filter(complete.cases(.))

  cat(sprintf("  Analytic frame: %d rating-level rows.\n", nrow(dat)))

  # Center release_year at median to improve numerical stability
  dat$release_year_c <- dat$release_year - stats::median(dat$release_year, na.rm = TRUE)

  # Variable map: instrument name → column name in dat
  z_col_map <- list(
    num_genres   = list(col = "num_genres",    rxu = c( 0.0,  0.8)),
    release_year = list(col = "release_year_c", rxu = c(-0.8,  0.0)),
    is_sequel    = list(col = "is_sequel",      rxu = c(-0.8,  0.0))
  )

  Y        <- "rating"
  X        <- "log_pop"
  controls <- "user_mean_rating"

  res_list <- vector("list", length(instruments))

  for (j in seq_along(instruments)) {
    instr <- instruments[j]
    spec  <- z_col_map[[instr]]
    if (is.null(spec)) {
      warning("Unknown MovieLens instrument: ", instr); next
    }

    Z_col <- spec$col; rxu <- spec$rxu

    # Subsample for movie-feature instruments (release_year, is_sequel)
    dat_j <- dat
    if (instr %in% c("release_year","is_sequel") &&
        !is.null(max_n_movie_features) &&
        nrow(dat) > max_n_movie_features) {
      cat(sprintf("  [%s] Subsampling to %d rows (seed = %d).\n",
                  instr, max_n_movie_features, seed))
      dat_j <- dat[sample(nrow(dat), max_n_movie_features), ]
    }

    cat(sprintf("\n  [%s] n = %d | rxu in [%.1f, %.1f]\n",
                instr, nrow(dat_j), rxu[1], rxu[2]))

    # Rename Z column for .residualise_and_project
    dat_j$Z_instr <- dat_j[[Z_col]]

    rp <- .residualise_and_project(dat_j, Y, X, "Z_instr", controls)
    res_list[[j]] <- .run_cr(rp, i = i_start + j - 1L, alpha = alpha,
                             rxu = rxu, instrument = instr,
                             domain_label = "canonical")
  }

  res <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
  rownames(res) <- NULL
  res
}

run_movielens_item_level <- function(
    instruments   = c("num_genres", "release_year", "is_sequel"),
    alpha         = 0.05,
    i_start       = 1L,
    ml_dir        = NULL,
    min_ratings   = 1L,       # keep movies with at least this many ratings
    compute_ratio = TRUE,     # first-stage F, 2SLS, DML and R
    seed          = 42L,
    tau_R         = 1.5
) {
  cat("\n=== MovieLens 100K Analysis (item level) ===\n")

  ml      <- .load_ml100k(ml_dir)
  ratings <- ml$ratings
  items   <- ml$items

  # Rater leniency: each user's mean rating over all their ratings
  user_ctrl <- ratings %>%
    group_by(user_id) %>%
    summarise(user_mean_rating = mean(rating), .groups = "drop")

  # One row per movie
  movies <- ratings %>%
    left_join(user_ctrl, by = "user_id") %>%
    group_by(movie_id) %>%
    summarise(mean_rating       = mean(rating),
              n_ratings         = n(),
              audience_leniency = mean(user_mean_rating),
              .groups = "drop") %>%
    mutate(log_pop = log(1 + n_ratings)) %>%
    left_join(items[, c("movie_id", "num_genres", "release_year", "is_sequel")],
              by = "movie_id") %>%
    filter(n_ratings >= min_ratings) %>%
    as.data.frame()

  # Center release_year at its median over movies (numerical stability)
  movies$release_year_c <- movies$release_year -
    stats::median(movies$release_year, na.rm = TRUE)

  cat(sprintf("  Item-level frame: %d movies (min_ratings = %d)\n",
              nrow(movies), min_ratings))

  z_col_map <- list(
    num_genres   = list(col = "num_genres",     rxu = c( 0.0,  0.8)),
    release_year = list(col = "release_year_c", rxu = c(-0.8,  0.0)),
    is_sequel    = list(col = "is_sequel",      rxu = c(-0.8,  0.0))
  )

  Y        <- "mean_rating"
  X        <- "log_pop"
  controls <- "audience_leniency"

  res_list <- vector("list", length(instruments))

  for (j in seq_along(instruments)) {
    instr <- instruments[j]
    spec  <- z_col_map[[instr]]
    if (is.null(spec)) {
      warning("Unknown MovieLens instrument: ", instr); next
    }

    d <- movies
    d$Z_instr <- d[[spec$col]]

    cat(sprintf("\n  [%s] rxu in [%.1f, %.1f]\n", instr, spec$rxu[1], spec$rxu[2]))

    rp  <- .residualise_and_project(d, Y, X, "Z_instr", controls)
    out <- .run_cr(rp, i = i_start + j - 1L, alpha = alpha,
                   rxu = spec$rxu, instrument = instr,
                   domain_label = "canonical")

    if (!is.null(out) && compute_ratio) {
      est   <- .iv_dml_ratio(d, Y, X, "Z_instr", controls, seed = seed + j)
      # Stage 2 of Algorithm 1: the bootstrap rule (p_hat_B >= tau_V) when
      # it ran; the single-sample check is kept in DIAL_single.
      # Stage 3: the sign gate (estimates of opposite sign give "B (sign)")
      # and the symmetric ratio rule (B when |log R| >= log tau_R).
      valid_1 <- grepl("✓", out$Zero_in_CI, fixed = TRUE)
      valid   <- if (!is.null(out$Valid_tau2) && !is.na(out$Valid_tau2))
                   out$Valid_tau2 == "Valid" else valid_1
      .verdict <- function(v)
        if (is.na(est$F) || est$F < 10) "C (weak)"
        else if (!v) "C (invalid)"
        else if (is.na(est$beta_dml) || is.na(est$R)) "?"
        else if (sign(est$beta_iv) != sign(est$beta_dml)) "B (sign)"
        else if (abs(log(est$R)) >= log(tau_R)) "B" else "A"
      out$F_first     <- round(est$F, 1)
      out$beta_2sls   <- round(est$beta_iv, 4)
      out$beta_DML    <- round(est$beta_dml, 4)
      out$R           <- round(est$R, 2)
      out$DIAL        <- .verdict(valid)
      out$DIAL_single <- .verdict(valid_1)
      cat(sprintf("  [%s] F = %.1f  beta_2SLS = %.4f  beta_DML = %.4f  R = %.2f  p_hat_B = %s  -> %s\n",
                  instr, est$F, est$beta_iv, est$beta_dml, est$R,
                  if (is.null(out$p_hat_B)) "n/a" else sprintf("%.3f", out$p_hat_B),
                  out$DIAL))
    }
    res_list[[j]] <- out
  }

  res <- dplyr::bind_rows(res_list[!vapply(res_list, is.null, logical(1))])
  rownames(res) <- NULL
  res
}

#================================================================
# APPLICATION 2 — PAPKE (1995) 401(k) PARTICIPATION
#
# Structural story
#   Y  = prate  — participation rate in 401(k) plan (0–1)
#   X  = mrate  — employer match rate (dollars matched per dollar
#                 contributed); endogenous because firms that set
#                 high match rates are also more likely to actively
#                 promote plan participation through other means
#   Z1 = sole   — indicator: firm offers no other retirement plan;
#        [paper name: sole_plan] workers in sole-plan firms have
#        stronger incentive to contribute, but sole-plan status is
#        driven by firm governance, not worker preferences
#   Z2 = age    — age of the 401(k) plan in years;
#        [paper name: plan_age] older plans have had longer to build
#        participation norms, but plan age is set at inception
#   Controls = ltotemp — log total employment (firm size)
#
# Sign prior: positive half-line
#   Unobserved firm-quality variable U (e.g., management quality) is
#   positively correlated with X (better-managed firms offer higher
#   match) and positively correlated with Y (participation). So
#   rho_DU > 0. Canonical (0, 0.3) reflects that match rates are
#   moderately but not fully driven by management quality; the
#   within-sign wide run uses (0, 0.8).
#
# Paper results (tab:empirical_full, tab:domain_sens)
#   sole_plan: canonical [0.101, 0.964], wide [0.101, 1.000]  p_hat_B 0.284, 0.258  C
#   plan_age:  canonical [0.663, 1.000], wide [0.685, 1.000]  p_hat_B 0.000, 0.000  C
#   (upper endpoints 1.082, 1.021 and 1.023 before clipping)
#   Both verdicts are robust to within-sign widening.
#================================================================

run_papke_analysis_fixed <- function(alpha    = 0.05,
                                     i_start  = 1L,
                                     rxu_can  = c(0.0, 0.3),
                                     rxu_wide = c(0.0, 0.8),
                                     run_wide = TRUE) {
  cat("\n=== Papke (1995) 401(k) Analysis ===\n")
  cat("  Instruments: sole_plan (= sole), plan_age (= age)\n")
  cat(sprintf("  Canonical D: c(%.1f, %.1f); within-sign wide D: c(%.1f, %.1f)\n\n",
              rxu_can[1], rxu_can[2], rxu_wide[1], rxu_wide[2]))

  # ---- Load data -------------------------------------------------------
  if (!requireNamespace("wooldridge", quietly = TRUE)) {
    message("  Installing wooldridge package...")
    utils::install.packages("wooldridge",
                            repos = "https://cloud.r-project.org/",
                            quiet = TRUE)
  }
  dat <- wooldridge::k401k

  # Defensive column check
  required_cols <- c("prate", "mrate", "sole", "age", "ltotemp")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("wooldridge::k401k missing expected columns: ",
         paste(missing_cols, collapse = ", "),
         "\nInspect names(wooldridge::k401k) and update required_cols.")
  }

  cat(sprintf("  wooldridge::k401k loaded: %d rows\n", nrow(dat)))

  Y        <- "prate"
  X        <- "mrate"
  controls <- "ltotemp"

  # Instrument specs: paper name → data column
  instr_spec <- list(
    sole_plan = "sole",
    plan_age  = "age"
  )

  res_list <- list()

  for (j in seq_along(instr_spec)) {
    instr_name <- names(instr_spec)[j]
    Z_col      <- instr_spec[[j]]
    i_j        <- i_start + j - 1L

    dat$Z_instr <- dat[[Z_col]]
    rp <- .residualise_and_project(dat, Y, X, "Z_instr", controls)

    cat(sprintf("  [%s] data col = '%s' | canonical c(%.1f, %.1f)\n",
                instr_name, Z_col, rxu_can[1], rxu_can[2]))
    res_list[[length(res_list) + 1L]] <-
      .run_cr(rp, i = i_j, alpha = alpha, rxu = rxu_can,
              instrument = instr_name, domain_label = "canonical")

    if (run_wide) {
      cat(sprintf("  [%s] within-sign wide c(%.1f, %.1f)\n",
                  instr_name, rxu_wide[1], rxu_wide[2]))
      res_list[[length(res_list) + 1L]] <-
        .run_cr(rp, i = i_j, alpha = alpha, rxu = rxu_wide,
                instrument = instr_name, domain_label = "wide_canonical")
    }
  }

  res <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
  rownames(res) <- NULL
  res
}

#================================================================
# APPLICATION 3 — NHANES FISH CONSUMPTION
#
# Structural story
#   Y  = log_mercury  — log total blood mercury (ug/L); the health
#                       outcome of interest; mercury bioaccumulates
#                       through fish consumption
#   X  = fish_freq    — total fish servings over 30 days; the
#                       endogenous treatment; individuals who eat
#                       more fish are also likely to eat more
#                       healthy foods generally, confounding the
#                       mercury–health relationship
#   Z1 = income       — poverty-income ratio; wealthier individuals
#                       eat more fish (first stage) but income also
#                       correlates with other healthy-eating behaviours
#                       (mercury-independent channels → exclusion
#                       restriction concern)
#   Z2 = education    — ordinal 1–5; same story as income
#   Controls = age, female, bmi
#     nhanes.fish has no BMI column, so bmi is the age-band x sex cell
#     median imputed by 01_nhanes_build.R (a deterministic function of
#     age band and sex).
#
# Sign prior: positive half-line, canonical (0, 0.8)
#   Income and education correlate positively with unobserved
#   healthy-eating habits (U), which also affect mercury exposure (e.g.,
#   via supplement use, fruit/vegetable intake). So rho_DU > 0. The
#   canonical (0, 0.8) is already the widest admissible positive
#   subinterval, so there is no within-sign wide run.
#
# Paper results (tab:empirical_full, MCUB sets)
#   income:    [0.511, 1.000]  p_hat_B 0.000  C  (upper endpoint 1.053 before clipping)
#   education: [0.157, 1.000]  p_hat_B 0.062  C  (upper endpoint 1.015 before clipping)
#
# DATA LOADING ORDER
#   1. If data_url points to a .csv (built by 01_nhanes_build.R): read.csv
#   2. Else if CrossScreening is installed: use nhanes.fish directly
#   3. Else: error with instructions
#================================================================

run_nhanes_analysis_fixed <- function(
    data_url = "nhanes_mercury.csv",    # path to CSV from 01_nhanes_build.R
    alpha    = 0.05,
    i_start  = 1L
) {
  cat("\n=== NHANES Fish Consumption Analysis ===\n")
  cat("  Instruments: income, education\n")
  cat("  rxu_range: c(0, 0.8) for both (canonical = widest)\n\n")

  # ---- Load data -------------------------------------------------------
  # Priority 1: CSV from 01_nhanes_build.R
  dat <- tryCatch({
    if (file.exists(data_url)) {
      cat(sprintf("  Loading from CSV: %s\n", data_url))
      d <- utils::read.csv(data_url, stringsAsFactors = FALSE)
      cat(sprintf("  %d rows loaded.\n", nrow(d)))
      d
    } else {
      stop("CSV not found")
    }
  }, error = function(e) {

    # Priority 2: CrossScreening package
    cat("  CSV not found; trying CrossScreening::nhanes.fish...\n")
    if (!requireNamespace("CrossScreening", quietly = TRUE)) {
      stop(
        "Data not found. Either:\n",
        "  (a) Run 01_nhanes_build.R first to create nhanes_mercury.csv, OR\n",
        "  (b) Install CrossScreening: install.packages('CrossScreening')\n",
        "      then pass data_url = NULL."
      )
    }
    # Load from CrossScreening
    env_tmp <- new.env()
    utils::data("nhanes.fish", package = "CrossScreening", envir = env_tmp)
    raw <- env_tmp$nhanes.fish
    if (!is.data.frame(raw)) raw <- as.data.frame(raw)

    # Require the columns 01_nhanes_build.R would have produced
    needed <- c("o.LBXTHG", "fish", "income", "education", "age", "gender")
    missing <- setdiff(needed, names(raw))
    if (length(missing) > 0)
      stop("nhanes.fish missing columns: ", paste(missing, collapse = ", "))

    bmi_col <- grep("bmi", tolower(names(raw)), value = TRUE, fixed = TRUE)
    if (length(bmi_col) > 0) {
      raw$bmi <- raw[[bmi_col[1]]]
    } else {
      # Age-sex BMI imputation (same lookup as 01_nhanes_build.R)
      bmi_lookup <- function(age, female) {
        band <- cut(age, breaks = c(-Inf,20,40,60,Inf),
                    labels = c("u20","20_39","40_59","60p"), right = FALSE)
        lookup <- list(u20 = c(M=24.0, F=24.5), `20_39` = c(M=27.4, F=27.8),
                       `40_59` = c(M=29.0, F=29.4), `60p` = c(M=28.5, F=28.9))
        vapply(seq_along(age), function(i) {
          b <- as.character(band[i])
          if (is.na(b)) NA_real_
          else lookup[[b]][if (female[i] == 1L) "F" else "M"]
        }, numeric(1))
      }
      raw$bmi <- bmi_lookup(raw$age, as.integer(raw$gender == 2L))
    }

    data.frame(
      log_mercury = log(pmax(raw$o.LBXTHG, 0.1)),
      fish_freq   = raw$fish,
      income      = raw$income,
      education   = raw$education,
      age         = raw$age,
      female      = as.integer(raw$gender == 2L),
      bmi         = raw$bmi,
      stringsAsFactors = FALSE
    )
  })

  # Defensive column check
  required_cols <- c("log_mercury","fish_freq","income","education",
                     "age","female","bmi")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0)
    stop("NHANES data missing columns: ", paste(missing_cols, collapse = ", "))

  Y        <- "log_mercury"
  X        <- "fish_freq"
  controls <- c("age", "female", "bmi")

  instr_spec <- list(
    income    = list(rxu = c(0.0, 0.8)),
    education = list(rxu = c(0.0, 0.8))
  )

  res_list <- vector("list", length(instr_spec))

  for (j in seq_along(instr_spec)) {
    instr_name <- names(instr_spec)[j]
    rxu        <- instr_spec[[instr_name]]$rxu

    cat(sprintf("  [%s] rxu in [%.1f, %.1f]\n", instr_name, rxu[1], rxu[2]))

    dat$Z_instr <- dat[[instr_name]]
    rp <- .residualise_and_project(dat, Y, X, "Z_instr", controls)
    res_list[[j]] <- .run_cr(rp, i = i_start + j - 1L, alpha = alpha,
                             rxu = rxu, instrument = instr_name,
                             domain_label = "canonical")
  }

  res <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
  rownames(res) <- NULL
  res
}

#================================================================
# APPLICATION 4 — CARD (1995) RETURNS TO EDUCATION
#
# Structural story
#   Y  = lwage    — log hourly wage
#   X  = educ     — years of schooling; endogenous because ability
#                   affects both education decisions and wages
#   Z  = fatheduc — father's years of schooling; affects child's
#        [paper name: fatheduc] schooling through credit constraints,
#        peer networks, and parental investment (first stage) but
#        Card argues father's education does not directly cause
#        the child's wage beyond its effect on the child's schooling
#   Controls = exper, expersq, black, smsa, south
#
# Sign prior: positive half-line
#   Unobserved ability U is positively correlated with years of
#   schooling (higher ability → more education) and positively
#   correlated with wages (higher ability → higher productivity).
#   So rho_DU > 0. Canonical (0, 0.4); within-sign wide (0, 0.8).
#
# Paper results (tab:empirical_full, tab:domain_sens)
#   fatheduc: canonical [-0.071, 0.584], wide [-0.051, 0.904]
#   Both sets contain zero, but only 81% and 80% of the resampled sets do
#   (p_hat_B = 0.814 and 0.798 < tau_V = 0.95). DIAL therefore returns C
#   in both within-sign domains, so the verdict is robust.
#
# DATA
#   wooldridge::card — 3,010 rows. fatheduc is missing for 690 of them,
#   so the complete-case analytic sample is n = 2,320. Report this n.
#
# CROSS-SIGN DIAGNOSTIC (optional, off by default)
#   run_cross_sign = TRUE adds a c(-0.5, 0.5) run. Extending D across the
#   sign boundary is NOT a within-sign robustness check (cor:domain, part
#   c): it re-admits rho_DU = 0. The row is labelled
#   "cross_sign_sensitivity" and is excluded from the summary table.
#================================================================

run_card_analysis_fixed <- function(alpha          = 0.05,
                                    i_start        = 1L,
                                    rxu_can        = c(0.0, 0.4),
                                    rxu_wide       = c(0.0, 0.8),
                                    run_cross_sign = FALSE) {
  cat("\n=== Card (1995) Returns to Education Analysis ===\n")
  cat("  Instrument: fatheduc\n")
  cat(sprintf("  Canonical D: c(%.1f, %.1f); within-sign wide D: c(%.1f, %.1f)\n\n",
              rxu_can[1], rxu_can[2], rxu_wide[1], rxu_wide[2]))

  # ---- Load data -------------------------------------------------------
  if (!requireNamespace("wooldridge", quietly = TRUE)) {
    message("  Installing wooldridge package...")
    utils::install.packages("wooldridge",
                            repos = "https://cloud.r-project.org/",
                            quiet = TRUE)
  }
  dat <- wooldridge::card

  required_cols <- c("lwage","educ","fatheduc","exper","expersq",
                     "black","smsa","south")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("wooldridge::card missing columns: ",
         paste(missing_cols, collapse = ", "),
         "\nInspect names(wooldridge::card).")
  }

  n_raw       <- nrow(dat)
  n_complete  <- sum(complete.cases(dat[, required_cols]))
  cat(sprintf("  wooldridge::card: %d total rows, %d complete on key vars.\n",
              n_raw, n_complete))

  Y        <- "lwage"
  X        <- "educ"
  Z        <- "fatheduc"
  controls <- c("exper", "expersq", "black", "smsa", "south")

  dat$Z_instr <- dat[[Z]]
  rp <- .residualise_and_project(dat, Y, X, "Z_instr", controls)

  # Canonical and wide runs share the seed index i_start.
  cat(sprintf("  [fatheduc] canonical c(%.1f, %.1f)\n", rxu_can[1], rxu_can[2]))
  out_c <- .run_cr(rp, i = i_start, alpha = alpha, rxu = rxu_can,
                   instrument = "fatheduc", domain_label = "canonical")

  cat(sprintf("  [fatheduc] within-sign wide c(%.1f, %.1f)\n", rxu_wide[1], rxu_wide[2]))
  out_w <- .run_cr(rp, i = i_start, alpha = alpha, rxu = rxu_wide,
                   instrument = "fatheduc", domain_label = "wide_canonical")

  out_list <- list(out_c, out_w)

  if (run_cross_sign) {
    cat("  [fatheduc] cross-sign diagnostic c(-0.5, 0.5) (NOT a robustness check)\n")
    out_s <- .run_cr(rp, i = i_start + 99L, alpha = alpha, rxu = c(-0.5, 0.5),
                     instrument = "fatheduc",
                     domain_label = "cross_sign_sensitivity")
    out_list[[3L]] <- out_s
  }

  res <- do.call(rbind, out_list[!vapply(out_list, is.null, logical(1))])
  rownames(res) <- NULL

  # Print the within-sign robustness verdict
  if (!is.null(out_c) && !is.null(out_w)) {
    .v <- function(o) if (!is.null(o$Valid_tau2) && !is.na(o$Valid_tau2)) o$Valid_tau2
                      else ifelse(grepl("✓", o$Zero_in_CI), "Valid", "Invalid")
    v_can  <- .v(out_c)
    v_wide <- .v(out_w)
    cat(sprintf("\n  Within-sign robustness: canonical = %s, wide = %s -> %s\n",
                v_can, v_wide,
                if (v_can == v_wide) "robust" else "VERDICT FLIPS"))
  }

  res
}

#================================================================
# RUN ALL PHASE 2 APPLICATIONS
#================================================================

run_all_phase2_applications <- function(
    nhanes_csv      = "nhanes_mercury.csv",
    alpha           = 0.05,
    verbose         = TRUE,
    movielens_level = c("item", "rating"),   # item = the analysis the paper reports
    ml_dir          = NULL                   # unzipped ml-100k folder, or NULL to download
) {
  movielens_level <- match.arg(movielens_level)
  cat("\n", strrep("=", 68), "\n", sep = "")
  cat("DIAL Phase 2 — All four application families\n")
  cat(strrep("=", 68), "\n\n", sep = "")

  all_results <- list()

  safe_run <- function(label, expr) {
    cat(sprintf("\n%s\n", strrep("-", 60)))
    res <- tryCatch(expr, error = function(e) {
      message("  ERROR in ", label, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(res) && nrow(res) > 0) res$Study <- label
    res
  }

  all_results$movielens <- safe_run(
    "MovieLens 100K",
    if (movielens_level == "item")
      run_movielens_item_level(alpha = alpha, ml_dir = ml_dir)
    else
      run_movielens_analysis_fixed(alpha = alpha, ml_dir = ml_dir)
  )
  all_results$papke     <- safe_run("Papke 401(k)",
                                    run_papke_analysis_fixed(alpha = alpha))
  all_results$nhanes    <- safe_run("NHANES Fish",
                                    run_nhanes_analysis_fixed(data_url = nhanes_csv,
                                                              alpha = alpha))
  all_results$card      <- safe_run("Card Returns to Educ.",
                                    run_card_analysis_fixed(alpha = alpha))

  # Keep canonical and within-sign wide rows; drop cross-sign diagnostics.
  all_results <- lapply(all_results, function(r) {
    if (is.null(r)) return(NULL)
    r[r$domain_label %in% c("canonical", "wide_canonical"), , drop = FALSE]
  })
  all_results <- all_results[!vapply(all_results, is.null, logical(1))]

  if (!length(all_results)) {
    warning("All Phase 2 analyses failed.")
    return(NULL)
  }

  combined <- dplyr::bind_rows(all_results)
  rownames(combined) <- NULL

  # ---- Summary table ---------------------------------------------------
  cat("\n\n", strrep("=", 68), "\n", sep = "")
  cat("PHASE 2 SUMMARY (canonical and within-sign wide domains)\n")
  cat(strrep("=", 68), "\n\n", sep = "")

  # Determine verdict column (handles both naming conventions)
  zi_col <- if ("Zero_in_CI" %in% names(combined)) "Zero_in_CI"
  else if ("Zero_MCUB" %in% names(combined)) "Zero_MCUB"
  else NA_character_

  ci_col <- if ("CI_Bei" %in% names(combined)) "CI_Bei"
  else if ("CI_MCUB" %in% names(combined)) "CI_MCUB"
  else NA_character_

  if (!is.na(zi_col) && !is.na(ci_col)) {
    has_rule <- "Valid_tau2" %in% names(combined)
    for (i in seq_len(nrow(combined))) {
      r       <- combined[i, ]
      # Verdict from the bootstrap rule when it ran, else the single-sample check
      ok <- if (has_rule && !is.na(r$Valid_tau2)) r$Valid_tau2 == "Valid"
            else grepl("✓", r[[zi_col]], fixed = TRUE)
      verdict <- if (ok) "VALID (A/B)" else "INVALID (C)"
      cat(sprintf("  %-14s %-15s %-12s n=%-6d %s  p_hat_B=%s  %s\n",
                  r$Instrument,
                  r$domain_label,
                  r$rxu_range,
                  as.integer(r$n_sample),
                  r[[ci_col]],
                  if (has_rule) sprintf("%.3f", r$p_hat_B) else "n/a",
                  verdict))
    }
  } else {
    print(combined)
  }

  # ---- Comparison with the paper ---------------------------------------
  cat("\n\nPaper targets (MCUB sets; tab:empirical_full, tab:domain_sens):\n")
  targets <- data.frame(
    instrument   = c("num_genres", "release_year", "is_sequel",
                     "sole_plan", "sole_plan", "plan_age", "plan_age",
                     "income", "education", "fatheduc", "fatheduc"),
    domain_label = c("canonical", "canonical", "canonical",
                     "canonical", "wide_canonical", "canonical", "wide_canonical",
                     "canonical", "canonical", "canonical", "wide_canonical"),
    # MovieLens targets are the item-level sets (movielens_level = "item").
    # Endpoints before clipping to [-1, 1]; raw-data bootstrap, B = 500.
    lo = c(-0.400, -1.042, -1.072, 0.101, 0.101, 0.663, 0.685,
            0.511,  0.157, -0.071, -0.051),
    hi = c( 0.695, -0.244, -0.524, 0.964, 1.082, 1.021, 1.023,
            1.053,  1.015,  0.584,  0.904),
    p_hat_B = c(1.000, 0.016, 0.000, 0.284, 0.258, 0.000, 0.000,
                0.000, 0.062, 0.814, 0.798),
    verdict = c("Valid", rep("Invalid", 10)),
    stringsAsFactors = FALSE
  )
  for (k in seq_len(nrow(targets))) {
    tgt <- targets[k, ]
    cat(sprintf("  %-14s %-15s paper: [%+.3f, %+.3f]  p_hat_B=%.3f  (%s)\n",
                tgt$instrument, tgt$domain_label, tgt$lo, tgt$hi,
                tgt$p_hat_B, tgt$verdict))
  }

  cat("\n")
  tryCatch({
    write.csv(combined, "phase2_results.csv", row.names = FALSE)
    cat("Results written to: phase2_results.csv\n")
  }, error = function(e) invisible(NULL))

  invisible(combined)
}

#================================================================
# QUICK-START REFERENCE
#
#   DIAL_SKIP_AUTORUN <- TRUE
#   source("../foundation/Fin_sim3_clean.R")
#   source("../foundation/Fin_Empirical5_clean.R")
#   source("03_phase2_empirical.R")
#
#   # Individual applications (ml_dir = unzipped ml-100k folder, or NULL
#   # to download from GroupLens):
#   res_ml     <- run_movielens_item_level(ml_dir = NULL)      # used by the paper
#   res_ml_r   <- run_movielens_analysis_fixed(ml_dir = NULL)  # rating level, comparison only
#   res_papke  <- run_papke_analysis_fixed()          # canonical + wide
#   res_nhanes <- run_nhanes_analysis_fixed("nhanes_mercury.csv")
#   res_card   <- run_card_analysis_fixed()           # canonical + wide
#
#   # Full Phase 2 table:
#   res_p2 <- run_all_phase2_applications()
#
#   # Combine Phase 1 + Phase 2 into one table:
#   p1 <- run_phase1_all()                # from 02_phase1_empirical.R
#   p2 <- run_all_phase2_applications()
#================================================================
