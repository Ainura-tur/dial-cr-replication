# ==============================================================================
# dial_specs.R  -  the 13 empirical instruments: data loaders and settings
# ==============================================================================
#
# Shared by 09_implied_endogeneity.R, 11_diagnostics.R and 12_rstar_gate.R,
# so the three scripts analyse the same rows under the same domains. Sourcing this
# file defines LOADERS, SPECS and dial_load_spec_data(); it runs nothing.
#
# Needs dial_data_path() from foundation/Fin_Empirical5_clean.R, and, for
# MovieLens, .load_ml100k() from 03_phase2_empirical.R.
#
# LOCAL DATA FILES (not in the repository)
#   MovieLens 100K: a folder ml-100k holding u.data and u.item
#   NHANES:         nhanes_mercury.csv, written by 01_nhanes_build.R
# Each is looked for, in this order, at
#   1. the path in the environment variable DIAL_ML_DIR or DIAL_NHANES_CSV
#   2. ML_DIR or NHANES_CSV, if the calling script set them
#   3. the scripts folder, 4. the working directory, 5. <repo>/data
# so the scripts find them whatever folder R starts in. Without the MovieLens
# folder, 03's reader downloads MovieLens from GroupLens; without the NHANES
# file, the two NHANES rows are NA. dial_data_report() lists what was found.
# Papke and Card come from the wooldridge package. Without it, the loaders
# read k401k.csv and card.csv through dial_data_path(); set
# options(dial.data_dir = <folder>) to point them at another folder.
# ==============================================================================

.dial_here <- if (exists(".script_dir")) .script_dir else getwd()
.dial_root <- if (exists(".repo_root")) .repo_root else normalizePath(file.path(.dial_here, ".."))
.dial_first <- function(cands, ok) {
  for (p in cands) if (length(p) == 1L && !is.na(p) && nzchar(p) && ok(p))
    return(normalizePath(p, winslash = "/"))
  NULL
}
.dial_is_ml  <- function(p) file.exists(file.path(p, "u.data")) && file.exists(file.path(p, "u.item"))
ML_DIR <- .dial_first(c(Sys.getenv("DIAL_ML_DIR", ""),
                        if (exists("ML_DIR") && is.character(ML_DIR)) ML_DIR else "",
                        file.path(.dial_here, "ml-100k"), "ml-100k",
                        file.path(.dial_root, "data", "ml-100k")), .dial_is_ml)
.dial_nh <- .dial_first(c(Sys.getenv("DIAL_NHANES_CSV", ""),
                          if (exists("NHANES_CSV") && is.character(NHANES_CSV)) NHANES_CSV else "",
                          file.path(.dial_here, "nhanes_mercury.csv"), "nhanes_mercury.csv",
                          file.path(.dial_root, "data", "nhanes_mercury.csv")), file.exists)
NHANES_CSV <- if (is.null(.dial_nh)) "nhanes_mercury.csv" else .dial_nh

# Which data sources this machine has; prints one line each and returns a
# named logical vector (TRUE = found).
dial_data_report <- function() {
  dta <- c(ditella = "Di_Tella_Schargrodsky_2013.dta", galiani = "Crime.dta",
           banerjee = "yld_sett_aug03.dta", burde = "afgan.dta")
  has_dta <- vapply(dta, function(f) isTRUE(tryCatch(file.exists(dial_data_path(f)),
                                                     error = function(e) FALSE)), logical(1))
  has_wd <- requireNamespace("wooldridge", quietly = TRUE)
  has_csv <- function(f) isTRUE(tryCatch(file.exists(dial_data_path(f)), error = function(e) FALSE))
  ok <- c(has_dta, movielens = !is.null(ML_DIR), nhanes = file.exists(NHANES_CSV),
          papke = has_wd || has_csv("k401k.csv"), card = has_wd || has_csv("card.csv"))
  src <- c(dta, movielens = if (is.null(ML_DIR)) "ml-100k/ not found (03 will try to download it)" else ML_DIR,
           nhanes = if (file.exists(NHANES_CSV)) NHANES_CSV else "nhanes_mercury.csv not found",
           papke = if (has_wd) "wooldridge package" else "k401k.csv",
           card  = if (has_wd) "wooldridge package" else "card.csv")
  cat("Data sources:\n")
  for (k in names(ok)) cat(sprintf("  %-9s %-5s %s\n", k, if (ok[[k]]) "ok" else "MISS", src[[k]]))
  if (!all(ok)) cat("  Rows whose data are missing come out as NA. See the header of dial_specs.R.\n")
  invisible(ok)
}

.dial_unlabel <- function(d) {
  d <- as.data.frame(d)
  for (v in names(d))
    if (inherits(d[[v]], "haven_labelled")) d[[v]] <- as.vector(unclass(d[[v]]))
  d
}
.dial_read_dta <- function(file) .dial_unlabel(haven::read_dta(dial_data_path(file)))
.dial_wooldridge <- function(name, csv) {
  if (requireNamespace("wooldridge", quietly = TRUE)) {
    e <- new.env()
    utils::data(list = name, package = "wooldridge", envir = e)
    return(as.data.frame(get(name, envir = e)))
  }
  utils::read.csv(dial_data_path(csv), stringsAsFactors = FALSE)
}

LOADERS <- list(
  ditella = function() {
    d <- .dial_read_dta("Di_Tella_Schargrodsky_2013.dta")
    d[which(d$offendersPerJudge > 9), , drop = FALSE]
  },
  galiani = function() {
    d <- .dial_read_dta("Crime.dta")
    d[which(d$cohort > 1957 & d$cohort < 1963), , drop = FALSE]
  },
  banerjee = function() {
    d <- .dial_read_dta("yld_sett_aug03.dta")
    d[which(d$year >= 1965 & d$phwht <= 1), , drop = FALSE]
  },
  burde = function() .dial_read_dta("afgan.dta"),
  movielens = function() {
    if (!exists(".load_ml100k", mode = "function"))
      stop("03_phase2_empirical.R not loaded (it provides .load_ml100k)")
    ml <- .load_ml100k(ML_DIR)
    user_ctrl <- ml$ratings %>% dplyr::group_by(user_id) %>%
      dplyr::summarise(user_mean_rating = mean(rating), .groups = "drop")
    movies <- ml$ratings %>% dplyr::left_join(user_ctrl, by = "user_id") %>%
      dplyr::group_by(movie_id) %>%
      dplyr::summarise(mean_rating = mean(rating), n_ratings = dplyr::n(),
                       audience_leniency = mean(user_mean_rating), .groups = "drop") %>%
      dplyr::mutate(log_pop = log(1 + n_ratings)) %>%
      dplyr::left_join(ml$items[, c("movie_id", "num_genres", "release_year", "is_sequel")],
                       by = "movie_id") %>%
      as.data.frame()
    movies$release_year_c <- movies$release_year -
      stats::median(movies$release_year, na.rm = TRUE)
    movies
  },
  papke = function() .dial_wooldridge("k401k", "k401k.csv"),
  nhanes = function() {
    f <- NHANES_CSV
    if (!file.exists(f) && exists(".script_dir")) f <- file.path(.script_dir, NHANES_CSV)
    if (!file.exists(f)) stop("nhanes_mercury.csv not found; run 01_nhanes_build.R")
    utils::read.csv(f, stringsAsFactors = FALSE)
  },
  card = function() .dial_wooldridge("card", "card.csv")
)

C_DT <- c("mostSeriousCrime", "age", "ageSquared", "argentine",
          "numberPreviousImprisonments", "judicialDistrict", "yearOfImprisonment")
C_BU <- c("headchild", "age", "yrsvill", "farsi", "tajik", "farmers", "agehead",
          "educhead", "nhh", "land", "sheep", "distschool", "chagcharan")
C_BA <- c("alt", "totrain", "so_black", "so_red", "so_all", "lat", "coastal",
          "brule1", "year")

# valid = passes both gates of Algorithm 1 in tab:empirical_full: F >= 10
# and the r* gate of foundation/rstar_gate.R (results/rstar_gate.csv, from
# 12_rstar_gate.R); seed = DML seed;
# phase = the output of 02 or 03 that holds the instrument's MCUB set.
SPECS <- list(
  list(study = "Di Tella", instrument = "judgeAlreadyUsedEM", data = "ditella",
       Y = "recidivism", X = "electronicMonitoring", Z = "judgeAlreadyUsedEM",
       controls = C_DT, D = c(0.15, 0.25), valid = TRUE, seed = 1101L, phase = 1L),
  list(study = "Di Tella", instrument = "percJudgeSentToEM", data = "ditella",
       Y = "recidivism", X = "electronicMonitoring", Z = "percJudgeSentToEM",
       controls = C_DT, D = c(0, 0.25), valid = TRUE, seed = 1102L, phase = 1L),
  list(study = "MovieLens", instrument = "num_genres", data = "movielens",
       Y = "mean_rating", X = "log_pop", Z = "num_genres",
       controls = "audience_leniency", D = c(0, 0.8), valid = TRUE, seed = 43L, phase = 2L),
  list(study = "MovieLens", instrument = "release_year", data = "movielens",
       Y = "mean_rating", X = "log_pop", Z = "release_year_c",
       controls = "audience_leniency", D = c(-0.8, 0), valid = FALSE, seed = 44L, phase = 2L),
  list(study = "MovieLens", instrument = "is_sequel", data = "movielens",
       Y = "mean_rating", X = "log_pop", Z = "is_sequel",
       controls = "audience_leniency", D = c(-0.8, 0), valid = FALSE, seed = 45L, phase = 2L),
  list(study = "Galiani", instrument = "highnumber", data = "galiani",
       Y = "crimerate", X = "sm", Z = "highnumber",
       controls = c("cohort", "draftnumber", "navy"), D = c(-0.2, 0),
       valid = TRUE, seed = 1103L, phase = 1L),
  list(study = "Papke", instrument = "sole_plan", data = "papke",
       Y = "prate", X = "mrate", Z = "sole", controls = "ltotemp",
       D = c(0, 0.3), valid = FALSE, seed = 1104L, phase = 2L),
  list(study = "Papke", instrument = "plan_age", data = "papke",
       Y = "prate", X = "mrate", Z = "age", controls = "ltotemp",
       D = c(0, 0.3), valid = FALSE, seed = 1105L, phase = 2L),
  list(study = "NHANES", instrument = "income", data = "nhanes",
       Y = "log_mercury", X = "fish_freq", Z = "income",
       controls = c("age", "female", "bmi"), D = c(0, 0.8), valid = FALSE,
       seed = 1106L, phase = 2L),
  list(study = "NHANES", instrument = "education", data = "nhanes",
       Y = "log_mercury", X = "fish_freq", Z = "education",
       controls = c("age", "female", "bmi"), D = c(0, 0.8), valid = FALSE,
       seed = 1107L, phase = 2L),
  list(study = "Card", instrument = "fatheduc", data = "card",
       Y = "lwage", X = "educ", Z = "fatheduc",
       controls = c("exper", "expersq", "black", "smsa", "south"),
       D = c(0, 0.4), valid = TRUE, seed = 1108L, phase = 2L),
  list(study = "Burde", instrument = "buildschool", data = "burde",
       Y = "testscore", X = "enrolled", Z = "buildschool", controls = C_BU,
       D = c(0, 0.4), valid = FALSE, seed = 1109L, phase = 1L),
  list(study = "Banerjee", instrument = "instru", data = "banerjee",
       Y = "phwht", X = "p_nland", Z = "instru", controls = C_BA,
       D = c(-0.6, -0.4), valid = TRUE, seed = 1110L, phase = 1L)
)

# Load (and cache) the data of a specification; NULL with a message if the
# data are not available on this machine.
.dial_data_cache <- new.env()
dial_load_spec_data <- function(s) {
  if (!is.null(.dial_data_cache[[s$data]])) return(.dial_data_cache[[s$data]])
  d <- tryCatch(LOADERS[[s$data]](), error = function(e) {
    cat("    DATA NOT AVAILABLE:", conditionMessage(e), "\n"); NULL })
  if (!is.null(d)) assign(s$data, d, envir = .dial_data_cache)
  d
}
