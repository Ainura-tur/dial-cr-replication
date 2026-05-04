# ==============================================================================
# nhanes_build.R  -  Build the NHANES analytic frame for DIAL Table 4
# ==============================================================================
#
# PURPOSE
#   Reproduces the n = 1,107 NHANES analytic extract used for the
#   "income" and "education" instrument rows of Table 4 (NHANES Fish
#   Consumption, Zhao et al. 2018 extract).
#
#   This is the script cited by name in the paper:
#       "the NHANES analytic frame is constructed by nhanes_build.R in
#        the supplementary material, which reads data(nhanes.fish) from
#        the CrossScreening R package (the same 1,107-row extract
#        documented in Rosenbaum 2014, Zhao et al. 2018)."
#
# OUTPUT
#   nhanes_mercury.csv  -- 1,107-row analytic frame with columns:
#       log_mercury, fish_freq, income, education, age, female, bmi
#
# DOWNSTREAM CONSUMERS
#   DIAL_NeurIPS_fin_empirical5_phase2a.R::run_nhanes_analysis_fixed()
#   reads this CSV from the working directory.
#
# REQUIRES
#   CrossScreening R package (Zhao et al. 2018; CRAN)
#
# RUN
#   Rscript nhanes_build.R
# ==============================================================================

if (!requireNamespace("CrossScreening", quietly = TRUE)) {
  message("Installing CrossScreening...")
  install.packages("CrossScreening", repos = "https://cloud.r-project.org/")
}

# ---- Load nhanes.fish ------------------------------------------------------

env_tmp <- new.env()
utils::data("nhanes.fish", package = "CrossScreening", envir = env_tmp)
raw <- env_tmp$nhanes.fish
if (!is.data.frame(raw)) raw <- as.data.frame(raw)

cat(sprintf("CrossScreening::nhanes.fish loaded: %d rows, %d cols\n",
            nrow(raw), ncol(raw)))

# Defensive column check -- column names match Zhao et al. 2018 extract
required <- c("o.LBXTHG", "fish", "income", "education", "age", "gender")
missing  <- setdiff(required, names(raw))
if (length(missing) > 0) {
  stop("nhanes.fish is missing expected columns: ",
       paste(missing, collapse = ", "), "\n",
       "Available columns: ", paste(names(raw), collapse = ", "))
}

# ---- Build analytic frame --------------------------------------------------

# Outcome: log total blood mercury (LBXTHG, in ug/L).  We floor at 0.1
# before logging because LBXTHG = 0 occurs in a handful of rows (below
# the assay's lower limit of detection).
log_mercury <- log(pmax(raw$o.LBXTHG, 0.1))

# Treatment: fish meals in last 30 days (12-level ordinal in nhanes.fish)
fish_freq <- raw$fish

# Instruments
income    <- raw$income       # 11-level ordinal (PIR brackets)
education <- raw$education    # 5-level ordinal

# Demographics / controls
age       <- raw$age
female    <- as.integer(raw$gender == 2L)   # NHANES coding: 1 = M, 2 = F

# ---- BMI: use direct column if available, else age x sex impute -----------

bmi_col <- intersect(names(raw),
                     c("bmi", "BMI", "BMXBMI", "o.BMXBMI"))
if (length(bmi_col) > 0) {
  bmi <- raw[[bmi_col[1]]]
  cat(sprintf("BMI column found: '%s' (n missing = %d)\n",
              bmi_col[1], sum(is.na(bmi))))
} else {
  cat("No direct BMI column. Using age x sex cell-median imputation.\n")
  bmi <- NA_real_ * seq_len(nrow(raw))
}

# Age x sex cell-median imputation (replaces NA where present)
# Cell medians from CDC NHANES 2003-2004 reference, the wave from which
# nhanes.fish was extracted.
bmi_lookup <- function(age_v, female_v) {
  band <- cut(age_v, breaks = c(-Inf, 20, 40, 60, Inf),
              labels = c("u20", "20_39", "40_59", "60p"),
              right = FALSE)
  lookup <- list(
    "u20"   = c(M = 24.0, F = 24.5),
    "20_39" = c(M = 27.4, F = 27.8),
    "40_59" = c(M = 29.0, F = 29.4),
    "60p"   = c(M = 28.5, F = 28.9)
  )
  vapply(seq_along(age_v), function(i) {
    b <- as.character(band[i])
    if (is.na(b)) return(NA_real_)
    lookup[[b]][if (isTRUE(female_v[i] == 1L)) "F" else "M"]
  }, numeric(1))
}

needs_impute <- is.na(bmi)
if (any(needs_impute)) {
  bmi[needs_impute] <- bmi_lookup(age[needs_impute], female[needs_impute])
  cat(sprintf("Imputed BMI for %d rows.\n", sum(needs_impute)))
}

# ---- Assemble and write ---------------------------------------------------

analytic <- data.frame(
  log_mercury = log_mercury,
  fish_freq   = fish_freq,
  income      = income,
  education   = education,
  age         = age,
  female      = female,
  bmi         = bmi,
  stringsAsFactors = FALSE
)

# Drop rows with any remaining NA on the analytic columns
n_pre  <- nrow(analytic)
analytic <- analytic[stats::complete.cases(analytic), ]
n_post <- nrow(analytic)
cat(sprintf("Analytic frame: %d rows -> %d rows after complete-case filter.\n",
            n_pre, n_post))

if (n_post != 1107L) {
  warning(sprintf(paste(
    "Analytic frame has %d rows, expected 1,107 (Zhao et al. 2018).",
    "If this is a CrossScreening version skew, downstream Table 4",
    "rows for income/education will be slightly off.", sep = "\n  "),
    n_post))
}

out_path <- "nhanes_mercury.csv"
utils::write.csv(analytic, out_path, row.names = FALSE)
cat(sprintf("Wrote %d rows x %d cols to %s\n",
            nrow(analytic), ncol(analytic), out_path))

cat("\nSummary statistics:\n")
print(summary(analytic))

cat("\nDONE.\n")
