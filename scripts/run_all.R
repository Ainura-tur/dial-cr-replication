# ==============================================================================
# run_all.R  -  run the replication scripts in order, each in a fresh R
# ==============================================================================
#
# PURPOSE
#   Runs the numbered scripts in the order that the paper's tables and
#   figures need them. Each script runs in its own Rscript process, started
#   in the scripts folder, so no objects, options or parallel workers carry
#   over from one script to the next. Each script's console output goes to
#   results/logs/<script>.log. The runner checks the setup before it starts,
#   stops at the first script that fails, and prints the end of that log.
#   Run 00_setup.R once before the first run.
#
# STEPS
#   PRESET chooses the steps; STEPS overrides it with a list of step numbers.
#   The steps always run in the order of ORDER below (01, 02, 03, 04, 12, 09,
#   11, 06, 13, 14, 05, 07, 10, 08), whatever order STEPS lists them in.
#     "rerun" (default): the steps that use the validity gate or DoubleML.
#       They read the outputs of 01-04, which must already exist.
#         12  validity gate for the 13 instruments     under a minute
#         09  r* and the intervals for R (DoubleML)    about 20 minutes
#         11  diagnostics and labels (DoubleML)        about 20 minutes
#         06  structural simulation, rows SIM_ROWS     about 30 minutes
#         13  Monte Carlo of the gate                  about 1-1.5 hours
#         14  HSIC and the gate over repeated samples  about 12-17 minutes
#         08  Figures 1-4                              under 2 minutes
#     "paper": every script behind the paper, from the start: 01, 02, 03 and
#       04, then the steps above, with all seven rows of 06. About 6-7 hours.
#   Neither preset runs 05, 07 and 10, the Monte Carlo studies of the
#   MCUB-share rule (app:mc), because no table or figure uses their output.
#   List them in STEPS to run them (about 3, 2.7 and 8.5 hours).
#
# DATA
#   02, 09, 11 and 12 read the four .dta files in data/ at the repository
#   root (the folder next to scripts/ and foundation/). The runner stops
#   before it starts if one of them is missing, because those scripts would
#   otherwise drop the Di Tella, Galiani, Burde and Banerjee instruments.
#
# DOUBLEML
#   06, 09 and 11 must use DoubleML. Before it starts, the runner checks that
#   DoubleML, mlr3, mlr3learners and ranger load. It sets
#   DIAL_REQUIRE_DOUBLEML=1, so that 09 and 11 stop instead of falling back,
#   and it stops if the log of 06 does not report "DML engine: DoubleML".
#
# OUTPUT
#   results/logs/<script>.log for each step. On the console: the files that
#   each step wrote, and at the end the numbers that the paper quotes from
#   06, 11 and 12.
#
# RUN (from the scripts folder)
#   RStudio:  open run_all.R, check the settings below, choose Session >
#             Set Working Directory > To Source File Location, then
#             source("run_all.R"). To follow a long step, open its log.
#   Terminal: Rscript run_all.R             (the PRESET below)
#             Rscript run_all.R paper
#             Rscript run_all.R 12 13 08    (only these steps)
#             Rscript run_all.R dry         (check the setup, run nothing)
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))

# Namespace protection (avoid MASS/stats masking dplyr verbs)
select    <- dplyr::select;    filter    <- dplyr::filter
mutate    <- dplyr::mutate;    slice     <- dplyr::slice
recode    <- dplyr::recode;    rename    <- dplyr::rename
summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange   <- dplyr::arrange;   count     <- dplyr::count
lag       <- dplyr::lag

# ==============================================================================
# SETTINGS
# ==============================================================================

PRESET     <- "rerun"         # "rerun" or "paper" (see STEPS above)
STEPS      <- NULL            # e.g. c("12", "13", "08") to run only these
SIM_ROWS   <- "A,Rev,RevHid"  # rows of 06 (all seven under "paper"; "" = all)
GATE_PARTS <- "A,B,C,D"       # parts of 13; "A,C,D" skips the power curves
ML_DIR     <- "ml-100k"       # MovieLens folder for 03 (downloaded if absent)
DRY_RUN    <- FALSE           # TRUE: check the setup and list the steps only

# ==============================================================================
# STEP TABLE
# ==============================================================================

ORDER   <- c("01", "02", "03", "04", "12", "09", "11", "06", "13", "14", "05", "07", "10", "08")
PRESETS <- list(rerun = c("12", "09", "11", "06", "13", "14", "08"),
                paper = c("01", "02", "03", "04", "06", "09", "11", "12", "13", "14", "08"))
MINUTES <- c("01" = "1", "02" = "60", "03" = "60", "04" = "45", "05" = "180",
             "06" = "30-60", "07" = "160", "09" = "20", "10" = "510", "11" = "20",
             "12" = "1", "13" = "60-90", "14" = "12-17", "08" = "2")
DML_STEPS <- c("06", "09", "11")
# Files a step reads, named by the step that writes them
NEEDS <- list(
  "11" = c("scripts/phase1_results.csv" = "02", "scripts/phase2_results.csv" = "03",
           "results/implied_endogeneity.csv" = "09"),
  "12" = c("scripts/phase1_results.csv" = "02", "scripts/phase2_results.csv" = "03"),
  "08" = c("results/rstar_gate.csv" = "12", "results/implied_endogeneity.csv" = "09",
           "results/rstar_gate_power.csv" = "13"))
FOUNDATION_FILES <- c("Fin_sim3_clean.R", "Fin_Empirical5_clean.R", "bootstrap_rule.R",
                      "jp_dgp.R", "rstar_gate.R", "dml_fallback.R")
# The bundled data of the Phase 1 studies, read from data/ by 02, 09, 11, 12
DATA_FILES <- c("Di_Tella_Schargrodsky_2013.dta", "Crime.dta", "yld_sett_aug03.dta",
                "afgan.dta")
DATA_STEPS <- c("02", "09", "11", "12")

# ==============================================================================
# PATHS
# ==============================================================================

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    f    <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) && nchar(f)) dirname(normalizePath(f)) else getwd()
  }
)
if (!length(.script_dir) || is.na(.script_dir[1])) .script_dir <- getwd()
SCRIPTS_DIR <- normalizePath(.script_dir, winslash = "/")
REPO_ROOT   <- normalizePath(file.path(SCRIPTS_DIR, ".."), winslash = "/")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
LOG_DIR     <- file.path(RESULTS_DIR, "logs")
RSCRIPT     <- file.path(R.home("bin"),
                         if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

# ==============================================================================
# HELPERS
# ==============================================================================

# The script of a step: the one file in the scripts folder named <step>_*.R
find_script <- function(id) {
  f <- list.files(SCRIPTS_DIR, pattern = paste0("^", id, "_.*\\.R$"))
  if (length(f) != 1L)
    stop(sprintf("step %s: expected one file named %s_*.R in %s, found %s", id, id,
                 SCRIPTS_DIR, if (length(f)) paste(f, collapse = ", ") else "none"),
         call. = FALSE)
  f
}

# Environment variables for a step
step_env <- function(id, preset) {
  e <- c(DIAL_REQUIRE_DOUBLEML = "1")
  if (id == "06") e[["DIAL_SIM_ROWS"]] <- if (identical(preset, "paper")) "" else SIM_ROWS
  if (id == "13") e[["DIAL_GATE_PARTS"]] <- GATE_PARTS
  if (id == "10") e[["DIAL_EDGE_CELL"]] <- "0"
  e
}

# 02 and 03 define functions; this file loads the foundation and runs them
wrapper <- function(id, script) {
  call <- if (id == "02") "run_phase1_all()" else
    sprintf("run_all_phase2_applications(ml_dir = %s)",
            if (dir.exists(file.path(SCRIPTS_DIR, ML_DIR))) deparse(ML_DIR) else "NULL")
  f <- file.path(tempdir(), sprintf("run_step_%s.R", id))
  writeLines(c("DIAL_SKIP_AUTORUN <- TRUE",
               'source("../foundation/Fin_sim3_clean.R")',
               'source("../foundation/Fin_Empirical5_clean.R")',
               sprintf('source("%s")', script),
               sprintf("invisible(%s)", call)), f)
  f
}

# Path relative to the repository root
rel <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  ifelse(startsWith(p, REPO_ROOT), substring(p, nchar(REPO_ROOT) + 2L), p)
}

# Modification times of the files in results/, scripts/ and
# scripts/paper/figures; comparing two snapshots shows what a step wrote
snapshot <- function() {
  dirs <- c(RESULTS_DIR, SCRIPTS_DIR, file.path(SCRIPTS_DIR, "paper", "figures"))
  f <- unlist(lapply(dirs[dir.exists(dirs)], list.files, full.names = TRUE))
  f <- f[!dir.exists(f)]
  stats::setNames(as.numeric(file.mtime(f)), rel(f))
}
written_between <- function(before, after) {
  old <- before[names(after)]
  sort(names(after)[is.na(old) | after > old])
}

# The CSV in results/ that has a given column (output names can change)
csv_with_column <- function(col) {
  for (f in list.files(RESULTS_DIR, pattern = "\\.csv$", full.names = TRUE)) {
    h <- tryCatch(names(utils::read.csv(f, nrows = 1L, check.names = FALSE)),
                  error = function(e) character(0))
    if (col %in% h) return(f)
  }
  NULL
}

nz <- function(x) if (is.null(x) || !length(x)) NA else x

# ==============================================================================
# PREFLIGHT
# ==============================================================================

preflight <- function(steps, preset) {
  problems <- character(0); notes <- character(0)
  if (!file.exists(RSCRIPT)) problems <- c(problems, paste("Rscript not found at", RSCRIPT))
  miss <- FOUNDATION_FILES[!file.exists(file.path(REPO_ROOT, "foundation", FOUNDATION_FILES))]
  if (length(miss))
    problems <- c(problems, paste("missing in foundation/:", paste(miss, collapse = ", ")))
  if (any(DATA_STEPS %in% steps)) {
    miss <- DATA_FILES[!file.exists(file.path(REPO_ROOT, "data", DATA_FILES))]
    if (length(miss))
      problems <- c(problems, sprintf(
        "missing in %s/data: %s. Copy the data folder of the repository there; without it, steps %s drop the Di Tella, Galiani, Burde and Banerjee instruments",
        REPO_ROOT, paste(miss, collapse = ", "),
        paste(intersect(DATA_STEPS, steps), collapse = ", ")))
  }
  for (id in steps) {
    s <- tryCatch(find_script(id), error = function(e) conditionMessage(e))
    if (!grepl("\\.R$", s)) problems <- c(problems, s)
    for (k in names(NEEDS[[id]])) {
      by <- NEEDS[[id]][[k]]
      if (!(by %in% steps) && !file.exists(file.path(REPO_ROOT, k)))
        problems <- c(problems, sprintf("step %s reads %s, which step %s writes; run %s too",
                                        id, k, by, by))
    }
  }
  if ("08" %in% steps && !file.exists(file.path(RESULTS_DIR, "rstar_gate_power.csv")) &&
      !("13" %in% steps && grepl("B", GATE_PARTS)))
    notes <- c(notes, "08 will skip Figure 3: results/rstar_gate_power.csv (13, part B) is missing")
  if (any(c("03", "09", "11", "12") %in% steps)) {
    if (!("01" %in% steps) && !file.exists(file.path(SCRIPTS_DIR, "nhanes_mercury.csv")))
      notes <- c(notes, "scripts/nhanes_mercury.csv is missing (01 writes it): the NHANES rows will be NA")
    if (!dir.exists(file.path(SCRIPTS_DIR, ML_DIR)))
      notes <- c(notes, sprintf("no scripts/%s folder: MovieLens will be downloaded from GroupLens", ML_DIR))
  }
  if (any(DML_STEPS %in% steps)) {
    pk <- c("DoubleML", "mlr3", "mlr3learners", "ranger", "data.table")
    bad <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
    if (length(bad))
      problems <- c(problems, paste("06, 09 and 11 need DoubleML, but these packages do not load:",
                                    paste(bad, collapse = ", ")))
  }
  cat("================================================================================\n")
  cat("run_all.R\n")
  cat("================================================================================\n")
  cat("Scripts folder:", SCRIPTS_DIR, "\n")
  cat("R:             ", RSCRIPT, "\n")
  cat("Steps:\n")
  for (id in steps) {
    s <- tryCatch(find_script(id), error = function(e) "(not found)")
    env <- step_env(id, preset)
    env <- env[names(env) != "DIAL_REQUIRE_DOUBLEML" | id %in% DML_STEPS]
    cat(sprintf("  %s  %-30s about %s min%s\n", id, s, MINUTES[[id]],
                if (length(env)) paste0("  [", paste(sprintf("%s=%s", names(env), env),
                                                     collapse = " "), "]") else ""))
  }
  if (length(notes)) cat(paste0("Note: ", notes), sep = "\n")
  if (length(problems)) {
    cat(paste0("PROBLEM: ", problems), sep = "\n")
    stop("fix the problems above before running", call. = FALSE)
  }
  cat("Setup OK.\n\n")
}

# ==============================================================================
# RUN ONE STEP
# ==============================================================================

run_step <- function(id, preset) {
  script <- find_script(id)
  log <- file.path(LOG_DIR, sub("\\.R$", ".log", script))
  env <- step_env(id, preset)
  old <- Sys.getenv(names(env), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(env))
  on.exit(for (k in names(old)) {
    if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k))
  }, add = TRUE)
  target <- if (id %in% c("02", "03")) wrapper(id, script) else script
  before <- snapshot()
  t0 <- Sys.time()
  cat(sprintf("[%s] step %s: %s (about %s min); log: %s\n", format(t0, "%H:%M"), id, script,
              MINUTES[[id]], rel(log)))
  st <- system2(RSCRIPT, shQuote(target), stdout = log, stderr = log)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  lines <- readLines(log, warn = FALSE)
  if (!identical(as.integer(st), 0L)) {
    cat(sprintf("    FAILED after %.1f min (exit status %s). End of the log:\n", mins, st))
    cat(paste0("      ", utils::tail(lines, 20)), sep = "\n")
    stop(sprintf("step %s (%s) failed; see %s", id, script, rel(log)), call. = FALSE)
  }
  if (id == "06" && !any(grepl("DML engine: DoubleML", lines, fixed = TRUE)))
    stop("06 did not use DoubleML (see the line 'DML engine' in ", rel(log),
         "), so its DML estimates cannot go in the paper", call. = FALSE)
  wrote <- written_between(before, snapshot())
  cat(sprintf("    done in %.1f min; wrote %s\n", mins,
              if (length(wrote)) paste(wrote, collapse = ", ") else "no files in results/ or scripts/"))
  data.frame(step = id, script = script, minutes = round(mins, 1), stringsAsFactors = FALSE)
}

# ==============================================================================
# REPORT
# ==============================================================================

report <- function(done, preset) {
  cat("\n================================================================================\n")
  cat("Summary\n")
  cat("================================================================================\n")
  print(done, row.names = FALSE)
  cat(sprintf("Total: %.1f min\n", sum(done$minutes)))

  if ("12" %in% done$step && file.exists(file.path(RESULTS_DIR, "rstar_gate.csv"))) {
    g <- utils::read.csv(file.path(RESULTS_DIR, "rstar_gate.csv"), stringsAsFactors = FALSE)
    cat(sprintf("\nGate (12): %d pass, %d fail, %d weak (the paper has 6, 5 and 2)\n",
                sum(g$verdict == "pass", na.rm = TRUE), sum(g$verdict == "C (invalid)", na.rm = TRUE),
                sum(g$verdict == "C (weak)", na.rm = TRUE)))
  }

  if ("11" %in% done$step) {
    f <- csv_with_column("label3")
    if (!is.null(f)) {
      d <- utils::read.csv(f, stringsAsFactors = FALSE)
      d <- d[d$valid %in% TRUE & !is.na(d$R), ]
      cat(sprintf("\nInstruments that pass the gate (11, DML engine %s):\n",
                  paste(unique(stats::na.omit(d$dml_engine)), collapse = " + ")))
      for (i in seq_len(nrow(d)))
        cat(sprintf("  %-20s R = %.2f [%.2f, %.2f]  B share %.3f  splits %.2f to %.2f  DML %.4f (SE %.4f)  label %s\n",
                    d$instrument[i], d$R[i], d$R_lo[i], d$R_hi[i], d$share_B[i],
                    d$R_split_min[i], d$R_split_max[i], d$b_DML_fit[i], d$se_DML[i], d$label3[i]))
    }
  }

  rds <- file.path(RESULTS_DIR, "sim_dial_v7_results.rds")
  if ("06" %in% done$step && file.exists(rds)) {
    res  <- readRDS(rds)
    rows <- if (identical(preset, "paper") || !nzchar(SIM_ROWS)) names(res) else
      intersect(trimws(strsplit(SIM_ROWS, ",")[[1]]), names(res))
    cat("\nStructural simulation (06; the LaTeX rows are in its log):\n")
    for (k in rows) {
      r <- res[[k]]
      cat(sprintf("  %-26s DML %.3f  IV %.3f  R %.2f  %-12s strict share %s  engine %s\n",
                  nz(r$label), nz(r$beta_dml), nz(r$beta_iv), nz(r$ratio), nz(r$dial),
                  if (is.na(nz(r$p_hat_strict))) "NA" else sprintf("%.3f", r$p_hat_strict),
                  nz(r$dml_engine)))
    }
  }

  if ("08" %in% done$step)
    cat("\nCopy the four PNGs in scripts/paper/figures next to the .tex.\n")
  cat("\nDONE.\n")
}

# ==============================================================================
# MAIN
# ==============================================================================

main <- function() {
  args   <- commandArgs(trailingOnly = TRUE)
  dry    <- DRY_RUN || "dry" %in% args
  args   <- setdiff(args, "dry")
  preset <- PRESET
  steps  <- STEPS
  if (length(args) == 1L && args %in% names(PRESETS)) preset <- args
  else if (length(args)) steps <- args
  if (is.null(steps)) {
    if (!(preset %in% names(PRESETS)))
      stop("PRESET must be one of: ", paste(names(PRESETS), collapse = ", "), call. = FALSE)
    steps <- PRESETS[[preset]]
  } else preset <- "custom"
  steps <- sprintf("%02d", as.integer(steps))
  bad <- setdiff(steps, ORDER)
  if (length(bad))
    stop("unknown steps: ", paste(bad, collapse = ", "), " (known: ",
         paste(ORDER, collapse = ", "), ")", call. = FALSE)
  steps <- ORDER[ORDER %in% steps]

  old_wd <- setwd(SCRIPTS_DIR)
  on.exit(setwd(old_wd), add = TRUE)
  dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
  preflight(steps, preset)
  if (dry) { cat("Dry run: nothing was run.\n"); return(invisible(NULL)) }

  done <- list()
  for (id in steps) done[[id]] <- run_step(id, preset)
  report(do.call(rbind, done), preset)
}

main()
