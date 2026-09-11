rbf_analysis_path <- function() {
  if (!requireNamespace("rprojroot", quietly = TRUE)) {
    stop("Package 'rprojroot' is required.")
  }
  
  root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
  
  cond1 <- file.exists(file.path(getwd(), ".Rprofile"))
  cond2 <- basename(dirname(getwd())) == "data" && 
    normalizePath(dirname(dirname(getwd()))) == normalizePath(root)
    
  if (cond1 && cond2) {
    return(getwd())
  }
  
  return(root)
}

rbf_new <- function(name) {
  if (nchar(trimws(name)) == 0) {
    stop("name cannot be empty.")
  }
  if (grepl("\\s", name)) {
    stop("name cannot contain whitespace.")
  }
  if (grepl("[/\\\\]", name)) {
    stop("name cannot contain path separators.")
  }
  
  root <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
  target <- file.path(root, "data", name)
  
  if (dir.exists(target)) {
    stop(paste0("Analysis '", name, "' already exists at ", target))
  }
  
  dir.create(target, recursive = TRUE)
  
  writeLines('source(file.path("..", "..", "R", "source_all.R"))', file.path(target, ".Rprofile"))
  writeLines("{}", file.path(target, "wf_context.json"))
  writeLines(paste0("# Analysis: ", name, "\nCreated: ", Sys.Date(), "\n"), file.path(target, "README.md"))
  
  cat(paste0("Analysis '", name, "' created at data/", name, "/. Open that folder as your working directory, then call init_workflow().\n"))
}

guide <- function(wf) {
  if (!is.null(wf$stage) && wf$stage == "exit") {
    phase <- "Phase: Exit"
    expect <- "You have exited the Bayesian path."
    next_step <- "See exit.yaml for the log."
  } else if (isTRUE(wf$loo_complete)) {
    phase <- "Phase 7: Reporting"
    expect <- "Your analysis is ready to render."
    next_step <- "Open templates/bayesflow_report.qmd"
  } else if (isTRUE(wf$ppc_complete)) {
    phase <- "Phase 6: Model Comparison"
    expect <- "You will see a LOO comparison table."
    next_step <- "Run: wf <- run_phase6(wf, fit)"
  } else if (isTRUE(wf$diagnostics$passed)) {
    phase <- "Phase 5: Posterior Predictive Checks"
    expect <- "You will see PPC density overlay plots."
    next_step <- "Open templates/phase5_ppc.qmd"
  } else if (isFALSE(wf$diagnostics$passed) && !isTRUE(wf$diagnostics$acknowledged)) {
    phase <- "Phase 4: Diagnostics (blocked)"
    expect <- "Diagnostics failed and have not been acknowledged."
    next_step <- "Run: wf <- wf$diagnose()"
  } else if (isFALSE(wf$diagnostics$passed) && isTRUE(wf$diagnostics$acknowledged)) {
    phase <- "Phase 4: Diagnostics (acknowledged)"
    expect <- "Diagnostics failed but have been reviewed."
    next_step <- "Open templates/phase5_ppc.qmd (results carry a caveat)"
  } else if (!is.null(wf$fit_timestamp) && !is.na(wf$fit_timestamp)) {
    phase <- "Phase 4: Diagnostics"
    expect <- "You will see Rhat, ESS, and divergence summaries."
    next_step <- "Open templates/phase4_diagnostics.qmd"
  } else if (!is.null(wf$priors_objects)) {
    phase <- "Phase 3: Model Fitting"
    expect <- "You will see sampling progress and a fit summary."
    next_step <- "Open templates/phase3_fit.qmd"
  } else if (!is.null(wf$declared_goal)) {
    phase <- "Phase 2: Prior Specification"
    expect <- "You will see prior predictive draws."
    next_step <- "Open templates/phase2_priors.qmd"
  } else {
    phase <- "Phase 1: Goal Declaration and Data Inspection"
    expect <- "You will see exploratory plots and off-ramp options."
    next_step <- "Open templates/phase1_exploration.qmd"
  }
  
  cat(paste0("Current phase : ", phase, "\n"))
  cat(paste0("What to expect: ", expect, "\n"))
  cat(paste0("Next step     : ", next_step, "\n"))
  
  invisible(wf)
}

rbf_install <- function(dry_run = FALSE) {
  res <- c(
    r_version = FALSE,
    toolchain = FALSE,
    renv = FALSE,
    renv_restore = FALSE,
    cmdstan_present = FALSE,
    cmdstan_install = FALSE,
    smoke_test = FALSE
  )
  
  if (as.numeric(R.Version()$major) >= 4 && as.numeric(R.Version()$minor) >= 3) {
    cat("✓ Step 1: R version >= 4.3\n")
    res["r_version"] <- TRUE
  } else {
    cat("✗ Step 1: R version >= 4.3\n  → Install R 4.3 or later from https://r-project.org\n")
  }
  
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cat("✗ Step 2: C++ toolchain\n  → Install cmdstanr first — see docs/user-guide/01-installation.md\n")
    stop()
  }
  
  tc_check <- tryCatch({
    cmdstanr::check_cmdstan_toolchain(fix = FALSE)
    TRUE
  }, error = function(e) FALSE)
  
  if (tc_check) {
    cat("✓ Step 2: C++ toolchain\n")
    res["toolchain"] <- TRUE
  } else {
    cat("✗ Step 2: C++ toolchain\n  → Install cmdstanr first — see docs/user-guide/01-installation.md\n")
    stop()
  }
  
  if (requireNamespace("renv", quietly = TRUE)) {
    cat("✓ Step 3: renv installed\n")
    res["renv"] <- TRUE
  } else {
    if (!dry_run) {
      tryCatch({
        install.packages("renv")
      }, error = function(e) NULL)
    }
    if (requireNamespace("renv", quietly = TRUE)) {
      cat("✓ Step 3: renv installed\n")
      res["renv"] <- TRUE
    } else {
      cat("✗ Step 3: renv\n  → Run: install.packages('renv')\n")
    }
  }
  
  if (!dry_run) {
    rest_check <- tryCatch({
      renv::restore(prompt = FALSE)
      TRUE
    }, error = function(e) {
      cat(paste0("  [renv error] ", conditionMessage(e), "\n"))
      FALSE
    })
    if (rest_check) {
      cat("✓ Step 4: renv::restore()\n")
      res["renv_restore"] <- TRUE
    } else {
      cat("✗ Step 4: renv::restore()\n  → See docs/user-guide/01-installation.md\n")
    }
  } else {
    cat("✗ Step 4: renv::restore()\n  → See docs/user-guide/01-installation.md\n")
  }
  
  cs_pres <- tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) FALSE)
  
  if (cs_pres) {
    cat("✓ Step 5: CmdStan present\n")
    res["cmdstan_present"] <- TRUE
  } else {
    cat("✗ Step 5: CmdStan not found\n")
  }
  
  if (res["cmdstan_present"]) {
    cat("✓ Step 6: CmdStan already present — skipped\n")
    res["cmdstan_install"] <- TRUE
  } else {
    if (!dry_run) {
      cs_inst <- tryCatch({
        cmdstanr::install_cmdstan()
        TRUE
      }, error = function(e) FALSE)
      if (cs_inst) {
        cat("✓ Step 6: CmdStan installed\n")
        res["cmdstan_install"] <- TRUE
      } else {
        cat("✗ Step 6: CmdStan install failed\n  → Run cmdstanr::install_cmdstan() manually\n")
      }
    } else {
      cat("✗ Step 6: CmdStan install failed\n  → Run cmdstanr::install_cmdstan() manually\n")
    }
  }
  
  if (!dry_run) {
    tryCatch({
      source("examples/stan_demo.R")
      cat("✓ Step 7: Stan smoke test\n")
      res["smoke_test"] <- TRUE
    }, error = function(e) {
      cat(paste0("✗ Step 7: Stan smoke test failed\n  → ", e$message, "\n"))
    })
  } else {
    cat("✗ Step 7: Stan smoke test failed\n  → Skipped due to dry run\n")
  }
  
  cat(sprintf("RBayesflow environment: %d/7 checks passed.\n", sum(res)))
  invisible(res)
}
