# RBayesflow — Installation

Follow these steps in order. Do not skip ahead — each step depends
on the previous one.

---

## Step 1: Install R ≥ 4.3

Download and install R 4.3 or later from
<https://r-project.org>.

On Windows, also install **RTools 4.5** (the C++ toolchain R needs
to compile Stan models):
<https://cran.r-project.org/bin/windows/Rtools/>

After installing RTools, restart R before continuing.

---

## Step 2: Get RBayesflow

Clone the repository:

```r
# In a terminal (not the R console):
# git clone https://github.com/your-org/RBayesflow.git
# cd RBayesflow
```

Or download the ZIP from GitHub and unzip it. Open the resulting
`RBayesflow/` folder as your working directory in RStudio or
Positron.

---

## Step 3: Install cmdstanr (one manual step)

This is the one step `rbf_install()` cannot do for you. Run this
in a fresh R session from the RBayesflow project root:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

Wait for the installation to complete before continuing.

---

## Step 4: Run rbf_install()

Source `R/install.R` and run the environment checker:

```r
source("R/install.R")
rbf_install()
```

The function runs seven checks and prints `✓` or `✗` for each:

| Step | What it checks |
|------|----------------|
| 1 | R version ≥ 4.3 |
| 2 | C++ toolchain (RTools on Windows) |
| 3 | renv installed |
| 4 | renv::restore() — installs all pinned packages |
| 5 | CmdStan present |
| 6 | CmdStan install (if absent) |
| 7 | Stan smoke test — compiles and samples a minimal model |

If any step prints `✗`, the remediation instruction is printed
directly below it. Follow the instruction, then re-run
`rbf_install()`. If Step 2 fails, fix the toolchain before
continuing — subsequent steps cannot succeed without it.

Full remediation guidance is in this document (Steps 1–3 above)
and in the inline output of `rbf_install()`.

---

## Step 5: Verify

From the RBayesflow project root, run:

```r
source("R/source_all.R")
wf <- init_workflow(mode = "learn")
guide(wf)
```

Expected output from `guide(wf)`:

````
Current phase : Phase 1: Goal Declaration and Data Inspection
What to expect: You will see exploratory plots and off-ramp options.
Next step     : Open templates/phase1_exploration.qmd
````

If you see this, installation is complete. Proceed to
`docs/user-guide/02-posit-assistant-setup.md` to configure
Posit Assistant, or to `docs/user-guide/03-starting-an-analysis.md`
to begin your first analysis.
