# ADR-001: Workflow-as-Project-Folder vs. R Package

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** SDD §8 (Constraints)

---

## Context

RBayesflow must be deliverable to students and working scientists. Two distribution strategies are available:

**Option A — R package (submitted to CRAN or installed from GitHub)**  
- Standard R distribution mechanism
- Requires package structure: NAMESPACE, DESCRIPTION, `R/`, `man/`, tests via `testthat`
- CRAN submission triggers ongoing maintenance burden (policy compliance, reverse-dependency checks, R version compatibility)
- Requires exporting functions as a public API, making breaking changes more expensive
- Grants `install.packages()` or `remotes::install_github()` convenience

**Option B — Project folder (R scripts + Quarto templates + configuration)**  
- Clone or download the folder; `source()` the initialization script
- No package infrastructure; no NAMESPACE; no CRAN
- Designed to be read and modified, not treated as a black box
- Perfectly suited to a workflow whose purpose is to be understood, not just executed
- The Antigravity/Spec Kit model: the deliverable is the project, not the package

---

## Decision

**Option B: project folder.**

RBayesflow is delivered as a structured directory (`RBayesflow/`) containing R scripts, Quarto templates, and a workflow initialization script. Users clone or copy the folder and source the setup file.

---

## Rationale

1. **No new statistical code means no package need.** A package is the right form when you are distributing reusable functions with a stable API. RBayesflow provides sequencing, display logic, and templates — all of which benefit from being read and modified, not abstracted.
2. **Transparency is a design goal.** The workflow explicitly aims to keep statistical evidence visible. Wrapping it in a package obscures the plumbing that students are supposed to understand.
3. **CRAN maintenance is out of scope.** The brief explicitly excludes CRAN submission. Package overhead (R CMD check, dependency pinning, NEWS.md conventions) would consume development capacity for no user benefit.
4. **Installation is trivial for the target audience.** Students and working scientists who can install Stan can clone a Git repo or unzip a folder.

---

## Consequences

- **Positive:** Zero CRAN maintenance burden. Users can inspect and modify every function. Versioning is just Git tags.
- **Positive:** Quarto templates ship as first-class files, not bundled `system.file()` calls.
- **Negative:** No `install.packages()` convenience. Mitigated by a one-line setup script in the README.
- **Negative:** No automatic dependency resolution. Mitigated by a documented `renv.lock` or `DESCRIPTION`-style dependency list that users install manually.
- **Watch:** If RBayesflow gains significant adoption and a stable API emerges organically, packaging may become worthwhile. Track in a future ADR.
