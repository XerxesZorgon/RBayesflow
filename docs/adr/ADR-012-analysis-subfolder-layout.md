# ADR-012 — Analysis Subfolder Layout and `export_context()` Path

**Status:** Accepted
**Date:** 2026-09-10
**Owner:** John Peach
**Supersedes:** —
**Superseded by:** —

---

## Context

RBayesflow v0.1.0 supports only one analysis at a time per RBayesflow installation. `export_context()` writes `wf_context.json` to the project root, so a second analysis run in the same RBayesflow clone would overwrite the first analysis's context file. Users want to maintain multiple concurrent analyses (for example, one per dataset, or one per experiment) under a single shared RBayesflow install without cloning the whole project folder each time.

The v0.2.0 guidance layer introduces `rbf_new("name")`, which needs a well-defined place to put per-analysis artifacts. `wf_context.json` is the immediate driver, but the same folder is the natural home for the analysis's Quarto notebooks, saved `brmsfit` objects, and any generated report.

## Options considered

**Option A — one RBayesflow clone per analysis.** Each analysis is its own copy of the whole project. Rejected: redundant R source files that must be kept in sync manually; violates the "one install, many analyses" mental model that `rbf_install()` is designed to support; wastes disk on every clone.

**Option B — analysis subfolders under `data/`, context file written per-analysis.** `rbf_new("name")` creates `data/<name>/`, and `export_context()` writes `wf_context.json` to whichever folder the user's working directory is in. Selected.

**Option C — a database or registry of analyses.** A top-level `analyses/` registry file that maps names to state. Rejected: out of scope for v0.2.0; adds a new file format and coordination burden with no benefit at the single-user scale RBayesflow targets.

## Decision

**Option B.** `rbf_new("name")` creates `data/<name>/` and populates it with:

- `.Rprofile` — sources the shared workflow code from the project root via `source(file.path("..", "..", "R", "source_all.R"))`
- `wf_context.json` — initialised as `{}`
- `README.md` — the analysis name and creation date

The analysis subfolder is self-contained from the user's perspective: they open it as a working directory in RStudio or Positron, and everything runs from there. The R source files are not duplicated; the `.Rprofile` sources them from the project root.

`export_context()` gains a new default for its `path` argument: `file.path(rbf_analysis_path(), "wf_context.json")`. The helper `rbf_analysis_path()` returns `getwd()` when called from inside a `data/<name>/` folder (detected by an `.Rprofile` in the working directory and that directory being an immediate child of a `data/` folder under the RBayesflow project root, as located by `rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))`). Otherwise it returns the project root.

## Consequences

**Positive:**

- Multiple analyses coexist without conflict. A user can run three concurrent analyses under `data/pilot/`, `data/main_study/`, and `data/replication/`, each with its own `wf_context.json`.
- The `.Rprofile` mechanism keeps analysis subfolders self-contained without duplicating R source files. There is one canonical copy of every workflow function in `R/`.
- Backward compatibility is preserved: callers of `export_context(wf, path = "some/explicit/path.json")` are unaffected. The change is in the default only.
- Posit Assistant's context file is always adjacent to the analysis it describes, so the "read `wf_context.json` in the current folder" instruction in the phase templates works uniformly whether the user opens the analysis folder or the project root.

**Negative:**

- A user with a strongly configured global `.Rprofile` may see the generated one silently override or supplement their global settings. Mitigation: documented in `docs/user-guide/03-starting-an-analysis.md`, and the generated `.Rprofile` contains only the single `source()` call so it does not disturb unrelated global settings.
- The `rbf_analysis_path()` detection depends on `rprojroot` correctly identifying the project root; a `DESCRIPTION` file must exist at the project root (it does, per PHASE-0).

**Neutral:**

- `data/` is now a functional folder in the project layout rather than a scratch directory. Existing files under `data/` (if any) are unaffected; the change is additive.

---

*End of ADR-012.*
