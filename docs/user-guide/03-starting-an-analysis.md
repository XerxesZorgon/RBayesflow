---
{
  "id": "file_hka1y2ke",
  "filetype": "document",
  "filename": "03-starting-an-analysis",
  "created_at": "2026-09-24T19:08:24.728Z",
  "updated_at": "2026-09-24T19:08:24.728Z",
  "meta": {
    "location": "/",
    "tags": [],
    "categories": [],
    "description": "",
    "source": "markdown"
  }
}
---
# RBayesflow — Starting an Analysis

---

## Create an Analysis Subfolder

Each analysis lives in its own subfolder under `data/`. From the RBayesflow project root, run:

```r
source("R/source_all.R")
rbf_new("my_analysis")
```

Replace `"my_analysis"` with a short name for your analysis — no spaces, use underscores. `rbf_new()` creates `data/my_analysis/` containing:

- `.Rprofile` — sources the shared RBayesflow R code automatically when you open this folder as a working directory
- `wf_context.json` — the Posit Assistant context file (starts empty)
- `README.md` — the analysis name and creation date

---

## Open the Analysis Folder as Your Working Directory

In **RStudio**: Session → Set Working Directory → Choose Directory → select `data/my_analysis/`.

In **Positron**: File → Open Folder → select `data/my_analysis/`.

Once you set the working directory, the `.Rprofile` in the folder sources `R/source_all.R` automatically, loading all RBayesflow functions.

---

## Initialise the Workflow

Call `init_workflow()` with your chosen mode:

```r
wf <- init_workflow(mode = "learn")
```

| Mode | Use when |
|------|----------|
| `"learn"` | You are new to Bayesian analysis or want guided feedback at every step |
| `"practice"` | You use Bayesian methods occasionally and want standard output with reproducibility safeguards |
| `"expert"` | You already own a Bayesian workflow and want audit-bundle access with no guardrails |

You can also set the stage:

```r
wf <- init_workflow(mode = "practice", stage = "confirm")
```

The default stage is `"explore"`. Use `"confirm"` when you are running a pre-registered or publication-ready analysis that requires full prior justification and sensitivity checks.

---

## Find Out Where You Are: guide(wf)

At any point in a session, call:

```r
guide(wf)
```

This reads your current `wf_state` and prints three lines:

````
Current phase : <phase name>
What to expect: <one sentence on what you will see>
Next step     : <exact function call or template to open>
````

Call it whenever you are unsure what to do next.

---

## A Note on .Rprofile

The `.Rprofile` in each analysis subfolder contains a single line that sources `R/source_all.R` from the RBayesflow project root. This is what loads all workflow functions automatically.

If you have a global `~/.Rprofile` that sets options or loads packages, the subfolder `.Rprofile` will run after it on session start — the two do not conflict unless your global `.Rprofile` also sources RBayesflow files. If you see unexpected behaviour, check your global `.Rprofile` for conflicts and remove or reorder the relevant lines.
