# ADR-003: Quarto vs. RMarkdown for Report Templates

**Status:** Accepted  
**Date:** 2026-09-08  
**Deciders:** John Peach  
**Linked specs:** SDD §5 Phase 7 (Reporting), SDD §8 (Constraints)

---

## Context

RBayesflow's Phase 7 (Reporting) requires document templates that users fill in as they move through the workflow. Two document systems are viable in the R/RStudio ecosystem:

**Option A — RMarkdown (`.Rmd`)**  
- Mature, widely used in the R community
- Native RStudio support via knitr
- HTML, PDF, and Word output formats
- Development has slowed significantly since Quarto's introduction (Posit has shifted investment to Quarto)
- Uses R-specific chunk options (`{r, echo=TRUE, ...}`)

**Option B — Quarto (`.qmd`)**  
- Posit's current-generation scientific publishing system
- Language-agnostic (R, Python, Julia, Observable JS)
- Active development and Posit's primary investment
- Unified chunk option syntax (`#| echo: true`)
- Better cross-reference system (equations, figures, tables)
- Native support for Posit Connect, GitHub Pages, and other publishing targets
- RStudio ≥ 2022.07 ships with Quarto built in
- Supports `params:` for parameterized reports (useful for the workflow's mode/stage injection)

---

## Decision

**Option B: Quarto.**

All report templates are `.qmd` files. Users require RStudio ≥ 2024.04 with the bundled Quarto version (≥ 1.4).

---

## Rationale

1. **Active development and longevity.** Posit's team has explicitly redirected investment from RMarkdown to Quarto. New features (better cross-references, better figure layout, `revealjs` presentations) land in Quarto first; RMarkdown receives maintenance only.
2. **Parameterized reports.** Quarto's `params:` YAML header allows injecting `mode`, `stage`, and `wf_state` summary fields at render time. This maps cleanly to RBayesflow's mode/stage architecture.
3. **Cross-reference quality.** The Gelman et al. workflow produces numbered equations and model-comparison tables that benefit from Quarto's superior cross-reference system.
4. **Posit Assistant alignment.** Posit Assistant is aware of Quarto document structure in RStudio; it can read `.qmd` chunk contexts better than `.Rmd` in newer versions.
5. **Target audience compatibility.** Students using recent RStudio (≥ 2024.04) have Quarto available without a separate install step.

---

## Consequences

- **Positive:** Templates benefit from all current Quarto features; no maintenance cliff.
- **Positive:** `params:` injection allows mode-aware rendering without runtime R logic in the template.
- **Positive:** Quarto's cross-reference system handles numbered equations cleanly (relevant for model comparison tables and diagnostic summaries).
- **Negative:** Users on very old RStudio (< 2022.07) need a separate Quarto install. Documented in README prerequisites.
- **Negative:** `.qmd` syntax differs slightly from `.Rmd`; users migrating existing RMarkdown reports must convert. Not a concern for new users.
- **Not affected:** The workflow R scripts (`.R` files) are format-agnostic; this decision applies only to report templates.
