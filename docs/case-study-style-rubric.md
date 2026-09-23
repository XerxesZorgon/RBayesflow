# RBayesflow Case Study Style Rubric

*For articles accompanying the Bayesian Workflow case study series.*
*Derived from the Chapter 17 (Sleep Study) reference article.*

Apply the /humanizer skill after writing the document.

---

## 1. Document Structure

Every case study article follows this skeleton, in this order. Each section is separated by a horizontal rule (`---`).

```
# Chapter N Case Study: [Short Descriptive Title]

*[Attribution line]*
*[Implementation line]*

---

## The Question
## The Data
## What Is RBayesflow?  (first article only; link to it in subsequent articles)
## Setting Up: Phase 1
## [Model or Analysis Sections]
## How the [Concept] Evolved / Summary Table
## How RBayesflow Guided the Analysis
## Extensions for the Student
---
## Glossary
```

The H1 title uses the form `Chapter N Case Study: [Topic]`. All body sections use H2. Sub-sections within a model section (The Model, Choosing the Priors, etc.) use H3. Nothing deeper than H3 is used.

---

## 2. Opening Hook

The very first line after the section heading `## The Question` is a single, short rhetorical question in its own paragraph. It must be answerable in plain English without any notation. It anchors the entire article around a human decision, not a mathematical one.

Examples from Ch. 17:
- "How do you decide what to believe before you see the data?"

The hook is followed by two or three paragraphs that:
1. Introduce the core concept being illustrated.
2. Name the dataset and the case study goal.
3. Preview the arc: "We build N models... each one exposes..."

---

## 3. Attribution and Italics Block

The two lines immediately below the H1 title are italicised attribution lines:

```markdown
*Based on Gelman, Vehtari et al., Bayesian Workflow (2026), Chapter N.*
*Implemented using RBayesflow and brms in R.*
```

No citation keys. No footnotes. Link to the book if a stable URL is available.

---

## 4. Figures

Insert figures using the bare Markdown image syntax with a `figs/` relative path and an SVG filename:

```markdown
![Figure N.M](./figs/Fig-N.M.svg)
```

The figure reference always appears immediately after the paragraph that introduces it. The paragraph before the figure ends with a forward reference: "This is Figure N.M." or "...as shown in Figure N.M." The figure image tag is on its own line, preceded and followed by a blank line. Figure numbers should match numbers given in the folder containing the figures.

After the image tag, write one to three sentences of interpretation in plain prose (no sub-heading). This interpretation is mandatory — do not leave a figure without commentary.

**The interpretation must describe what the figure actually shows, not what it should show.** Read the rendered figure before writing the interpretation. When the visual result differs from the expected result (for example, a posterior that is shifted relative to a true-value marker, or a K-M curve that reaches zero faster than expected), the interpretation must name the discrepancy and explain it. Copying an expected result from the book or from code comments without checking the figure is a specific error to avoid.

**Each interpretation must identify at least one concrete visual feature** the reader can locate in the figure: a gap between a posterior slab and a dashed reference line, the rate of descent of a K-M curve, the width of a credible interval, the overlap between two distributions. Confirmation statements of the form "the posteriors sit on the true values" are not sufficient on their own — name the specific visual evidence that supports that claim, or report the discrepancy if the evidence does not support it.

Number figures sequentially as `N.M` where `N` is the chapter number and `M` increments from 1. Never skip a number.

Do not use `CenteredFigure` MDX components. Do not include captions inside the image tag.

---

## 5. Code Blocks

Every code block uses the `r` language tag:

````markdown
```r
code here
```
````

Code blocks serve two purposes:

**Specification blocks** show the exact R calls used — model formulas, prior definitions, RBayesflow phase calls. These are complete and runnable.

**Console output blocks** show what RBayesflow printed when the code was run. Use plain (untagged) fenced blocks for console output to visually distinguish it from runnable code:

````markdown
```
!!! DIAGNOSTIC FAILURE !!! Failed criteria: ...
```
````

Every code block is immediately preceded by a prose sentence that explains what it does and why. Code never appears without context.

---

## 6. Mathematics

All non-trivial formulas use KaTeX: inline math in `$...$`, display math in `$$...$$`. Display equations are used for model definitions and named distributions. Inline math is used for parameter names within prose. Display math should have the double dollar signs on separate lines from the math formulas. If more than one line of math is used, use \begin{aligned}, \end{aligned} and insert & before equal signs, or whatever symbol is appropriate.

Immediately after every display equation, define every symbol that appears in it. Use a bulleted list if there are three or more symbols; use inline prose if there are one or two. Never leave a symbol undefined.

Exponential distributions are parameterised by rate $\lambda$ with mean $1/\lambda$ — state this explicitly the first time the distribution appears, because the rate / scale ambiguity trips readers.

---

## 7. Prose Style

**Sentence length.** Vary sentence length deliberately. Short declarative sentences anchor key claims. Longer sentences develop reasoning. Never use more than three long sentences consecutively.

**Active voice.** Prefer active constructions. "The model assumes..." not "It is assumed by the model..." Write in first-person plural when describing choices made in the analysis: "We drop the first two days..." "We centre the time axis..."

**No em-dashes, or double dashes (--)** Replace every em-dash with a comma, a semicolon, a colon, or a restructured sentence. This is a hard constraint.

Avoid colons mid-sentence.

**No jargon without definition.** Every technical term is defined on its first use. The definition is in the same sentence or the next sentence; it is not deferred to the glossary. The glossary entry then deepens the definition or adds a link.

**Annotation of design choices.** Every modelling choice — every prior, every structural decision — is accompanied by a "because" sentence. "We centre the prior at zero because we want the data to determine the direction, not the prior." The reader should never encounter a number (250, 0.02, LKJ(1)) without an explanation of where it came from.

**Annotation of simulation recovery.** Every figure that compares a posterior against a simulated true value requires a "because" sentence for any visible offset. "The posterior peaks slightly left of the true value because the Beta(1,10) prior exerts downward pressure; the bias is small relative to the posterior width." The reader should not be left to wonder whether a shifted posterior is a model failure or expected behaviour.

**Tables for comparisons.** Use Markdown tables when comparing more than two versions of the same thing (prior evolution across models, LOO comparison, diagnostic results). Tables use the pipe format with a header row and alignment dashes. Every table has an immediately preceding sentence that introduces it and explains what to look for.

You are a first-class author tasked with writing the text to sound indistinguishable from a skilled human writer. Apply the following strict adjustments:

1. Burstiness (Sentence Variety): Human writing naturally shifts rhythms. Mix some shorter sentences with longer, more complex ones. Break up uniform, mid-length paragraphs.
2. Vocabulary Scrubber: Completely remove common AI catchphrases. Do not use words like "delve," "testament," "transformative," "beacon," "seamlessly," "moreover," "foster," or "tapestry."
3. Natural Transitions: Eliminate rigid essay-style transitions (e.g., "In conclusion," "Furthermore," "It is important to remember"). Use conversational flow instead.
4. Active Voice & Contractions: Convert passive phrasing to active. Use natural contractions (e.g., "don't" instead of "do not", "it's" instead of "it is") to relax the tone.
5. Casual Professionalism: Write like a subject-matter expert explaining a concept to a peer over coffee — informed, clear, but completely unpretentious.

---

## 8. Section: "How the [X] Evolved"

This section appears near the end of the main body, before "How RBayesflow Guided the Analysis". It presents a summary table showing how the key concept (priors, models, diagnostics) accumulated across the analysis. Columns for the prior-evolution table are: Model, New parameter, Prior, Reasoning. The Reasoning column uses plain English phrases, not notation alone.

---

## 9. Section: "How RBayesflow Guided the Analysis"

This section is a numbered list, one item per RBayesflow phase that was exercised. Each item uses the pattern:

```
N. **Phase name (Phase N):** One sentence saying what the phase did.
   One sentence saying what the output was.
```

This section is short — five to eight list items. It is a checkpoint that connects the modelling narrative back to the workflow structure.

When a chapter uses `cmdstanr` directly rather than brms, the Phase 3 item must say so and name the reason: "Because the chapter's likelihoods cannot be expressed through brms formula interfaces, Phase 3 used `cmdstanr` directly via a thin `cstan()` helper." When LOO-CV was not run, the Phase 6 item must state why and name the qualitative comparison method used instead.

---

## 10. Section: "Extensions for the Student"

A bulleted list of four to six items. Each item:
- Names the topic in bold.
- Names the figure or package relevant to that extension.
- States one concrete action: "Try it on Model 1 to see..." "Replace `gaussian()` with `student()`..."
- Does not include code.

The items are ordered from lower to higher difficulty.

---

## 11. Glossary

The glossary appears as the final section with the H2 heading `## Glossary`. It is separated from the preceding section by `---`.

**Entry format.** Each entry starts with the term in bold. If an external authoritative link exists (Wikipedia, Stan docs, package homepage), hyperlink the term itself. If the term is an acronym, expand it in the first sentence.

```markdown
**[Rhat](https://mc-stan.org/rstan/reference/Rhat.html) (R-hat):** A
diagnostic number...
```

**Content.** Each entry gives:
1. A plain-English definition (one to two sentences).
2. The typical value or threshold relevant to this case study, when applicable (e.g., "values above 1.01 suggest a problem").
3. A connection to how the term appeared in the article, when this adds clarity.

**Coverage.** Include every bolded term from the body, every acronym expanded on first use, every distribution name, every RBayesflow-specific term, and every concept that a reader without a statistics background might not know. Err on the side of more entries rather than fewer.

**Order.** Alphabetical by the first word of the term, ignoring "the" and "a". Backticks and special characters are ignored for sort order (so `` `wf_state` `` sorts under W). Acronyms sort by their spelled-out form — E-BFMI sorts under E, not under the expansion.

**Do not include.** Plain R function names (e.g., `brm()`, `loo_compare()`), dataset names, or author names.

---

## 12. Links

Use hyperlinks generously for:
- Package names (`brms`, `lme4`, `bayesplot`, `loo`, `priorsense`) linked to their homepage or CRAN page.
- Named statistical concepts linked to Wikipedia or authoritative explainers.
- Diagnostic terms linked to the Stan documentation.
- Datasets linked to their source paper or repository.

Do not use bare URLs. Always use `[anchor text](url)` format. Do not use footnote-style links.

No citation keys (e.g., `[@gelmanBayesianWorkflow2020]`). No reference section at the end. All attribution is done through hyperlinks or the attribution italics block at the top.

---

## 13. What to Omit

The following items from the Wild Peaches MDX publication format do NOT appear in the case study articles:

- YAML front matter (no `---` block at the top of the file)
- `CenteredFigure` MDX components
- `[@citekey]` citation keys
- References / Bibliography section
- NotebookLM link
- `hero:` image specification
- `software:` tag list
- `tags:` or `keywords:` fields
- Line breaks in paragraphs. All standard text paragraphs and list items should be formatted as single continuous lines, while preserving the structural formatting for headers, code blocks, math equations, and tables.

---

## 14. Quick Checklist Before Submission

- [ ] H1 follows `Chapter N Case Study: [Title]` form
- [ ] Attribution lines are present and italicised
- [ ] Opening hook is a single rhetorical question in its own paragraph
- [ ] Every figure uses `![Figure N.M](figs/Fig-N.M.svg)` syntax
- [ ] Every figure interpretation describes what the figure **actually shows**, verified by reading the rendered output
- [ ] Every figure interpretation names at least one concrete visual feature the reader can locate
- [ ] Any visible discrepancy between a posterior and a true-value marker is named and explained, not silently passed over
- [ ] Every code block has a preceding context sentence
- [ ] Every symbol in every display equation is defined
- [ ] No em-dashes anywhere in the document
- [ ] No citation keys anywhere in the document
- [ ] No YAML front matter
- [ ] Every bolded term appears in the glossary
- [ ] Glossary entries are in alphabetical order (backticks and acronyms handled per §11)
- [ ] Every modelling choice has a "because" sentence
- [ ] Every simulation recovery figure has a "because" sentence for any visible offset from the true value
- [ ] "Extensions for the Student" has four to six items, ordered easy to hard
- [ ] "How RBayesflow Guided the Analysis" is a numbered list, phases only
- [ ] Phase 3 item names `cmdstanr` directly when brms was not used, with a reason
- [ ] Phase 6 item states why LOO was skipped and names the qualitative comparison method, when LOO was not run
