---
{
  "id": "file_iy0x07fg",
  "filetype": "document",
  "filename": "02-posit-assistant-setup",
  "created_at": "2026-09-24T19:08:26.716Z",
  "updated_at": "2026-09-24T19:08:26.716Z",
  "meta": {
    "location": "/",
    "tags": [],
    "categories": [],
    "description": "",
    "source": "markdown"
  }
}
---
# RBayesflow — Posit Assistant Setup

Posit Assistant is the AI layer in RBayesflow. It reads `wf_context.json` — a JSON snapshot of your current workflow state written by `export_context(wf)` — and uses it to give context-aware guidance: interpreting diagnostics, explaining Rhat and ESS values, pointing to relevant sections of Gelman et al. (2020), and suggesting next steps in plain language.

Posit Assistant **never runs code for you**. All suggestions appear as non-executable previews; you decide whether to run them.

RBayesflow supports Posit Assistant in both RStudio and Positron. Follow the track for your IDE below.

---

## Getting an OpenRouter API Key

OpenRouter provides pay-as-you-go access to many language models through a single API, so you can pick the model that fits your budget and reasoning needs.

1. Go to <https://openrouter.ai> and create an account.
2. In the dashboard, click **Keys** → **Create Key**.
3. Give the key a name (e.g. `rbayesflow`), click **Create**, and copy the key. Store it somewhere safe — you will not see it again.
4. Add credit to your account (Settings → Credits). A few dollars is enough for many sessions with Gemini 3.8 Flash; more is prudent for regular Claude Opus 5 use.

**Recommended models:**

| Model | ID | Input / Output ($/M tokens) | Context |
|---|---|---|---|
| Claude Opus 5 (primary) | `anthropic/claude-opus-5` | $5.00 / $25.00 | 1,000,000 |
| Gemini 3.8 Flash (low-cost) | `google/gemini-3.8-flash` | $0.75 / $3.75 | 1,048,576 |

Claude Opus 5 is Anthropic's flagship reasoning model and the recommended default for interpreting diagnostics and generating Socratic prompts in learn mode. Gemini 3.8 Flash costs about one-seventh as much per token with a comparable context window — a reasonable substitute when cost matters more than the last increment of reasoning quality. Prices are pay-per-use through OpenRouter, listed as of September 2026; check the OpenRouter dashboard for current rates. If a listed model is retired, pick any OpenRouter model whose context window comfortably exceeds the size of `wf_context.json` (typically well under 10 000 tokens).

---

## Option A — RStudio (≥ 2026.04, Posit Assistant 0.7.7+)

1. Open RStudio. If the Posit Assistant panel is not visible, go to **Help → Copilot** (or the Assistant toolbar button) and follow the prompt to enable it.
2. Click the **gear icon** in the Posit Assistant panel.
3. Select **Providers** → **Add Provider**.
4. Choose **OpenRouter** from the provider list.
5. Paste your API key into the **API Key** field.
6. Select your model (e.g. `anthropic/claude-opus-5` or `google/gemini-3.8-flash`).
7. Click **Done**.

Posit Assistant is now configured. It will be available in every RStudio session.

---

## Option B — Positron (current release)

1. Open the Command Palette (`Ctrl+Shift+P` / `Cmd+Shift+P`).
2. Type **Positron Assistant: Configure Language Model Providers** and select it.
3. Click **Add Model Provider** → **Custom Provider**.
4. Set **Base URL** to: `https://openrouter.ai/api/v1`
5. Paste your API key into the **API Key** field.
6. Enter the model name (e.g. `anthropic/claude-opus-5` or `google/gemini-3.8-flash`).
7. Click **Done**.

---

## Using Posit Assistant with RBayesflow

At the start of each session, refresh the context file so Posit Assistant has your current workflow state:

```r
export_context(wf)
```

This writes `wf_context.json` to your analysis subfolder (`data/<name>/`). Posit Assistant reads it automatically when you are working in that folder.

To ask for guidance, open the Posit Assistant chat panel and paste a prompt such as:

````
Read wf_context.json in the current folder. I am in learn mode, explore stage. My diagnostics just ran. What should I check next?
````

---

## What Posit Assistant Will and Won't Do

**Will:**

- Interpret Rhat, ESS, divergence counts, and BFMI values in plain language
- Explain what a diagnostic failure means for your specific model family
- Suggest next steps based on your current `wf_state`
- Point to the relevant section of Gelman et al. (2020)
- Generate code suggestions as non-executable previews for you to review

**Won't:**

- Run code automatically
- Display coefficient summaries if `wf$diagnostics$passed` is not `TRUE`
- Make modelling decisions for you
- Access the internet or any data source other than `wf_context.json`
