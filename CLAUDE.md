# CLAUDE.md — masgenomics-docs

> **This file is local-only. It is gitignored and must never be staged, committed, or pushed.**
> If you see it appear in `git status`, the .gitignore entry has regressed — fix the .gitignore, do not commit the file.

## ⚠️ CRITICAL workflow rule — build locally, get APPROVE, push at the END

**Every edit to this site follows this loop. No exceptions, no per-stage commits.**

1. Implement file edits for the current stage.
2. **Build the full local preview** with `quarto render` so the rendered site sits at `/Users/aguswibowo/Library/CloudStorage/OneDrive-JamesCookUniversity/JCU/Research/Genomics/Analysis/genomic_prediction/masgenomics-docs/_site/`. Verify 0 warnings, 0 errors across all changed pages.
3. **Verify sibling-package cleanliness** (`git status --porcelain` in `../masbayes/`, `../masreml/`, `../maspipeline/` must each be empty).
4. **Surface the changed pages to the user** — list the modified / new / deleted files with their `_site/<path>.html` targets. Do NOT run any `git` write commands.
5. **Wait for the user's explicit `APPROVE`** in PLAN.md (or chat). If they request changes, fold them in and rebuild the preview before re-surfacing — do NOT commit intermediate states.
6. Flip the stage's `Status` to FINISHED in `PLAN.md`. Move to the next stage.
7. The **single consolidated commit + push** happens only at the very end of the entire plan, after the user gives the final APPROVE — that is `PLAN.md → Stage Z`. Until then nothing is pushed.

If you ever feel the urge to `git add` between stages, **stop**. The site is built at `_site/`; the user reviews; only Stage Z at the end pushes.

## ⚠️ Always purge deployment history after every push

**Every successful push triggers a new GitHub Pages deployment record. Always clean up the older records, leaving only the latest one (or none if the next workflow is still in_progress).** The user expects this as a routine post-push step.

After `git push` completes and the new Publish workflow has at least been queued:

```bash
# 1. List current deployment records
gh api "repos/bowo1698/masgenomics-docs/deployments?per_page=100" --paginate \
  --jq '.[] | "\(.id) sha=\(.sha[0:7]) created=\(.created_at)"' \
  > /tmp/deps.txt
cat /tmp/deps.txt

# 2. Keep the most recent ID (line 1). If the user explicitly asks to delete
#    "termasuk yang terakhir" / "including the last one", pass tail -n +1 instead
#    of tail -n +2. The in-progress workflow does NOT yet have a deployment
#    record, so deleting "the last one" is safe — its eventual record will be
#    created fresh.
tail -n +2 /tmp/deps.txt | awk '{print $1}' | while read id; do
  gh api -X POST "repos/bowo1698/masgenomics-docs/deployments/$id/statuses" -f state=inactive --silent
  gh api -X DELETE "repos/bowo1698/masgenomics-docs/deployments/$id" --silent
  echo "deleted: $id"
done
rm -f /tmp/deps.txt
```

**Why this matters:** GitHub Pages keeps the deployment history visible on the repo home page (`Environments → github-pages`) and on each commit's "Deployments" pane. The live site is served from the `gh-pages` branch (independent of these records), so deleting the records is purely cosmetic — the deployed HTML stays online. Leaving stale records around clutters the UI and can confuse readers about which build is current.

**Hard invariants:**

- **Never touch `refs/heads/gh-pages`.** The records live in the GitHub API; the branch is what serves traffic.
- **Never delete a deployment for the in-progress workflow.** If `status=in_progress` for the most recent Publish run, its deployment record has not yet been created — the running build itself is safe. Just keep an eye on the cleanup loop so it does not race a freshly-created record.
- **Cleanup is post-push, not pre-push.** Do not delete records before the push; the new push needs to land first so the new record (when created) becomes the survivor.

Aliases the user uses interchangeably for this task: "hapus riwayat deployment", "remove deployment history", "purge deployments", "kosongkan environments". Treat all of them as the same routine and run the loop above.

## CI / `_freeze` cache pitfalls (learned the hard way)

Three issues bit hard during the May 15 marker-qtl-congruency rollout. Avoid them on future tutorial pages.

1. **Local Quarto version must match the version pinned in `.github/workflows/publish.yml`.** If local Quarto is newer, the freeze cache hash stored in `_freeze/<page>/execute-results/html.json` is in a format the CI Quarto cannot validate, so CI falls through to re-executing the chunk — which fails because `masbayes` / `masreml` are not installed in CI. Pin both sides. The workflow currently uses Quarto `1.9.37`; bump local to match before regenerating freeze caches.
2. **`actions/cache@v4` for `_freeze` is harmful and has been removed.** The cache key `Linux-freeze-${{ hashFiles('**/*.qmd') }}` shifts every time any qmd changes, so the `restore-keys` fallback brings back a *stale* `_freeze/` directory from an earlier run. That stale directory overwrites the committed cache for the page that just changed, and CI re-executes its chunks. The freeze cache is shipped via git (see § "What this repo is") — that is the single source of truth. **Do not re-add an `actions/cache` step for `_freeze`.**
3. **`DT::datatable` keeps drifting on CI.** `DT` ships an htmlwidget binding under `site_libs/datatables-binding-<version>/`, plus transitive bits like `crosstalk-1.2.1/js/crosstalk.min.js`. The path versions baked into the freeze cache must match the versions CI emits, or the widget div renders but the binding JS hits 404 and the table never initialises. The site convention — and the only reliably portable choice for now — is **`knitr::kable`** for tabular results in `tasks/`. Every other `tasks/*.qmd` already uses it. Reserve DT for `input-data/*.qmd` where the page is small enough to absorb the drift, and only after re-rendering against the CI-installed `DT` version.

## What this repo is

A **Quarto-based GitHub Pages documentation site** that unifies the docs for the masgenomics suite:

- [`masreml`](https://github.com/bowo1698/masreml) — REML mixed models / BLUP / GBLUP / GWABLUP / GWAS / cross-validation / threshold models (R + Rust via extendr).
- [`masbayes`](https://github.com/bowo1698/masbayes) — Bayesian genomic prediction (BayesA, BayesR; R + Rust via extendr).
- [`maspipeline`](https://github.com/bowo1698/maspipeline) — Rust CLI preprocessing pipeline (phasing → haploblock discovery → microhaplotype genotyping; produces the inputs that masreml / masbayes consume).

All three packages share Rust backends and target SNP and microhaplotype genotypes. This site is the single public reference for the three.

- **Live URL:** https://bowo1698.github.io/masgenomics-docs/
- **GitHub:** https://github.com/bowo1698/masgenomics-docs
- **Default branch:** `main`. Rendered HTML lives on `gh-pages` (force-orphan, deployed by CI).
- **Latest tagged release:** `v0.1.0` (2026-05-11)
- **Local path:** sibling of `masbayes/`, `masreml/`, `maspipeline/` under the workspace root.

## Navbar structure (current, post-revisions)

```
Home · Installation · Theory · Task · Input Data · Reference ▼ · Internals
                                                      ├─ maspipeline
                                                      ├─ masreml
                                                      └─ masbayes
```

- Demo Data tab **removed** (content folded into Input Data).
- Choosing tab **removed** entirely.
- Installation tab **added** (single page: Rust + cargo setup, then `install.packages(github tar.gz)` for both R packages).

## Theory section structure (current, multi-page redesign 2026-05-13)

Four top-level theory pages, two of which have been redesigned into **multi-page sidebar structures** with sub-pages. Pages use a mix of overview + method-specific sub-pages with code and references consolidated at the bottom of each sub-page.

### 1. `theory/biallelic-and-multiallelic.qmd`
Single page. Data-layer: biallelic vs multi-allelic biology, raw data forms for SNP + MH (Da 2015 $\mathbf{W}_{\alpha h}$ coding + worked example), matrix for models, important functions.

### 2. `theory/mixed-model-blup-family.qmd` — overview + 4 sub-pages
Sidebar entry "Mixed Model & BLUP Family" expands to:
- **`blup-reml.qmd`** — REML algorithms: why ML is biased, restricted log-likelihood, EM / NR / AI-REML derivation + pseudocode, algorithm comparison table
- **`blup-pblup.qmd`** — PBLUP: K = A, pedigree kinship, IBD intuition
- **`blup-gblup.qmd`** — GBLUP: SNP-GBLUP (VanRaden 2008, 3 formulations, allele freq callout, rare-allele weighting) + MH-GBLUP (Da 2015, A_h = T_{αh}T_{αh}', LD erosion robustness)
- **`blup-gwablup.qmd`** — GWABLUP: 5-step Meuwissen 2024 pipeline, EMMAX rationale, inflation bias callout, graceful degradation note

Overview page (`mixed-model-blup-family.qmd`) contains: TL;DR (4-method table), shared assumptions, Henderson MME + biological intuition, navigation table to sub-pages, REML algorithm table, When to Use table.

### 3. `theory/bayesian-alphabet.qmd` — overview + 3 sub-pages
Sidebar entry "Bayesian Alphabet" expands to:
- **`bayes-mcmc-em.qmd`** — MCMC & EM Algorithms: Gibbs sampler (stationary distribution, Fundamental Theorem, mixing time, spectral gap), EM (KL-decomposition / ELBO, Jensen monotonicity, MAP extension), full math + pseudocode, MCMC vs EM comparison table
- **`bayes-bayesa.qmd`** — BayesA: scaled inv-χ² prior, Student-t marginal, S hyperparameter derivation, SNP vs MH W matrix differences + biological implications
- **`bayes-bayesr.qmd`** — BayesR: 4-component mixture prior, Dirichlet on π, PIP diagnostics + thresholds, SNP vs MH W matrix + biological implications (PIP locality, π̂₀ comparison, stability under selection)

Overview page (`bayesian-alphabet.qmd`) contains: TL;DR (BayesA vs BayesR vs SNP/MH vs MCMC/EM), shared linear model, shared assumptions, navigation table, When to Use table.

### 4. `theory/continuous-and-binary-trait.qmd`
Single page. Trait-type page: (1) Traits in quantitative genetics, (2) How traits are calculated (continuous vs binary Laplace in masreml vs Albert-Chib in masbayes), then Variance & h², When to use, See it in code, References.

**References policy:** primary publications only. Wibowo Chapter 1 used as source for biological prose in MH-GBLUP and BayesA/BayesR sub-pages — **never listed in any References block** because thesis is not yet published. Internal workspace PDFs are private working sources, not in the public References block.

## API dependencies (READ THIS BEFORE EDITING ANY TUTORIAL OR INPUT-DATA PAGE)

This site is downstream of three sibling repos. Their **public APIs are the source of truth** for everything documented here — schemas, tutorials, reference pages, demo data:

| Sibling repo | What this site uses from it |
|---|---|
| `../masreml/` | exported R functions (REML, BLUP, GBLUP, GWAS, cross-validation, threshold models), `inst/extdata/demo_data{,_small}.rds`, Rust source under `src/rust/src/` |
| `../masbayes/` | exported R functions (BayesA, BayesR, `construct_wah_matrix`, etc.), `inst/extdata/demo_data{,_small}.rds`, Rust source under `src/rust/src/` |
| `../maspipeline/` | CLI tools (`convert-to-vcf`, `convert-from-vcf`, `haplotype-hybrid`) plus its `README.md` (split across `reference/maspipeline/{index,data-preparation,phasing,genotype-to-haplotype,microhaplotype-discovery}.qmd` — overview keeps General info + Installation + footer meta verbatim; 4 stage pages carry the relevant README section plus a Goal / Pipeline-position intro callout per stage) |

Rules whenever you touch sibling-package code:

- **Sibling working trees must stay pristine.** Never edit files in `../masreml/`, `../masbayes/`, or `../maspipeline/` from this repo. If you need to *consume* anything from them at build time (e.g. regenerate reference docs, sync the maspipeline README into the reference page), do it via a copy/temp pattern. After any such operation, verify `git status --porcelain` in all three sibling repos shows zero changes.
- **Demo data lives upstream.** `demo-data/main/demo_data.rds` and `demo-data/toy/demo_data_small.rds` are byte-identical to the files shipped under each sibling's `inst/extdata/`. They are produced by `../tools/make_demo_data.R` (which depends on `library(masbayes)`). masgenomics-docs only consumes them; canonical regeneration happens upstream.
- **maspipeline reference is split into overview + 4 stage pages.** `reference/maspipeline/index.qmd` keeps the README's General info, Installation, and footer meta sections verbatim (LICENSE link patched to the GitHub URL). The four stage pages (`data-preparation.qmd`, `phasing.qmd`, `genotype-to-haplotype.qmd`, `microhaplotype-discovery.qmd`) each carry the corresponding README section verbatim plus a small Goal + Pipeline-position intro callout. If the upstream README changes substantially, refresh the affected stage page(s) plus index.qmd manually rather than divining the intent.
- **Track API breakage**. Any change to a `masreml::` or `masbayes::` exported symbol that appears in this site can silently invalidate a tutorial. The planned `_scripts/check-api.R` (PLAN.md Stage A) detects that; until it lands, rely on manual `quarto render` after sibling-package updates.

### Authorised sibling-repo edits (exceptions log)

The "sibling working trees stay pristine" rule has been overridden once, with explicit user authorisation:

- **2026-05-12** — `masreml/R/load_data.R` + `man/load_data.Rd` (commit `2ba45c3`) and `masbayes/R/load_data.R` + `man/load_data.Rd` (commit `e0696fe`): synced `load_data()` roxygen with the actual `demo_data{,_small}.rds` dimensions (pre-fix the docstrings claimed 200×400 SNP / 50 QTL; actual is 200×100 / 10 QTL). Both commits authored by `bowo1698 <aguswibowo1698@gmail.com>`, no Claude trailer; only `R/load_data.R` and `man/load_data.Rd` modified per repo.

If a future task requires touching a sibling repo from a masgenomics-docs session, append a new dated entry here.

## Canonical references (in the parent workspace, NOT in this repo)

- `../docs/plans/2026-05-10-masgenomics-docs-design.md` — design doc + Implementation Status snapshot. Read first when picking work up cold.
- `../docs/plans/2026-05-13-blup-family-theory-redesign.md` — brainstorm design for BLUP family multi-page redesign (Stages B1–B2).
- `../docs/plans/2026-05-13-bayesian-alphabet-redesign.md` — brainstorm design for Bayesian Alphabet multi-page redesign (Stage B3).
- `./PLAN.md` (this repo, local-only) — **current** plan.

## Hard rules (do not break)

### Repo hygiene

1. **Never stage, commit, or push `CLAUDE.md` or `PLAN.md`** (this file and the plan file). Both are in `.gitignore`. `REVISION_PLAN.md` stays in `.gitignore` too even though the file is currently absent — guards against future revision plans being accidentally committed.
2. **Never commit anything under `../docs/`** (workspace planning artifacts) into this repo. `docs/` is also in this repo's `.gitignore`.
3. **Never touch `../masbayes/`, `../masreml/`, or `../maspipeline/` working trees** from this repo. For altdoc-style regeneration of reference pages, use the `rsync`-to-temp pattern in `_scripts/render-reference.R`. Verify `git status --porcelain` in all three sibling repos before AND after any cross-repo read.
4. **Targeted `git add <file>` only**, never `git add .` or `git add -A`. Verify the staged set with `git diff --cached --name-only` before commit. (Reminder: there are no per-stage commits in this plan — only Stage Z at the very end. See CRITICAL workflow rule above.)

### Commit authorship — no Claude attribution, ever

- **Author of every commit on `main` must be `bowo1698 <aguswibowo1698@gmail.com>` only.**
- **Never add `Co-Authored-By: Claude …` trailers** to commit messages. Never use `--author=Claude…`. Do not include "Generated with Claude Code" or any similar credit lines.
- This is enforced retroactively: 5 earlier commits had Claude trailers and were rewritten via `git filter-branch --msg-filter` and force-pushed. A local backup branch `backup-with-claude-trailers` preserves the pre-rewrite SHAs (not pushed).
- If you ever accidentally include a trailer, rewrite history before pushing.
- PR/Issue descriptions, comments, release notes follow the same rule — no Claude credit.

### Subagent dispatch

If you ever dispatch subagents for work on this repo, **the subagent prompt must explicitly forbid Claude trailers, forbid touching CLAUDE.md / PLAN.md / docs/ / sibling-package repos, and forbid per-stage `git` writes**. Without that, the agent will default to bad behaviour and you will have to redo or rewrite history.

## Project-specific gotchas

1. **GitHub handle is `bowo1698`**, not `aguswibowo` (the macOS username). All URLs, `install.packages(...)`-style strings, CITATION entries use `bowo1698`. Local filesystem paths use `aguswibowo`.

2. **Lua shortcode location.** `{{< rust-src <pkg> <path> >}}` only resolves when the Lua file is wired up as a **Quarto extension**:
   - `_extensions/rust-src/_extension.yml` declares `contributes: shortcodes:`
   - `_extensions/rust-src/rust-src.lua` implements it
   - `_quarto.yml` has **no** `filters:` line for it.

3. **altdoc output path.** altdoc 0.7.2 writes `.qmd` files to `<pkg>/_quarto/man/`, not `<pkg>/docs/reference/`. The script also must `setwd(pkg_path)`. Generated files have no YAML frontmatter — the script prepends one before copying into `reference/<pkg>/`.

4. **Theme is single, not dual.** `_quarto.yml` uses `[superhero, assets/custom.scss]` only. No light fallback, no navbar toggle. Inline `<code>` is amber `#ffc107` on rgba(255,255,255,0.08) via custom SCSS.

5. **Interactive table widget is DT, not reactable.** Reactable's runtime JS clashes with the dark theme.

6. **Demo data layout (post-revision)** matches the upstream sibling `inst/extdata/` byte-for-byte. The canonical artifact is a single `.rds` per size:
   ```
   demo-data/
   ├── generate.R                  # adapted from ../tools/make_demo_data.R
   ├── main/                       # large config (n=200, n_snp=100, n_blocks=50, n_qtl=10, h2=0.5)
   │   ├── demo_data.rds                          # CANONICAL — byte-identical to ../masbayes/inst/extdata/demo_data.rds
   │   ├── pheno_continuous_G1.csv                # derived projections
   │   ├── pheno_binary_G1.csv
   │   ├── geno_G1.csv
   │   ├── snp_map.txt
   │   ├── microhaplotype_coordinates.csv
   │   ├── mh_genotypes/hap_geno_<chr>            # 4-col-per-block MicrohapsSel format, 5 files
   │   └── qtl/{qtl_snp_ids,qtl_mh_blocks,effects_snp,effects_mh}.rds
   ├── toy/                        # small config (n=100, n_snp=50, n_blocks=25, n_qtl=5, h2=0.5)
   ├── main.zip                    # bundles every file under main/
   └── toy.zip
   ```
   `set.seed(42)` (genetic stream) + `set.seed(41)` (split) at the top of `generate.R` guarantees bit-identical reruns. CI does **not** regenerate; everything is a static committed asset.

7. **CI does NOT install masreml / masbayes / altdoc.** Reference `.qmd` are pre-committed. **Task tutorial chunks evaluate locally and ship outputs via the `_freeze/tasks/` cache** (Quarto `freeze: auto` is set globally in `_quarto.yml`); CI renders from cache without re-evaluating. The author must `quarto render` locally with `masreml`, `masbayes`, and `CMplot` installed before committing — the `_freeze/tasks/<page>/figure-html/*.png` and `execute-results/html.json` files must be committed alongside any qmd source change. Reference pages and other non-task qmds that would still call the packages use `eval: false`. CI installs `tidyverse, DT, gt, knitr, rmarkdown, data.table, downlit, xml2`. `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: 'true'` is set for the June 2026 default switch.

8. **Phenotype scenarios shipped: 2** (continuous + binary).

## Rust source paths (verified 2026-05-12)

Used by `{{< rust-src <pkg> <path> >}}` in `internals/*.qmd`. Verified on disk:

| Package | Files at `<pkg>/src/rust/src/` |
|---|---|
| `masreml` | `lib.rs`, dirs: `gwas/`, `matrix/`, `reml/`, `solver/`, `utils/` |
| `masbayes` | `lib.rs`, `bayesa.rs`, `bayesa_em.rs`, `bayesr.rs`, `bayesr_em.rs`, `matrix.rs`, `types.rs`, `utils.rs` |

The current `internals/masreml-rust.qmd` still references placeholder paths (`gblup/mod.rs`, `gwas/mod.rs`) that need fixing in PLAN.md Stage I. `gblup/` does NOT exist (code is split across `reml/`, `solver/`, `matrix/`). `gwas/` exists as a directory but verify `mod.rs` inside before referencing.

## Common dev commands

```bash
cd /Users/aguswibowo/Library/CloudStorage/OneDrive-JamesCookUniversity/JCU/Research/Genomics/Analysis/genomic_prediction/masgenomics-docs

# Live preview while editing (auto-reload on save)
quarto preview

# One-shot render (uses _freeze cache)
quarto render

# Force clean rebuild
rm -rf _site _freeze && quarto render

# Re-pull reference docs from latest masreml + masbayes roxygen
# (manual — CI does not run this; package repos are NOT modified)
Rscript _scripts/render-reference.R

# Re-generate demo data (manual; rarely needed)
Rscript demo-data/generate.R main
Rscript demo-data/generate.R toy

# Publish path (ONLY at Stage Z, after user APPROVE in PLAN.md):
git add <specific-files>           # NEVER git add . / -A
git commit -m "type: short message"  # NO Co-Authored-By trailers
git push
```

## Open work (current state)

Site structure + demo data + input-data + theory (full multi-page redesign including REML, BLUP sub-pages, Bayesian sub-pages) + Home + Installation + maspipeline reference (overview + 4 stage pages) + Task tutorials (GP/GWAS each split into overview + masreml + masbayes sub-pages, with live-eval chunks and `_freeze/tasks/` cache populated) are all done locally (not yet pushed). Remaining work is in `internals/`.

| Area | Files | Status |
|---|---|---|
| Task tutorials | `tasks/{index, genomic-prediction/{index,masreml,masbayes}, gwas/{index,masreml,masbayes}}.qmd` | **Implemented** — live-eval chunks render summary tables, GEBV histograms, CV per-fold tables, CMplot Manhattan / PIP / WPPA plots, QTL recovery sanity checks. Output cached to `_freeze/tasks/` (commit alongside qmd) |
| Internals | `internals/{masreml-rust,masbayes-rust}.qmd` | PLAN.md Stage I — module map per file, FFI boundary |
| API tracker | `_scripts/check-api.R` + baseline | PLAN.md Stage A — drift detection |
| Final push | every locally-changed file since `dbb8ac8` | PLAN.md Stage Z — single consolidated commit + push, only after the user's final APPROVE |

## Repo layout (current)

```
masgenomics-docs/
├── _quarto.yml                          # navbar: Home/Installation/Theory/Task/Input Data/Reference▼/Internals
├── _extensions/rust-src/                # Quarto extension (NOT _scripts/shortcodes/)
│   ├── _extension.yml
│   └── rust-src.lua
├── _scripts/render-reference.R          # altdoc-on-temp-copy
├── assets/custom.scss                   # brand color, code-copy, DT dark overrides, inline-code
├── index.qmd                            # Home — Welcome, Why Rust?, Where to go next, Citation
├── installation/index.qmd               # Rust+cargo + masbayes/masreml install
├── theory/
│   ├── index.qmd
│   ├── biallelic-and-multiallelic.qmd
│   ├── mixed-model-blup-family.qmd      # BLUP overview (TL;DR, assumptions, MME, When to Use)
│   ├── blup-reml.qmd                    # REML algorithms (EM, NR, AI-REML)
│   ├── blup-pblup.qmd                   # PBLUP sub-page
│   ├── blup-gblup.qmd                   # GBLUP sub-page (SNP + MH)
│   ├── blup-gwablup.qmd                 # GWABLUP sub-page
│   ├── bayesian-alphabet.qmd            # Bayesian overview (TL;DR, shared model, When to Use)
│   ├── bayes-mcmc-em.qmd                # MCMC & EM algorithms sub-page
│   ├── bayes-bayesa.qmd                 # BayesA sub-page (SNP + MH)
│   ├── bayes-bayesr.qmd                 # BayesR sub-page (SNP + MH, PIP)
│   └── continuous-and-binary-trait.qmd
├── input-data/                          # demo-data presentation merged in here
│   ├── index.qmd                        # bundle entry-point (load_data + ZIPs)
│   ├── snp.qmd
│   ├── microhaplotypes.qmd
│   ├── phenotype.qmd                    # incl. QTL ground truth
│   └── pedigree.qmd
├── demo-data/                           # asset tree only — no rendered page
│   ├── generate.R
│   ├── main/                            # canonical .rds + derived flat exports + per-chr MH + qtl/
│   ├── toy/
│   ├── main.zip
│   └── toy.zip
├── tasks/                               # live-eval tutorials, freeze cache committed
│   ├── index.qmd                        # top-level Task overview
│   ├── genomic-prediction/
│   │   ├── index.qmd                    # GP overview (what / why / prerequisites / GP vs GWAS)
│   │   ├── masreml.qmd                  # GBLUP + CV + binary, SNP + MH
│   │   └── masbayes.qmd                 # BayesA + BayesR (MCMC + EM), manual CV, binary, SNP + MH
│   └── gwas/
│       ├── index.qmd                    # GWAS overview
│       ├── masreml.qmd                  # run_gwas + GWABLUP + QTL recovery (mirrors examples/04_gwas.R)
│       └── masbayes.qmd                 # BayesR PIP / WPPA (mirrors masbayes/examples/04_gwas.R)
├── reference/
│   ├── maspipeline/                     # overview + 4 stage pages (split from README)
│   │   ├── index.qmd                    # General info + Installation + footer meta
│   │   ├── data-preparation.qmd
│   │   ├── phasing.qmd
│   │   ├── genotype-to-haplotype.qmd
│   │   └── microhaplotype-discovery.qmd
│   ├── masreml/<index + 16 .qmd>        # auto-generated via _scripts/render-reference.R, committed
│   └── masbayes/<index + 9 .qmd>        # auto-generated, committed (run_bayesa_em/_mcmc and run_bayesr_em/_mcmc intentionally not surfaced — umbrella run_bayesa/run_bayesr only)
├── internals/<3 .qmd>                   # module-map TODO — PLAN.md Stage I
├── .github/workflows/publish.yml
├── README.md
├── CITATION.bib
├── .gitignore                           # CLAUDE.md, PLAN.md, REVISION_PLAN.md, docs/ ignored
├── CLAUDE.md                            # ← THIS FILE — local only
└── PLAN.md                              # local only
```
