# masgenomics-docs

Documentation site for the [masbayes](https://github.com/bowo1698/masbayes) and
[masreml](https://github.com/bowo1698/masreml) R/Rust packages for genomic
prediction with SNP and microhaplotype data.

**Live site:** https://bowo1698.github.io/masgenomics-docs/

## Local development

Requires Quarto ≥1.5 and R ≥4.4.

```bash
# Install R deps
Rscript -e 'install.packages(c("altdoc","reactable","gt","tidyverse"))'

# Live preview (auto-reload)
quarto preview

# One-shot render to ./_site/
quarto render

# Regenerate auto reference pages from package roxygen
Rscript _scripts/render-reference.R
```

## Repo layout

See `docs/plans/2026-05-10-masgenomics-docs-design.md` in the parent project
for the full design rationale.

## License

MIT (content) + package licenses for code examples.
