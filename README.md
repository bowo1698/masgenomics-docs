# masgenomics-docs

Documentation site for the **masgenomics** suite, companion R packages with
shared Rust backends for genomic prediction on SNP and microhaplotype data.

- [**masbayes**](https://github.com/bowo1698/masbayes): Bayesian genomic
  prediction (BayesA, BayesR; MCMC and EM samplers).
- [**masreml**](https://github.com/bowo1698/masreml): REML mixed models,
  BLUP / GBLUP / GWABLUP, GWAS, cross-validation, threshold models for
  binary traits.
- [**maspipeline**](https://github.com/bowo1698/maspipeline): Rust CLI
  for phasing and microhaplotype discovery; produces inputs for the two
  R packages above.

Both R packages support **SNP** and **microhaplotype** genotypes natively
and handle **continuous** and **binary** traits through a single unified
API.

**Live site:** https://bowo1698.github.io/masgenomics-docs/

The site covers theory (mixed models, Bayesian alphabet), task tutorials
(genomic prediction, GWAS), input data formats with bundled demo, function
reference, and Rust internals.

## Citation

See [`CITATION.bib`](CITATION.bib).

## License

MIT — see [`LICENSE`](LICENSE). Code examples follow their respective
package licenses.
