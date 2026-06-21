# masgenomics-docs

Documentation site for the **masgenomics** suite, companion R packages with
shared Rust backends for genomic prediction on SNP and microhaplotype data.

- [**maspipeline**](https://github.com/bowo1698/maspipeline): Rust CLI
  for phasing and microhaplotype discovery; produces inputs for the two
  R packages above.
- [**masreml**](https://github.com/bowo1698/masreml): REML-BLUP, GWAS, and GWABLUP for
  biallelic SNP and multi-allelic markers
- [**masbayes**](https://github.com/bowo1698/masbayes): Bayesian genomic prediction for
  biallelic SNP and multi-allelic markers

Both R packages support **SNP** and **microhaplotype** genotypes natively
and handle **continuous** and **binary** traits through a single unified
API.

**Live site:** https://bowo1698.github.io/masgenomics-docs/

The site covers theory (mixed models, Bayesian alphabet), task tutorials
(genomic prediction, GWAS), input data formats with bundled demo, function
reference, and Rust internals..

## Citation

See [`CITATION.bib`](CITATION.bib).

## License

MIT — see [`LICENSE`](LICENSE). Code examples follow their respective
package licenses.
