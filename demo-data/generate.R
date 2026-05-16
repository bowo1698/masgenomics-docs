#!/usr/bin/env Rscript
# demo-data/generate.R
#
# Regenerate the masgenomics-docs demo dataset locally. The canonical
# upstream copy lives in `../masbayes/inst/extdata/demo_data{,_small}.rds`
# and `../masreml/inst/extdata/demo_data{,_small}.rds` (byte-identical
# between packages, produced by `../tools/make_demo_data.R`). This script
# is the masgenomics-docs-local counterpart of that generator: same
# numerical logic, same seeds, but output is routed to
# `demo-data/<size>/<canonical-filename>` rather than into either
# sibling-package `inst/extdata/`.
#
# Mapping between this script's CLI sizes and the upstream config labels:
#
#   Rscript demo-data/generate.R main   # large config -> demo-data/main/demo_data.rds
#   Rscript demo-data/generate.R toy    # small config -> demo-data/toy/demo_data_small.rds
#
# The "main" / "toy" labels are the masgenomics-docs reading-audience
# convention; "large" / "small" are the upstream config keys. The byte
# layout of the .rds is identical to upstream, so a freshly regenerated
# `demo-data/main/demo_data.rds` should match the upstream copy under
# `masbayes/inst/extdata/demo_data.rds` byte-for-byte.
#
# Dependency: `library(masbayes)` is required at generation time -- the
# QTL@MH simulation uses `masbayes::construct_wah_matrix()` to build the
# W_alpha-h matrix. Install masbayes locally before running this script
# (CI does NOT run generate.R; the .rds files are committed assets).

suppressPackageStartupMessages({
  library(masbayes)
})

# Same two configurations as ../tools/make_demo_data.R. Kept verbatim so
# `main` / `toy` outputs here remain byte-identical to the upstream
# `inst/extdata/demo_data{,_small}.rds`. If upstream changes, sync the
# values here in lockstep.
CONFIGS <- list(
  main = list(
    n_sires            = 10L,
    n_dams             = 10L,
    n_offspring_family = 20L,
    n_snp              = 500L,
    n_blocks           = 250L,
    n_snp_per_block    = 2L,
    n_qtl              = 20L,
    h2                 = 0.5,
    seed               = 42L,
    split_seed         = 41L,
    n_train_per_family = 16L,
    n_test_per_family  = 4L,
    out_file           = "demo_data.rds"
  ),
  toy = list(
    n_sires            = 10L,
    n_dams             = 10L,
    n_offspring_family = 10L,
    n_snp              = 50L,
    n_blocks           = 25L,
    n_snp_per_block    = 2L,
    n_qtl              = 5L,
    h2                 = 0.5,
    seed               = 42L,
    split_seed         = 41L,
    n_train_per_family = 8L,
    n_test_per_family  = 2L,
    out_file           = "demo_data_small.rds"
  )
)

# Base-3 polynomial encoder. Same convention used in
# `examples/03_marker_QTL_congruency_theory.R` upstream -- base 3 (not 2)
# leaves room for multi-allelic SNPs even though current binary inputs
# yield a dense subset of integer ids.
encode_hap <- function(mat) {
  apply(mat, 1L, function(x) sum(x * 3L ^ (seq_along(x) - 1L)))
}

# Derived MicrohapsSel-pipeline-compatible per-chromosome MH TSVs.
#
# Why this exists separately from the canonical .rds: the .rds carries
# block-level integer-encoded MH alleles in d$mh; MicrohapsSel-style
# tooling consumes a per-chromosome tab-delimited file with the
# duplicated 4-column-per-block layout. The duplication encodes
# "haplotype value applies to both SNPs in the block", which only makes
# unambiguous sense when n_snp_per_block == 2 (each polynomial digit is
# in {1, 2} and decodes uniquely from the encoded id). Both demo configs
# satisfy this; the stopifnot below guards against future changes.
#
# Output is restricted to demo-data/<size>/mh_genotypes/. Sibling-
# package inst/extdata/ is NEVER written from this script.
write_per_chr_mh_files <- function(demo_data, out_dir, n_snp_per_block) {
  stopifnot(n_snp_per_block == 2L)
  mh_dir <- file.path(out_dir, "mh_genotypes")
  if (!dir.exists(mh_dir)) dir.create(mh_dir, recursive = TRUE)

  ind_ids <- rownames(demo_data$mh)
  n_ind   <- length(ind_ids)
  chrs    <- sort(unique(demo_data$map_mh$chr))

  for (chr in chrs) {
    blocks_global <- which(demo_data$map_mh$chr == chr)
    n_blocks_chr  <- length(blocks_global)
    n_cols        <- 1L + 4L * n_blocks_chr

    # Header uses per-chr LOCAL block numbering (1..n_blocks_chr).
    # Each block contributes 4 columns: (snp1_h1, snp2_h1, snp1_h2,
    # snp2_h2). The duplicated `hap_<i>_1, hap_<i>_1` names mirror the
    # MicrohapsSel convention "haplotype i column-1 applies to SNPs 1
    # and 2 of the block".
    headers <- character(n_cols)
    headers[1L] <- "ID"
    for (b_local in seq_len(n_blocks_chr)) {
      base_col <- 1L + 4L * (b_local - 1L) + 1L
      hdr <- sprintf("hap_%d", b_local)
      headers[base_col]      <- paste0(hdr, "_1")
      headers[base_col + 1L] <- paste0(hdr, "_1")
      headers[base_col + 2L] <- paste0(hdr, "_2")
      headers[base_col + 3L] <- paste0(hdr, "_2")
    }

    mat <- matrix("", nrow = n_ind, ncol = n_cols)
    mat[, 1L] <- ind_ids
    for (b_local in seq_len(n_blocks_chr)) {
      b_global <- blocks_global[b_local]
      e_h1 <- demo_data$mh[, 2L * b_global - 1L]
      e_h2 <- demo_data$mh[, 2L * b_global]

      # Decode encoded = snp1 + 3*snp2 with snp1, snp2 in {1, 2}.
      # snp2 holds the higher-order base-3 digit; recover it first.
      s2_h1 <- (e_h1 - 1L) %/% 3L
      s1_h1 <- e_h1 - 3L * s2_h1
      s2_h2 <- (e_h2 - 1L) %/% 3L
      s1_h2 <- e_h2 - 3L * s2_h2

      base_col <- 1L + 4L * (b_local - 1L) + 1L
      mat[, base_col]      <- as.character(s1_h1)
      mat[, base_col + 1L] <- as.character(s2_h1)
      mat[, base_col + 2L] <- as.character(s1_h2)
      mat[, base_col + 3L] <- as.character(s2_h2)
    }

    # writeLines instead of write.table because data.frame and matrix
    # column-name machinery in R deduplicates names ("hap_1_1" ->
    # "hap_1_1.1"). The MicrohapsSel reader expects raw duplicates.
    out_file <- file.path(mh_dir, sprintf("hap_geno_%d", chr))
    writeLines(
      c(paste(headers,                     collapse = "\t"),
        apply(mat, 1L, paste,              collapse = "\t")),
      out_file
    )
  }
}

# Derived flat exports of the canonical .rds. Same one-direction rule
# as write_per_chr_mh_files: these files are projections of demo_data
# slots and live ONLY in demo-data/<size>/. They are convenient for
# users who want plain CSV/TXT inputs (matching the legacy
# MicrohapsSel-pipeline layout) without having to readRDS() the
# canonical bundle themselves.
write_derived_flat_files <- function(demo_data, out_dir) {

  # ── Phenotypes split by trait kind ──────────────────────────────────────
  # Continuous: SNP- and MH-architecture y plus the true breeding values
  # used to simulate each. Binary: same but with the median-thresholded
  # y. tbv columns are repeated in both files so each phenotype CSV is
  # self-contained for downstream regression / classification tooling.
  pheno_cont <- data.frame(
    id             = demo_data$pheno$id,
    sex            = demo_data$pheno$sex,
    y_cont_qtl_snp = demo_data$pheno$y_cont_qtl_snp,
    y_cont_qtl_mh  = demo_data$pheno$y_cont_qtl_mh,
    tbv_qtl_snp    = demo_data$pheno$tbv_qtl_snp,
    tbv_qtl_mh     = demo_data$pheno$tbv_qtl_mh,
    stringsAsFactors = FALSE
  )
  pheno_bin <- data.frame(
    id            = demo_data$pheno$id,
    sex           = demo_data$pheno$sex,
    y_bin_qtl_snp = demo_data$pheno$y_bin_qtl_snp,
    y_bin_qtl_mh  = demo_data$pheno$y_bin_qtl_mh,
    tbv_qtl_snp   = demo_data$pheno$tbv_qtl_snp,
    tbv_qtl_mh    = demo_data$pheno$tbv_qtl_mh,
    stringsAsFactors = FALSE
  )
  write.csv(pheno_cont,
            file.path(out_dir, "pheno_continuous_G1.csv"),
            row.names = FALSE)
  write.csv(pheno_bin,
            file.path(out_dir, "pheno_binary_G1.csv"),
            row.names = FALSE)

  # ── SNP genotype dosage matrix as wide CSV ──────────────────────────────
  geno_df <- data.frame(id = rownames(demo_data$snp),
                        demo_data$snp,
                        stringsAsFactors = FALSE,
                        check.names = FALSE)
  write.csv(geno_df, file.path(out_dir, "geno_G1.csv"), row.names = FALSE)

  # ── SNP map (tab-delimited so it can be read with read.table default) ──
  write.table(demo_data$map_snp,
              file.path(out_dir, "snp_map.txt"),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)

  # ── Microhaplotype coordinates (verbatim d$map_mh, schema matches the
  # maspipeline microhaplotype_coordinates.csv) ───────────────────────────
  write.csv(demo_data$map_mh,
            file.path(out_dir, "microhaplotype_coordinates.csv"),
            row.names = FALSE)

  # ── QTL info: one RDS per slot, mirroring the legacy `qtl/` directory.
  # snp_ids and mh_blocks are the human-readable identifiers (SNP names
  # and block_ids) rather than numeric indices, so downstream QTL
  # recovery analysis can join against d$map_snp / d$map_mh without an
  # index lookup.
  qtl_dir <- file.path(out_dir, "qtl")
  if (!dir.exists(qtl_dir)) dir.create(qtl_dir, recursive = TRUE)

  snp_ids <- colnames(demo_data$snp)[demo_data$qtl$snp_idx]

  # d$qtl$mh_idx indexes columns of W_ah, which holds (n_alleles_b - 1)
  # cols per block (one allele is dummy-encoded as reference). The
  # block label of W_ah column k is reconstructed by repeating each
  # block_id (n_alleles_b - 1) times in d$allele_freq order.
  af        <- demo_data$allele_freq
  blk       <- unique(af$haplotype)
  counts    <- as.integer(table(factor(af$haplotype, levels = blk))) - 1L
  wah_blk   <- rep(blk, times = counts)
  mh_blocks <- wah_blk[demo_data$qtl$mh_idx]

  saveRDS(snp_ids,                   file.path(qtl_dir, "qtl_snp_ids.rds"))
  saveRDS(mh_blocks,                 file.path(qtl_dir, "qtl_mh_blocks.rds"))
  saveRDS(demo_data$qtl$effects_snp, file.path(qtl_dir, "effects_snp.rds"))
  saveRDS(demo_data$qtl$effects_mh,  file.path(qtl_dir, "effects_mh.rds"))
}

generate_demo <- function(cfg) {

  N_SIRES                <- cfg$n_sires
  N_DAMS                 <- cfg$n_dams
  N_OFFSPRING_PER_FAMILY <- cfg$n_offspring_family
  N                      <- N_SIRES * N_OFFSPRING_PER_FAMILY
  N_SNP                  <- cfg$n_snp
  N_BLOCKS               <- cfg$n_blocks
  N_SNP_PER_BLOCK        <- cfg$n_snp_per_block
  N_QTL                  <- cfg$n_qtl
  H2                     <- cfg$h2
  SEED                   <- cfg$seed
  N_TRAIN_PER_FAMILY     <- cfg$n_train_per_family
  N_TEST_PER_FAMILY      <- cfg$n_test_per_family

  IND_NAMES  <- sprintf("IND%03d", seq_len(N))
  SIRE_NAMES <- sprintf("SIRE%02d", seq_len(N_SIRES))
  DAM_NAMES  <- sprintf("DAM%02d",  seq_len(N_DAMS))
  SNP_NAMES  <- sprintf("SNP%03d",  seq_len(N_SNP))
  FAM_NAMES  <- sprintf("fam_%02d", seq_len(N_SIRES))

  stopifnot(N_BLOCKS * N_SNP_PER_BLOCK == N_SNP)
  stopifnot(N_TRAIN_PER_FAMILY + N_TEST_PER_FAMILY == N_OFFSPRING_PER_FAMILY)

  set.seed(SEED)

  # ── 1. Founder haplotypes ──────────────────────────────────────────────────
  # Family structure is engineered into the GRM by having each offspring
  # inherit one strand per locus from its parents' strands. Free recombination
  # per locus is biologically extreme but produces enough genetic diversity
  # that all full-sibs per family are distinct while still sharing the
  # parental allele pools that drive related-pair G.
  maf      <- runif(N_SNP, 0.1, 0.5)
  sire_h1  <- matrix(0L, N_SIRES, N_SNP)
  sire_h2  <- matrix(0L, N_SIRES, N_SNP)
  dam_h1   <- matrix(0L, N_DAMS,  N_SNP)
  dam_h2   <- matrix(0L, N_DAMS,  N_SNP)
  for (j in seq_len(N_SNP)) {
    sire_h1[, j] <- as.integer(rbinom(N_SIRES, 1L, maf[j]))
    sire_h2[, j] <- as.integer(rbinom(N_SIRES, 1L, maf[j]))
    dam_h1 [, j] <- as.integer(rbinom(N_DAMS,  1L, maf[j]))
    dam_h2 [, j] <- as.integer(rbinom(N_DAMS,  1L, maf[j]))
  }

  # ── 2. Offspring haplotypes via inheritance ────────────────────────────────
  hap1      <- matrix(0L, N, N_SNP)
  hap2      <- matrix(0L, N, N_SNP)
  snp       <- matrix(0L, N, N_SNP)
  family_id <- character(N)
  sire_of   <- character(N)
  dam_of    <- character(N)

  idx <- 1L
  for (k in seq_len(N_SIRES)) {
    for (m in seq_len(N_OFFSPRING_PER_FAMILY)) {
      s_pick <- as.integer(rbinom(N_SNP, 1L, 0.5))
      d_pick <- as.integer(rbinom(N_SNP, 1L, 0.5))
      hap1[idx, ] <- ifelse(s_pick == 0L, sire_h1[k, ], sire_h2[k, ])
      hap2[idx, ] <- ifelse(d_pick == 0L, dam_h1 [k, ], dam_h2 [k, ])
      snp [idx, ] <- hap1[idx, ] + hap2[idx, ]
      family_id[idx] <- FAM_NAMES[k]
      sire_of  [idx] <- SIRE_NAMES[k]
      dam_of   [idx] <- DAM_NAMES [k]
      idx <- idx + 1L
    }
  }
  rownames(snp) <- IND_NAMES
  colnames(snp) <- SNP_NAMES

  # ── 3. Hap-block matrix consumable by both packages ────────────────────────
  hap_all <- matrix(0L, N, N_SNP * 2L)
  for (j in seq_len(N_SNP)) {
    hap_all[, 2L * j - 1L] <- hap1[, j] + 1L
    hap_all[, 2L * j     ] <- hap2[, j] + 1L
  }

  HAP_COLS_PER_BLOCK <- N_SNP_PER_BLOCK * 2L
  hap_reordered <- matrix(0L, N, N_SNP * 2L)
  col_out <- 1L
  for (b in seq_len(N_BLOCKS)) {
    snp_range <- ((b - 1L) * N_SNP_PER_BLOCK + 1L):(b * N_SNP_PER_BLOCK)
    for (j in snp_range) {
      hap_reordered[, col_out]      <- hap_all[, 2L * j - 1L]
      hap_reordered[, col_out + 1L] <- hap_all[, 2L * j     ]
      col_out <- col_out + 2L
    }
  }

  hap_block        <- matrix(0L, N, N_BLOCKS * 2L)
  allele_freq_list <- list(haplotype = character(),
                           allele    = integer(),
                           freq      = numeric())

  for (b in seq_len(N_BLOCKS)) {
    cols    <- ((b - 1L) * HAP_COLS_PER_BLOCK + 1L):(b * HAP_COLS_PER_BLOCK)
    hap_sub <- hap_reordered[, cols, drop = FALSE]
    h1_id   <- encode_hap(hap_sub[, seq(1L, HAP_COLS_PER_BLOCK, by = 2L),
                                  drop = FALSE])
    h2_id   <- encode_hap(hap_sub[, seq(2L, HAP_COLS_PER_BLOCK, by = 2L),
                                  drop = FALSE])
    hap_block[, 2L * b - 1L] <- h1_id
    hap_block[, 2L * b     ] <- h2_id
    tbl     <- table(c(h1_id, h2_id))
    freqs   <- as.numeric(tbl) / sum(tbl)
    alleles <- as.integer(names(tbl))
    allele_freq_list$haplotype <- c(allele_freq_list$haplotype,
                                    rep(paste0("block_", b), length(alleles)))
    allele_freq_list$allele    <- c(allele_freq_list$allele, alleles)
    allele_freq_list$freq      <- c(allele_freq_list$freq,   freqs)
  }
  storage.mode(hap_block) <- "integer"
  rownames(hap_block) <- IND_NAMES

  block_id <- paste0("block_", rep(seq_len(N_BLOCKS), each = 2L))
  attr(hap_block, "block_id") <- block_id

  # ── 4. W_mh for QTL@MH simulation ──────────────────────────────────────────
  wah  <- construct_wah_matrix(hap_block, block_id, allele_freq_list, NULL, TRUE)
  W_mh <- wah$W_ah

  # ── 5. QTL positions and effects ───────────────────────────────────────────
  qtl_snp_idx <- sort(sample.int(N_SNP,      N_QTL))
  qtl_mh_idx  <- sort(sample.int(ncol(W_mh), N_QTL))

  raw_snp     <- rnorm(N_QTL)
  raw_mh      <- rnorm(N_QTL)
  effects_snp <- raw_snp / sqrt(sum(raw_snp ^ 2))
  effects_mh  <- raw_mh  / sqrt(sum(raw_mh  ^ 2))

  beta_snp <- rep(0, N_SNP)
  beta_snp[qtl_snp_idx] <- effects_snp
  beta_mh  <- rep(0, ncol(W_mh))
  beta_mh [qtl_mh_idx]  <- effects_mh

  # ── 6. TBV and standardise ─────────────────────────────────────────────────
  snp_centred <- scale(snp, center = TRUE, scale = FALSE)
  tbv_snp_raw <- as.numeric(snp_centred %*% beta_snp)
  tbv_mh_raw  <- as.numeric(W_mh        %*% beta_mh)

  standardise <- function(x) (x - mean(x)) / sd(x)
  tbv_qtl_snp <- standardise(tbv_snp_raw)
  tbv_qtl_mh  <- standardise(tbv_mh_raw)

  # ── 7. Fixed-effect covariate (sex) ────────────────────────────────────────
  sex_codes <- rep(c(0L, 1L), length.out = N)
  sex_codes <- sample(sex_codes)
  sex_effect_target <- 0.5

  # ── 8. Simulate phenotypes ─────────────────────────────────────────────────
  sigma2_e <- (1 - H2) / H2

  simulate_y <- function(tbv) {
    e        <- rnorm(N, 0, sqrt(sigma2_e))
    y_base   <- tbv + e
    sex_beta <- sex_effect_target * sd(y_base)
    y_cont   <- y_base + sex_codes * sex_beta
    list(
      y_cont   = y_cont,
      y_bin    = as.integer(y_cont > median(y_cont)),
      sex_beta = sex_beta
    )
  }

  set.seed(SEED + 1L)
  sim_snp <- simulate_y(tbv_qtl_snp)
  set.seed(SEED + 2L)
  sim_mh  <- simulate_y(tbv_qtl_mh)

  # ── 9. Within-family train/test split ──────────────────────────────────────
  # Split seed chosen empirically during R3 tuning: with n_test small, the
  # Pearson r between GEBV and TBV has high sampling SD. The fixed seed
  # yields a split whose r_test_g is close to the underlying CV-stable
  # accuracy for both packages -- a representative draw rather than a
  # lucky tail.
  set.seed(cfg$split_seed)
  train_idx <- integer()
  test_idx  <- integer()
  for (k in seq_len(N_SIRES)) {
    fam_members  <- which(family_id == FAM_NAMES[k])
    picked_test  <- sample(fam_members, N_TEST_PER_FAMILY)
    picked_train <- setdiff(fam_members, picked_test)
    train_idx    <- c(train_idx, picked_train)
    test_idx     <- c(test_idx,  picked_test)
  }
  train_idx <- sort(train_idx)
  test_idx  <- sort(test_idx)

  # ── 9.5. Physical-position maps (synthetic, deterministic) ────────────────
  # PURE FUNCTION OF SIZES — no RNG. Placed AFTER the last set.seed() /
  # rnorm() / sample() call (Step 9) so the random stream that produced
  # `snp`, `hap_block`, `qtl_*`, phenotypes, and the train/test split is
  # byte-identical to upstream. Adding this block earlier would disturb
  # downstream RNG draws and break byte-equality with
  # `inst/extdata/demo_data{,_small}.rds`.
  #
  # Layout: N_CHR synthetic chromosomes, evenly partitioning the SNP set;
  # each MH block spans the (2i-1, 2i) SNP pair (already the existing
  # generator's relation since N_SNP == N_BLOCKS * N_SNP_PER_BLOCK with
  # N_SNP_PER_BLOCK == 2). The maspipeline `microhaplotype_coordinates.csv`
  # schema (block_id, chr, start_pos, end_pos, n_snps) is reproduced
  # verbatim so production pipelines can consume `d$map_mh` without
  # translation.
  N_CHR          <- 5L
  SNP_SPACING_BP <- 100000L
  CHR_BASE_BP    <- 1000000L

  stopifnot(
    N_SNP %% N_CHR == 0L,
    N_BLOCKS %% N_CHR == 0L,
    N_SNP == N_BLOCKS * N_SNP_PER_BLOCK
  )

  snps_per_chr <- N_SNP %/% N_CHR
  snp_chr      <- rep(seq_len(N_CHR), each = snps_per_chr)
  snp_pos      <- rep(seq(CHR_BASE_BP, by = SNP_SPACING_BP,
                          length.out = snps_per_chr), N_CHR)
  map_snp <- data.frame(
    SNP   = SNP_NAMES,
    CHROM = as.integer(snp_chr),
    POS   = as.integer(snp_pos),
    stringsAsFactors = FALSE
  )

  block_first <- seq(1L, by = N_SNP_PER_BLOCK, length.out = N_BLOCKS)
  block_last  <- seq(N_SNP_PER_BLOCK, by = N_SNP_PER_BLOCK,
                     length.out = N_BLOCKS)
  # Each block must live entirely on one chromosome. The N_SNP %% N_CHR == 0
  # guard above plus equal partitioning makes this true by construction;
  # the assertion guards against future size changes that violate it.
  stopifnot(snp_chr[block_first] == snp_chr[block_last])

  map_mh <- data.frame(
    block_id  = paste0("block_", seq_len(N_BLOCKS)),
    chr       = as.integer(snp_chr[block_first]),
    start_pos = as.integer(snp_pos[block_first]),
    end_pos   = as.integer(snp_pos[block_last]),
    n_snps    = as.integer(rep(N_SNP_PER_BLOCK, N_BLOCKS)),
    stringsAsFactors = FALSE
  )

  # ── 10. Pedigree ───────────────────────────────────────────────────────────
  pedigree <- data.frame(
    id   = c(SIRE_NAMES, DAM_NAMES, IND_NAMES),
    sire = c(rep(NA_character_, N_SIRES),
             rep(NA_character_, N_DAMS),
             sire_of),
    dam  = c(rep(NA_character_, N_SIRES),
             rep(NA_character_, N_DAMS),
             dam_of),
    stringsAsFactors = FALSE
  )

  # ── 11. Phenotype data frame ───────────────────────────────────────────────
  pheno <- data.frame(
    id             = IND_NAMES,
    sex            = factor(c("F", "M")[sex_codes + 1L], levels = c("F", "M")),
    y_cont_qtl_snp = sim_snp$y_cont,
    y_cont_qtl_mh  = sim_mh$y_cont,
    y_bin_qtl_snp  = sim_snp$y_bin,
    y_bin_qtl_mh   = sim_mh$y_bin,
    tbv_qtl_snp    = tbv_qtl_snp,
    tbv_qtl_mh     = tbv_qtl_mh,
    stringsAsFactors = FALSE
  )

  # ── 12. Assemble ───────────────────────────────────────────────────────────
  # allele_freq is the table required by construct_wah_matrix() when no
  # reference_structure is supplied. Bundling it here lets users call
  # construct_wah_matrix(d$mh, attr(d$mh, "block_id"), d$allele_freq)
  # directly without recomputing.
  demo_data <- list(
    snp         = snp,
    mh          = hap_block,
    allele_freq = allele_freq_list,
    pheno       = pheno,
    pedigree    = pedigree,
    qtl         = list(
      snp_idx     = qtl_snp_idx,
      mh_idx      = qtl_mh_idx,
      effects_snp = effects_snp,
      effects_mh  = effects_mh
    ),
    meta = list(
      n               = N,
      n_snp           = N_SNP,
      n_blocks        = N_BLOCKS,
      n_snp_per_block = N_SNP_PER_BLOCK,
      n_qtl           = N_QTL,
      n_families      = N_SIRES,
      n_per_family    = N_OFFSPRING_PER_FAMILY,
      h2_target       = H2,
      sex_beta_snp    = sim_snp$sex_beta,
      sex_beta_mh     = sim_mh$sex_beta,
      seed            = SEED,
      split_seed      = cfg$split_seed,
      size            = cfg$out_file
    ),
    family_id = family_id,
    train_idx = train_idx,
    test_idx  = test_idx,
    map_snp   = map_snp,
    map_mh    = map_mh
  )

  # ── 13. Sanity checks ──────────────────────────────────────────────────────
  resid_snp       <- sim_snp$y_cont - tbv_qtl_snp - sex_codes * sim_snp$sex_beta
  resid_mh        <- sim_mh$y_cont  - tbv_qtl_mh  - sex_codes * sim_mh$sex_beta
  realised_h2_snp <- var(tbv_qtl_snp) / (var(tbv_qtl_snp) + var(resid_snp))
  realised_h2_mh  <- var(tbv_qtl_mh)  / (var(tbv_qtl_mh)  + var(resid_mh))

  train_fam_counts <- table(family_id[train_idx])
  test_fam_counts  <- table(family_id[test_idx])
  stopifnot(all(train_fam_counts == N_TRAIN_PER_FAMILY),
            all(test_fam_counts  == N_TEST_PER_FAMILY))

  sibs_in_train <- vapply(test_idx, function(i) {
    sum(family_id[train_idx] == family_id[i])
  }, integer(1))
  stopifnot(all(sibs_in_train == N_TRAIN_PER_FAMILY))

  test_wah <- construct_wah_matrix(
    hap_block, block_id, allele_freq_list, NULL, TRUE
  )
  stopifnot(is.numeric(test_wah$W_ah), nrow(test_wah$W_ah) == N)
  if (requireNamespace("masreml", quietly = TRUE)) {
    g_mh_check <- masreml::build_G_mh(mh_list = hap_block, ids = IND_NAMES)
    stopifnot(nrow(g_mh_check) == N, ncol(g_mh_check) == N)
  }

  list(
    demo_data       = demo_data,
    realised_h2_snp = realised_h2_snp,
    realised_h2_mh  = realised_h2_mh,
    prevalence_snp  = mean(sim_snp$y_bin),
    prevalence_mh   = mean(sim_mh$y_bin),
    n_wmh_columns   = ncol(W_mh)
  )
}

# ─── Run requested size(s) ──────────────────────────────────────────────────
# Output is routed only into demo-data/<size>/<canonical-filename>.
# Sibling-package inst/extdata/ is NEVER written from this script — that
# is governed by ../tools/make_demo_data.R upstream.

args <- commandArgs(trailingOnly = TRUE)
sizes <- if (length(args) == 0L) names(CONFIGS) else args
unknown <- setdiff(sizes, names(CONFIGS))
if (length(unknown) > 0L) {
  stop(sprintf("Unknown size(s): %s. Valid: %s",
               paste(unknown, collapse = ", "),
               paste(names(CONFIGS), collapse = ", ")))
}

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1L)$ofile)),
  error = function(e) "demo-data"
)
out_root <- if (basename(script_dir) == "demo-data") {
  script_dir
} else {
  file.path(script_dir, "demo-data")
}

for (size_label in sizes) {
  cfg <- CONFIGS[[size_label]]
  out_dir <- file.path(out_root, size_label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat(sprintf(">>> Generating size='%s' -> %s/%s\n",
              size_label, out_dir, cfg$out_file))
  res <- generate_demo(cfg)

  saveRDS(res$demo_data, file.path(out_dir, cfg$out_file))

  # Derived per-chromosome MH TSVs alongside the canonical .rds.
  # Only into demo-data/<size>/mh_genotypes/ — never into sibling-package
  # inst/extdata/.
  write_per_chr_mh_files(res$demo_data, out_dir, cfg$n_snp_per_block)

  # Derived flat phenotype / genotype / map / QTL files, also only
  # under demo-data/<size>/. Pure projections of the .rds; no RNG.
  write_derived_flat_files(res$demo_data, out_dir)

  cat(sprintf("  n=%d  n_snp=%d  n_blocks=%d  n_qtl=%d\n",
              cfg$n_sires * cfg$n_offspring_family,
              cfg$n_snp, cfg$n_blocks, cfg$n_qtl))
  cat(sprintf("  W_mh ncol         : %d\n", res$n_wmh_columns))
  cat(sprintf("  h2 realised       : SNP=%.3f  MH=%.3f\n",
              res$realised_h2_snp, res$realised_h2_mh))
  cat(sprintf("  prevalence        : SNP=%.3f  MH=%.3f\n",
              res$prevalence_snp, res$prevalence_mh))
  cat(sprintf("  split             : %d train / %d test (within-family %d/%d)\n",
              length(res$demo_data$train_idx),
              length(res$demo_data$test_idx),
              cfg$n_train_per_family, cfg$n_test_per_family))
  cat(sprintf("  written           : %s\n",
              file.path(out_dir, cfg$out_file)))
  cat("\n")
}
