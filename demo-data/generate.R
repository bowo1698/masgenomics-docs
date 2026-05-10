#!/usr/bin/env Rscript
# generate.R — parameterised demo-data generator for masgenomics-docs
#
# Usage:
#   Rscript demo-data/generate.R [main|toy]   # defaults to "main"
#
# Outputs: demo-data/<size>/genotypes/... and demo-data/<size>/{genotype,
#          microhap,phenotype,pedigree}.csv

suppressPackageStartupMessages(library(tidyverse))

# ── 0.  CLI argument ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
size <- if (length(args) >= 1) args[1] else "main"
if (!size %in% c("main", "toy"))
  stop("Argument must be 'main' or 'toy'. Got: ", size)

# ── 1.  Fixed seed (must come first, before any RNG call) ────────────────────
set.seed(42)

# ── 2.  Size parameters ───────────────────────────────────────────────────────
if (size == "main") {
  N_IND         <- 1000L
  N_CHR         <- 10L
  N_SNP_PER_CHR <- 500L
} else {                  # toy
  N_IND         <- 200L
  N_CHR         <- 5L
  N_SNP_PER_CHR <- 100L
}
IND_NAMES <- sprintf("IND%04d", seq_len(N_IND))

# ── 3.  Output paths ──────────────────────────────────────────────────────────
# Locate this script's directory robustly under Rscript
get_script_dir <- function() {
  # commandArgs gives "--file=<path>" when called via Rscript
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])),
                         mustWork = FALSE))
  }
  # Fallback: assume working directory contains demo-data/
  normalizePath(file.path(getwd(), "demo-data"), mustWork = FALSE)
}

BASE_ROOT <- get_script_dir()    # demo-data/
OUT_DIR   <- file.path(BASE_ROOT, size)   # demo-data/<size>/
BASE      <- file.path(OUT_DIR, "genotypes")  # demo-data/<size>/genotypes/

dirs <- c(
  file.path(BASE, "haploG0"),
  file.path(BASE, "mh_genotypes_G0"),
  file.path(BASE, "all_genotypes_G0"),
  file.path(BASE, "mh_info_ld_haploblock_G0", "stats"),
  file.path(BASE, "all_info_haploblock_G0",   "stats")
)
walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ── 4.  SNP positions & haplotype matrices ───────────────────────────────────
allele_pairs <- list(c("A","T"), c("G","C"), c("T","A"), c("C","G"))

chr_snps  <- map(seq_len(N_CHR), ~sort(sample(1000:500000, N_SNP_PER_CHR)))
chr_haplo <- map(seq_len(N_CHR),
                 ~matrix(sample(0:1, N_IND * 2L * N_SNP_PER_CHR, replace = TRUE),
                          nrow = N_IND * 2L, ncol = N_SNP_PER_CHR))

# ── 5.  VCF ──────────────────────────────────────────────────────────────────
vcf_header <- c(
  "##fileformat=VCFv4.1", "##source=dummy",
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
          IND_NAMES), collapse = "\t"))

vcf_rows <- do.call(rbind, lapply(seq_len(N_CHR), function(c) {
  do.call(rbind, lapply(seq_len(N_SNP_PER_CHR), function(si) {
    pos <- chr_snps[[c]][si]
    ap  <- allele_pairs[[(si - 1L) %% 4L + 1L]]
    gts <- sapply(seq_len(N_IND), function(ind)
      paste(chr_haplo[[c]][(ind - 1L)*2L + 1L, si],
            chr_haplo[[c]][(ind - 1L)*2L + 2L, si], sep = "|"))
    c(as.character(c), pos,
      sprintf("ssa%02d_%d_%s_%s", c, pos, ap[1], ap[2]),
      ap[1], ap[2], ".", "PASS", ".", "GT", gts)
  }))
}))

vcf_data_lines <- apply(vcf_rows, 1, paste, collapse = "\t")
writeLines(c(vcf_header, vcf_data_lines),
           file.path(BASE, "salmon_60k_pairs125bp.vcf"))

# ── 6.  haploG0 per chr ───────────────────────────────────────────────────────
walk(seq_len(N_CHR), function(c) {
  lines <- map_chr(seq_len(N_IND), function(i) {
    h1 <- paste(chr_haplo[[c]][(i - 1L)*2L + 1L, ], collapse = " ")
    h2 <- paste(chr_haplo[[c]][(i - 1L)*2L + 2L, ], collapse = " ")
    paste(IND_NAMES[i], h1, h2, sep = "\t")
  })
  writeLines(lines, file.path(BASE, "haploG0", paste0("chr", c)))
})

# ── 7.  Block-maker functions (verbatim from original) ───────────────────────
make_mh_blocks <- function(chr_snps, max_span = 125) {
  global_blk <- 1
  all_blocks <- list()
  for (c in seq_along(chr_snps)) {
    pos <- chr_snps[[c]]
    i   <- 1
    while (i <= length(pos)) {
      j <- i
      while (j < length(pos) && (pos[j + 1] - pos[i]) <= max_span &&
             (j - i) < 3) j <- j + 1
      snp_idx <- i:j
      all_blocks[[length(all_blocks) + 1]] <- list(
        blk = global_blk, chr = c, snp_idx = snp_idx, pos = pos[snp_idx])
      global_blk <- global_blk + 1
      i <- j + 1
    }
  }
  all_blocks
}

make_haplo_blocks <- function(chr_snps, window = 100000) {
  global_blk <- 1
  all_blocks <- list()
  for (c in seq_along(chr_snps)) {
    pos       <- chr_snps[[c]]
    win_start <- pos[1]
    i         <- 1
    while (i <= length(pos)) {
      win_end <- win_start + window
      snp_idx <- which(pos >= win_start & pos < win_end)
      if (length(snp_idx) > 0) {
        all_blocks[[length(all_blocks) + 1]] <- list(
          blk = global_blk, chr = c, snp_idx = snp_idx, pos = pos[snp_idx])
        global_blk <- global_blk + 1
      }
      win_start <- win_end
      i         <- max(snp_idx) + 1
      if (length(snp_idx) == 0) break
    }
  }
  all_blocks
}

mh_blocks    <- make_mh_blocks(chr_snps)
haplo_blocks <- make_haplo_blocks(chr_snps)

# ── 8.  Write block files + genotypes (verbatim from original) ───────────────
write_block_files <- function(blocks, info_dir, geno_dir) {
  by_chr <- split(blocks, map_int(blocks, ~.x$chr))
  walk(names(by_chr), function(c) {
    blk_list <- by_chr[[c]]
    ci       <- as.integer(c)

    writeLines(
      map_chr(blk_list,
              ~paste0("blk", .x$blk, "\t ",
                      paste(.x$snp_idx - 1, collapse = " "))),
      file.path(info_dir, paste0("hap_block_", c)))
    writeLines(
      map_chr(blk_list,
              ~paste0("blk", .x$blk, "\t",
                      paste(.x$snp_idx - 1, collapse = "\t"))),
      file.path(info_dir, paste0("hap_block_info_", c)))

    hdr <- c("ID", unlist(map(seq_along(blk_list), function(local_idx) c(
      paste0("hap_", local_idx, "_1"), paste0("hap_", local_idx, "_1"),
      paste0("hap_", local_idx, "_2"), paste0("hap_", local_idx, "_2")))))

    rows <- map_chr(seq_len(N_IND), function(i) {
      vals <- unlist(map(blk_list, function(b) {
        h1 <- chr_haplo[[ci]][(i - 1L)*2L + 1L, b$snp_idx]
        h2 <- chr_haplo[[ci]][(i - 1L)*2L + 2L, b$snp_idx]
        a1 <- strtoi(paste(h1, collapse = ""), base = 2L) + 1L
        a2 <- strtoi(paste(h2, collapse = ""), base = 2L) + 1L
        c(a1, a1, a2, a2)
      }))
      paste(c(IND_NAMES[i], vals), collapse = "\t")
    })
    writeLines(c(paste(hdr, collapse = "\t"), rows),
               file.path(geno_dir, paste0("hap_geno_", c)))
  })
}

write_block_files(mh_blocks,
                  file.path(BASE, "mh_info_ld_haploblock_G0"),
                  file.path(BASE, "mh_genotypes_G0"))
write_block_files(haplo_blocks,
                  file.path(BASE, "all_info_haploblock_G0"),
                  file.path(BASE, "all_genotypes_G0"))

# ── 9.  Stats CSVs (verbatim from original) ──────────────────────────────────
write_stats <- function(blocks, stats_dir) {
  coord <- map_dfr(blocks, ~tibble(
    block_id         = paste0("blk", .x$blk),
    chr              = .x$chr,
    start_pos        = .x$pos[1],
    end_pos          = .x$pos[length(.x$pos)],
    n_snps           = length(.x$snp_idx),
    physical_span_bp = .x$pos[length(.x$pos)] - .x$pos[1],
    pic              = round(runif(1, 0.05, 0.5), 4),
    NA_col           = NA))

  snp_det <- map_dfr(blocks, function(b) {
    map_dfr(seq_along(b$snp_idx), function(k) {
      ap <- allele_pairs[[(b$snp_idx[k] - 1L) %% 4L + 1L]]
      tibble(block_id  = paste0("blk", b$blk),
             chr       = b$chr,
             snp_index = b$snp_idx[k] - 1L,
             snp_id    = sprintf("ssa%02d_%d_%s_%s",
                                 b$chr, b$pos[k], ap[1], ap[2]),
             position  = b$pos[k])
    })
  })

  write_csv(coord,   file.path(stats_dir, "microhaplotype_coordinates.csv"))
  write_csv(snp_det, file.path(stats_dir, "snp_selection_detailed.csv"))
}

write_stats(mh_blocks,    file.path(BASE, "mh_info_ld_haploblock_G0", "stats"))
write_stats(haplo_blocks, file.path(BASE, "all_info_haploblock_G0",   "stats"))

# ── 10.  Wide-format CSVs for website ────────────────────────────────────────

## 10a. genotype.csv  (dosage 0/1/2 = h1 + h2) --------------------------------
n_snp_total  <- N_CHR * N_SNP_PER_CHR
snp_col_names <- unlist(lapply(seq_len(N_CHR), function(c)
  sprintf("snp_%02d_%04d", c, seq_len(N_SNP_PER_CHR))))

geno_mat <- do.call(cbind, lapply(seq_len(N_CHR), function(c) {
  h_mat <- chr_haplo[[c]]                 # rows: ind*2, cols: snps
  # haplotype 1 rows: 1,3,5,...  haplotype 2 rows: 2,4,6,...
  h1_rows <- seq(1, N_IND * 2L, by = 2L)
  h2_rows <- seq(2, N_IND * 2L, by = 2L)
  h_mat[h1_rows, ] + h_mat[h2_rows, ]    # dosage
}))                                       # N_IND x n_snp_total

geno_df <- as.data.frame(geno_mat)
colnames(geno_df) <- snp_col_names
geno_df <- cbind(data.frame(id = IND_NAMES, stringsAsFactors = FALSE), geno_df)
write.csv(geno_df, file.path(OUT_DIR, "genotype.csv"), row.names = FALSE)

## 10b. microhap.csv  (diploid MH genotype strings like "3/5") ---------------
n_mh   <- length(mh_blocks)
mh_ids <- sprintf("mh_%03d", seq_len(n_mh))

mh_mat <- matrix("", nrow = N_IND, ncol = n_mh)
for (bi in seq_len(n_mh)) {
  b  <- mh_blocks[[bi]]
  ci <- b$chr
  for (i in seq_len(N_IND)) {
    h1  <- chr_haplo[[ci]][(i - 1L)*2L + 1L, b$snp_idx]
    h2  <- chr_haplo[[ci]][(i - 1L)*2L + 2L, b$snp_idx]
    a1  <- strtoi(paste(h1, collapse = ""), base = 2L) + 1L
    a2  <- strtoi(paste(h2, collapse = ""), base = 2L) + 1L
    mh_mat[i, bi] <- paste(a1, a2, sep = "/")
  }
}

mh_df <- as.data.frame(mh_mat, stringsAsFactors = FALSE)
colnames(mh_df) <- mh_ids
mh_df <- cbind(data.frame(id = IND_NAMES, stringsAsFactors = FALSE), mh_df)
write.csv(mh_df, file.path(OUT_DIR, "microhap.csv"), row.names = FALSE)

## 10c. phenotype.csv ----------------------------------------------------------
n_qtl    <- max(20L, round(0.01 * n_snp_total))
qtl_idx  <- sort(sample(seq_len(n_snp_total), n_qtl, replace = FALSE))
qtl_eff  <- rnorm(n_qtl, 0, 1)

# dosage matrix already built as geno_mat
bv       <- scale(geno_mat[, qtl_idx] %*% qtl_eff)[, 1]
h2       <- 0.4
sigma_e  <- sqrt((1 - h2) / h2)
y_cont   <- bv + rnorm(N_IND, 0, sigma_e)
y_binary <- as.integer(y_cont > median(y_cont))

pheno_df <- data.frame(
  id           = IND_NAMES,
  sex          = sample(c("M", "F"), N_IND, replace = TRUE),
  age          = sample(1:5,        N_IND, replace = TRUE),
  y_continuous = round(y_cont, 4),
  y_binary     = y_binary,
  stringsAsFactors = FALSE)
write.csv(pheno_df, file.path(OUT_DIR, "phenotype.csv"), row.names = FALSE)

## 10d. pedigree.csv -----------------------------------------------------------
n_founders <- max(50L, round(0.1 * N_IND))
generation <- pmin(4L, ceiling(seq_len(N_IND) / (N_IND / 4)))

sire_col <- character(N_IND)
dam_col  <- character(N_IND)
sire_col[seq_len(n_founders)] <- NA_character_
dam_col[seq_len(n_founders)]  <- NA_character_
for (i in (n_founders + 1L):N_IND) {
  sire_col[i] <- IND_NAMES[sample(seq_len(i - 1L), 1L)]
  dam_col[i]  <- IND_NAMES[sample(seq_len(i - 1L), 1L)]
}

ped_df <- data.frame(
  id         = IND_NAMES,
  sire       = sire_col,
  dam        = dam_col,
  generation = generation,
  stringsAsFactors = FALSE)
write.csv(ped_df, file.path(OUT_DIR, "pedigree.csv"), row.names = FALSE)

# ── 11.  Summary ─────────────────────────────────────────────────────────────
cat(sprintf(
  "\n=== generate.R [%s] complete ===\n  Individuals : %d\n  Chromosomes : %d\n  SNPs/chr    : %d\n  Total SNPs  : %d\n  MH blocks   : %d\n  Haplo blocks: %d\n  n_qtl       : %d\n  n_founders  : %d\n  Output dir  : %s\n  Files written:\n    genotype.csv   (%d rows x %d cols)\n    microhap.csv   (%d rows x %d cols)\n    phenotype.csv  (%d rows x 5 cols)\n    pedigree.csv   (%d rows x 4 cols)\n",
  size, N_IND, N_CHR, N_SNP_PER_CHR,
  n_snp_total,
  length(mh_blocks), length(haplo_blocks),
  n_qtl, n_founders,
  OUT_DIR,
  nrow(geno_df),  ncol(geno_df),
  nrow(mh_df),    ncol(mh_df),
  nrow(pheno_df),
  nrow(ped_df)))
