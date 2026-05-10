#!/usr/bin/env Rscript
# Demo dataset generator for masgenomics-docs.
# Mirrors the file layout consumed by MicrohapsSel/simulation/model_fitting.
# Pre-generated offline; CSVs are committed as static assets so every visitor
# sees identical numbers.

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
size <- if (length(args) >= 1) args[1] else "main"
stopifnot(size %in% c("main", "toy"))

CONFIGS <- list(
  main = list(N_IND = 1000, N_CHR = 10, N_SNP_PER_CHR = 500, N_QTL = 50),
  toy  = list(N_IND = 200,  N_CHR = 5,  N_SNP_PER_CHR = 100, N_QTL = 20)
)
cfg <- CONFIGS[[size]]

set.seed(42)

N_IND          <- cfg$N_IND
N_CHR          <- cfg$N_CHR
N_SNP_PER_CHR  <- cfg$N_SNP_PER_CHR
N_QTL          <- cfg$N_QTL
IND_NAMES      <- sprintf("IND%05d", 1:N_IND)
BASE           <- file.path("demo-data", size)

dirs <- c(
  BASE,
  file.path(BASE, "mh_genotypes"),
  file.path(BASE, "all_genotypes"),
  file.path(BASE, "mh_info_ld_haploblock_G0/stats"),
  file.path(BASE, "all_info_haploblock_G0/stats"),
  file.path(BASE, "qtl")
)
walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ── 1. SNP positions and haplotypes per chromosome (same as dummy_data.R) ──
allele_pairs <- list(c("A","T"), c("G","C"), c("T","A"), c("C","G"))
chr_snps  <- map(1:N_CHR, ~sort(sample(1000:5000000, N_SNP_PER_CHR)))
chr_haplo <- map(1:N_CHR, ~matrix(sample(0:1, N_IND * 2 * N_SNP_PER_CHR, replace = TRUE),
                                  nrow = N_IND * 2, ncol = N_SNP_PER_CHR))

# ── 2. VCF ──
vcf_header <- c(
  "##fileformat=VCFv4.1",
  "##source=masgenomics-docs-demo",
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", IND_NAMES),
        collapse = "\t")
)
vcf_rows <- do.call(rbind, lapply(1:N_CHR, function(c) {
  do.call(rbind, lapply(1:N_SNP_PER_CHR, function(si) {
    pos <- chr_snps[[c]][si]
    ap  <- allele_pairs[[(si - 1) %% 4 + 1]]
    gts <- sapply(1:N_IND, function(ind)
      paste(chr_haplo[[c]][(ind - 1) * 2 + 1, si],
            chr_haplo[[c]][(ind - 1) * 2 + 2, si], sep = "|"))
    c(as.character(c), pos, sprintf("ssa%02d_%d_%s_%s", c, pos, ap[1], ap[2]),
      ap[1], ap[2], ".", "PASS", ".", "GT", gts)
  }))
}))
vcf_data_lines <- apply(vcf_rows, 1, paste, collapse = "\t")
writeLines(c(vcf_header, vcf_data_lines), file.path(BASE, "salmon_60k.vcf"))

# ── 3. SNP IDs and per-individual dosage matrix (for geno_G1.csv) ──
snp_id <- character(N_CHR * N_SNP_PER_CHR)
geno_dosage <- matrix(0L, nrow = N_IND, ncol = N_CHR * N_SNP_PER_CHR)
col_offset <- 0L
for (c in 1:N_CHR) {
  for (si in 1:N_SNP_PER_CHR) {
    pos <- chr_snps[[c]][si]
    ap  <- allele_pairs[[(si - 1) %% 4 + 1]]
    snp_id[col_offset + si] <- sprintf("ssa%02d_%d_%s_%s", c, pos, ap[1], ap[2])
    geno_dosage[, col_offset + si] <-
      chr_haplo[[c]][seq(1, 2 * N_IND, by = 2), si] +
      chr_haplo[[c]][seq(2, 2 * N_IND, by = 2), si]
  }
  col_offset <- col_offset + N_SNP_PER_CHR
}
colnames(geno_dosage) <- snp_id

# ── 4. geno_G1.csv ──
geno_df <- tibble(
  individual_id   = IND_NAMES,
  generation      = 1L,
  population_type = "reference"
) %>%
  bind_cols(as_tibble(geno_dosage))
write_csv(geno_df, file.path(BASE, "geno_G1.csv"))

# ── 5. MH blocks and haplotype blocks (same logic as dummy_data.R) ──
make_mh_blocks <- function(chr_snps, max_span = 125) {
  global_blk <- 1; all_blocks <- list()
  for (c in seq_along(chr_snps)) {
    pos <- chr_snps[[c]]; i <- 1
    while (i <= length(pos)) {
      j <- i
      while (j < length(pos) && (pos[j + 1] - pos[i]) <= max_span && (j - i) < 3) j <- j + 1
      snp_idx <- i:j
      all_blocks[[length(all_blocks) + 1]] <- list(blk = global_blk, chr = c,
                                                   snp_idx = snp_idx, pos = pos[snp_idx])
      global_blk <- global_blk + 1
      i <- j + 1
    }
  }
  all_blocks
}
make_haplo_blocks <- function(chr_snps, window = 100000) {
  global_blk <- 1; all_blocks <- list()
  for (c in seq_along(chr_snps)) {
    pos <- chr_snps[[c]]
    win_start <- pos[1]
    repeat {
      win_end <- win_start + window
      snp_idx <- which(pos >= win_start & pos < win_end)
      if (length(snp_idx) > 0) {
        all_blocks[[length(all_blocks) + 1]] <- list(blk = global_blk, chr = c,
                                                     snp_idx = snp_idx, pos = pos[snp_idx])
        global_blk <- global_blk + 1
      }
      win_start <- win_end
      if (win_start > pos[length(pos)]) break
    }
  }
  all_blocks
}
mh_blocks    <- make_mh_blocks(chr_snps)
haplo_blocks <- make_haplo_blocks(chr_snps)

# ── 6. Block files (4-col-per-block format expected by load_hap_geno) ──
write_block_files <- function(blocks, info_dir, geno_dir) {
  by_chr <- split(blocks, map_int(blocks, ~ .x$chr))
  walk(names(by_chr), function(c) {
    blk_list <- by_chr[[c]]
    ci <- as.integer(c)

    writeLines(map_chr(blk_list,
                       ~ paste0("blk", .x$blk, "\t ", paste(.x$snp_idx - 1, collapse = " "))),
               file.path(info_dir, paste0("hap_block_", c)))
    writeLines(map_chr(blk_list,
                       ~ paste0("blk", .x$blk, "\t", paste(.x$snp_idx - 1, collapse = "\t"))),
               file.path(info_dir, paste0("hap_block_info_", c)))

    hdr <- c("ID", unlist(map(seq_along(blk_list), function(local_idx) c(
      paste0("hap_", local_idx, "_1"), paste0("hap_", local_idx, "_1"),
      paste0("hap_", local_idx, "_2"), paste0("hap_", local_idx, "_2")))))

    rows <- map_chr(1:N_IND, function(i) {
      vals <- unlist(map(blk_list, function(b) {
        h1 <- chr_haplo[[ci]][(i - 1) * 2 + 1, b$snp_idx]
        h2 <- chr_haplo[[ci]][(i - 1) * 2 + 2, b$snp_idx]
        a1 <- strtoi(paste(h1, collapse = ""), base = 2L) + 1
        a2 <- strtoi(paste(h2, collapse = ""), base = 2L) + 1
        c(a1, a1, a2, a2)
      }))
      paste(c(IND_NAMES[i], vals), collapse = "\t")
    })
    writeLines(c(paste(hdr, collapse = "\t"), rows),
               file.path(geno_dir, paste0("hap_geno_", c)))
  })
}
write_block_files(mh_blocks,    file.path(BASE, "mh_info_ld_haploblock_G0"),  file.path(BASE, "mh_genotypes"))
write_block_files(haplo_blocks, file.path(BASE, "all_info_haploblock_G0"),    file.path(BASE, "all_genotypes"))

# ── 7. Stats CSVs (coordinates + snp selection detail) ──
write_stats <- function(blocks, stats_dir) {
  coord <- map_dfr(blocks, ~ tibble(
    block_id         = paste0("blk", .x$blk),
    chr              = .x$chr,
    start_pos        = .x$pos[1],
    end_pos          = .x$pos[length(.x$pos)],
    n_snps           = length(.x$snp_idx),
    physical_span_bp = .x$pos[length(.x$pos)] - .x$pos[1],
    pic              = round(runif(1, 0.05, 0.5), 4),
    NA_col           = NA
  ))
  snp_det <- map_dfr(blocks, function(b) {
    map_dfr(seq_along(b$snp_idx), function(k) {
      ap <- allele_pairs[[(b$snp_idx[k] - 1) %% 4 + 1]]
      tibble(
        block_id  = paste0("blk", b$blk),
        chr       = b$chr,
        snp_index = b$snp_idx[k] - 1,
        snp_id    = sprintf("ssa%02d_%d_%s_%s", b$chr, b$pos[k], ap[1], ap[2]),
        position  = b$pos[k]
      )
    })
  })
  write_csv(coord,   file.path(stats_dir, "microhaplotype_coordinates.csv"))
  write_csv(snp_det, file.path(stats_dir, "snp_selection_detailed.csv"))
}
write_stats(mh_blocks,    file.path(BASE, "mh_info_ld_haploblock_G0/stats"))
write_stats(haplo_blocks, file.path(BASE, "all_info_haploblock_G0/stats"))

# ── 8. Pedigree ──
n_founders   <- max(50, round(0.1 * N_IND))
generation   <- pmin(4, ceiling((1:N_IND) / (N_IND / 4)))
sire         <- character(N_IND)
dam          <- character(N_IND)
family_id    <- character(N_IND)

# Founders: NA parents, each its own family
sire[1:n_founders]      <- NA_character_
dam[1:n_founders]       <- NA_character_
family_id[1:n_founders] <- sprintf("F%03d", 1:n_founders)

for (i in seq.int(n_founders + 1, N_IND)) {
  parent_pool <- 1:(i - 1)
  s_idx <- sample(parent_pool, 1)
  d_idx <- sample(parent_pool, 1)
  sire[i]      <- IND_NAMES[s_idx]
  dam[i]       <- IND_NAMES[d_idx]
  family_id[i] <- family_id[s_idx]                     # inherit from sire
}

ped_df <- tibble(
  individual_id   = IND_NAMES,
  family_id       = family_id,
  population_type = "reference",
  sire            = sire,
  dam             = dam,
  generation      = generation
)
write_csv(ped_df, file.path(BASE, "pedigree_G1.csv"))

# ── 9. QTL info ──
qtl_idx        <- sample.int(ncol(geno_dosage), N_QTL)
qtl_eff        <- rnorm(N_QTL, 0, 1)
qtl_snp_ids    <- snp_id[qtl_idx]
qtl_mh_blocks  <- sample(map_chr(mh_blocks,    ~ paste0("blk", .x$blk)),
                         min(20, length(mh_blocks)))
qtl_hap_blocks <- sample(map_chr(haplo_blocks, ~ paste0("blk", .x$blk)),
                         min(10, length(haplo_blocks)))

snp_meta <- tibble(
  snp_id = snp_id,
  chr    = rep(1:N_CHR, each = N_SNP_PER_CHR),
  pos    = unlist(chr_snps)
)

saveRDS(qtl_snp_ids,    file.path(BASE, "qtl/qtl_snp_ids.rds"))
saveRDS(qtl_mh_blocks,  file.path(BASE, "qtl/qtl_mh_blocks.rds"))
saveRDS(qtl_hap_blocks, file.path(BASE, "qtl/qtl_hap_blocks.rds"))
saveRDS(snp_meta,       file.path(BASE, "qtl/snp_meta.rds"))

# ── 10. SNP map ──
snp_map <- snp_meta %>% rename(position = pos)
write.table(snp_map, file.path(BASE, "snp_map.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ── 11. Phenotypes (continuous + binary) ──
geno_centered <- scale(geno_dosage[, qtl_idx], center = TRUE, scale = FALSE)
tbv           <- as.numeric(geno_centered %*% qtl_eff)
tbv           <- as.numeric(scale(tbv)[, 1])    # standardise to var=1
h2            <- 0.4
sigma_e       <- sqrt((1 - h2) / h2)
y_cont        <- tbv + rnorm(N_IND, 0, sigma_e)

write_csv(tibble(
  individual_id = IND_NAMES,
  phenotype     = round(y_cont, 4),
  tbv           = round(tbv, 4)
), file.path(BASE, "pheno_continuous_G1.csv"))

y_bin <- as.integer(y_cont > median(y_cont))
write_csv(tibble(
  individual_id = IND_NAMES,
  phenotype     = y_bin,
  tbv           = round(tbv, 4)
), file.path(BASE, "pheno_binary_G1.csv"))

# ── 12. Summary ──
cat(sprintf("\n[generate.R] %s done. Files in %s:\n", size, BASE))
cat(sprintf("  geno_G1.csv             : %d ind x %d SNPs\n",     N_IND, length(snp_id)))
cat(sprintf("  mh_genotypes/           : %d MH blocks across %d chrs\n", length(mh_blocks),    N_CHR))
cat(sprintf("  all_genotypes/          : %d hap blocks across %d chrs\n", length(haplo_blocks), N_CHR))
cat(sprintf("  pedigree_G1.csv         : %d individuals (%d founders)\n", N_IND, n_founders))
cat(sprintf("  pheno_continuous_G1.csv : h2=%.2f\n", h2))
cat(sprintf("  pheno_binary_G1.csv     : threshold @ median (prevalence ~50%%)\n"))
cat(sprintf("  qtl/                    : %d SNP QTL, %d MH QTL, %d hap QTL\n",
            length(qtl_snp_ids), length(qtl_mh_blocks), length(qtl_hap_blocks)))
cat(sprintf("  snp_map.txt             : %d rows\n", nrow(snp_map)))
cat("\n")
