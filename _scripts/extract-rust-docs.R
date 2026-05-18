#!/usr/bin/env Rscript
# Extract `//!` module docstrings from Rust source files in the sibling
# masbayes/ and masreml/ packages, and emit one .qmd snippet per .rs file
# under internals/_extracted/<package>/. The Quarto pages
# internals/masbayes-rust.qmd and internals/masreml-rust.qmd pull these
# snippets in via {{< include >}} directives.
#
# The script is read-only with respect to masbayes/ and masreml/: it never
# writes back to those source trees, only reads .rs files.
#
# Usage (from masgenomics-docs/ root):
#   Rscript _scripts/extract-rust-docs.R
#
# Or via the Quarto pre-render hook in _quarto.yml.

# ---- Path resolution -------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg))))
  }
  normalizePath(".")
}

script_dir <- get_script_dir()
repo_root  <- normalizePath(file.path(script_dir, ".."))
parent_dir <- normalizePath(file.path(repo_root, ".."))

packages <- list(
  masbayes = file.path(parent_dir, "masbayes", "src", "rust", "src"),
  masreml  = file.path(parent_dir, "masreml",  "src", "rust", "src")
)

out_root <- file.path(repo_root, "internals", "_extracted")

# ---- Helpers ---------------------------------------------------------------

# Find consecutive `//!` lines at the top of `rs_path` and return their
# stripped contents. Returns character(0) if no module docstring is found.
extract_module_doc <- function(rs_path) {
  lines <- tryCatch(
    readLines(rs_path, warn = FALSE),
    error = function(e) character()
  )
  if (length(lines) == 0) {
    return(character())
  }

  trimmed <- trimws(lines)
  is_doc <- startsWith(trimmed, "//!")
  if (!any(is_doc)) {
    return(character())
  }

  # First `//!` line in the file (skipping any leading regular comments,
  # blank lines, `use` declarations, attributes, etc.).
  start <- which(is_doc)[1]

  # Walk forward, collecting consecutive `//!` lines.
  end <- start
  for (i in seq(start, length(lines))) {
    if (startsWith(trimws(lines[i]), "//!")) {
      end <- i
    } else {
      break
    }
  }

  doc <- lines[start:end]
  # Strip the `//! ` (or `//!`) prefix, preserving any further indent.
  vapply(
    doc,
    function(l) sub("^\\s*//!\\s?", "", l, perl = TRUE),
    character(1),
    USE.NAMES = FALSE
  )
}

# List relative `.rs` paths under `rust_src_dir`, sorted.
list_rust_files <- function(rust_src_dir) {
  files <- list.files(
    rust_src_dir, pattern = "\\.rs$",
    recursive = TRUE, full.names = FALSE
  )
  sort(files)
}

# "gwas/emmax.rs" -> "gwas-emmax"   (used as snippet filename + heading id)
snippet_key <- function(rel_path) {
  noext <- sub("\\.rs$", "", rel_path)
  gsub("/", "-", noext)
}

# Adjust markdown headings inside an extracted docstring so they nest
# correctly under the snippet's `### file.rs` header (H3) in the parent
# Quarto page:
#   1. Drop the first H1 (`# Title`) if present — it duplicates the snippet
#      header that we prepend separately.
#   2. Shift remaining headings down by `shift` levels (## -> #####, etc.),
#      preserving relative hierarchy.
shift_headings <- function(lines, shift = 2L, drop_first_h1 = TRUE) {
  if (length(lines) == 0) return(lines)

  if (drop_first_h1) {
    nonblank <- which(nzchar(trimws(lines)))
    if (length(nonblank) >= 1) {
      first_idx <- nonblank[1]
      # Match `# Foo` but not `## Foo`
      if (grepl("^#\\s+[^#]", lines[first_idx])) {
        lines <- lines[-first_idx]
      }
    }
  }

  if (shift > 0) {
    pad <- strrep("#", shift)
    lines <- vapply(
      lines,
      function(l) if (grepl("^#+\\s", l)) paste0(pad, l) else l,
      character(1),
      USE.NAMES = FALSE
    )
  }

  lines
}

# Render a single Quarto snippet for one Rust file.
write_snippet <- function(snippet_path, pkg, rel_path, doc_lines) {
  key <- snippet_key(rel_path)
  header <- sprintf("### `%s` {#sec-%s-%s}", rel_path, pkg, key)
  link   <- sprintf("{{< rust-src %s %s >}}", pkg, rel_path)

  body <- if (length(doc_lines) == 0) {
    "*(no module-level documentation)*"
  } else {
    shifted <- shift_headings(doc_lines)
    paste(shifted, collapse = "\n")
  }

  content <- paste(c(header, "", link, "", body, ""), collapse = "\n")
  writeLines(content, snippet_path)
}

# ---- Main ------------------------------------------------------------------

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

total <- 0L
for (pkg in names(packages)) {
  src_dir <- packages[[pkg]]
  if (!dir.exists(src_dir)) {
    warning(sprintf(
      "Skipping %s: source directory not found at %s",
      pkg, src_dir
    ))
    next
  }

  pkg_out_dir <- file.path(out_root, pkg)
  dir.create(pkg_out_dir, showWarnings = FALSE, recursive = TRUE)

  rs_files <- list_rust_files(src_dir)
  for (rel_path in rs_files) {
    rs_full   <- file.path(src_dir, rel_path)
    doc_lines <- extract_module_doc(rs_full)
    snippet   <- file.path(pkg_out_dir, paste0(snippet_key(rel_path), ".qmd"))
    write_snippet(snippet, pkg, rel_path, doc_lines)
  }

  message(sprintf(
    "  [%s] %d file(s) -> %s",
    pkg, length(rs_files), pkg_out_dir
  ))
  total <- total + length(rs_files)
}

message(sprintf("Extracted %d Rust source snippet(s) total.", total))
