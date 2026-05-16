#!/usr/bin/env Rscript
# Regenerate Quarto reference pages from the roxygen man pages of
# masreml and masbayes via altdoc. Runs on a temp copy of each package
# so the package repos are NEVER touched.
#
# EMPIRICAL FINDINGS (altdoc 0.7.2):
#   - setup_docs() requires the working dir to be set to the package root (use setwd).
#   - render_docs() produces .qmd files in <pkg>/_quarto/man/ (NOT docs/reference/).
#   - The .qmd files have NO YAML frontmatter; they start with ## headers.
#   - We add a minimal YAML frontmatter when copying so they render as standalone docs.

suppressPackageStartupMessages({
  library(altdoc)
})

# Get the directory this script lives in (works under Rscript)
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
  masreml  = file.path(parent_dir, "masreml"),
  masbayes = file.path(parent_dir, "masbayes")
)

excludes <- c(".git", "target", "src/rust/target", ".DS_Store",
              "_freeze", "_site", ".quarto", ".Rproj.user", "docs",
              "_quarto", "altdoc")

tmp_root <- tempfile("masgenomics-altdoc-")
dir.create(tmp_root)
on.exit(unlink(tmp_root, recursive = TRUE), add = TRUE)

# Helper: add YAML frontmatter to a .qmd fragment from altdoc.
# altdoc produces headerless markdown starting with ## Title {.unnumbered}.
# We extract the title from that first line and add a proper YAML block.
# Also normalises any `eval=TRUE` chunk option to `eval=FALSE` — CI does not
# install masbayes/masreml, so reference example chunks must never execute.
add_frontmatter <- function(src_file, dst_file) {
  lines <- readLines(src_file, warn = FALSE)
  # Find the title: first line matching ^## ... {.unnumbered}
  title_line <- grep("^## .+ \\{\\.unnumbered\\}", lines, value = TRUE)[1]
  title <- if (!is.na(title_line)) {
    # Strip the ## prefix and {.unnumbered} suffix
    sub("^## (.+) \\{\\.unnumbered\\}$", "\\1", title_line)
  } else {
    tools::file_path_sans_ext(basename(src_file))
  }
  # Escape any double-quotes in title for YAML
  title_yaml <- gsub('"', '\\\\"', title)
  # Force eval=FALSE on every chunk header — CI lacks the sibling packages,
  # so altdoc's default eval=TRUE (from unwrapped @examples) would break the
  # build. Matches both `eval=TRUE` and `eval = TRUE` (any whitespace).
  is_chunk_header <- grepl("^```\\{[rR][ ,}]", lines)
  lines[is_chunk_header] <- gsub("eval\\s*=\\s*TRUE",
                                 "eval=FALSE",
                                 lines[is_chunk_header])
  frontmatter <- c(
    "---",
    sprintf('title: "%s"', title_yaml),
    "---",
    ""
  )
  writeLines(c(frontmatter, lines), dst_file)
}

# Build a reference index.qmd for a package listing all exported functions
make_index_qmd <- function(pkg_name, qmd_files, out_path) {
  fn_names <- tools::file_path_sans_ext(basename(qmd_files))
  fn_names <- fn_names[!grepl(paste0("^", pkg_name, "-package$"), fn_names, ignore.case = TRUE)]
  links <- sprintf("- [%s](%s.qmd)", fn_names, fn_names)
  content <- c(
    "---",
    sprintf('title: "%s Reference"', pkg_name),
    "---",
    "",
    sprintf("Reference documentation for the `%s` package.", pkg_name),
    "",
    "## Functions",
    "",
    links
  )
  writeLines(content, out_path)
}

for (pkg_name in names(packages)) {
  pkg_path <- packages[[pkg_name]]
  if (!dir.exists(pkg_path)) stop(sprintf("Package source not found: %s", pkg_path))

  tmp_pkg <- file.path(tmp_root, pkg_name)
  dir.create(tmp_pkg)

  message(sprintf("\n[render-reference] ── %s ──────────────────────────────", pkg_name))
  message(sprintf("[render-reference] copy %s -> %s", pkg_path, tmp_pkg))

  exclude_args <- paste(paste0("--exclude=", shQuote(excludes)), collapse = " ")
  cmd <- sprintf("rsync -a %s %s/ %s/", exclude_args, shQuote(pkg_path), shQuote(tmp_pkg))
  ret <- system(cmd)
  if (ret != 0) stop("rsync failed for ", pkg_name)
  message(sprintf("[render-reference] rsync complete. Contents: %s",
                  paste(list.files(tmp_pkg), collapse = ", ")))

  # altdoc requires working directory to be the package root
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_pkg)

  message(sprintf("[render-reference] altdoc::setup_docs(%s)", tmp_pkg))
  altdoc::setup_docs(tool = "quarto_website", overwrite = TRUE)

  message(sprintf("[render-reference] altdoc::render_docs(%s)", tmp_pkg))
  altdoc::render_docs(parallel = FALSE, verbose = FALSE)

  setwd(old_wd)

  # Locate the .qmd source files. altdoc 0.7.2 places them in _quarto/man/.
  candidates <- c(
    file.path(tmp_pkg, "_quarto", "man"),
    file.path(tmp_pkg, "docs", "reference"),
    file.path(tmp_pkg, "docs", "man")
  )
  src <- candidates[dir.exists(candidates)][1]
  if (is.na(src)) {
    message("DEBUG: full tree of ", tmp_pkg)
    print(list.files(tmp_pkg, recursive = TRUE, all.files = TRUE))
    stop(sprintf("Could not find altdoc .qmd output for %s", pkg_name))
  }
  message(sprintf("[render-reference] found qmd source at: %s", src))

  qmd_files <- list.files(src, pattern = "\\.qmd$", full.names = TRUE)
  message(sprintf("[render-reference] found %d .qmd files", length(qmd_files)))

  out_dir <- file.path(repo_root, "reference", pkg_name)
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)

  # Copy each .qmd with added frontmatter
  for (f in qmd_files) {
    dst <- file.path(out_dir, basename(f))
    add_frontmatter(f, dst)
  }

  # Generate an index.qmd listing all functions
  index_path <- file.path(out_dir, "index.qmd")
  make_index_qmd(pkg_name, qmd_files, index_path)

  n_copied <- length(list.files(out_dir, pattern = "\\.qmd$"))
  message(sprintf("[render-reference] wrote %d .qmd files (incl index) to %s",
                  n_copied, out_dir))
}

message("\n[render-reference] done.")
