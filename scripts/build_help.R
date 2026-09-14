# Run from the package root after editing roxygen comments or the vignette.
# Generates the browsable manual and executes the tutorial using local code.
stopifnot(file.exists("DESCRIPTION"), file.exists("vignettes/alrd-introduction.Rmd"))
roxygen2::roxygenise()
pkgload::load_all(quiet = TRUE)
destination <- file.path(getwd(), "output", "help")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
# Use R's bundled math renderer so local help also works without internet.
katex_source <- file.path(R.home("doc"), "html", "katex")
stopifnot(dir.exists(katex_source))
file.copy(katex_source, destination, recursive = TRUE, overwrite = TRUE)
topics <- c("ALRD-package", "fit_alrd", "fit_lrd", "estimate_ps",
            "working_covariance", "predict_ite", "simulate_alrd_data", "print.alrd_fit")
links <- setNames(paste0(topics, ".html"), topics)
links <- c(links, ALRD = "ALRD-package.html", print.lrd_fit = "print.alrd_fit.html")
for (topic in topics) {
  tools::Rd2HTML(file.path("man", paste0(topic, ".Rd")),
                 out = file.path(destination, paste0(topic, ".html")),
                 package = "ALRD", Links = links, stylesheet = "R.css",
                 texmath = "katex")
  page <- file.path(destination, paste0(topic, ".html"))
  html <- readLines(page, warn = FALSE)
  html <- gsub("https://cdn[.]jsdelivr[.]net/npm/katex@[^/]+/dist/katex[.]min[.]css",
               "katex/katex.css", html)
  html <- gsub("https://cdn[.]jsdelivr[.]net/npm/katex@[^/]+/dist/katex[.]min[.]js",
               "katex/katex.js", html)
  writeLines(html, page)
}
file.copy(file.path(R.home("doc"), "html", "R.css"), destination, overwrite = TRUE)
rmarkdown::render("vignettes/alrd-introduction.Rmd",
                   output_file = "alrd-user-guide.html", output_dir = destination,
                   envir = new.env(parent = globalenv()), quiet = TRUE)
# Render Markdown companions with native MathML, including offline formulas.
companions <- c("README.md", "docs/package_usage_guide.md",
                "docs/weekly_method_check.md", "docs/code_review_and_collaborator_tasks.md")
for (source in companions) {
  page <- file.path(destination, paste0(tools::file_path_sans_ext(basename(source)), ".html"))
  rmarkdown::pandoc_convert(normalizePath(source), from = "markdown+tex_math_dollars", to = "html5",
    output = page,
    options = c("--standalone", "--mathml", "--css=R.css",
                paste0("--metadata=title:", tools::file_path_sans_ext(basename(source)))))
  # Keep project-document links working from the generated help directory.
  html <- readLines(page, warn = FALSE)
  matches <- gregexpr('href="[^"#][^"]*"', html)
  regmatches(html, matches) <- lapply(regmatches(html, matches), function(attrs) {
    vapply(attrs, function(attr) {
      target <- sub('^href="(.*)"$', '\\1', attr)
      if (grepl("^([[:alpha:]][[:alnum:]+.-]*:|/)", target) || target == "R.css") return(attr)
      resolved <- file.path(dirname(source), target)
      if (sub("^\\./", "", resolved) %in% companions) {
        target <- paste0(tools::file_path_sans_ext(basename(target)), ".html")
      } else {
        target <- paste0("../../", resolved)
      }
      paste0('href="', target, '"')
    }, character(1), USE.NAMES = FALSE)
  })
  writeLines(html, page)
}
labels <- c("Package overview", "fit_alrd(): Augmented learning",
            "fit_lrd(): Direct learning", "estimate_ps(): Propensity scores",
            "working_covariance(): Working correlation matrices",
            "predict_ite(): Treatment-effect prediction",
            "simulate_alrd_data(): Example data", "print(): Fitted-object display")
items <- paste0('<li><a href="', topics, '.html">', labels, '</a></li>')
writeLines(c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
             '<title>ALRD documentation</title><link rel="stylesheet" href="R.css">',
             '</head><body><main><h1>ALRD documentation</h1>',
             '<p>Development version 0.0.0.9000. Continuous repeated outcomes with binary baseline treatment.</p>',
             '<p><a href="alrd-user-guide.html">Start here: executable user guide</a></p>',
             '<p>The guide covers installation, data requirements, randomized and observational workflows,',
             'lasso selection, prediction, output interpretation and troubleshooting.</p>',
             '<h2>Function reference</h2><ul>', items, '</ul>',
             '<h2>Project documents</h2><ul>',
             '<li><a href="README.html">Package README</a></li>',
             '<li><a href="package_usage_guide.html">Companion usage guide</a></li>',
             '<li><a href="weekly_method_check.html">Weekly report</a></li>',
             '<li><a href="code_review_and_collaborator_tasks.html">Code review and collaborator tasks</a></li></ul>',
             '<p>Read the validation limits in each fitting and prediction topic before substantive use.</p>',
             '</main></body></html>'), file.path(destination, "index.html"))
file.copy(file.path(destination, "index.html"),
           file.path(destination, "00Index.html"), overwrite = TRUE)
cat("Documentation generated at", file.path(destination, "index.html"), "\n")
