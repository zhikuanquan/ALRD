# Run from the package root after editing roxygen comments or the vignette.
# Generates the browsable manual and executes the tutorial using local code.
stopifnot(file.exists("DESCRIPTION"), file.exists("vignettes/alrd-introduction.Rmd"))
roxygen2::roxygenise()
pkgload::load_all(quiet = TRUE)
destination <- file.path(getwd(), "output", "help")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
topics <- c("ALRD-package", "fit_alrd", "fit_lrd", "estimate_ps",
            "working_covariance", "predict_ite", "simulate_alrd_data", "print.alrd_fit")
links <- setNames(paste0(topics, ".html"), topics)
links <- c(links, ALRD = "ALRD-package.html", print.lrd_fit = "print.alrd_fit.html")
for (topic in topics) {
  tools::Rd2HTML(file.path("man", paste0(topic, ".Rd")),
                 out = file.path(destination, paste0(topic, ".html")),
                 package = "ALRD", Links = links, stylesheet = "R.css")
}
file.copy(file.path(R.home("doc"), "html", "R.css"), destination, overwrite = TRUE)
rmarkdown::render("vignettes/alrd-introduction.Rmd",
                   output_file = "alrd-user-guide.html", output_dir = destination,
                   envir = new.env(parent = globalenv()), quiet = TRUE)
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
             '<p>Read the validation limits in each fitting and prediction topic before substantive use.</p>',
             '</main></body></html>'), file.path(destination, "index.html"))
file.copy(file.path(destination, "index.html"),
           file.path(destination, "00Index.html"), overwrite = TRUE)
cat("Documentation generated at", file.path(destination, "index.html"), "\n")
