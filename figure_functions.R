# shared figure settings, labels, and save functions

# inputs:
# none

# outputs:
# figure settings and label functions loaded by config.R

# author: Sidney Bush
# date: 2026-02-13

# shared font sizes and export settings
FIG_BASE_SIZE <- 18
FIG_AXIS_TEXT_SIZE <- 16
FIG_AXIS_TITLE_SIZE <- 18
FIG_STRIP_TEXT_SIZE <- 16
FIG_ANNOT_TEXT_SIZE <- 5
FIG_TILE_TEXT_SIZE <- 6
FIG_POINT_SIZE_SMALL <- 1.5
FIG_POINT_SIZE_LARGE <- 3.0
FIG_WIDTH_SCALE <- 1.35
FIG_HEIGHT_SCALE <- 1.35
FIG_PREVIEW_DPI <- 300
FIG_PRODUCTION_DPI <- 300

# shared label spacing settings
FIG_LABEL_CLIP <- "off"
FIG_LABEL_PLOT_MARGIN_PT <- 18

# site colors used in figures
SITE_COLORS <- c(
  "WS09" = "#882255",
  "WS10" = "#AA4499",
  "WS01" = "#CC6677",
  "Look" = "#DDCC77",
  "WS02" = "#999933",
  "WS03" = "#117733",
  "WS06" = "#44AA99",
  "WS07" = "#88CCEE",
  "WS08" = "#6699CC",
  "Mack" = "#332288"
)

# plot theme
# individual figure scripts add any figure specific styling
theme_pub <- function(base_size = FIG_BASE_SIZE) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(hjust = 0),
      strip.text.y = ggplot2::element_text(hjust = 0)
    )
}

# save one figure as PNG, PDF, and TIFF
save_figure_set <- function(
  plot_obj,
  file_base,
  width,
  height,
  png_dir = MS_FIG_MAIN_DIR,
  pdf_dir = MS_FIG_MAIN_PDF_DIR,
  tiff_dir = MS_FIG_MAIN_TIFF_DIR,
  bg = "white",
  preview_dpi = FIG_PREVIEW_DPI,
  production_dpi = FIG_PRODUCTION_DPI
) {
  for (d in c(png_dir, pdf_dir, tiff_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  output_files <- c(
    png = file.path(png_dir, paste0(file_base, ".png")),
    pdf = file.path(pdf_dir, paste0(file_base, ".pdf")),
    tiff = file.path(tiff_dir, paste0(file_base, ".tiff"))
  )

  ggplot2::ggsave(
    output_files[["png"]],
    plot_obj,
    width = width,
    height = height,
    bg = bg,
    dpi = preview_dpi
  )
  ggplot2::ggsave(
    output_files[["pdf"]],
    plot_obj,
    width = width,
    height = height,
    bg = bg
  )
  ggplot2::ggsave(
    output_files[["tiff"]],
    plot_obj,
    width = width,
    height = height,
    bg = bg,
    dpi = production_dpi,
    compression = "lzw"
  )

  invisible(output_files)
}
