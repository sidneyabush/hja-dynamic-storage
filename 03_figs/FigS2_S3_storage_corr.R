# Supporting information storage metric correlation figures
# inputs: outputs/master/*.csv
# outputs: ms_materials/supp/FigS2_dynamic_storage_corr.* and FigS3_mobile_storage_corr.*

librarian::shelf(dplyr, readr, tidyr, ggplot2, scales, cran_repo = "https://cloud.r-project.org")

rm(list = ls())
source("config.R")

supp_dir <- MS_FIG_SUPP_DIR
supp_pdf_dir <- MS_FIG_SUPP_PDF_DIR
supp_tiff_dir <- MS_FIG_SUPP_TIFF_DIR
for (d in c(supp_dir, supp_pdf_dir, supp_tiff_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

annual_file <- file.path(OUTPUT_DIR, "master", MASTER_ANNUAL_FILE)
site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
fdc_wy_file <- file.path(OUT_MET_DYNAMIC_DIR, "fdc_slopes_wy.csv")

for (f in c(annual_file, site_file)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}

annual <- read_csv(annual_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

if (file.exists(fdc_wy_file)) {
  fdc_wy <- read_csv(fdc_wy_file, show_col_types = FALSE) %>%
    mutate(
      site = standardize_site_code(as.character(site)),
      year = suppressWarnings(as.numeric(
        if ("WaterYear" %in% names(.)) WaterYear else if ("year" %in% names(.)) year else NA_real_
      )),
      FDC_annual = suppressWarnings(as.numeric(
        if ("Slope" %in% names(.)) Slope else if ("fdc_slope" %in% names(.)) fdc_slope else if ("FDC" %in% names(.)) FDC else NA_real_
      ))
    ) %>%
    filter(site %in% SITE_ORDER_HYDROMETRIC, is.finite(year)) %>%
    select(site, year, FDC_annual) %>%
    distinct(site, year, .keep_all = TRUE)

  annual <- annual %>%
    mutate(site_chr = as.character(site)) %>%
    left_join(fdc_wy, by = c("site_chr" = "site", "year" = "year")) %>%
    mutate(
      FDC = dplyr::coalesce(FDC_annual, FDC),
      site = factor(site_chr, levels = SITE_ORDER_HYDROMETRIC)
    ) %>%
    select(-site_chr, -FDC_annual)
}

site_summary <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC) %>%
  mutate(site = factor(site, levels = SITE_ORDER_HYDROMETRIC))

make_corr_matrix <- function(data_df, metric_map, metric_labels, parse_metric_labels = FALSE) {
  corr_input <- data_df %>%
    select(any_of(unname(metric_map))) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

  if (ncol(corr_input) < 2) stop("At least two metrics are needed for a correlation plot.")

  colnames(corr_input) <- names(metric_map)
  corr_mat <- suppressWarnings(cor(corr_input, use = "pairwise.complete.obs"))

  idx <- which(lower.tri(corr_mat), arr.ind = TRUE)
  if (nrow(idx) == 0) stop("At least two metrics are needed for a correlation plot.")

  x_levels <- names(metric_map)[-length(metric_map)]
  y_levels <- names(metric_map)[-1]
  corr_pvals <- apply(idx, 1, function(rc) {
    x <- corr_input[[rc[2]]]
    y <- corr_input[[rc[1]]]
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 3) return(NA_real_)
    suppressWarnings(
      tryCatch(cor.test(x[keep], y[keep], method = "pearson")$p.value, error = function(e) NA_real_)
    )
  })

  corr_df <- tibble(
    row_metric = rownames(corr_mat)[idx[, 1]],
    col_metric = colnames(corr_mat)[idx[, 2]],
    r = corr_mat[idx],
    p_value = as.numeric(corr_pvals)
  ) %>%
    mutate(
      row_metric = factor(row_metric, levels = y_levels),
      col_metric = factor(col_metric, levels = x_levels),
      label = ifelse(is.finite(r), ifelse(abs(r) < 0.01, "|r|<0.01", sprintf("%.2f", r)), ""),
      sig_fontface = ifelse(is.finite(p_value) & p_value < 0.05, "bold", "plain")
    )

  ggplot(corr_df, aes(x = col_metric, y = row_metric, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = label, fontface = sig_fontface),
      size = FIG_TILE_TEXT_SIZE * 0.78,
      family = "sans",
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "firebrick3",
      mid = "white",
      high = "dodgerblue3",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = number_format(accuracy = 0.1),
      na.value = "grey85",
      name = "Pearson's r"
    ) +
    scale_x_discrete(
      labels = function(x) {
        label_text <- metric_labels[x]
        label_text[is.na(label_text)] <- x[is.na(label_text)]
        if (isTRUE(parse_metric_labels)) {
          parse(text = unname(label_text))
        } else {
          unname(label_text)
        }
      },
      limits = x_levels,
      drop = FALSE
    ) +
    scale_y_discrete(
      labels = function(x) {
        label_text <- metric_labels[x]
        label_text[is.na(label_text)] <- x[is.na(label_text)]
        if (isTRUE(parse_metric_labels)) {
          parse(text = unname(label_text))
        } else {
          unname(label_text)
        }
      },
      limits = y_levels,
      drop = FALSE
    ) +
    coord_fixed(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.35),
      axis.line.y.left = element_line(color = "black", linewidth = 0.35),
      axis.line.x.top = element_blank(),
      axis.line.y.right = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.x.bottom = element_line(color = "black", linewidth = 0.35),
      axis.ticks.y.left = element_line(color = "black", linewidth = 0.35),
      axis.ticks.x.top = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE - 1, family = "sans"),
      axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE - 1, family = "sans"),
      axis.ticks.length = grid::unit(4.5, "pt"),
      legend.position = "right",
      legend.title = element_text(size = FIG_AXIS_TEXT_SIZE - 1, family = "sans"),
      legend.text = element_text(size = FIG_AXIS_TEXT_SIZE - 3, family = "sans"),
      legend.key.height = grid::unit(18, "pt"),
      legend.key.width = grid::unit(12, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(8, 6, 8, 6)
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(90, "pt"),
        barwidth = grid::unit(12, "pt")
      )
    )
}

dynamic_metric_map <- c(
  RBI = "RBI",
  RCS = "RCS",
  FDC = "FDC",
  SD = "SD",
  WB = "WB"
)

mobile_metric_map <- c(
  BF = "BF_mean",
  DR = "DR",
  Fyw = "Fyw",
  MTT = "MTT"
)

p_dynamic <- make_corr_matrix(
  annual,
  dynamic_metric_map,
  metric_labels = c(
    RBI = "RBI",
    RCS = "RCS",
    FDC = "FDC",
    SD = "SD",
    WB = "WB"
  )
)

p_mobile <- make_corr_matrix(
  site_summary,
  mobile_metric_map,
  metric_labels = c(
    BF = "plain(BF)",
    DR = "plain(DR)",
    Fyw = "plain(F)[plain(yw)]",
    MTT = "plain(MTT)"
  ),
  parse_metric_labels = TRUE
)

corr_fig_width <- 5.7 * FIG_WIDTH_SCALE
corr_fig_height <- 5.0 * FIG_HEIGHT_SCALE

ggsave(
  file.path(supp_dir, "FigS2_dynamic_storage_corr.png"),
  p_dynamic,
  width = corr_fig_width,
  height = corr_fig_height,
  bg = "white",
  dpi = FIG_PREVIEW_DPI
)
ggsave(
  file.path(supp_pdf_dir, "FigS2_dynamic_storage_corr.pdf"),
  p_dynamic,
  width = corr_fig_width,
  height = corr_fig_height,
  bg = "white"
)
ggsave(
  file.path(supp_tiff_dir, "FigS2_dynamic_storage_corr.tiff"),
  p_dynamic,
  width = corr_fig_width,
  height = corr_fig_height,
  bg = "white",
  dpi = FIG_PRODUCTION_DPI,
  compression = "lzw"
)

ggsave(
  file.path(supp_dir, "FigS3_mobile_storage_corr.png"),
  p_mobile,
  width = corr_fig_width,
  height = corr_fig_height,
  bg = "white",
  dpi = FIG_PREVIEW_DPI
)
ggsave(
  file.path(supp_pdf_dir, "FigS3_mobile_storage_corr.pdf"),
  p_mobile,
  width = corr_fig_width,
  height = corr_fig_height,
  bg = "white"
)
ggsave(
  file.path(supp_tiff_dir, "FigS3_mobile_storage_corr.tiff"),
  p_mobile,
  width = corr_fig_width,
  height = corr_fig_height,
  bg = "white",
  dpi = FIG_PRODUCTION_DPI,
  compression = "lzw"
)
