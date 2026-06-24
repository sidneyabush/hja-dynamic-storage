# Figure S4 dynamic and mobile storage correlation matrix

# inputs:
# outputs/master/master_site.csv

# outputs:
# figs_tables_pub/supp/FigS4_dynamic_mobile_corr.*

# author: Sidney Bush
# date: 2026-02-13

librarian::shelf(dplyr, readr, tidyr, ggplot2, cran_repo = "https://cloud.r-project.org")

rm(list = ls())
source("config.R")

safe_ggsave <- function(filename, plot_obj, width, height, dpi = NULL, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    {
      if (is.null(dpi)) {
        ggplot2::ggsave(filename, plot_obj, width = width, height = height, bg = "white", ...)
      } else {
        ggplot2::ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, bg = "white", ...)
      }
      TRUE
    },
    error = function(e) {
      warning("Failed to save plot: ", filename, " (", conditionMessage(e), ")")
      FALSE
    }
  )
}

supp_dir <- MS_FIG_SUPP_DIR
supp_pdf_dir <- MS_FIG_SUPP_PDF_DIR
supp_tiff_dir <- MS_FIG_SUPP_TIFF_DIR
for (d in c(supp_dir, supp_pdf_dir, supp_tiff_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
site_df <- read_csv(site_file, show_col_types = FALSE) %>%
  filter(site %in% SITE_ORDER_HYDROMETRIC)

dynamic_map <- c(
  RBI = "RBI_mean",
  RCS = "RCS_mean",
  FDC = "FDC_mean",
  SD = "SD_mean",
  WB = "WB_mean"
)
mobile_map <- c(
  BF = "BF_mean",
  DR = "DR",
  Fyw = "Fyw",
  MTT = "MTT"
)

dynamic_map <- dynamic_map[dynamic_map %in% names(site_df)]
mobile_map <- mobile_map[mobile_map %in% names(site_df)]

# stop if the master site table no longer has dynamic or mobile metrics
if (length(dynamic_map) == 0 || length(mobile_map) == 0) {
  stop("Needed dynamic/mobile metric columns are missing from master_site.csv")
}

dynamic_axis_labels_plotmath <- c(
  RBI = "plain(RBI)",
  RCS = "plain(RCS)",
  FDC = "plain(FDC)",
  SD = "plain(SD)",
  WB = "plain(WB)"
)
mobile_axis_labels_plotmath <- c(
  BF = "plain(BF)",
  DR = "plain(DR)",
  Fyw = "plain(F)[yw]",
  MTT = "plain(MTT)"
)

corr_df <- tidyr::expand_grid(
  Mobile = names(mobile_map),
  Dynamic = names(dynamic_map)
) %>%
  mutate(
    mobile_col = unname(mobile_map[Mobile]),
    dynamic_col = unname(dynamic_map[Dynamic]),
    r = mapply(
      function(m_col, d_col) {
        x <- suppressWarnings(as.numeric(site_df[[m_col]]))
        y <- suppressWarnings(as.numeric(site_df[[d_col]]))
        suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
      },
      mobile_col,
      dynamic_col
    ),
    p = mapply(
      function(m_col, d_col) {
        x <- suppressWarnings(as.numeric(site_df[[m_col]]))
        y <- suppressWarnings(as.numeric(site_df[[d_col]]))
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 3) {
          return(NA_real_)
        }
        suppressWarnings(cor.test(x[ok], y[ok], method = "pearson")$p.value)
      },
      mobile_col,
      dynamic_col
    )
  ) %>%
  mutate(
    Mobile = factor(Mobile, levels = rev(names(mobile_map))),
    Dynamic = factor(Dynamic, levels = names(dynamic_map)),
    label = ifelse(
      is.finite(r),
      ifelse(abs(r) < 0.05, "|r| < 0.05", sprintf("%.2f", r)),
      ""
    ),
    sig_fontface = ifelse(!is.na(p) & p < 0.05, "bold", "plain")
  )

p <- ggplot(corr_df, aes(x = Dynamic, y = Mobile, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, fontface = sig_fontface), size = FIG_TILE_TEXT_SIZE) +
  scale_x_discrete(
    labels = function(x) {
      parsed <- dynamic_axis_labels_plotmath[x]
      parsed[is.na(parsed)] <- x[is.na(parsed)]
      parse(text = unname(parsed))
    }
  ) +
  scale_y_discrete(
    labels = function(x) {
      parsed <- mobile_axis_labels_plotmath[x]
      parsed[is.na(parsed)] <- x[is.na(parsed)]
      parse(text = unname(parsed))
    }
  ) +
  scale_fill_gradient2(
    low = "firebrick3",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",
    name = "Pearson's r"
  ) +
  labs(x = "Dynamic Storage Metrics", y = "Mobile Storage Metrics") +
  theme_pub() +
  theme(
    axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = FIG_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1)
  )

invisible(safe_ggsave(
  file.path(supp_dir, "FigS4_dynamic_mobile_corr.png"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.8 * FIG_HEIGHT_SCALE,
  dpi = FIG_PREVIEW_DPI
))

invisible(safe_ggsave(
  file.path(supp_pdf_dir, "FigS4_dynamic_mobile_corr.pdf"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.8 * FIG_HEIGHT_SCALE
))
invisible(safe_ggsave(
  file.path(supp_tiff_dir, "FigS4_dynamic_mobile_corr.tiff"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.8 * FIG_HEIGHT_SCALE,
  dpi = FIG_PRODUCTION_DPI,
  compression = "lzw"
))
