# figure 6 dynamic mobile storage correlation matrix

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

rm(list = ls())
source("config.R")

main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
for (d in c(main_dir, main_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

site_file <- file.path(OUTPUT_DIR, "master", MASTER_SITE_FILE)
if (!file.exists(site_file)) {
  stop("Missing required file: ", site_file)
}

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
  CHS = "CHS_mean",
  DR = "DR",
  Fyw = "Fyw",
  MTT1 = "MTT1",
  MTT2 = "MTT2"
)

dynamic_map <- dynamic_map[dynamic_map %in% names(site_df)]
mobile_map <- mobile_map[mobile_map %in% names(site_df)]

if (length(dynamic_map) == 0 || length(mobile_map) == 0) {
  stop("Required dynamic/mobile metric columns are missing from master_site.csv")
}

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
    )
  ) %>%
  mutate(
    Mobile = factor(Mobile, levels = rev(names(mobile_map))),
    Dynamic = factor(Dynamic, levels = names(dynamic_map)),
    label = ifelse(is.finite(r), sprintf("%.2f", r), "")
  )

p <- ggplot(corr_df, aes(x = Dynamic, y = Mobile, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = FIG_TILE_TEXT_SIZE) +
  scale_fill_gradient2(
    low = "firebrick3",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",
    name = "Pearson's r"
  ) +
  labs(x = "Dynamic", y = "Mobile") +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1)
  )

ggsave(
  file.path(main_dir, "Fig6_dynamic_mobile_corr.png"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.8 * FIG_HEIGHT_SCALE,
  dpi = 300
)

ggsave(
  file.path(main_pdf_dir, "Fig6_dynamic_mobile_corr.pdf"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.8 * FIG_HEIGHT_SCALE
)
