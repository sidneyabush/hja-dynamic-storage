# figure 8 pooled eco response mlr beta coefficients

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

rm(list = ls())
source("config.R")

output_dir <- OUT_MODELS_STORAGE_ECO_RESPONSE_MLR_DIR
main_dir <- MS_FIG_MAIN_DIR
main_pdf_dir <- MS_FIG_MAIN_PDF_DIR
for (d in c(main_dir, main_pdf_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

results_file <- file.path(output_dir, "storage_eco_response_mlr_results.csv")
if (!file.exists(results_file)) {
  stop("Missing required file: ", results_file)
}

norm_response <- function(x) {
  out <- as.character(x)
  out[out == "Q_7Q5"] <- "Q7Q5"
  out[out == "T_7DMax"] <- "T7DMax"
  out
}

norm_predictor <- function(x) {
  out <- as.character(x)
  out[out %in% c("P_WetSeason", "Pwetseason", "Pwetseason")] <- "Pws"
  out
}

response_order <- c("Q7Q5", "T7DMax")
predictor_order <- c("Pws", "RBI", "RCS", "FDC", "SD", "WB", "CHS")

coef_df <- read_csv(results_file, show_col_types = FALSE) %>%
  mutate(
    Response = norm_response(Response),
    Predictor = norm_predictor(Predictor),
    Beta_Std = suppressWarnings(as.numeric(Beta_Std)),
    p_value = suppressWarnings(as.numeric(p_value))
  ) %>%
  filter(Response %in% response_order, Predictor %in% predictor_order)

if (nrow(coef_df) == 0) {
  stop("No rows found for Q7Q5/T7DMax in eco model results")
}

plot_df <- expand_grid(
  Response = response_order,
  Predictor = predictor_order
) %>%
  left_join(
    coef_df %>%
      group_by(Response, Predictor) %>%
      slice(1) %>%
      ungroup(),
    by = c("Response", "Predictor")
  ) %>%
  mutate(
    Response = factor(Response, levels = response_order),
    Predictor = factor(Predictor, levels = rev(predictor_order)),
    sig = ifelse(is.finite(p_value) & p_value <= 0.05, "*", ""),
    label = ifelse(is.finite(Beta_Std), paste0(sprintf("%.2f", Beta_Std), sig), "")
  )

predictor_axis_labels <- function(vals) {
  txt <- as.character(vals)
  txt[txt == "Pws"] <- "P[ws]"
  parse(text = txt)
}

p <- ggplot(plot_df, aes(x = Response, y = Predictor, fill = Beta_Std)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = FIG_TILE_TEXT_SIZE) +
  scale_fill_gradient2(
    low = "firebrick3",
    mid = "white",
    high = "dodgerblue3",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = expression(italic(beta))
  ) +
  scale_y_discrete(labels = predictor_axis_labels) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
    axis.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.title = element_text(size = FIG_AXIS_TITLE_SIZE + 1),
    legend.text = element_text(size = FIG_AXIS_TEXT_SIZE + 1)
  )

ggsave(
  file.path(main_dir, "Fig8_eco_mlr_beta.png"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.4 * FIG_HEIGHT_SCALE,
  dpi = 300
)

ggsave(
  file.path(main_pdf_dir, "Fig8_eco_mlr_beta.pdf"),
  p,
  width = 7.2 * FIG_WIDTH_SCALE,
  height = 5.4 * FIG_HEIGHT_SCALE
)
