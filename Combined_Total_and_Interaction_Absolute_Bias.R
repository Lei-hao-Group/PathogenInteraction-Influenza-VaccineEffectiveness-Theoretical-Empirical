# ================================================================
# Combined code: quantify total absolute bias and additional absolute
# bias caused by pathogen interaction
#
# Bias scale used in this script:
#   Total absolute bias = (VE_estimated - VE_true) * 100
#   Extra interaction bias = (VE_estimated(theta) - VE_estimated(theta = 1)) * 100
# Values are shown in percentage points (pp), not relative bias (%).
# ================================================================

# Optional: set your own working directory
setwd("C:/Users/wghost/Desktop/TND/热图代码")

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

# ================================================================
# 1. Parameter grid
# ================================================================
param_grid <- expand.grid(
  theta = seq(0, 2, length.out = 500),
  lambda_flu = seq(0.01, 0.3, length.out = 500),
  ratio_XO = c(0.5, 1.0, 2.0),
  VE = c(0.3, 0.5, 0.7),
  KEEP.OUT.ATTRS = FALSE
)

# ================================================================
# 2. Bias calculation
#    A. Total absolute bias:
#       VEe - VE, VE1 - VE, VE2 - VE
#    B. Extra absolute bias caused by pathogen interaction:
#       VE(theta) - VE(theta = 1)
# ================================================================
calculate_bias <- function(df) {
  df %>%
    filter(theta * lambda_flu < 1) %>%
    mutate(
      A = 1 - VE,

      # Denominators under observed theta
      D_O = 1 - lambda_flu * A,
      D_X = 1 - theta * lambda_flu * A,
      D_mix = ratio_XO * D_X + D_O,

      # Denominators under no pathogen interaction (theta = 1)
      D_X_no_interaction = D_O,
      D_mix_no_interaction = (ratio_XO + 1) * D_O,

      # ------------------------------------------------------------
      # Estimated ORs under different control definitions
      # ------------------------------------------------------------
      # VEe: standard influenza-negative controls, mixed X and O
      OR_Mixed = A * (
        ratio_XO * (1 - theta * lambda_flu) + (1 - lambda_flu)
      ) / D_mix,

      # VE1: influenza-negative but pathogen X-positive controls
      OR_X = A * (1 - theta * lambda_flu) / D_X,

      # VE2: influenza-negative and pathogen X-negative controls
      OR_O = A * (1 - lambda_flu) / D_O,

      # No-interaction reference estimates (theta = 1)
      OR_Mixed_no_interaction = A * (
        ratio_XO * (1 - lambda_flu) + (1 - lambda_flu)
      ) / D_mix_no_interaction,

      OR_X_no_interaction = A * (1 - lambda_flu) / D_X_no_interaction,

      # ------------------------------------------------------------
      # VE estimates
      # ------------------------------------------------------------
      VEe_est = 1 - OR_Mixed,
      VE1_est = 1 - OR_X,
      VE2_est = 1 - OR_O,

      VEe_no_interaction = 1 - OR_Mixed_no_interaction,
      VE1_no_interaction = 1 - OR_X_no_interaction,

      # ------------------------------------------------------------
      # A. Total absolute bias, in percentage points
      # ------------------------------------------------------------
      Total_Bias_Mixed = (VEe_est - VE) * 100,
      Total_Bias_X = (VE1_est - VE) * 100,
      Total_Bias_O = (VE2_est - VE) * 100,

      # Equivalent closed-form formulas:
      # Total_Bias_X     = VE * A * theta * lambda_flu / D_X * 100
      # Total_Bias_O     = VE * A * lambda_flu / D_O * 100
      # Total_Bias_Mixed = VE * A * lambda_flu * (ratio_XO * theta + 1) / D_mix * 100

      # ------------------------------------------------------------
      # B. Extra absolute bias caused by pathogen interaction,
      #    in percentage points
      # ------------------------------------------------------------
      Extra_Bias_Mixed = (VEe_est - VEe_no_interaction) * 100,
      Extra_Bias_X = (VE1_est - VE1_no_interaction) * 100,
      Extra_Bias_O = 0,

      # Equivalent closed-form formulas:
      # Extra_Bias_X     = A * VE * lambda_flu * (theta - 1) / D_O / D_X * 100
      # Extra_Bias_Mixed = A * VE * ratio_XO * lambda_flu * (theta - 1) / D_O / D_mix * 100

      # Facet labels
      VE_label = factor(
        paste0("VE = ", scales::percent(VE, accuracy = 1)),
        levels = c("VE = 30%", "VE = 50%", "VE = 70%")
      ),
      ratio_label = factor(
        paste0("lambda[X]/lambda[O] == ", ratio_XO),
        levels = c(
          "lambda[X]/lambda[O] == 0.5",
          "lambda[X]/lambda[O] == 1",
          "lambda[X]/lambda[O] == 2"
        )
      )
    )
}

plot_data <- calculate_bias(param_grid)
plot_data_base <- plot_data %>% filter(ratio_XO == 1.0)

# ================================================================
# 3. Quantification tables
# ================================================================
quant_long <- bind_rows(
  plot_data %>%
    transmute(
      measure = "Total absolute bias: mixed endpoint VEe",
      VE, ratio_XO = as.character(ratio_XO), theta, lambda_flu,
      bias_pp = Total_Bias_Mixed
    ),
  plot_data_base %>%
    transmute(
      measure = "Total absolute bias: pathogen X VE1",
      VE, ratio_XO = "not applicable", theta, lambda_flu,
      bias_pp = Total_Bias_X
    ),
  plot_data_base %>%
    transmute(
      measure = "Total absolute bias: pathogen O VE2",
      VE, ratio_XO = "not applicable", theta, lambda_flu,
      bias_pp = Total_Bias_O
    ),
  plot_data %>%
    transmute(
      measure = "Extra interaction absolute bias: mixed endpoint VEe",
      VE, ratio_XO = as.character(ratio_XO), theta, lambda_flu,
      bias_pp = Extra_Bias_Mixed
    ),
  plot_data_base %>%
    transmute(
      measure = "Extra interaction absolute bias: pathogen X VE1",
      VE, ratio_XO = "not applicable", theta, lambda_flu,
      bias_pp = Extra_Bias_X
    )
)

bias_summary <- quant_long %>%
  group_by(measure, VE, ratio_XO) %>%
  summarise(
    min_bias_pp = min(bias_pp, na.rm = TRUE),
    p25_bias_pp = quantile(bias_pp, 0.25, na.rm = TRUE),
    median_bias_pp = median(bias_pp, na.rm = TRUE),
    p75_bias_pp = quantile(bias_pp, 0.75, na.rm = TRUE),
    max_bias_pp = max(bias_pp, na.rm = TRUE),
    .groups = "drop"
  )

# ================================================================
# 4. Plot settings
# ================================================================
# Absolute bias is expressed as percentage points (pp).
# In this parameter space, total absolute bias is approximately 0 to 22 pp.
total_limits <- c(0, 25)
total_breaks <- seq(0, 25, 5)

extra_limits <- c(-10, 15)
extra_breaks <- seq(-10, 15, 5)

fill_total <- scale_fill_gradient(
  low = "white",
  high = "#E74C3C",
  limits = total_limits,
  breaks = total_breaks,
  oob = scales::squish,
  name = "Total Absolute\nBias (%)",
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(5.2, "cm"),
    barwidth = unit(1.0, "cm"),
    ticks = TRUE
  )
)

fill_extra <- scale_fill_gradient2(
  low = "#2E86C1",
  mid = "white",
  high = "#E74C3C",
  limits = extra_limits,
  breaks = extra_breaks,
  oob = scales::squish,
  name = "Pathogen Interaction Extra\nAbsolute Bias (%)",
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(5.2, "cm"),
    barwidth = unit(1.0, "cm"),
    ticks = TRUE
  )
)

my_theme <- theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55),

    axis.text = element_text(size = 20, color = "gray25"),
    axis.title = element_text(size = 20, color = "black"),

    strip.background = element_rect(fill = "#f2f2f2", color = "black", linewidth = 0.55),
    strip.text.x = element_text(
      face = "bold", size = 20,
      margin = margin(t = 6, b = 6, l = 8, r = 8)
    ),
    strip.text.y.right = element_text(
      face = "bold",
      size = 20,
      margin = margin(l = 8, r = 8)
    ),

    panel.spacing.x = unit(1.25, "cm"),
    panel.spacing.y = unit(0.75, "cm"),

    plot.title = element_text(size = 24, hjust = 0, face = "bold"),
    plot.margin = margin(5, 5, 5, 5),

    legend.position = "right",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

# Helper function for heatmaps
make_heatmap <- function(data, bias_var, fill_scale, title_text,
                         facet_type = c("grid", "wrap")) {
  facet_type <- match.arg(facet_type)

  p <- ggplot(
    data,
    aes(
      x = theta,
      y = lambda_flu,
      fill = .data[[bias_var]],
      z = .data[[bias_var]]
    )
  ) +
    geom_tile() +
    geom_contour(color = "gray80", alpha = 0.6, bins = 8, linewidth = 0.25) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.85) +
    fill_scale +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 0.5)) +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
    labs(title = title_text, x = expression(theta), y = expression(lambda[F])) +
    my_theme

  if (facet_type == "grid") {
    p <- p + facet_grid(
      VE_label ~ ratio_label,
      labeller = labeller(ratio_label = label_parsed)
    )
  } else {
    p <- p + facet_wrap(~ VE_label, nrow = 1)
  }

  p
}

# ================================================================
# 5. Figure A: total absolute bias
# ================================================================
p_total_mixed <- make_heatmap(
  plot_data,
  bias_var = "Total_Bias_Mixed",
  fill_scale = fill_total,
  title_text = "A)",
  facet_type = "grid"
)

p_total_x <- make_heatmap(
  plot_data_base,
  bias_var = "Total_Bias_X",
  fill_scale = fill_total,
  title_text = "B)",
  facet_type = "wrap"
)

p_total_o <- make_heatmap(
  plot_data_base,
  bias_var = "Total_Bias_O",
  fill_scale = fill_total,
  title_text = "C)",
  facet_type = "wrap"
)

fig_total <- p_total_mixed / p_total_x / p_total_o +
  plot_layout(
    ncol = 1,
    heights = c(3.1, 1.25, 1.25),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 6)
  )

# ================================================================
# 6. Figure B: extra absolute bias due to pathogen interaction
# ================================================================
p_extra_mixed <- make_heatmap(
  plot_data,
  bias_var = "Extra_Bias_Mixed",
  fill_scale = fill_extra,
  title_text = "A)",
  facet_type = "grid"
)

p_extra_x <- make_heatmap(
  plot_data_base,
  bias_var = "Extra_Bias_X",
  fill_scale = fill_extra,
  title_text = "B)",
  facet_type = "wrap"
)

fig_extra <- p_extra_mixed / p_extra_x +
  plot_layout(
    ncol = 1,
    heights = c(3.1, 1.25),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 6)
  )

# ================================================================
# 7. Optional full combined figure
#    The total-bias and extra-bias panels use different color scales,
#    so the two legends are intentionally kept separate.
# ================================================================
fig_combined <- wrap_elements(full = fig_total) / wrap_elements(full = fig_extra) +
  plot_layout(
    ncol = 1,
    heights = c(15.0, 11.0)
  )

# ================================================================
# 8. Save outputs
# ================================================================
ggsave(
  "Figure_Total_Absolute_Bias_0603.png",
  plot = fig_total,
  width = 16,
  height = 15,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure_Interaction_Extra_Absolute_Bias_0603.png",
  plot = fig_extra,
  width = 16,
  height = 15,
  dpi = 300,
  bg = "white"
)

# Optional: save combined figure
# ggsave(
#   "Figure_Combined_Total_and_Interaction_Absolute_Bias_0603.png",
#   plot = fig_combined,
#   width = 16,
#   height = 26,
#   dpi = 300,
#   bg = "white"
# )
