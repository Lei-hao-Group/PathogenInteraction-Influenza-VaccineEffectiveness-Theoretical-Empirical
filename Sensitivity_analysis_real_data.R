setwd("C:/Users/wghost/Desktop/TND/logistic回归（年龄连续变量输入）")

# ==============================================================================
# 1. Environment setup and package loading
# ==============================================================================
library(tidyverse)
library(lubridate)
library(broom)
library(forestplot)
library(grid)
library(flextable)
library(officer)
library(table1)
library(webshot2)
library(ggrepel)
library(scales)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(coda)
library(ggdist)
library(ggh4x)

windowsFonts(Times = windowsFont("Times New Roman"))
windowsFonts(SimHei = windowsFont("SimHei"))

# ==============================================================================
# 2. Data reading and preprocessing
# ==============================================================================
if(!file.exists("All_data.csv")) stop("Error: File All_data.csv not found!")

raw_data <- read.csv("All_data.csv", stringsAsFactors = FALSE, encoding = "UTF-8")

clean_data <- raw_data %>%
  distinct() %>%
  mutate(Date = as.Date(Date),
         Age_months = as.numeric(Age_months),
         ID_CARD = toupper(ID_CARD)) %>%
  filter(Age_months >= 6 & Age_months <= 204) %>%
  group_by(ID_CARD) %>%
  arrange(Date) %>%
  slice(1) %>%
  ungroup()

# ==============================================================================
# 3. Vaccination status definition
# ==============================================================================
vac_records <- clean_data %>%
  select(ID_CARD, Date, starts_with("流感")) %>%
  pivot_longer(
    cols = matches("^流感"),
    names_to = c("Vac_Index", ".value"),
    names_pattern = "流感(\\d+)(.*)"
  ) %>%
  mutate(Vac_Date_Parsed = as.Date(接种时间)) %>%
  filter(!is.na(Vac_Date_Parsed)) %>%
  filter(Vac_Date_Parsed >= (Date - 730) & Vac_Date_Parsed <= (Date - 14)) %>%
  arrange(ID_CARD, Date, desc(Vac_Date_Parsed)) %>%
  distinct(ID_CARD, Date, .keep_all = TRUE) %>%
  transmute(
    ID_CARD, Date, Vaccinated = 1,
    Effective_Vac_Date = 接种时间, 
    Effective_Vac_Name = 制品名称
  )

final_df <- clean_data %>%
  left_join(vac_records, by = c("ID_CARD", "Date")) %>%
  mutate(Vaccinated = replace_na(Vaccinated, 0)) %>%
  relocate(Effective_Vac_Date, Effective_Vac_Name, .after = Vaccinated)

final_df <- final_df %>%
  mutate(
    Vac_Stratum = case_when(
      Vaccinated == 0 ~ "Unvaccinated",
      Vaccinated == 1 & !str_detect(Effective_Vac_Name, "四价") & str_detect(Effective_Vac_Name, "三价|亚单位|全病毒") ~ "IIV3 (Trivalent)",
      Vaccinated == 1 & str_detect(Effective_Vac_Name, "四价") ~ "IIV4 (Quadrivalent)",
      Vaccinated == 1 & str_detect(Effective_Vac_Name, "鼻喷|减毒") ~ "LAIV3 (Nasal Spray)",
      Vaccinated == 1 ~ "Other/Unknown"
    )
  )

# ==============================================================================
# 4. Variable definition
# ==============================================================================
final_df <- final_df %>%
  mutate(
    Case_AnyFlu = if_else(InfA == 1 | InfB == 1, 1, 0),
    Other_Pathogens_Count = HPIV + HADV + Ch + HRV + HMPV + HCOV + HRSV + Mp + Boca
  )

# --- Define control group markers ---
final_df <- final_df %>%
  mutate(
    Is_Control_Overall = Case_AnyFlu == 0,
    
    # Positive control: Influenza negative + target pathogen positive (other pathogens not restricted)
    Is_Control_PanPos  = Case_AnyFlu == 0 & Other_Pathogens_Count > 0,
    Is_Control_MpPos   = Case_AnyFlu == 0 & Mp == 1,
    Is_Control_HRVPos  = Case_AnyFlu == 0 & HRV == 1,
    Is_Control_HRSVPos = Case_AnyFlu == 0 & HRSV == 1,
    Is_Control_HADVPos = Case_AnyFlu == 0 & HADV == 1,
    Is_Control_HPIVPos = Case_AnyFlu == 0 & HPIV == 1,
    
    # Negative control: Influenza negative + target pathogen negative (other pathogens not restricted)
    Is_Control_PanNeg  = Case_AnyFlu == 0 & Other_Pathogens_Count == 0, 
    Is_Control_MpNeg   = Case_AnyFlu == 0 & Mp == 0,
    Is_Control_HRVNeg  = Case_AnyFlu == 0 & HRV == 0,
    Is_Control_HRSVNeg = Case_AnyFlu == 0 & HRSV == 0,
    Is_Control_HADVNeg = Case_AnyFlu == 0 & HADV == 0,
    Is_Control_HPIVNeg = Case_AnyFlu == 0 & HPIV == 0,
    
    Gender = as.factor(Gender), 
    Type = as.factor(Type), 
    Month = as.factor(format(Date, "%m")),
    Sample_Group = as.factor(case_when(Sample %in% c("咽拭子", "鼻咽拭子") ~ "Upper", TRUE ~ "Lower")),
    Age_Subgroup = case_when(Age_months >= 6 & Age_months <= 35 ~ "6-35 months", 
                             Age_months > 35 & Age_months <= 60 ~ "3-5 years", 
                             TRUE ~ "6-17 years")
  )

# ==============================================================================
# 5. Statistical modeling functions
# ==============================================================================
run_logistic_model <- function(data, ctrl_col, ctrl_label, subgroup_label, 
                               adjust_age = TRUE, adjust_gender = TRUE, adjust_type = TRUE,
                               target_vac_group = "All") {
  
  subset_df <- data
  if(target_vac_group != "All") {
    subset_df <- subset_df %>% filter(Vac_Stratum == "Unvaccinated" | Vac_Stratum == target_vac_group)
  }
  
  subset_df <- subset_df %>%
    filter(Case_AnyFlu == 1 | !!sym(ctrl_col) == TRUE) %>%
    mutate(y = Case_AnyFlu)
  
  n_case_total <- sum(subset_df$y == 1)
  n_case_vac   <- sum(subset_df$y == 1 & subset_df$Vaccinated == 1)
  n_ctrl_total <- sum(subset_df$y == 0)
  n_ctrl_vac   <- sum(subset_df$y == 0 & subset_df$Vaccinated == 1)
  
  if(n_case_vac < 5 | n_ctrl_vac < 5) return(NULL)
  
  case_str <- sprintf("%d/%d", n_case_vac, n_case_total)
  ctrl_str <- sprintf("%d/%d", n_ctrl_vac, n_ctrl_total)
  
  tryCatch({
    # --- 1. Unadjusted model ---
    model_unadj <- glm(y ~ Vaccinated, data = subset_df, family = binomial())
    res_unadj <- tidy(model_unadj, conf.int = TRUE) %>% filter(term == "Vaccinated")
    ve_unadj <- (1 - exp(res_unadj$estimate)) * 100
    unadj_text <- sprintf("%.1f (%.1f, %.1f)", ve_unadj, (1-exp(res_unadj$conf.high))*100, (1-exp(res_unadj$conf.low))*100)
    
    # --- 2. Adjusted model ---
    f_str <- "y ~ Vaccinated"
    if(adjust_gender && length(unique(subset_df$Gender)) > 1) f_str <- paste(f_str, "+ Gender")
    if(length(unique(subset_df$Month)) > 1) f_str <- paste(f_str, "+ Month")
    if(length(unique(subset_df$Sample_Group)) > 1) f_str <- paste(f_str, "+ Sample_Group")
    if(adjust_type && length(unique(subset_df$Type)) > 1) f_str <- paste(f_str, "+ Type")
    if(adjust_age && length(unique(subset_df$Age_months)) > 1) f_str <- paste(f_str, "+  Age_months")
    
    model_adj <- glm(as.formula(f_str), data = subset_df, family = binomial())
    res_adj <- tidy(model_adj, conf.int = TRUE) %>% filter(term == "Vaccinated")
    ve_adj <- (1 - exp(res_adj$estimate)) * 100
    adj_text <- sprintf("%.1f (%.1f, %.1f)", ve_adj, (1-exp(res_adj$conf.high))*100, (1-exp(res_adj$conf.low))*100)
    
    return(list(
      Subgroup = subgroup_label, Control_Type = ctrl_label,
      Case_Info = case_str, Ctrl_Info = ctrl_str,
      Unadj_VE_Text = unadj_text, Adj_VE_Text = adj_text,
      VE = ve_adj, Lower = (1-exp(res_adj$conf.high))*100, Upper = (1-exp(res_adj$conf.low))*100,
      Is_Summary = FALSE
    ))
  }, error = function(e) return(NULL))
}

# ============================================================================================
# 6. Perform the analysis cycle (distinguishing between positive and negative control groups).
# ============================================================================================

pos_control_list <- list(
  list(col = "Is_Control_Overall", label = "Influenza Negative (Flu-)"),
  list(col = "Is_Control_MpPos",   label = "MP Positive (Flu-, MP+)"),
  list(col = "Is_Control_HRVPos",  label = "HRV Positive (Flu-, HRV+)"),
  list(col = "Is_Control_HRSVPos", label = "HRSV Positive (Flu-, HRSV+)"),
  list(col = "Is_Control_HADVPos", label = "HAdV Positive (Flu-, HAdV+)"),
  list(col = "Is_Control_HPIVPos", label = "HPIV Positive (Flu-, HPIV+)"), 
  list(col = "Is_Control_PanPos",  label = "All Non-Influenza Pathogens Positive (Flu-, ORV+)")
)

neg_control_list <- list(
  list(col = "Is_Control_Overall", label = "Influenza Negative (Flu-)"),
  list(col = "Is_Control_MpNeg",   label = "MP Negative (Flu-, MP-)"),
  list(col = "Is_Control_HRVNeg",  label = "HRV Negative (Flu-, HRV-)"),
  list(col = "Is_Control_HRSVNeg", label = "HRSV Negative (Flu-, HRSV-)"),
  list(col = "Is_Control_HADVNeg", label = "HAdV Negative (Flu-, HAdV-)"),
  list(col = "Is_Control_HPIVNeg", label = "HPIV Negative (Flu-, HPIV-)"), 
  list(col = "Is_Control_PanNeg",  label = "All Non-Influenza Pathogens Negative (Pan-)")
)

run_analysis_block <- function(ctrl_list) {
  res_all <- list()
  
  # PART 1: Overall
  for(ctrl in ctrl_list) {
    res <- run_logistic_model(final_df, ctrl$col, ctrl$label, "Overall")
    if(!is.null(res)) { res$Sub_Cat <- "Overall"; res_all[[length(res_all)+1]] <- res }
  }
  
  # PART 2: Age
  age_groups <- list(list(code="6-35 months", label="6-35 months"), 
                     list(code="3-5 years", label="3-5 years"), 
                     list(code="6-17 years", label="6-17 years"))
  for(age in age_groups) {
    sub_df <- final_df %>% filter(Age_Subgroup == age$code)
    for(ctrl in ctrl_list) {
      res <- run_logistic_model(sub_df, ctrl$col, ctrl$label, age$label, adjust_age=FALSE)
      if(!is.null(res)) { res$Sub_Cat <- "Age"; res_all[[length(res_all)+1]] <- res }
    }
  }
  
  # PART 3: Gender
  gender_groups <- list(list(code="男", label="Male"), list(code="女", label="Female"))
  for(g in gender_groups) {
    sub_df <- final_df %>% filter(Gender == g$code)
    for(ctrl in ctrl_list) {
      res <- run_logistic_model(sub_df, ctrl$col, ctrl$label, g$label, adjust_gender=FALSE)
      if(!is.null(res)) { res$Sub_Cat <- "Gender"; res_all[[length(res_all)+1]] <- res }
    }
  }
  
  # PART 4: Type
  type_groups <- list(list(code="门诊", label="Outpatient"), list(code="住院", label="Inpatient"))
  for(t in type_groups) {
    sub_df <- final_df %>% filter(str_detect(Type, t$code))
    for(ctrl in ctrl_list) {
      res <- run_logistic_model(sub_df, ctrl$col, ctrl$label, t$label, adjust_type=FALSE)
      if(!is.null(res)) { res$Sub_Cat <- "Type"; res_all[[length(res_all)+1]] <- res }
    }
  }
  
  # PART 5: Vaccine Type
  target_types <- list(
    list(code = "IIV3 (Trivalent)", label = "Trivalent (IIV3)"),
    list(code = "IIV4 (Quadrivalent)", label = "Quadrivalent (IIV4)"),
    list(code = "LAIV3 (Nasal Spray)", label = "Nasal Spray (LAIV3)")
  )
  for(vac in target_types) {
    if(sum(final_df$Vac_Stratum == vac$code) >= 10) {
      for(ctrl in ctrl_list) {
        res <- run_logistic_model(final_df, ctrl$col, ctrl$label, vac$label, target_vac_group = vac$code)
        if(!is.null(res)) { res$Sub_Cat <- "Vaccine"; res_all[[length(res_all)+1]] <- res }
      }
    }
  }
  
  return(bind_rows(res_all))
}

# Perform positive and negative control group analyses separately.
df_pos <- run_analysis_block(pos_control_list)
df_neg <- run_analysis_block(neg_control_list)

# =================================================================================
# 7. Draw a detailed forest plot (split by positive/negative and overall/subgroup).
# =================================================================================

pos_levels <- sapply(pos_control_list, function(x) x$label)
neg_levels <- sapply(neg_control_list, function(x) x$label)
subgroup_levels <- c("Overall", "6-35 months", "3-5 years", "6-17 years", 
                     "Male", "Female", "Outpatient", "Inpatient", 
                     "Trivalent (IIV3)", "Quadrivalent (IIV4)", "Nasal Spray (LAIV3)")

df_pos <- df_pos %>% mutate(
  Control_Type_Factor = factor(Control_Type, levels = pos_levels),
  Subgroup = factor(Subgroup, levels = subgroup_levels)
)
df_neg <- df_neg %>% mutate(
  Control_Type_Factor = factor(Control_Type, levels = neg_levels),
  Subgroup = factor(Subgroup, levels = subgroup_levels)
)

draw_forest_detailed <- function(data, filename) {
  
  header_df <- data.frame(
    Subgroup = "Group", 
    Control_Type = "Control Type", 
    Case_Info = "Case (Vac/Total)",    
    Ctrl_Info = "Control (Vac/Total)", 
    Unadj_VE_Text = "Unadjusted VE (%) (95% CI)",
    Adj_VE_Text = "Adjusted VE (%) (95% CI)",
    Is_Summary = TRUE, 
    VE = NA, Lower = NA, Upper = NA,
    stringsAsFactors = FALSE
  )
  
  if(nrow(data) == 0) return()
  
  plot_data <- bind_rows(header_df, data %>% arrange(Subgroup, Control_Type_Factor))
  
  plot_data <- plot_data %>%
    mutate(
      Subgroup      = replace_na(as.character(Subgroup), ""),
      Control_Type  = replace_na(as.character(Control_Type), ""),
      Case_Info     = replace_na(Case_Info, ""),
      Ctrl_Info     = replace_na(Ctrl_Info, ""),
      Unadj_VE_Text = replace_na(Unadj_VE_Text, ""),
      Adj_VE_Text   = replace_na(Adj_VE_Text, "")
    )
  
  plot_data$VE    <- as.numeric(plot_data$VE)
  plot_data$Lower <- as.numeric(plot_data$Lower)
  plot_data$Upper <- as.numeric(plot_data$Upper)
  
  plot_data$Subgroup_Display <- plot_data$Subgroup
  plot_data$Subgroup_Display[duplicated(plot_data$Subgroup)] <- ""
  
  row_height_px <- 130
  base_height_px <- 400
  png_height <- base_height_px + (nrow(plot_data) * row_height_px)
  
  png(filename, width = 7750, height = png_height, res = 300)
  
  my_lines <- list("2" = gpar(lwd = 2, col = "black"))
  
  p <- forestplot(
    labeltext = list(
      plot_data$Subgroup_Display, 
      plot_data$Control_Type, 
      plot_data$Case_Info, 
      plot_data$Ctrl_Info, 
      plot_data$Unadj_VE_Text, 
      plot_data$Adj_VE_Text     
    ),
    mean  = plot_data$VE, 
    lower = plot_data$Lower, 
    upper = plot_data$Upper,
    graph.pos = 7, 
    graphwidth = unit(0.1, "npc"), 
    colgap = unit(4.5, "mm"),
    align = c("l", "l", "c", "c", "c", "c"),
    zero = 0,
    is.summary = plot_data$Is_Summary,
    boxsize = 0.15, 
    hrzl_lines = my_lines,
    col = fpColors(box = "black", lines = "black", zero = "gray60", summary = "black"),
    lwd.zero = 1, lwd.ci = 2.2, ci.vertices = TRUE,
    txt_gp = fpTxtGp(
      label   = gpar(fontfamily = "Times", cex = 1.1, fontsize = 20), 
      ticks   = gpar(fontfamily = "Times", cex = 1.0, fontsize = 16),
      xlab    = gpar(fontfamily = "Times", cex = 1.0, fontsize = 18, fontface="bold"),
      summary = gpar(fontfamily = "Times", cex = 1.1, fontsize = 20, fontface = "bold")
    ),
    xlab = "Adjusted VE (%)",
    xticks = seq(-50, 100, 50),
    clip = c(-50, 100)) %>%
    fp_set_zebra_style(val = "#F5F5F5")
  
  print(p)
  dev.off()
  print(paste("Saved Image:", filename))
}

# ==============================================================================
# 8. First, create individual forest images, then piece them together.
# ==============================================================================

# --- 8.1 First, generate single images for positive, negative, overall, and subgroups. ---
plot_categories <- list(
  list(cat = "Overall", name = "Overall_Analysis"),
  list(cat = "Age",     name = "Age_Subgroups"),
  list(cat = "Gender",  name = "Gender_Subgroups"),
  list(cat = "Type",    name = "Type_Subgroups"),
  list(cat = "Vaccine", name = "Vaccine_Subgroups")
)

for (pcat in plot_categories) {
  
  # Positive control group series diagram
  d_pos <- df_pos %>% filter(Sub_Cat == pcat$cat)
  if (nrow(d_pos) > 0) {
    draw_forest_detailed(
      d_pos,
      paste0("PosControl_", pcat$name, "_SA.png")
    )
  }
  
  # Negative control group series diagram
  d_neg <- df_neg %>% filter(Sub_Cat == pcat$cat)
  if (nrow(d_neg) > 0) {
    draw_forest_detailed(
      d_neg,
      paste0("NegControl_", pcat$name, "_SA.png")
    )
  }
}


# ====================================================================================================
# 9. The jigsaw puzzle function reads a PNG image, adds A) B) C) D) to the top left corner, and then stitches them vertically.
# ====================================================================================================

make_png_panel <- function(file, tag, tag_size = 7) {
  
  img <- png::readPNG(file, native = TRUE)
  
  ggplot() +
    
    annotation_custom(
      grob = grid::rasterGrob(img, interpolate = TRUE),
      xmin = 0, xmax = 1,
      ymin = 0, ymax = 0.94
    ) +
    annotate(
      "text",
      x = 0.05, y = 0.985,
      label = tag,
      hjust = 0, vjust = 1,
      family = "Times",
      fontface = "bold",
      size = tag_size
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}


get_png_aspect <- function(file) {
  img <- png::readPNG(file, native = TRUE)
  # dim(img)[1] = height, dim(img)[2] = width
  dim(img)[1] / dim(img)[2] / 0.94
}


combine_png_panels <- function(files,
                               output_file,
                               tags = NULL,
                               width = 16,
                               dpi = 300,
                               tag_size = 7) {
  
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    warning("No files found for: ", output_file)
    return(invisible(NULL))
  }
  
  if (is.null(tags)) {
    tags <- paste0(LETTERS[seq_along(files)], ")")
  } else {
    tags <- tags[seq_along(files)]
  }
  
  panel_list <- purrr::map2(
    files,
    tags,
    ~ make_png_panel(.x, .y, tag_size = tag_size)
  )
  
  panel_heights <- purrr::map_dbl(files, get_png_aspect)
  
  combined_plot <- patchwork::wrap_plots(
    panel_list,
    ncol = 1,
    heights = panel_heights
  )
  
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = width,
    height = width * sum(panel_heights),
    units = "in",
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
  
  message("Saved combined figure: ", output_file)
}


# ==============================================================================
# 10. Overall puzzle: A) Positive control group; B) Negative control group
# ==============================================================================

combine_png_panels(
  files = c(
    "PosControl_Overall_Analysis_SA.png",
    "NegControl_Overall_Analysis_SA.png"
  ),
  output_file = "Figure_Overall_Positive_Negative_Controls_SA.png",
  tags = c("A)", "B)"),
  width = 16,
  dpi = 300,
  tag_size = 8
)


# ======================================================================================================
# 11. Subgroup puzzle pieces: Positive group puzzled separately, negative group puzzled separately.
# ======================================================================================================

# --- 11.1 Positive group：Age, Gender, Type, Vaccine ---
combine_png_panels(
  files = c(
    "PosControl_Age_Subgroups_SA.png",
    "PosControl_Vaccine_Subgroups_SA.png", 
    "PosControl_Gender_Subgroups_SA.png",
    "PosControl_Type_Subgroups_SA.png"
  ),
  output_file = "Figure_Subgroups_Positive_Controls_SA.png",
  tags = c("A)", "B)", "C)", "D)"),
  width = 26,
  dpi = 300,
  tag_size = 12
)


# --- 11.2 Negative group：Age, Gender, Type, Vaccine ---
combine_png_panels(
  files = c(
    "NegControl_Age_Subgroups_SA.png",
    "NegControl_Vaccine_Subgroups_SA.png", 
    "NegControl_Gender_Subgroups_SA.png",
    "NegControl_Type_Subgroups_SA.png"
  ),
  output_file = "Figure_Subgroups_Negative_Controls_SA.png",
  tags = c("A)", "B)", "C)", "D)"),
  width = 26,
  dpi = 300,
  tag_size = 12
)