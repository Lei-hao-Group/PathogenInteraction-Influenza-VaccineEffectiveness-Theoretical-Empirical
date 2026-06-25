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
library(patchwork)
library(magick)
library(cowplot)

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
  dplyr::select(ID_CARD, Date, tidyselect::starts_with("流感")) %>%
  pivot_longer(
    cols = tidyselect::matches("^流感"),
    names_to = c("Vac_Index", ".value"),
    names_pattern = "流感(\\d+)(.*)"
  ) %>%
  mutate(Vac_Date_Parsed = as.Date(接种时间)) %>%
  filter(!is.na(Vac_Date_Parsed)) %>%
  filter(Vac_Date_Parsed >= (Date - 365) & Vac_Date_Parsed <= (Date - 14)) %>%
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

# ================================================================================================
# 6. Perform the analysis cycle (distinguishing between positive and negative control groups).
# ================================================================================================

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
      paste0("PosControl_", pcat$name, ".png")
    )
  }
  
  # Negative control group series diagram
  d_neg <- df_neg %>% filter(Sub_Cat == pcat$cat)
  if (nrow(d_neg) > 0) {
    draw_forest_detailed(
      d_neg,
      paste0("NegControl_", pcat$name, ".png")
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
    "PosControl_Overall_Analysis.png",
    "NegControl_Overall_Analysis.png"
  ),
  output_file = "Figure_Overall_Positive_Negative_Controls.png",
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
    "PosControl_Age_Subgroups.png",
    "PosControl_Vaccine_Subgroups.png", 
    "PosControl_Gender_Subgroups.png",
    "PosControl_Type_Subgroups.png"
  ),
  output_file = "Figure_Subgroups_Positive_Controls.png",
  tags = c("A)", "B)", "C)", "D)"),
  width = 26,
  dpi = 300,
  tag_size = 12
)


# --- 11.2 Negative group：Age, Gender, Type, Vaccine ---
combine_png_panels(
  files = c(
    "NegControl_Age_Subgroups.png",
    "NegControl_Vaccine_Subgroups.png", 
    "NegControl_Gender_Subgroups.png",
    "NegControl_Type_Subgroups.png"
  ),
  output_file = "Figure_Subgroups_Negative_Controls.png",
  tags = c("A)", "B)", "C)", "D)"),
  width = 26,
  dpi = 300,
  tag_size = 12
)

# ==============================================================================
# 8. Descriptive statistical analysis (Table 2) -> Export to Word
# ==============================================================================
table_data <- final_df %>%
  mutate(
    Flu_Status = factor(if_else(Case_AnyFlu == 1, "Influenza Positive", "Influenza Negative"),
                        levels = c("Influenza Negative", "Influenza Positive")),
    Gender = factor(Gender, levels = c("男", "女"), labels = c("Male", "Female")), 
    Age_Group = factor(Age_Subgroup, levels = c("6-35 months", "3-5 years", "6-17 years")),
    Vaccination_Status = factor(Vaccinated, levels = c(0, 1), labels = c("Unvaccinated", "Vaccinated")),
    Vaccine_Type = factor(Vac_Stratum,
                          levels = c("Unvaccinated", "IIV3 (Trivalent)", "IIV4 (Quadrivalent)", "LAIV3 (Nasal Spray)"),
                          labels = c("Unvaccinated", "IIV3", "IIV4", "LAIV3")),
    Mp_Label   = factor(Mp, levels = c(0, 1), labels = c("Negative", "Positive")),
    HRV_Label  = factor(HRV, levels = c(0, 1), labels = c("Negative", "Positive")),
    HRSV_Label = factor(HRSV, levels = c(0, 1), labels = c("Negative", "Positive")),
    HPIV_Label = factor(HPIV, levels = c(0, 1), labels = c("Negative", "Positive")),
    HADV_Label = factor(HADV, levels = c(0, 1), labels = c("Negative", "Positive")), 
    ORV_Label  = factor(ifelse(Other_Pathogens_Count > 0, 1, 0), levels = c(0, 1), labels = c("Negative", "Positive"))
  )

label(table_data$Age_months)         <- "Age, Months"
label(table_data$Gender)             <- "Gender"
label(table_data$Age_Group)          <- "Age Group"
label(table_data$Vaccination_Status) <- "Vaccination History"
label(table_data$Vaccine_Type)       <- "Vaccine Type"
label(table_data$Mp_Label)           <- "MP"
label(table_data$HRV_Label)          <- "HRV"
label(table_data$HRSV_Label)         <- "HRSV"
label(table_data$HPIV_Label)         <- "HPIV"
label(table_data$HADV_Label)         <- "HAdV"
label(table_data$ORV_Label)          <- "All Non-Influenza Pathogens (ORV)"

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("", "Median (IQR), months"=sprintf("%s (%s, %s)", MEDIAN, Q1, Q3)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.2f%%)", FREQ, PCT))))
}

table1_res <- table1(~ Gender + Age_months + Age_Group +
                       Vaccination_Status + Vaccine_Type +
                       Mp_Label + HRV_Label + HRSV_Label + HPIV_Label + HADV_Label + ORV_Label| Flu_Status,
                     data = table_data,
                     render.continuous = my.render.cont,
                     render.categorical = my.render.cat,
                     overall = c(left="Overall Population"))

ft <- t1flex(table1_res) %>%
  autofit() %>%  
  font(fontname = "Times New Roman", part = "all") %>% 
  fontsize(size = 10, part = "all") 

save_as_docx(ft, path = "Table1_Descriptive_Statistics.docx")

# ==============================================================================
# 9. Bayesian QBA / MCMC validation for the 0603 theoretical derivation
#    Paste this section after your original sections 1-8, after final_df has been created.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(MASS)
  library(coda)
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
  library(ggdist)
  library(scales)
  library(patchwork)
})

# Keep random-number generation reproducible across repeated runs.
try({
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}, silent = TRUE)
set.seed(2026)

base_qba_outdir <- "Bayesian_QBA_0603_V07_NUTS_prior_sensitivity_MCMC"
dir.create(base_qba_outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 9.1 Prior scenarios and base MCMC settings
# ------------------------------------------------------------------------------
# Scenario definitions follow the prior-sensitivity plan:
# Main, Weak, Flat VE, Wide theta, and Skeptical theta.
prior_scenarios <- tibble::tribble(
  ~scenario_id, ~Prior_Scenario,       ~VE_prior_a, ~VE_prior_b, ~lambda_prior_a, ~lambda_prior_b, ~logtheta_prior_mean, ~logtheta_prior_sd,
  1L,           "Main",                 6,           7,           13,              87,              0,                    0.50,
  2L,           "Weak",                 2,           2,            2,              18,              0,                    1.00,
  3L,           "Flat_VE",              1,           1,            2,              18,              0,                    1.00,
  4L,           "Wide_theta",           6,           7,           13,              87,              0,                    1.50,
  5L,           "Skeptical_theta",      6,           7,           13,              87,              0,                    0.25
) %>%
  dplyr::mutate(
    VE_prior_label = paste0("Beta(", VE_prior_a, ",", VE_prior_b, ")"),
    lambda_prior_label = paste0("Beta(", lambda_prior_a, ",", lambda_prior_b, ")"),
    logtheta_prior_label = paste0("Normal(0, ", logtheta_prior_sd, "^2)"),
    scenario_outdir = file.path(base_qba_outdir, Prior_Scenario)
  )

readr::write_csv(prior_scenarios, file.path(base_qba_outdir, "Prior_Scenarios.csv"))
print(prior_scenarios)

base_mcmc_settings <- list(
  # NUTS / cmdstanr settings.
  # The scientific model is kept the same as V03; only the sampler changes.
  sampler = "NUTS_cmdstanr",
  nuts_iter_warmup = 2000,
  nuts_iter_sampling = 4000,
  nuts_chains = 4,
  nuts_parallel_chains = 4,
  nuts_thin = 1,
  nuts_adapt_delta = 0.95,
  nuts_max_treedepth = 12,
  nuts_metric = "dense_e",
  nuts_seed_base = 2026000L,
  
  # Legacy fields retained so downstream output tables remain comparable to V03.
  n_iter = NA_integer_,
  burn   = NA_integer_,
  thin   = 1,
  n_chains = 4,
  
  # Priors; overwritten by each prior scenario.
  VE_prior_a = NA_real_,
  VE_prior_b = NA_real_,
  lambda_prior_a = NA_real_,
  lambda_prior_b = NA_real_,
  logtheta_prior_mean = 0,
  logtheta_prior_sd = NA_real_,
  
  # Main likelihood components.
  use_beta_1 = TRUE,
  use_beta_2 = TRUE,
  use_logR  = TRUE,
  
  # G-computation uncertainty.
  n_sim_gcomp = 1000,
  propagate_r_uncertainty = TRUE,
  
  # Seed for posterior predictive r_X uncertainty.
  r_seed = 2027001L,
  
  # Stan model filename. The script writes this file automatically if absent.
  stan_file = "qba_tnd_nuts_0603.stan"
)
# ------------------------------------------------------------------------------
# 9.2 Basic checks, pathogen definitions and helper functions
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 9.2 Basic checks and pathogen definitions
# ------------------------------------------------------------------------------
if (!exists("final_df")) {
  stop("final_df does not exist. Please run your original preprocessing code first.")
}

required_cols <- c("Case_AnyFlu", "Vaccinated", "Gender", "Month", "Sample_Group", "Type", "Age_months",
                   "Mp", "HRV", "HRSV", "HADV", "HPIV", "Other_Pathogens_Count",
                   "Is_Control_Overall", "Is_Control_MpPos", "Is_Control_HRVPos", "Is_Control_HRSVPos",
                   "Is_Control_HADVPos", "Is_Control_HPIVPos", "Is_Control_PanPos",
                   "Is_Control_MpNeg", "Is_Control_HRVNeg", "Is_Control_HRSVNeg",
                   "Is_Control_HADVNeg", "Is_Control_HPIVNeg", "Is_Control_PanNeg")
missing_cols <- setdiff(required_cols, names(final_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in final_df: ", paste(missing_cols, collapse = ", "))
}

final_df <- final_df %>%
  dplyr::mutate(
    Vaccinated = as.numeric(Vaccinated),
    ORV_Pos = as.integer(Other_Pathogens_Count > 0),
    Gender = as.factor(Gender),
    Type = as.factor(Type),
    Month = as.factor(Month),
    Sample_Group = as.factor(Sample_Group)
  )

pathogen_map <- tibble::tribble(
  ~Pathogen, ~x_col,    ~pos_ctrl,             ~neg_ctrl,
  "MP",      "Mp",     "Is_Control_MpPos",   "Is_Control_MpNeg",
  "HRV",     "HRV",    "Is_Control_HRVPos",  "Is_Control_HRVNeg",
  "HRSV",    "HRSV",   "Is_Control_HRSVPos", "Is_Control_HRSVNeg",
  "HAdV",    "HADV",   "Is_Control_HADVPos", "Is_Control_HADVNeg",
  "HPIV",    "HPIV",   "Is_Control_HPIVPos", "Is_Control_HPIVNeg",
  "ORV",     "ORV_Pos","Is_Control_PanPos",  "Is_Control_PanNeg"
)

pathogen_levels <- pathogen_map$Pathogen

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
clamp_p <- function(p) clamp(p, 1e-8, 1 - 1e-8)

# ------------------------------------------------------------------------------
# 9.3 Helper functions: adjusted TND VE and Flu-negative X-composition G-computation
# ------------------------------------------------------------------------------
make_adjusted_formula <- function(outcome, data) {
  rhs <- c("Vaccinated")
  if ("Gender" %in% names(data) && length(unique(na.omit(data$Gender))) > 1) rhs <- c(rhs, "Gender")
  if ("Month" %in% names(data) && length(unique(na.omit(data$Month))) > 1) rhs <- c(rhs, "Month")
  if ("Sample_Group" %in% names(data) && length(unique(na.omit(data$Sample_Group))) > 1) rhs <- c(rhs, "Sample_Group")
  if ("Type" %in% names(data) && length(unique(na.omit(data$Type))) > 1) rhs <- c(rhs, "Type")
  if ("Age_months" %in% names(data) && length(unique(na.omit(data$Age_months))) > 1) rhs <- c(rhs, "Age_months")
  stats::as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

fit_tnd_logOR <- function(data, ctrl_col, label = ctrl_col) {
  dat <- data %>%
    dplyr::filter(Case_AnyFlu == 1 | .data[[ctrl_col]] == TRUE) %>%
    dplyr::mutate(y = as.integer(Case_AnyFlu == 1)) %>%
    tidyr::drop_na(y, Vaccinated, Gender, Month, Sample_Group, Type, Age_months)
  
  if (sum(dat$y == 1 & dat$Vaccinated == 1) < 5 || sum(dat$y == 0 & dat$Vaccinated == 1) < 5) {
    return(tibble::tibble(Control = label, beta = NA_real_, se = NA_real_, VE_obs = NA_real_,
                          VE_obs_percent = NA_real_, n = nrow(dat), n_case = sum(dat$y == 1), n_ctrl = sum(dat$y == 0)))
  }
  
  f <- make_adjusted_formula("y", dat)
  fit <- stats::glm(f, data = dat, family = binomial())
  co <- summary(fit)$coefficients
  
  if (!"Vaccinated" %in% rownames(co)) {
    return(tibble::tibble(Control = label, beta = NA_real_, se = NA_real_, VE_obs = NA_real_,
                          VE_obs_percent = NA_real_, n = nrow(dat), n_case = sum(dat$y == 1), n_ctrl = sum(dat$y == 0)))
  }
  
  beta <- unname(co["Vaccinated", "Estimate"])
  se   <- unname(co["Vaccinated", "Std. Error"])
  
  tibble::tibble(
    Control = label,
    beta = beta,
    se = se,
    VE_obs = 1 - exp(beta),
    VE_obs_percent = 100 * (1 - exp(beta)),
    VE_obs_lower = 1 - exp(beta + 1.96 * se),
    VE_obs_upper = 1 - exp(beta - 1.96 * se),
    n = nrow(dat),
    n_case = sum(dat$y == 1),
    n_ctrl = sum(dat$y == 0),
    n_case_vac = sum(dat$y == 1 & dat$Vaccinated == 1),
    n_ctrl_vac = sum(dat$y == 0 & dat$Vaccinated == 1)
  )
}

fit_gcomp_logR <- function(data, x_col, n_sim = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  dat <- data %>%
    dplyr::filter(Case_AnyFlu == 0) %>%
    dplyr::mutate(x = as.integer(.data[[x_col]] == 1)) %>%
    tidyr::drop_na(x, Vaccinated, Gender, Month, Sample_Group, Type, Age_months)
  
  if (sum(dat$x == 1 & dat$Vaccinated == 1) < 5 || sum(dat$x == 1 & dat$Vaccinated == 0) < 5) {
    return(tibble::tibble(logR = NA_real_, se_logR = NA_real_, p_v = NA_real_, p_u = NA_real_,
                          log_odds_u = NA_real_, se_log_odds_u = NA_real_, n_flu_neg = nrow(dat),
                          n_x_pos = sum(dat$x == 1)))
  }
  
  f <- make_adjusted_formula("x", dat)
  fit <- stats::glm(f, data = dat, family = binomial())
  
  new_v <- dat; new_v$Vaccinated <- 1
  new_u <- dat; new_u$Vaccinated <- 0
  tt <- stats::delete.response(stats::terms(fit))
  Xv <- stats::model.matrix(tt, new_v)
  Xu <- stats::model.matrix(tt, new_u)
  
  beta_hat <- stats::coef(fit)
  vc <- stats::vcov(fit)
  
  calc_stats <- function(beta_vec) {
    pv_i <- stats::plogis(as.numeric(Xv %*% beta_vec))
    pu_i <- stats::plogis(as.numeric(Xu %*% beta_vec))
    pv <- clamp_p(mean(pv_i))
    pu <- clamp_p(mean(pu_i))
    odds_v <- pv / (1 - pv)
    odds_u <- pu / (1 - pu)
    c(logR = log(odds_v / odds_u), p_v = pv, p_u = pu, log_odds_u = log(odds_u))
  }
  
  point <- calc_stats(beta_hat)
  
  sim_stats <- tryCatch({
    draws <- MASS::mvrnorm(n = n_sim, mu = beta_hat, Sigma = vc)
    t(apply(draws, 1, calc_stats))
  }, error = function(e) NULL)
  
  if (is.null(sim_stats)) {
    se_logR <- abs(summary(fit)$coefficients["Vaccinated", "Std. Error"])
    se_log_odds_u <- 0.50
  } else {
    se_logR <- stats::sd(sim_stats[, "logR"], na.rm = TRUE)
    se_log_odds_u <- stats::sd(sim_stats[, "log_odds_u"], na.rm = TRUE)
  }
  
  tibble::tibble(
    logR = unname(point["logR"]),
    se_logR = max(se_logR, 1e-6),
    p_v = unname(point["p_v"]),
    p_u = unname(point["p_u"]),
    log_odds_u = unname(point["log_odds_u"]),
    se_log_odds_u = max(se_log_odds_u, 1e-6),
    n_flu_neg = nrow(dat),
    n_x_pos = sum(dat$x == 1)
  )
}

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 9.4 Build observed inputs from hospital data ONCE
# ------------------------------------------------------------------------------
# Observed inputs do not depend on prior distributions. They are generated once
# using fixed seeds for G-computation uncertainty, and reused for all scenarios.
beta_e_fit <- fit_tnd_logOR(final_df, "Is_Control_Overall", "VEe: Flu- overall")

obs_inputs <- pathogen_map %>%
  dplyr::mutate(seed_gcomp = 202600L + dplyr::row_number()) %>%
  purrr::pmap_dfr(function(Pathogen, x_col, pos_ctrl, neg_ctrl, seed_gcomp) {
    b1 <- fit_tnd_logOR(
      final_df,
      ctrl_col = pos_ctrl,
      label = paste0("VE1: Flu-, ", Pathogen, "+")
    )
    
    b2 <- fit_tnd_logOR(
      final_df,
      ctrl_col = neg_ctrl,
      label = paste0("VE2: Flu-, ", Pathogen, "-")
    )
    
    gr <- fit_gcomp_logR(
      final_df,
      x_col = x_col,
      n_sim = base_mcmc_settings$n_sim_gcomp,
      seed = as.integer(seed_gcomp)
    )
    
    tibble::tibble(
      Pathogen = Pathogen,
      x_col = x_col,
      pos_ctrl = pos_ctrl,
      neg_ctrl = neg_ctrl,
      seed_gcomp = seed_gcomp,
      beta_e = beta_e_fit$beta,
      se_e = beta_e_fit$se,
      VEe_obs = beta_e_fit$VE_obs,
      VEe_obs_lower = beta_e_fit$VE_obs_lower,
      VEe_obs_upper = beta_e_fit$VE_obs_upper,
      beta_1 = b1$beta,
      se_1 = b1$se,
      VE1_obs = b1$VE_obs,
      VE1_obs_lower = b1$VE_obs_lower,
      VE1_obs_upper = b1$VE_obs_upper,
      beta_2 = b2$beta,
      se_2 = b2$se,
      VE2_obs = b2$VE_obs,
      VE2_obs_lower = b2$VE_obs_lower,
      VE2_obs_upper = b2$VE_obs_upper,
      logR = gr$logR,
      se_logR = gr$se_logR,
      p_v = gr$p_v,
      p_u = gr$p_u,
      log_odds_u = gr$log_odds_u,
      se_log_odds_u = gr$se_log_odds_u,
      n_flu_neg = gr$n_flu_neg,
      n_x_pos = gr$n_x_pos
    )
  })

readr::write_csv(obs_inputs, file.path(base_qba_outdir, "QBA_Observed_Inputs.csv"))
print(obs_inputs)

# ------------------------------------------------------------------------------
# 9.5 0603 theoretical functions
# ------------------------------------------------------------------------------
theory_0603 <- function(VE, lambda, theta) {
  VE     <- as.numeric(unname(VE))[1]
  lambda <- as.numeric(unname(lambda))[1]
  theta  <- as.numeric(unname(theta))[1]
  
  A <- 1 - VE
  
  if (!is.finite(VE) || !is.finite(lambda) || !is.finite(theta)) return(NULL)
  if (VE <= 0 || VE >= 1 || lambda <= 0 || lambda >= 1 || theta <= 0) return(NULL)
  if ((1 - theta * lambda) <= 0) return(NULL)
  if ((1 - theta * lambda * A) <= 0) return(NULL)
  if ((1 - lambda * A) <= 0) return(NULL)
  
  OR1 <- A * (1 - theta * lambda) / (1 - theta * lambda * A)
  OR2 <- A * (1 - lambda) / (1 - lambda * A)
  RX  <- ((1 - theta * lambda * A) * (1 - lambda)) /
    ((1 - lambda * A) * (1 - theta * lambda))
  
  if (any(!is.finite(c(OR1, OR2, RX))) || any(c(OR1, OR2, RX) <= 0)) return(NULL)
  
  out <- c(
    log(OR1),
    log(OR2),
    log(RX),
    1 - OR1,
    1 - OR2,
    RX,
    theta * lambda
  )
  
  names(out) <- c(
    "beta_1",
    "beta_2",
    "logR",
    "VE1_theory",
    "VE2_theory",
    "R_theory",
    "phi"
  )
  
  out
}

add_predictions_long <- function(posterior_wide, obs_inputs, settings) {
  theta_cols <- paste0("theta_", pathogen_levels)
  logtheta_cols <- paste0("logtheta_", pathogen_levels)
  
  long_list <- vector("list", length(pathogen_levels))
  
  for (j in seq_along(pathogen_levels)) {
    p <- pathogen_levels[j]
    long_list[[j]] <- posterior_wide %>%
      dplyr::transmute(
        postwarmup_iter, chain, accept_stat_mean_over_chains,
        Pathogen = p,
        VE_true,
        lambda_F,
        log_theta = .data[[paste0("logtheta_", p)]],
        theta = .data[[paste0("theta_", p)]],
        phi = .data[[paste0("phi_", p)]]
      )
  }
  
  d <- dplyr::bind_rows(long_list) %>%
    dplyr::left_join(
      obs_inputs %>%
        dplyr::select(Pathogen, VEe_obs, beta_e, se_e, log_odds_u, se_log_odds_u),
      by = "Pathogen"
    )
  
  if (isTRUE(settings$propagate_r_uncertainty)) {
    set.seed(settings$r_seed)
    d <- d %>%
      dplyr::mutate(
        se_log_odds_u_use = dplyr::if_else(
          is.finite(se_log_odds_u) & se_log_odds_u > 0,
          se_log_odds_u,
          0
        ),
        log_odds_u_draw = stats::rnorm(
          n = dplyr::n(),
          mean = log_odds_u,
          sd = se_log_odds_u_use
        )
      )
  } else {
    d <- d %>% dplyr::mutate(log_odds_u_draw = log_odds_u)
  }
  
  d %>%
    dplyr::mutate(
      A = 1 - VE_true,
      OR1_theory = A * (1 - theta * lambda_F) / (1 - theta * lambda_F * A),
      OR2_theory = A * (1 - lambda_F) / (1 - lambda_F * A),
      R_theory = ((1 - theta * lambda_F * A) * (1 - lambda_F)) /
        ((1 - lambda_F * A) * (1 - theta * lambda_F)),
      VE1_theory = 1 - OR1_theory,
      VE2_theory = 1 - OR2_theory,
      beta_1_theory = log(OR1_theory),
      beta_2_theory = log(OR2_theory),
      logR_theory = log(R_theory),
      
      # r_X = lambda_X/lambda_O, derived from unvaccinated Flu-negative X odds.
      r_implied = exp(log_odds_u_draw) *
        (1 - lambda_F) / (1 - theta * lambda_F),
      
      valid_vee = is.finite(r_implied) & r_implied > 0 &
        (1 - theta * lambda_F) > 0 &
        (1 - theta * lambda_F * A) > 0 &
        (1 - lambda_F * A) > 0,
      
      ORe_theory = dplyr::if_else(
        valid_vee,
        A * (r_implied * (1 - theta * lambda_F) + (1 - lambda_F)) /
          (r_implied * (1 - theta * lambda_F * A) + (1 - lambda_F * A)),
        NA_real_
      ),
      beta_e_theory = log(ORe_theory),
      VEe_theory = 1 - ORe_theory
    )
}

# ------------------------------------------------------------------------------
# 9.6 Stan / NUTS implementation of the same 0603 joint model
# ------------------------------------------------------------------------------
# The model is equivalent to the V03 likelihood but uses Stan's NUTS sampler.
# To avoid the hard invalid region theta_X * lambda_F >= 1, the Stan model samples
# phi_X = theta_X * lambda_F directly, constrained to (0,1), then recovers
# theta_X = phi_X / lambda_F. The prior is still placed on log(theta_X).

stan_code_0603 <- '
// qba_tnd_nuts_0603.stan
// Stan/NUTS version of the 0603 Bayesian QBA model.

data {
  int<lower=1> P;

  vector[P] beta1_obs;
  vector[P] beta2_obs;
  vector[P] logR_obs;

  vector<lower=1e-12>[P] se_beta1;
  vector<lower=1e-12>[P] se_beta2;
  vector<lower=1e-12>[P] se_logR;

  int<lower=0, upper=1> use_beta1;
  int<lower=0, upper=1> use_beta2;
  int<lower=0, upper=1> use_logR;

  real<lower=0> VE_prior_a;
  real<lower=0> VE_prior_b;
  real<lower=0> lambda_prior_a;
  real<lower=0> lambda_prior_b;

  real logtheta_prior_mean;
  real<lower=0> logtheta_prior_sd;
}

parameters {
  real<lower=1e-9, upper=1 - 1e-9> VE_true;
  real<lower=1e-9, upper=1 - 1e-9> lambda_F;
  vector<lower=1e-9, upper=1 - 1e-9>[P] phi;
}

transformed parameters {
  vector[P] theta;
  vector[P] logtheta;
  vector[P] beta1_theory;
  vector[P] beta2_theory;
  vector[P] logR_theory;
  real A;

  A = 1 - VE_true;

  for (p in 1:P) {
    theta[p] = phi[p] / lambda_F;
    logtheta[p] = log(theta[p]);

    beta1_theory[p] = log(A) + log1m(phi[p]) - log1m(phi[p] * A);
    beta2_theory[p] = log(A) + log1m(lambda_F) - log1m(lambda_F * A);
    logR_theory[p] = log1m(phi[p] * A) + log1m(lambda_F)
                     - log1m(lambda_F * A) - log1m(phi[p]);
  }
}

model {
  VE_true ~ beta(VE_prior_a, VE_prior_b);
  lambda_F ~ beta(lambda_prior_a, lambda_prior_b);

  // Prior is on logtheta as in V03. Since logtheta = log(phi) - log(lambda_F),
  // d(logtheta)/d(phi) = 1/phi, so the Jacobian contribution is -log(phi).
  for (p in 1:P) {
    target += normal_lpdf(logtheta[p] | logtheta_prior_mean, logtheta_prior_sd);
    target += -log(phi[p]);
  }

  for (p in 1:P) {
    if (use_beta1 == 1) beta1_obs[p] ~ normal(beta1_theory[p], se_beta1[p]);
    if (use_beta2 == 1) beta2_obs[p] ~ normal(beta2_theory[p], se_beta2[p]);
    if (use_logR == 1) logR_obs[p] ~ normal(logR_theory[p], se_logR[p]);
  }
}

generated quantities {
  vector[P] VE1_theory;
  vector[P] VE2_theory;
  vector[P] R_theory;
  real A_gq;

  A_gq = 1 - VE_true;

  for (p in 1:P) {
    VE1_theory[p] = 1 - A_gq * (1 - phi[p]) / (1 - phi[p] * A_gq);
    VE2_theory[p] = 1 - A_gq * (1 - lambda_F) / (1 - lambda_F * A_gq);
    R_theory[p] = ((1 - phi[p] * A_gq) * (1 - lambda_F)) /
                  ((1 - lambda_F * A_gq) * (1 - phi[p]));
  }
}
'

write_stan_model_if_missing <- function(stan_file) {
  # Always write the model to avoid accidentally compiling an older Stan file
  # with the same name in the working directory.
  readr::write_lines(stan_code_0603, stan_file)
  invisible(stan_file)
}

# build_stan_data <- function(valid_obs, settings) {
#   list(
#     P = nrow(valid_obs),
#     beta1_obs = as.numeric(valid_obs$beta_1),
#     beta2_obs = as.numeric(valid_obs$beta_2),
#     logR_obs  = as.numeric(valid_obs$logR),
#     se_beta1 = pmax(as.numeric(valid_obs$se_1), 1e-12),
#     se_beta2 = pmax(as.numeric(valid_obs$se_2), 1e-12),
#     se_logR  = pmax(as.numeric(valid_obs$se_logR), 1e-12),
#     use_beta1 = as.integer(isTRUE(settings$use_beta_1)),
#     use_beta2 = as.integer(isTRUE(settings$use_beta_2)),
#     use_logR  = as.integer(isTRUE(settings$use_logR)),
#     VE_prior_a = as.numeric(settings$VE_prior_a),
#     VE_prior_b = as.numeric(settings$VE_prior_b),
#     lambda_prior_a = as.numeric(settings$lambda_prior_a),
#     lambda_prior_b = as.numeric(settings$lambda_prior_b),
#     logtheta_prior_mean = as.numeric(settings$logtheta_prior_mean),
#     logtheta_prior_sd = as.numeric(settings$logtheta_prior_sd)
#   )
# }

build_stan_data <- function(valid_obs, settings) {
  
  logR_obs_use <- if (isTRUE(settings$use_logR)) {
    as.numeric(valid_obs$logR)
  } else {
    rep(0, nrow(valid_obs))
  }
  
  se_logR_use <- if (isTRUE(settings$use_logR)) {
    pmax(as.numeric(valid_obs$se_logR), 1e-12)
  } else {
    rep(1, nrow(valid_obs))
  }
  
  list(
    P = nrow(valid_obs),
    beta1_obs = as.numeric(valid_obs$beta_1),
    beta2_obs = as.numeric(valid_obs$beta_2),
    logR_obs  = logR_obs_use,
    
    se_beta1 = pmax(as.numeric(valid_obs$se_1), 1e-12),
    se_beta2 = pmax(as.numeric(valid_obs$se_2), 1e-12),
    se_logR  = se_logR_use,
    
    use_beta1 = as.integer(isTRUE(settings$use_beta_1)),
    use_beta2 = as.integer(isTRUE(settings$use_beta_2)),
    use_logR  = as.integer(isTRUE(settings$use_logR)),
    
    VE_prior_a = as.numeric(settings$VE_prior_a),
    VE_prior_b = as.numeric(settings$VE_prior_b),
    lambda_prior_a = as.numeric(settings$lambda_prior_a),
    lambda_prior_b = as.numeric(settings$lambda_prior_b),
    logtheta_prior_mean = as.numeric(settings$logtheta_prior_mean),
    logtheta_prior_sd = as.numeric(settings$logtheta_prior_sd)
  )
}

make_nuts_init <- function(VE0, lambda0, theta0, obs_df) {
  VE0 <- clamp(VE0, 0.02, 0.95)
  lambda0 <- clamp(lambda0, 0.002, 0.60)
  theta0 <- rep(theta0, nrow(obs_df))
  theta0 <- clamp(theta0, 0.05, 5.0)
  phi0 <- clamp(theta0 * lambda0, 1e-5, 0.95)
  list(VE_true = VE0, lambda_F = lambda0, phi = as.numeric(phi0))
}

stan_draws_to_posterior_wide <- function(fit, valid_obs, scenario_name) {
  draws_df <- posterior::as_draws_df(fit$draws(
    variables = c("VE_true", "lambda_F", "theta", "logtheta", "phi")
  ))
  
  out <- tibble::tibble(
    postwarmup_iter = as.integer(draws_df$.iteration),
    chain = as.integer(draws_df$.chain),
    VE_true = as.numeric(draws_df$VE_true),
    lambda_F = as.numeric(draws_df$lambda_F)
  )
  
  for (j in seq_len(nrow(valid_obs))) {
    p <- as.character(valid_obs$Pathogen[j])
    out[[paste0("logtheta_", p)]] <- as.numeric(draws_df[[paste0("logtheta[", j, "]")]])
    out[[paste0("theta_", p)]] <- as.numeric(draws_df[[paste0("theta[", j, "]")]])
    out[[paste0("phi_", p)]] <- as.numeric(draws_df[[paste0("phi[", j, "]")]])
  }
  
  out %>%
    dplyr::mutate(
      Prior_Scenario = scenario_name
    )
}

summarise_nuts_diagnostics <- function(fit, valid_obs, scenario_name) {
  diag_vars <- c("VE_true", "lambda_F",
                 paste0("theta[", seq_len(nrow(valid_obs)), "]"),
                 paste0("logtheta[", seq_len(nrow(valid_obs)), "]"))
  diag_names <- c("VE_true", "lambda_F",
                  paste0("theta_", valid_obs$Pathogen),
                  paste0("logtheta_", valid_obs$Pathogen))
  
  ss <- posterior::summarise_draws(fit$draws(variables = diag_vars))
  ss <- as.data.frame(ss)
  ss$Parameter <- diag_names[match(ss$variable, diag_vars)]
  
  tibble::tibble(
    Prior_Scenario = scenario_name,
    Parameter = ss$Parameter,
    Rhat_rank_split = as.numeric(ss$rhat),
    ESS_bulk = as.numeric(ss$ess_bulk),
    ESS_tail = as.numeric(ss$ess_tail)
  )
}

summarise_nuts_sampler <- function(fit, scenario_name) {
  sd <- fit$sampler_diagnostics(format = "df")
  
  get_col <- function(df, choices) {
    hit <- intersect(choices, names(df))
    if (length(hit) == 0) return(rep(NA_real_, nrow(df)))
    as.numeric(df[[hit[1]]])
  }
  
  accept_stat <- get_col(sd, c("accept_stat__", "accept_stat"))
  divergent <- get_col(sd, c("divergent__", "divergent"))
  treedepth <- get_col(sd, c("treedepth__", "treedepth"))
  n_leapfrog <- get_col(sd, c("n_leapfrog__", "n_leapfrog"))
  energy <- get_col(sd, c("energy__", "energy"))
  
  ebfmi_by_chain <- tryCatch({
    sd %>%
      dplyr::group_by(.chain) %>%
      dplyr::summarise(
        E_BFMI = {
          e <- get_col(dplyr::cur_data(), c("energy__", "energy"))
          if (length(e) < 2 || stats::var(e, na.rm = TRUE) == 0) {
            NA_real_
          } else {
            mean(diff(e)^2, na.rm = TRUE) / stats::var(e, na.rm = TRUE)
          }
        },
        .groups = "drop"
      )
  }, error = function(e) tibble::tibble(.chain = sort(unique(sd$.chain)), E_BFMI = NA_real_))
  
  sd %>%
    dplyr::group_by(.chain) %>%
    dplyr::summarise(
      Prior_Scenario = scenario_name,
      Diagnostic_Type = "NUTS_sampler",
      accept_stat_mean = mean(accept_stat, na.rm = TRUE),
      n_divergent = sum(divergent, na.rm = TRUE),
      divergent_rate = mean(divergent, na.rm = TRUE),
      mean_treedepth = mean(treedepth, na.rm = TRUE),
      max_treedepth = max(treedepth, na.rm = TRUE),
      mean_n_leapfrog = mean(n_leapfrog, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(ebfmi_by_chain, by = ".chain") %>%
    dplyr::rename(chain = .chain)
}

run_nuts_joint <- function(valid_obs, settings, inits, scenario_id, scenario_name) {
  stan_file <- write_stan_model_if_missing(settings$stan_file)
  mod <- cmdstanr::cmdstan_model(stan_file)
  stan_data <- build_stan_data(valid_obs, settings)
  
  fit <- mod$sample(
    data = stan_data,
    chains = settings$nuts_chains,
    parallel_chains = settings$nuts_parallel_chains,
    iter_warmup = settings$nuts_iter_warmup,
    iter_sampling = settings$nuts_iter_sampling,
    thin = settings$nuts_thin,
    seed = as.integer(settings$nuts_seed_base + scenario_id),
    adapt_delta = settings$nuts_adapt_delta,
    max_treedepth = settings$nuts_max_treedepth,
    metric = settings$nuts_metric,
    init = inits,
    refresh = 500
  )
  
  posterior_wide <- stan_draws_to_posterior_wide(fit, valid_obs, scenario_name)
  acceptance_table <- summarise_nuts_sampler(fit, scenario_name)
  diagnostic_table <- summarise_nuts_diagnostics(fit, valid_obs, scenario_name)
  
  list(fit = fit, posterior_wide = posterior_wide, acceptance_table = acceptance_table,
       diagnostic_table = diagnostic_table)
}

# ------------------------------------------------------------------------------
# Helper: raincloud plot for posterior theta distributions
# Order within each pathogen: boxplot -> points -> density
# ------------------------------------------------------------------------------
make_density_components <- function(data,
                                    value_col,
                                    group_cols,
                                    x_base_col,
                                    width = 0.30,
                                    n = 512) {
  value_col <- rlang::ensym(value_col)
  x_base_col <- rlang::ensym(x_base_col)
  curve_df <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::group_modify(~{
      x <- .x %>% dplyr::pull(!!value_col)
      x <- x[is.finite(x)]
      if (length(x) < 5 || stats::sd(x, na.rm = TRUE) == 0) {
        return(tibble::tibble(x = numeric(0), y = numeric(0), x_base = numeric(0)))
      }
      den <- stats::density(x, n = n)
      dens_scaled <- den$y / max(den$y, na.rm = TRUE)
      x_base <- unique(.x %>% dplyr::pull(!!x_base_col))[1]
      tibble::tibble(
        y = den$x,
        x = x_base - dens_scaled * width,
        x_base = x_base
      )
    }) %>%
    dplyr::ungroup()
  poly_df <- curve_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::group_modify(~{
      d <- .x %>% dplyr::arrange(y)
      if (nrow(d) == 0) return(tibble::tibble(x = numeric(0), y = numeric(0)))
      tibble::tibble(
        y = c(d$y, rev(d$y)),
        x = c(d$x, rep(d$x_base[1], nrow(d)))
      )
    }) %>%
    dplyr::ungroup()
  edge_df <- curve_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      x_edge = max(x, na.rm = TRUE),
      y_min = min(y, na.rm = TRUE),
      y_max = max(y, na.rm = TRUE),
      .groups = "drop"
    )
  list(poly = poly_df, curve = curve_df, edge = edge_df)
}

summarise_interval <- function(data, group_cols, value_col) {
  value_col <- rlang::ensym(value_col)
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      mean = mean(!!value_col, na.rm = TRUE),
      median = median(!!value_col, na.rm = TRUE),
      q025 = stats::quantile(!!value_col, 0.025, na.rm = TRUE),
      q25 = stats::quantile(!!value_col, 0.25, na.rm = TRUE),
      q75 = stats::quantile(!!value_col, 0.75, na.rm = TRUE),
      q975 = stats::quantile(!!value_col, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

make_shared_parameter_plot <- function(posterior_wide,
                                       scenario_name = "Main",
                                       parameter = c("VE_true", "lambda_F"),
                                       panel_label = "A)",
                                       x_label = NULL,
                                       y_label = "VE posterior value (%)",
                                       fill_color = "#3B6FB6") {
  parameter <- match.arg(parameter)
  if (is.null(x_label)) {
    x_label <- if (parameter == "VE_true") "VE" else expression(lambda[F])
  }
  plot_df <- posterior_wide %>%
    dplyr::filter(Prior_Scenario == scenario_name) %>%
    dplyr::transmute(
      Parameter = parameter,
      value = .data[[parameter]] * 100,
      x_base = 1.00,
      x_box = 1.18,
      x_point = 1.30
    ) %>%
    dplyr::filter(is.finite(value))
  summary_df <- summarise_interval(plot_df, c("Parameter"), value) %>%
    dplyr::mutate(x_box = 1.18)
  dens <- make_density_components(plot_df, value, c("Parameter"), x_base, width = 0.36)
  set.seed(2026)
  point_df <- plot_df %>%
    dplyr::slice_sample(n = min(2500L, nrow(plot_df))) %>%
    dplyr::mutate(x_jitter = x_point + stats::runif(dplyr::n(), -0.035, 0.035))
  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = dens$poly,
      ggplot2::aes(x = x, y = y),
      fill = scales::alpha(fill_color, 0.28),
      color = NA
    ) +
    ggplot2::geom_path(
      data = dens$curve,
      ggplot2::aes(x = x, y = y),
      color = fill_color,
      linewidth = 0.90,
      lineend = "round",
      linejoin = "round"
    ) +
    ggplot2::geom_segment(
      data = dens$edge,
      ggplot2::aes(x = x_edge, xend = x_edge, y = y_min, yend = y_max),
      color = fill_color,
      linewidth = 0.90,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(x = x_box, xend = x_box, y = q025, yend = q975),
      color = fill_color,
      linewidth = 0.85
    ) +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(x = x_box - 0.045, xend = x_box + 0.045, y = q025, yend = q025),
      color = fill_color,
      linewidth = 0.85
    ) +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(x = x_box - 0.045, xend = x_box + 0.045, y = q975, yend = q975),
      color = fill_color,
      linewidth = 0.85
    ) +
    ggplot2::geom_rect(
      data = summary_df,
      ggplot2::aes(xmin = x_box - 0.07, xmax = x_box + 0.07, ymin = q25, ymax = q75),
      fill = scales::alpha(fill_color, 0.25),
      color = fill_color,
      linewidth = 0.85
    ) +
    ggplot2::geom_segment(
      data = summary_df,
      ggplot2::aes(x = x_box - 0.07, xend = x_box + 0.07, y = median, yend = median),
      color = fill_color,
      linewidth = 0.95
    ) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = x_jitter, y = value),
      color = scales::alpha(fill_color, 0.25),
      size = 0.55
    ) +
    ggplot2::geom_point(
      data = summary_df,
      ggplot2::aes(x = x_box, y = mean),
      shape = 23,
      fill = "white",
      color = "black",
      stroke = 0.8,
      size = 2.8
    ) +
    ggplot2::scale_x_continuous(limits = c(0.48, 1.42), breaks = NULL) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
    ggplot2::labs(x = x_label, y = y_label, title = panel_label) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      text = ggplot2::element_text(),
      plot.title = ggplot2::element_text(size = 28, face = "bold", hjust = 0, margin = ggplot2::margin(b = 0)),
      plot.title.position = "plot",
      axis.title.x = ggplot2::element_text(size = 18, margin = ggplot2::margin(t = 2)),
      axis.title.y = ggplot2::element_text(size = 18, margin = ggplot2::margin(r = 2)),
      axis.text.y = ggplot2::element_text(size = 18, color = "black"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 8, l = 0)
    )
}

make_theta_raincloud <- function(posterior_draws, pathogen_levels, qba_outdir = NULL, scenario_name = "Main") {
  theta_df <- posterior_draws %>%
    dplyr::filter(Prior_Scenario == scenario_name) %>%
    dplyr::mutate(
      Pathogen = factor(Pathogen, levels = pathogen_levels),
      x_id = as.numeric(Pathogen),
      x_base = x_id - 0.05,
      x_box = x_id + 0.12,
      x_point = x_id + 0.28
    ) %>%
    dplyr::filter(is.finite(theta), !is.na(Pathogen))
  theta_summary <- summarise_interval(theta_df, c("Pathogen", "x_box"), theta)
  set.seed(2026)
  theta_sample <- theta_df %>%
    dplyr::group_by(Pathogen) %>%
    dplyr::group_modify(~dplyr::slice_sample(.x, n = min(2500L, nrow(.x)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_jitter = x_point + stats::runif(dplyr::n(), -0.04, 0.04))
  dens <- make_density_components(theta_df, theta, c("Pathogen"), x_base, width = 0.34)
  theta_cols <- c(
    "MP" = "#1b4332",
    "HRV" = "#b8c000",
    "HRSV" = "#8c564b",
    "HAdV" = "#7e57c2",
    "HPIV" = "#00a6c8",
    "ORV" = "#43a047"
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.65) +
    ggplot2::geom_polygon(
      data = dens$poly,
      ggplot2::aes(x = x, y = y, group = Pathogen, fill = Pathogen),
      alpha = 0.28,
      color = NA
    ) +
    ggplot2::geom_path(
      data = dens$curve,
      ggplot2::aes(x = x, y = y, group = Pathogen, color = Pathogen),
      linewidth = 0.85,
      lineend = "round",
      linejoin = "round"
    ) +
    ggplot2::geom_segment(
      data = dens$edge,
      ggplot2::aes(x = x_edge, xend = x_edge, y = y_min, yend = y_max, color = Pathogen),
      linewidth = 0.85,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = theta_summary,
      ggplot2::aes(x = x_box, xend = x_box, y = q025, yend = q975, color = Pathogen),
      linewidth = 0.80
    ) +
    ggplot2::geom_segment(
      data = theta_summary,
      ggplot2::aes(x = x_box - 0.045, xend = x_box + 0.045, y = q025, yend = q025, color = Pathogen),
      linewidth = 0.80
    ) +
    ggplot2::geom_segment(
      data = theta_summary,
      ggplot2::aes(x = x_box - 0.045, xend = x_box + 0.045, y = q975, yend = q975, color = Pathogen),
      linewidth = 0.80
    ) +
    ggplot2::geom_rect(
      data = theta_summary,
      ggplot2::aes(xmin = x_box - 0.065, xmax = x_box + 0.065, ymin = q25, ymax = q75, fill = Pathogen, color = Pathogen),
      alpha = 0.35,
      linewidth = 0.80
    ) +
    ggplot2::geom_segment(
      data = theta_summary,
      ggplot2::aes(x = x_box - 0.065, xend = x_box + 0.065, y = median, yend = median, color = Pathogen),
      linewidth = 0.95
    ) +
    ggplot2::geom_point(
      data = theta_sample,
      ggplot2::aes(x = x_jitter, y = theta, color = Pathogen),
      alpha = 0.18,
      size = 0.45
    ) +
    ggplot2::geom_point(
      data = theta_summary,
      ggplot2::aes(x = x_box, y = mean),
      shape = 23,
      fill = "white",
      color = "black",
      stroke = 0.8,
      size = 2.8
    ) +
    ggplot2::scale_x_continuous(breaks = seq_along(pathogen_levels), labels = pathogen_levels, limits = c(0.40, length(pathogen_levels) + 0.75)) +
    ggplot2::scale_y_continuous(breaks = seq(0, ceiling(max(theta_df$theta, na.rm = TRUE)), by = 1)) +
    ggplot2::scale_color_manual(values = theta_cols) +
    ggplot2::scale_fill_manual(values = theta_cols) +
    ggplot2::labs(x = "Pathogen", y = expression(theta ~" posterior value"), title = "C)") +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      text = ggplot2::element_text(),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 28, face = "bold", hjust = 0, margin = ggplot2::margin(b = 0)),
      plot.title.position = "plot",
      axis.title.x = ggplot2::element_text(size = 18, margin = ggplot2::margin(t = 4)),
      axis.title.y = ggplot2::element_text(size = 18, margin = ggplot2::margin(r = 2)),
      axis.text.x = ggplot2::element_text(size = 18, color = "black"),
      axis.text.y = ggplot2::element_text(size = 18, color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    )
  if (!is.null(qba_outdir)) {
    ggplot2::ggsave(file.path(qba_outdir, "Joint_Posterior_theta_raincloud_by_pathogen.png"), p, width = 10, height = 6.5, dpi = 300)
  }
  return(p)
}

make_ve_validation_raincloud <- function(posterior_draws,
                                         obs_inputs,
                                         qba_outdir = NULL,
                                         scenario_name = "Main",
                                         pathogen_levels = c("MP", "HRV", "HRSV", "HAdV", "HPIV", "ORV")) {
  posterior_long <- posterior_draws %>%
    dplyr::filter(Prior_Scenario == scenario_name) %>%
    dplyr::select(Pathogen, VEe_theory, VE1_theory, VE2_theory) %>%
    tidyr::pivot_longer(cols = c(VEe_theory, VE1_theory, VE2_theory), names_to = "VE_Type", values_to = "VE_theory") %>%
    dplyr::mutate(
      VE_Type = dplyr::recode(VE_Type, VEe_theory = "VE[e]", VE1_theory = "VE[1]", VE2_theory = "VE[2]"),
      VE_Type = factor(VE_Type, levels = c("VE[e]", "VE[1]", "VE[2]")),
      Pathogen = factor(Pathogen, levels = pathogen_levels),
      x_id = as.numeric(Pathogen),
      x_base = x_id - 0.12,
      x_box = x_id + 0.03,
      VE_percent = 100 * VE_theory
    ) %>%
    dplyr::filter(is.finite(VE_percent), !is.na(Pathogen))
  theory_summary <- summarise_interval(posterior_long, c("VE_Type", "Pathogen", "x_box"), VE_percent)
  dens <- make_density_components(posterior_long, VE_percent, c("VE_Type", "Pathogen"), x_base, width = 0.30)
  obs_long <- dplyr::bind_rows(
    obs_inputs %>% dplyr::transmute(Pathogen, VE_Type = "VE[e]", VE_obs = 100 * VEe_obs, VE_lower = 100 * VEe_obs_lower, VE_upper = 100 * VEe_obs_upper),
    obs_inputs %>% dplyr::transmute(Pathogen, VE_Type = "VE[1]", VE_obs = 100 * VE1_obs, VE_lower = 100 * VE1_obs_lower, VE_upper = 100 * VE1_obs_upper),
    obs_inputs %>% dplyr::transmute(Pathogen, VE_Type = "VE[2]", VE_obs = 100 * VE2_obs, VE_lower = 100 * VE2_obs_lower, VE_upper = 100 * VE2_obs_upper)
  ) %>%
    dplyr::mutate(
      VE_Type = factor(VE_Type, levels = c("VE[e]", "VE[1]", "VE[2]")),
      Pathogen = factor(Pathogen, levels = pathogen_levels),
      x_id = as.numeric(Pathogen),
      x_obs = x_id + 0.25
    ) %>%
    dplyr::filter(is.finite(VE_obs), is.finite(VE_lower), is.finite(VE_upper))
  y_min <- 35
  y_max <- 75
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = dens$poly,
      ggplot2::aes(x = x, y = y, group = interaction(VE_Type, Pathogen), fill = "Theoretical distribution"),
      alpha = 0.55,
      color = NA
    ) +
    ggplot2::geom_path(
      data = dens$curve,
      ggplot2::aes(x = x, y = y, group = interaction(VE_Type, Pathogen)),
      color = "grey40",
      linewidth = 0.75,
      lineend = "round",
      linejoin = "round"
    ) +
    ggplot2::geom_segment(
      data = dens$edge,
      ggplot2::aes(x = x_edge, xend = x_edge, y = y_min, yend = y_max),
      color = "grey40",
      linewidth = 0.75,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = theory_summary,
      ggplot2::aes(x = x_box, xend = x_box, y = q025, yend = q975),
      linewidth = 0.80,
      color = "black"
    ) +
    ggplot2::geom_segment(
      data = theory_summary,
      ggplot2::aes(x = x_box - 0.045, xend = x_box + 0.045, y = q025, yend = q025),
      linewidth = 0.80,
      color = "black"
    ) +
    ggplot2::geom_segment(
      data = theory_summary,
      ggplot2::aes(x = x_box - 0.045, xend = x_box + 0.045, y = q975, yend = q975),
      linewidth = 0.80,
      color = "black"
    ) +
    ggplot2::geom_rect(
      data = theory_summary,
      ggplot2::aes(xmin = x_box - 0.070, xmax = x_box + 0.070, ymin = q25, ymax = q75),
      fill = "white",
      color = "black",
      linewidth = 0.80
    ) +
    ggplot2::geom_segment(
      data = theory_summary,
      ggplot2::aes(x = x_box - 0.070, xend = x_box + 0.070, y = median, yend = median),
      linewidth = 0.90,
      color = "black"
    ) +
    ggplot2::geom_point(
      data = theory_summary,
      ggplot2::aes(x = x_box, y = mean),
      shape = 23,
      fill = "white",
      color = "black",
      stroke = 0.75,
      size = 2.6
    ) +
    ggplot2::geom_errorbar(
      data = obs_long,
      ggplot2::aes(x = x_obs, ymin = VE_lower, ymax = VE_upper, color = "VE estimated from the logistic regression model"),
      width = 0.055,
      linewidth = 0.82
    ) +
    ggplot2::geom_point(
      data = obs_long,
      ggplot2::aes(x = x_obs, y = VE_obs, color = "VE estimated from the logistic regression model"),
      shape = 18,
      size = 2.8
    ) +
    ggplot2::facet_wrap(~VE_Type, ncol = 1, scales = "fixed", labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_continuous(breaks = seq_along(pathogen_levels), labels = pathogen_levels, limits = c(0.45, length(pathogen_levels) + 0.70), expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(breaks = seq(35, 75, by = 10), labels = scales::label_number(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = c("Theoretical distribution" = "#D9D4CC"), name = NULL, breaks = c("Theoretical distribution")) +
    ggplot2::scale_color_manual(values = c("VE estimated from the logistic regression model" = "#B2182B"), name = NULL, breaks = c("VE estimated from the logistic regression model")) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::labs(x = "Pathogen", y = "Vaccine effectiveness (%)", title = "D)") +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      text = ggplot2::element_text(),
      plot.title = ggplot2::element_text(size = 28, face = "bold", hjust = 0, margin = ggplot2::margin(b = 0)),
      plot.title.position = "plot",
      strip.text = ggplot2::element_text(size = 18, face = "plain"),
      axis.title = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 18, color = "black"),
      axis.text.y = ggplot2::element_text(size = 18, color = "black"),
      legend.position = "top",
      legend.justification = "left",
      legend.text = ggplot2::element_text(size = 18),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.85),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    )
  if (!is.null(qba_outdir)) {
    ggplot2::ggsave(file.path(qba_outdir, paste0("Validation_VE_posterior_vs_logistic_", scenario_name, ".png")), p, width = 10.5, height = 8.5, dpi = 300)
  }
  return(p)
}

make_combined_posterior_validation_plot <- function(posterior_wide,
                                                    posterior_draws,
                                                    obs_inputs,
                                                    qba_outdir,
                                                    scenario_name = "Main",
                                                    pathogen_levels = c("MP", "HRV", "HRSV", "HAdV", "HPIV", "ORV")) {
  p_shared_ve <- make_shared_parameter_plot(
    posterior_wide = posterior_wide,
    scenario_name = scenario_name,
    parameter = "VE_true",
    panel_label = NULL,
    x_label = "VE",
    y_label = "VE posterior value (%)",
    fill_color = "#3B6FB6"
  ) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 1, b = 2, l = 0)
    )
  p_shared_lambda <- make_shared_parameter_plot(
    posterior_wide = posterior_wide,
    scenario_name = scenario_name,
    parameter = "lambda_F",
    panel_label = NULL,
    x_label = expression(lambda[F]),
    y_label = expression(lambda[F]~" posterior value (%)"),
    fill_color = "#A64B3C"
  ) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 2, l = 1)
    )
  p_theta <- make_theta_raincloud(
    posterior_draws = posterior_draws,
    pathogen_levels = pathogen_levels,
    qba_outdir = qba_outdir,
    scenario_name = scenario_name
  ) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    )
  p_validation <- make_ve_validation_raincloud(
    posterior_draws = posterior_draws,
    obs_inputs = obs_inputs,
    qba_outdir = qba_outdir,
    scenario_name = scenario_name,
    pathogen_levels = pathogen_levels
  ) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    )
  make_panel <- function(plot,
                         label,
                         label_size = 18,
                         title_height = 0.28,
                         label_x = 0.02,
                         label_y = 0.55) {
    gridExtra::arrangeGrob(
      grid::textGrob(
        label,
        x = grid::unit(label_x, "npc"),
        y = grid::unit(label_y, "npc"),
        hjust = 0,
        vjust = 0.5,
        gp = grid::gpar(
          fontface = "bold",
          fontsize = label_size
        )
      ),
      plot,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(title_height, "in"),
        grid::unit(1, "null")
      )
    )
  }
  p_a <- make_panel(
    p_shared_ve,
    "A)",
    label_size = 18,
    title_height = 0.24,
    label_x = 0.02,
    label_y = 0.55
  )
  
  p_b <- make_panel(
    p_shared_lambda,
    "B)",
    label_size = 18,
    title_height = 0.24,
    label_x = 0.02,
    label_y = 0.55
  )
  
  p_c <- make_panel(
    p_theta,
    "C)",
    label_size = 18,
    title_height = 0.24,
    label_x = 0.02,
    label_y = 0.55
  )
  
  p_d <- make_panel(
    p_validation,
    "D)",
    label_size = 18,
    title_height = 0.24,
    label_x = 0.02,
    label_y = 0.55
  )
  
  p_top_left <- gridExtra::arrangeGrob(
    p_a,
    grid::nullGrob(), 
    p_b,
    ncol = 3,
    widths = c(1, 0.1, 1)
  )
  p_left <- gridExtra::arrangeGrob(
    p_top_left,
    grid::nullGrob(),
    p_c,
    ncol = 1,
    heights = c(0.74, 0.08, 1.48)
  )
  p_combined <- gridExtra::arrangeGrob(
    p_left,
    grid::nullGrob(), 
    p_d,
    ncol = 3,
    widths = c(1.55, 0.10, 2.55)
  )
  grDevices::png(
    filename = file.path(qba_outdir, paste0("Combined_Posterior_Parameters_and_VE_Validation_", scenario_name, ".png")),
    width = 17,
    height = 16,
    units = "in",
    res = 320
  )
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      x = 0.5,
      y = 0.5,
      width = 0.94,
      height = 0.94
    )
  )
  grid::grid.draw(p_combined)
  grid::popViewport()
  grDevices::dev.off()
  
  return(p_combined)
}

run_qba_prior_scenario <- function(scenario_row, obs_inputs) {
  scenario_row <- scenario_row[1, ]
  scenario_name <- as.character(scenario_row$Prior_Scenario[[1]])
  scenario_id <- as.integer(scenario_row$scenario_id[[1]])
  qba_outdir <- as.character(scenario_row$scenario_outdir[[1]])
  dir.create(qba_outdir, showWarnings = FALSE, recursive = TRUE)
  
  message("\n==============================")
  message("Running prior scenario: ", scenario_name)
  message("Output directory: ", normalizePath(qba_outdir, mustWork = FALSE))
  message("==============================")
  
  mcmc_settings <- base_mcmc_settings
  mcmc_settings$VE_prior_a <- scenario_row$VE_prior_a[[1]]
  mcmc_settings$VE_prior_b <- scenario_row$VE_prior_b[[1]]
  mcmc_settings$lambda_prior_a <- scenario_row$lambda_prior_a[[1]]
  mcmc_settings$lambda_prior_b <- scenario_row$lambda_prior_b[[1]]
  mcmc_settings$logtheta_prior_mean <- scenario_row$logtheta_prior_mean[[1]]
  mcmc_settings$logtheta_prior_sd <- scenario_row$logtheta_prior_sd[[1]]
  mcmc_settings$r_seed <- 2027000L + scenario_id
  
  readr::write_csv(scenario_row, file.path(qba_outdir, "QBA_Prior_Scenario_Setting.csv"))
  
  # Save actual MCMC settings used in this scenario for reproducibility.
  mcmc_run_settings <- tibble::tibble(
    Prior_Scenario = scenario_name,
    sampler = mcmc_settings$sampler,
    nuts_iter_warmup = mcmc_settings$nuts_iter_warmup,
    nuts_iter_sampling = mcmc_settings$nuts_iter_sampling,
    nuts_chains = mcmc_settings$nuts_chains,
    nuts_parallel_chains = mcmc_settings$nuts_parallel_chains,
    nuts_thin = mcmc_settings$nuts_thin,
    nuts_adapt_delta = mcmc_settings$nuts_adapt_delta,
    nuts_max_treedepth = mcmc_settings$nuts_max_treedepth,
    nuts_metric = mcmc_settings$nuts_metric,
    n_iter = mcmc_settings$n_iter,
    burn = mcmc_settings$burn,
    thin = mcmc_settings$thin,
    n_chains = mcmc_settings$n_chains,
    n_sim_gcomp = mcmc_settings$n_sim_gcomp,
    propagate_r_uncertainty = mcmc_settings$propagate_r_uncertainty,
    r_seed = mcmc_settings$r_seed
  )
  readr::write_csv(mcmc_run_settings, file.path(qba_outdir, "QBA_MCMC_Run_Settings.csv"))
  
  valid_obs <- obs_inputs %>%
    dplyr::filter(
      is.finite(beta_e), is.finite(se_e),
      is.finite(beta_1), is.finite(se_1),
      is.finite(beta_2), is.finite(se_2),
      is.finite(log_odds_u)
    )
  
  if (isTRUE(mcmc_settings$use_logR)) {
    valid_obs <- valid_obs %>%
      dplyr::filter(is.finite(logR), is.finite(se_logR))
  } else {
    valid_obs <- valid_obs %>%
      dplyr::mutate(
        logR = 0,
        se_logR = 1
      )
  }
  
  valid_obs <- valid_obs %>%
    dplyr::arrange(factor(Pathogen, levels = pathogen_levels))
  
  if (nrow(valid_obs) < 2) stop("Too few valid pathogens for joint MCMC.")
  
  # Starting values: shared VE from average VE2; lambda around fixed starting values;
  # theta varied by chain.
  VE2_start <- 1 - exp(stats::median(valid_obs$beta_2, na.rm = TRUE))
  VE_start <- ifelse(
    is.finite(VE2_start),
    clamp(VE2_start, 0.20, 0.80),
    mcmc_settings$VE_prior_a / (mcmc_settings$VE_prior_a + mcmc_settings$VE_prior_b)
  )
  
  inits <- list(
    make_nuts_init(VE_start, 0.08, 0.7, valid_obs),
    make_nuts_init(VE_start, 0.13, 1.0, valid_obs),
    make_nuts_init(VE_start, 0.20, 1.8, valid_obs),
    make_nuts_init(VE_start, 0.25, 2.2, valid_obs)
  )
  
  # If nuts_chains is changed, match the number of provided initial values.
  inits <- inits[seq_len(mcmc_settings$nuts_chains)]
  
  message("Running ", scenario_name, " using Stan/NUTS ...")
  nuts_run <- run_nuts_joint(
    valid_obs = valid_obs,
    settings = mcmc_settings,
    inits = inits,
    scenario_id = scenario_id,
    scenario_name = scenario_name
  )
  
  posterior_wide <- nuts_run$posterior_wide
  acceptance_table <- nuts_run$acceptance_table
  diagnostic_table <- nuts_run$diagnostic_table
  posterior_wide <- posterior_wide %>%
    dplyr::mutate(accept_stat_mean_over_chains = mean(acceptance_table$accept_stat_mean, na.rm = TRUE))
  
  posterior_draws <- add_predictions_long(posterior_wide, obs_inputs, mcmc_settings) %>%
    dplyr::mutate(Prior_Scenario = scenario_name)
  
  readr::write_csv(posterior_wide, file.path(qba_outdir, "QBA_Joint_Posterior_Draws_Wide.csv"))
  readr::write_csv(posterior_draws, file.path(qba_outdir, "QBA_Joint_Posterior_Draws_Long.csv"))
  readr::write_csv(acceptance_table, file.path(qba_outdir, "QBA_Joint_NUTS_Sampler_Diagnostics.csv"))
  
  # Posterior summaries.
  posterior_summary_shared <- posterior_wide %>%
    dplyr::summarise(
      Prior_Scenario = scenario_name,
      VE_prior = scenario_row$VE_prior_label[[1]],
      lambda_prior = scenario_row$lambda_prior_label[[1]],
      logtheta_prior = scenario_row$logtheta_prior_label[[1]],
      VE_true_mean = mean(VE_true, na.rm = TRUE),
      VE_true_median = median(VE_true, na.rm = TRUE),
      VE_true_CrI025 = stats::quantile(VE_true, 0.025, na.rm = TRUE),
      VE_true_CrI975 = stats::quantile(VE_true, 0.975, na.rm = TRUE),
      lambda_F_mean = mean(lambda_F, na.rm = TRUE),
      lambda_F_median = median(lambda_F, na.rm = TRUE),
      lambda_F_CrI025 = stats::quantile(lambda_F, 0.025, na.rm = TRUE),
      lambda_F_CrI975 = stats::quantile(lambda_F, 0.975, na.rm = TRUE),
      accept_stat_mean_over_chains = mean(accept_stat_mean_over_chains, na.rm = TRUE),
      .groups = "drop"
    )
  
  posterior_summary_pathogen <- posterior_draws %>%
    dplyr::group_by(Prior_Scenario, Pathogen) %>%
    dplyr::summarise(
      dplyr::across(
        c(log_theta, theta, phi, r_implied, VEe_theory, VE1_theory, VE2_theory),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          CrI025 = ~stats::quantile(.x, 0.025, na.rm = TRUE),
          CrI975 = ~stats::quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      scenario_row %>%
        dplyr::select(Prior_Scenario, VE_prior_label, lambda_prior_label, logtheta_prior_label),
      by = "Prior_Scenario"
    )
  
  readr::write_csv(posterior_summary_shared, file.path(qba_outdir, "QBA_Joint_Posterior_Summary_Shared.csv"))
  readr::write_csv(posterior_summary_pathogen, file.path(qba_outdir, "QBA_Joint_Posterior_Summary_Pathogen.csv"))
  
  # Parameter convergence diagnostics are generated directly from cmdstanr/posterior.
  # Rhat_rank_split uses rank-normalized split R-hat; ESS_bulk and ESS_tail are retained.
  readr::write_csv(diagnostic_table, file.path(qba_outdir, "QBA_Joint_Parameter_Convergence_Diagnostics.csv"))
  
  # Compact pathogen-level comparison table.
  posterior_for_compare <- posterior_draws %>%
    dplyr::group_by(Prior_Scenario, Pathogen) %>%
    dplyr::summarise(
      theta_median = median(theta, na.rm = TRUE),
      theta_CrI025 = stats::quantile(theta, 0.025, na.rm = TRUE),
      theta_CrI975 = stats::quantile(theta, 0.975, na.rm = TRUE),
      phi_median = median(phi, na.rm = TRUE),
      phi_CrI025 = stats::quantile(phi, 0.025, na.rm = TRUE),
      phi_CrI975 = stats::quantile(phi, 0.975, na.rm = TRUE),
      lambdaX_over_lambdaO_median = median(r_implied, na.rm = TRUE),
      lambdaX_over_lambdaO_CrI025 = stats::quantile(r_implied, 0.025, na.rm = TRUE),
      lambdaX_over_lambdaO_CrI975 = stats::quantile(r_implied, 0.975, na.rm = TRUE),
      VEe_theory_median = median(VEe_theory, na.rm = TRUE),
      VEe_theory_CrI025 = stats::quantile(VEe_theory, 0.025, na.rm = TRUE),
      VEe_theory_CrI975 = stats::quantile(VEe_theory, 0.975, na.rm = TRUE),
      VE1_theory_median = median(VE1_theory, na.rm = TRUE),
      VE1_theory_CrI025 = stats::quantile(VE1_theory, 0.025, na.rm = TRUE),
      VE1_theory_CrI975 = stats::quantile(VE1_theory, 0.975, na.rm = TRUE),
      VE2_theory_median = median(VE2_theory, na.rm = TRUE),
      VE2_theory_CrI025 = stats::quantile(VE2_theory, 0.025, na.rm = TRUE),
      VE2_theory_CrI975 = stats::quantile(VE2_theory, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      VE_true_median_shared = posterior_summary_shared$VE_true_median[1],
      VE_true_CrI025_shared = posterior_summary_shared$VE_true_CrI025[1],
      VE_true_CrI975_shared = posterior_summary_shared$VE_true_CrI975[1],
      lambda_F_median_shared = posterior_summary_shared$lambda_F_median[1],
      lambda_F_CrI025_shared = posterior_summary_shared$lambda_F_CrI025[1],
      lambda_F_CrI975_shared = posterior_summary_shared$lambda_F_CrI975[1]
    )
  
  compare_table <- obs_inputs %>%
    dplyr::select(Pathogen, VEe_obs, VEe_obs_lower, VEe_obs_upper,
                  VE1_obs, VE1_obs_lower, VE1_obs_upper,
                  VE2_obs, VE2_obs_lower, VE2_obs_upper,
                  n_flu_neg, n_x_pos) %>%
    dplyr::mutate(Prior_Scenario = scenario_name) %>%
    dplyr::left_join(posterior_for_compare, by = c("Prior_Scenario", "Pathogen")) %>%
    dplyr::mutate(
      VEe_obs_percent = 100 * VEe_obs,
      VE1_obs_percent = 100 * VE1_obs,
      VE2_obs_percent = 100 * VE2_obs,
      VEe_theory_percent = 100 * VEe_theory_median,
      VE1_theory_percent = 100 * VE1_theory_median,
      VE2_theory_percent = 100 * VE2_theory_median,
      VEe_diff_pp = 100 * (VEe_theory_median - VEe_obs),
      VE1_diff_pp = 100 * (VE1_theory_median - VE1_obs),
      VE2_diff_pp = 100 * (VE2_theory_median - VE2_obs),
      VEe_obs_point_in_theory_CrI = VEe_obs >= VEe_theory_CrI025 & VEe_obs <= VEe_theory_CrI975,
      VE1_obs_point_in_theory_CrI = VE1_obs >= VE1_theory_CrI025 & VE1_obs <= VE1_theory_CrI975,
      VE2_obs_point_in_theory_CrI = VE2_obs >= VE2_theory_CrI025 & VE2_obs <= VE2_theory_CrI975,
      Direction = dplyr::case_when(
        theta_median > 1 & theta_CrI025 > 1 ~ "Synergistic evidence: theta > 1",
        theta_median < 1 & theta_CrI975 < 1 ~ "Antagonistic evidence: theta < 1",
        TRUE ~ "Uncertain around theta = 1"
      ),
      VEe_obs_CI_overlap_theory_CrI = VEe_obs_upper >= VEe_theory_CrI025 & VEe_obs_lower <= VEe_theory_CrI975,
      VE1_obs_CI_overlap_theory_CrI = VE1_obs_upper >= VE1_theory_CrI025 & VE1_obs_lower <= VE1_theory_CrI975,
      VE2_obs_CI_overlap_theory_CrI = VE2_obs_upper >= VE2_theory_CrI025 & VE2_obs_lower <= VE2_theory_CrI975
    ) %>%
    dplyr::left_join(
      scenario_row %>% dplyr::select(Prior_Scenario, VE_prior_label, lambda_prior_label, logtheta_prior_label),
      by = "Prior_Scenario"
    )
  
  # Main validation plot: posterior predictive VE distributions vs adjusted logistic VE estimates.
  p_validation <- make_ve_validation_raincloud(
    posterior_draws = posterior_draws,
    obs_inputs = obs_inputs,
    qba_outdir = qba_outdir,
    scenario_name = scenario_name,
    pathogen_levels = pathogen_levels
  )
  
  readr::write_csv(compare_table, file.path(qba_outdir, "QBA_Joint_Theory_vs_Observed_Comparison.csv"))
  
  validation_long <- dplyr::bind_rows(
    compare_table %>%
      dplyr::transmute(
        Prior_Scenario, Pathogen,
        Estimand = "VEe: Flu- overall",
        Observed = 100 * VEe_obs,
        Observed_Lower = 100 * VEe_obs_lower,
        Observed_Upper = 100 * VEe_obs_upper,
        Theory_Median = 100 * VEe_theory_median,
        Theory_CrI025 = 100 * VEe_theory_CrI025,
        Theory_CrI975 = 100 * VEe_theory_CrI975,
        Difference_pp = VEe_diff_pp,
        Observed_Point_In_Theory_CrI = VEe_obs_point_in_theory_CrI,
        Observed_CI_Overlap_Theory_CrI = VEe_obs_CI_overlap_theory_CrI
      ),
    compare_table %>%
      dplyr::transmute(
        Prior_Scenario, Pathogen,
        Estimand = "VE1: Flu-, X+",
        Observed = 100 * VE1_obs,
        Observed_Lower = 100 * VE1_obs_lower,
        Observed_Upper = 100 * VE1_obs_upper,
        Theory_Median = 100 * VE1_theory_median,
        Theory_CrI025 = 100 * VE1_theory_CrI025,
        Theory_CrI975 = 100 * VE1_theory_CrI975,
        Difference_pp = VE1_diff_pp,
        Observed_Point_In_Theory_CrI = VE1_obs_point_in_theory_CrI,
        Observed_CI_Overlap_Theory_CrI = VE1_obs_CI_overlap_theory_CrI
      ),
    compare_table %>%
      dplyr::transmute(
        Prior_Scenario, Pathogen,
        Estimand = "VE2: Flu-, X-",
        Observed = 100 * VE2_obs,
        Observed_Lower = 100 * VE2_obs_lower,
        Observed_Upper = 100 * VE2_obs_upper,
        Theory_Median = 100 * VE2_theory_median,
        Theory_CrI025 = 100 * VE2_theory_CrI025,
        Theory_CrI975 = 100 * VE2_theory_CrI975,
        Difference_pp = VE2_diff_pp,
        Observed_Point_In_Theory_CrI = VE2_obs_point_in_theory_CrI,
        Observed_CI_Overlap_Theory_CrI = VE2_obs_CI_overlap_theory_CrI
      )
  ) %>%
    dplyr::mutate(
      Estimand = factor(Estimand, levels = c("VEe: Flu- overall", "VE1: Flu-, X+", "VE2: Flu-, X-"))
    )
  
  readr::write_csv(validation_long, file.path(qba_outdir, "QBA_Joint_VEe_VE1_VE2_Validation_Long.csv"))
  
  # Scenario-specific figures.
  p_theta_raincloud <- make_theta_raincloud(
    posterior_draws = posterior_draws,
    pathogen_levels = valid_obs$Pathogen,
    qba_outdir = qba_outdir,
    scenario_name = scenario_name
  )
  p_shared_ve <- make_shared_parameter_plot(
    posterior_wide = posterior_wide,
    scenario_name = scenario_name,
    parameter = "VE_true",
    panel_label = "A)",
    y_label = "VE posterior value (%)",
    fill_color = "#3B6FB6"
  )
  p_shared_lambda <- make_shared_parameter_plot(
    posterior_wide = posterior_wide,
    scenario_name = scenario_name,
    parameter = "lambda_F",
    panel_label = "B)",
    y_label = expression(lambda[F] ~" posterior value (%)"),
    fill_color = "#A64B3C"
  )
  p_shared <- patchwork::wrap_plots(p_shared_ve, p_shared_lambda, ncol = 2, widths = c(1, 1))
  ggplot2::ggsave(file.path(qba_outdir, "Joint_Posterior_shared_VE_lambda.png"), p_shared, width = 7.5, height = 4.2, dpi = 320)
  p_validation <- make_ve_validation_raincloud(
    posterior_draws = posterior_draws,
    obs_inputs = obs_inputs,
    qba_outdir = qba_outdir,
    scenario_name = scenario_name,
    pathogen_levels = valid_obs$Pathogen
  )
  ggplot2::ggsave(file.path(qba_outdir, "Joint_Unified_Validation_VEe_VE1_VE2.png"), p_validation, width = 12.5, height = 8.8, dpi = 320)
  p_combined <- make_combined_posterior_validation_plot(
    posterior_wide = posterior_wide,
    posterior_draws = posterior_draws,
    obs_inputs = obs_inputs,
    qba_outdir = qba_outdir,
    scenario_name = scenario_name,
    pathogen_levels = valid_obs$Pathogen
  )
  
  trace_cols <- c("VE_true", "lambda_F", paste0("theta_", valid_obs$Pathogen))
  
  trace_label_levels <- c(
    "VE",
    "\u03BBF",
    paste0("\u03B8-", valid_obs$Pathogen)
  )
  
  trace_long <- posterior_wide %>%
    dplyr::arrange(chain, postwarmup_iter) %>%
    dplyr::group_by(chain) %>%
    dplyr::mutate(draw_index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(draw_index, postwarmup_iter, chain, dplyr::all_of(trace_cols)) %>%
    tidyr::pivot_longer(
      cols = -c(draw_index, postwarmup_iter, chain),
      names_to = "Parameter",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      chain = factor(chain),
      Parameter_label = dplyr::case_when(
        Parameter == "VE_true" ~ "VE",
        Parameter == "lambda_F" ~ "\u03BBF",
        stringr::str_detect(Parameter, "^theta_") ~
          paste0("\u03B8-", stringr::str_remove(Parameter, "^theta_")),
        TRUE ~ Parameter
      ),
      Parameter_label = dplyr::case_when(
        Parameter == "VE_true" ~ "VE",
        Parameter == "lambda_F" ~ "lambda[F]",
        stringr::str_detect(Parameter, "^theta_") ~
          paste0("theta*'-", stringr::str_remove(Parameter, "^theta_"), "'"),
        TRUE ~ Parameter
      )
    )
  
  p_trace_theta <- ggplot2::ggplot(
    trace_long,
    ggplot2::aes(x = postwarmup_iter, y = Value, group = chain, color = chain)
  ) +
    ggplot2::geom_line(alpha = 0.6, linewidth = 0.25) +
    ggplot2::facet_wrap(
      ~ Parameter_label,
      scales = "free_y",
      ncol = 2,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::labs(
      x = "iteration",
      y = "Posterior draw",
      color = "Chain"
    ) +
    ggplot2::theme_bw()

  p_trace_theta_for_legend <- p_trace_theta +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 12),
      legend.text  = ggplot2::element_text(size = 11)
    )
  
  trace_common_legend <- cowplot::get_legend(p_trace_theta_for_legend)
  
  p_trace_theta_nolegend <- p_trace_theta +
    ggplot2::theme(
      legend.position = "none"
    )
  
  ggplot2::ggsave(
    file.path(qba_outdir, "Joint_MCMC_Trace_theta_shared_parameters_nolegend.png"),
    p_trace_theta_nolegend,
    width = 9,
    height = 7,
    dpi = 400
  )

  if (scenario_name == "Main") {
    ggplot2::ggsave(
      file.path(dirname(qba_outdir), "MCMC_Trace_Common_Legend.png"),
      cowplot::ggdraw(trace_common_legend),
      width = 4.5,
      height = 0.6,
      dpi = 400,
      bg = "white"
    )
  }

  ggplot2::ggsave(
    file.path(qba_outdir, "Joint_MCMC_Trace_theta_shared_parameters.png"),
    p_trace_theta,
    width = 10,
    height = 8,
    dpi = 400
  )
  
  list(
    shared = posterior_summary_shared,
    pathogen = posterior_summary_pathogen,
    compare = compare_table,
    validation = validation_long,
    diagnostics = diagnostic_table,
    acceptance = acceptance_table,
    run_settings = mcmc_run_settings
  )
}

# ------------------------------------------------------------------------------
# 9.10 Run all prior scenarios and create combined outputs
# ------------------------------------------------------------------------------
scenario_results <- vector("list", nrow(prior_scenarios))
for (ii in seq_len(nrow(prior_scenarios))) {
  scenario_results[[ii]] <- run_qba_prior_scenario(prior_scenarios[ii, ], obs_inputs)
}
names(scenario_results) <- prior_scenarios$Prior_Scenario

combined_shared <- dplyr::bind_rows(lapply(scenario_results, `[[`, "shared"))
combined_pathogen <- dplyr::bind_rows(lapply(scenario_results, `[[`, "pathogen"))
combined_compare <- dplyr::bind_rows(lapply(scenario_results, `[[`, "compare"))
combined_validation <- dplyr::bind_rows(lapply(scenario_results, `[[`, "validation"))
combined_diagnostics <- dplyr::bind_rows(lapply(scenario_results, `[[`, "diagnostics"))
combined_acceptance <- dplyr::bind_rows(lapply(scenario_results, `[[`, "acceptance"))
combined_run_settings <- dplyr::bind_rows(lapply(scenario_results, `[[`, "run_settings"))

add_mcmc_panel_label <- function(file,
                                 tag,
                                 scenario_name = NULL,
                                 show_scenario_name = TRUE,
                                 tag_height = 180,
                                 tag_x = 80,
                                 tag_y = 45,
                                 tag_size = 60) {
  
  img <- magick::image_read(file)
  info <- magick::image_info(img)
  
  img_with_margin <- magick::image_extent(
    img,
    geometry = paste0(info$width, "x", info$height + tag_height),
    gravity = "south",
    color = "white"
  )
  
  label_text <- if (show_scenario_name && !is.null(scenario_name)) {
    paste0(tag, " ", scenario_name)
  } else {
    tag
  }
  
  img_with_label <- magick::image_annotate(
    img_with_margin,
    text = label_text,
    gravity = "northwest",
    location = paste0("+", tag_x, "+", tag_y),
    size = tag_size,
    color = "black",
    weight = 700
  )
  
  return(img_with_label)
}


combine_mcmc_trace_panels_3col_common_legend <- function(prior_scenarios,
                                                         base_qba_outdir,
                                                         output_file = "Figure_MCMC_Trace_All_Prior_Scenarios_3col_CommonLegend.png",
                                                         show_scenario_name = TRUE,
                                                         tag_height = 180,
                                                         tag_x = 80,
                                                         tag_y = 45,
                                                         tag_size = 60) {
  
  trace_files <- file.path(
    prior_scenarios$scenario_outdir,
    "Joint_MCMC_Trace_theta_shared_parameters_nolegend.png"
  )
  
  legend_file <- file.path(
    base_qba_outdir,
    "MCMC_Trace_Common_Legend.png"
  )
  
  tags <- paste0(LETTERS[seq_along(trace_files)], ")")
  scenario_names <- prior_scenarios$Prior_Scenario
  
  keep <- file.exists(trace_files)
  
  if (!all(keep)) {
    warning(
      "Some no-legend trace plot files were not found: ",
      paste(trace_files[!keep], collapse = "; ")
    )
  }
  
  trace_files <- trace_files[keep]
  tags <- tags[keep]
  scenario_names <- scenario_names[keep]
  
  if (length(trace_files) == 0) {
    warning("No no-legend MCMC trace plot files found.")
    return(invisible(NULL))
  }
  
  panel_imgs <- purrr::pmap(
    list(trace_files, tags, scenario_names),
    function(file, tag, scenario_name) {
      add_mcmc_panel_label(
        file = file,
        tag = tag,
        scenario_name = scenario_name,
        show_scenario_name = show_scenario_name,
        tag_height = tag_height,
        tag_x = tag_x,
        tag_y = tag_y,
        tag_size = tag_size
      )
    }
  )
  
  panel_info <- purrr::map_dfr(
    panel_imgs,
    ~ as.data.frame(magick::image_info(.x))
  )
  
  target_width  <- max(panel_info$width)
  target_height <- max(panel_info$height)
  
  panel_imgs <- purrr::map(
    panel_imgs,
    ~ magick::image_extent(
      .x,
      geometry = paste0(target_width, "x", target_height),
      gravity = "center",
      color = "white"
    )
  )
  
  n_col <- 3
  n_blank <- (n_col - length(panel_imgs) %% n_col) %% n_col
  
  if (n_blank > 0) {
    blank_panels <- replicate(
      n_blank,
      magick::image_blank(
        width = target_width,
        height = target_height,
        color = "white"
      ),
      simplify = FALSE
    )
    panel_imgs <- c(panel_imgs, blank_panels)
  }
  
  row_imgs <- list()
  
  for (i in seq(1, length(panel_imgs), by = n_col)) {
    row_imgs[[length(row_imgs) + 1]] <- magick::image_append(
      do.call(c, panel_imgs[i:(i + n_col - 1)]),
      stack = FALSE
    )
  }
  
  combined_img <- magick::image_append(
    do.call(c, row_imgs),
    stack = TRUE
  )
  
  if (file.exists(legend_file)) {
    
    combined_info <- magick::image_info(combined_img)
    
    legend_img <- magick::image_read(legend_file)
    legend_info <- magick::image_info(legend_img)
    
    legend_canvas <- magick::image_blank(
      width = combined_info$width,
      height = legend_info$height + 120,
      color = "white"
    )
    
    legend_canvas <- magick::image_composite(
      legend_canvas,
      legend_img,
      gravity = "center"
    )
    
    final_img <- magick::image_append(
      c(combined_img, legend_canvas),
      stack = TRUE
    )
    
  } else {
    warning("Common legend file not found: ", legend_file)
    final_img <- combined_img
  }
  
  out_path <- file.path(base_qba_outdir, output_file)
  
  magick::image_write(
    final_img,
    path = out_path,
    format = "png",
    density = "600x600"
  )
  
  message("Saved combined MCMC trace figure with common legend: ", out_path)
}


combine_mcmc_trace_panels_3col_common_legend(
  prior_scenarios = prior_scenarios,
  base_qba_outdir = base_qba_outdir,
  output_file = "Figure_MCMC_Trace_All_Prior_Scenarios_3col_CommonLegend.png",
  show_scenario_name = FALSE,
  tag_height = 180,
  tag_x = 80,
  tag_y = 45,
  tag_size = 60
)

readr::write_csv(combined_shared, file.path(base_qba_outdir, "Sensitivity_Combined_Shared_Parameters.csv"))
readr::write_csv(combined_pathogen, file.path(base_qba_outdir, "Sensitivity_Combined_Pathogen_Summaries.csv"))
readr::write_csv(combined_compare, file.path(base_qba_outdir, "Sensitivity_Combined_Theory_vs_Observed_Comparison.csv"))
readr::write_csv(combined_validation, file.path(base_qba_outdir, "Sensitivity_Combined_VEe_VE1_VE2_Validation_Long.csv"))
readr::write_csv(combined_diagnostics, file.path(base_qba_outdir, "Sensitivity_Combined_Parameter_Convergence_Diagnostics.csv"))
readr::write_csv(combined_acceptance, file.path(base_qba_outdir, "Sensitivity_Combined_NUTS_Sampler_Diagnostics.csv"))
readr::write_csv(combined_run_settings, file.path(base_qba_outdir, "Sensitivity_Combined_MCMC_Run_Settings.csv"))

# Compact table for reporting prior sensitivity of theta direction and validation pattern.
sensitivity_direction_table <- combined_compare %>%
  dplyr::select(
    Prior_Scenario, VE_prior_label, lambda_prior_label, logtheta_prior_label,
    Pathogen,
    theta_median, theta_CrI025, theta_CrI975, Direction,
    # Empirical_R, R_theory_median, R_theory_CrI025, R_theory_CrI975, R_obs_point_in_theory_CrI,
    VE1_obs_percent, VE1_theory_percent, VE1_obs_point_in_theory_CrI, VE1_obs_CI_overlap_theory_CrI,
    VE2_obs_percent, VE2_theory_percent, VE2_obs_point_in_theory_CrI, VE2_obs_CI_overlap_theory_CrI,
    VEe_obs_percent, VEe_theory_percent, VEe_obs_point_in_theory_CrI, VEe_obs_CI_overlap_theory_CrI,
    VE_true_median_shared, lambda_F_median_shared
  ) %>%
  dplyr::arrange(factor(Prior_Scenario, levels = prior_scenarios$Prior_Scenario), Pathogen)

readr::write_csv(sensitivity_direction_table, file.path(base_qba_outdir, "Sensitivity_Theta_Direction_and_Validation_Summary.csv"))
print(sensitivity_direction_table)

# Optional combined plots across scenarios.
p_theta_sensitivity <- combined_compare %>%
  ggplot2::ggplot(ggplot2::aes(x = theta_median, y = Pathogen)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = theta_CrI025, xmax = theta_CrI975), height = 0.18) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~Prior_Scenario) +
  ggplot2::labs(
    x = expression(theta~"median and 95% CrI"),
    y = NULL,
    title = "Prior sensitivity analysis for pathogen-specific theta"
  ) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(base_qba_outdir, "Sensitivity_Theta_by_Prior_Scenario.png"), p_theta_sensitivity, width = 12, height = 7, dpi = 300)

p_shared_sensitivity <- combined_shared %>%
  dplyr::select(Prior_Scenario, VE_true_median, VE_true_CrI025, VE_true_CrI975,
                lambda_F_median, lambda_F_CrI025, lambda_F_CrI975) %>%
  tidyr::pivot_longer(
    cols = -Prior_Scenario,
    names_to = c("Parameter", ".value"),
    names_pattern = "(VE_true|lambda_F)_(median|CrI025|CrI975)"
  ) %>%
  dplyr::mutate(Parameter = dplyr::recode(Parameter, VE_true = "Shared true VE", lambda_F = "Shared lambda_F")) %>%
  ggplot2::ggplot(ggplot2::aes(x = median, y = Prior_Scenario)) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = CrI025, xmax = CrI975), height = 0.18) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~Parameter, scales = "free_x") +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = NULL, title = "Prior sensitivity analysis for shared parameters") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(base_qba_outdir, "Sensitivity_Shared_Parameters_by_Prior_Scenario.png"), p_shared_sensitivity, width = 10, height = 5, dpi = 300)

message("Prior sensitivity analysis completed. Combined outputs saved in: ", normalizePath(base_qba_outdir))

