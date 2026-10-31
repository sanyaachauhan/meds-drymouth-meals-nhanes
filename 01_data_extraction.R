# =====================================================================
# Project: Meds, Dry Mouth, and Meals (NHANES 2011-2018)
# Author: Sanya Chauhan | MSc Health Sciences, University of Lucerne
# 
# Purpose:
# Exploring if taking multiple dry-mouth medications (anticholinergic burden) 
# correlates with worse self-reported dental health and lower daily caloric intake. 
# Combines 4 NHANES cycles, applies CDC complex survey weights, runs GLMs, 
# and uses XGBoost + SHAP plots to check for non-linear tipping points.
# =====================================================================

library(tidyverse)
library(nhanesA)
library(survey)
library(xgboost)
library(caret)
library(shapviz)
library(Ckmeans.1d.dp)
library(ggplot2)

# ---------------------------------------------------------------------
# 1. Load Demographics (Cycles G to J)
# ---------------------------------------------------------------------

# NHANES cycles (2011-2018) to get a large sample size
demo_g <- nhanes("DEMO_G")
demo_h <- nhanes("DEMO_H")
demo_i <- nhanes("DEMO_I")
demo_j <- nhanes("DEMO_J")

clean_household_cols <- function(df) {
  cols_to_fix <- c("DMDHHSIZ", "DMDFMSIZ", "DMDHHSZA", "DMDHHSZB", "DMDHHSZC", "DMDHHSZD", "DMDHHSZE")
  for(col in cols_to_fix) {
    if(col %in% names(df)) {
      df[[col]] <- parse_number(as.character(df[[col]]))
    }
  }
  return(df)
}

demo_g <- clean_household_cols(demo_g)
demo_h <- clean_household_cols(demo_h)
demo_i <- clean_household_cols(demo_i)
demo_j <- clean_household_cols(demo_j)

df_demo <- bind_rows(demo_g, demo_h, demo_i, demo_j)


# ---------------------------------------------------------------------
# 2. Loading Oral Health Questionnaires
# ---------------------------------------------------------------------

ohq_g <- nhanes("OHQ_G")
ohq_h <- nhanes("OHQ_H")
ohq_i <- nhanes("OHQ_I")
ohq_j <- nhanes("OHQ_J")

# Standardizing columns to character before merging (except respondent ID / SEQN)
prep_for_bind <- function(df) {
  df %>% mutate(across(-any_of("SEQN"), as.character))
}

ohq_g <- prep_for_bind(ohq_g)
ohq_h <- prep_for_bind(ohq_h)
ohq_i <- prep_for_bind(ohq_i)
ohq_j <- prep_for_bind(ohq_j)

df_oral <- bind_rows(ohq_g, ohq_h, ohq_i, ohq_j)

# Joining oral health onto demographics by respondent ID
master_data <- left_join(df_demo, df_oral, by = "SEQN")


# ---------------------------------------------------------------------
# 3. Prescription Medications & ACB Scoring
# ---------------------------------------------------------------------

rx_g <- nhanes("RXQ_RX_G")
rx_h <- nhanes("RXQ_RX_H")
rx_i <- nhanes("RXQ_RX_I")
rx_j <- nhanes("RXQ_RX_J")

df_rx <- bind_rows(rx_g, rx_h, rx_i, rx_j)

# Score drugs based on Anticholinergic Cognitive Burden (ACB) / xerostomia risk:
# 3 = strong anticholinergics (high dry-mouth risk)
# 1 = moderate anticholinergics
# 0 = minimal/no anticholinergic effect
df_rx_scored <- df_rx %>%
  mutate(drug_clean = str_to_lower(RXDDRUG)) %>% 
  mutate(acb_points = case_when(
    str_detect(drug_clean, "oxybutynin|diphenhydramine|amitriptyline|darifenacin|solifenacin|tolterodine|imipramine|nortriptyline") ~ 3,
    str_detect(drug_clean, "citalopram|paroxetine|cetirizine|loratadine|furosemide|metoprolol|atenolol") ~ 1,
    TRUE ~ 0
  ))

# Calculate cumulative ACB score per person
df_acb_summary <- df_rx_scored %>%
  group_by(SEQN) %>%
  summarize(cumulative_ACB_score = sum(acb_points, na.rm = TRUE))

# Merge into main dataset; set missing scores to 0 (no reported meds)
master_data <- left_join(master_data, df_acb_summary, by = "SEQN") %>%
  mutate(cumulative_ACB_score = replace_na(cumulative_ACB_score, 0))


# ---------------------------------------------------------------------
# 4. Dentition & Dietary Data
# ---------------------------------------------------------------------

dent_g <- nhanes("OHXDEN_G")
dent_h <- nhanes("OHXDEN_H")
dent_i <- nhanes("OHXDEN_I")
dent_j <- nhanes("OHXDEN_J")

dent_g <- dent_g %>% mutate(across(-any_of("SEQN"), as.character))
dent_h <- dent_h %>% mutate(across(-any_of("SEQN"), as.character))
dent_i <- dent_i %>% mutate(across(-any_of("SEQN"), as.character))
dent_j <- dent_j %>% mutate(across(-any_of("SEQN"), as.character))
df_dent <- bind_rows(dent_g, dent_h, dent_i, dent_j)
master_data <- left_join(master_data, df_dent, by = "SEQN")

diet_g <- nhanes("DR1TOT_G")
diet_h <- nhanes("DR1TOT_H")
diet_i <- nhanes("DR1TOT_I")
diet_j <- nhanes("DR1TOT_J")

diet_g <- diet_g %>% mutate(across(-any_of("SEQN"), as.character))
diet_h <- diet_h %>% mutate(across(-any_of("SEQN"), as.character))
diet_i <- diet_i %>% mutate(across(-any_of("SEQN"), as.character))
diet_j <- diet_j %>% mutate(across(-any_of("SEQN"), as.character))
df_diet <- bind_rows(diet_g, diet_h, diet_i, diet_j)
master_data <- left_join(master_data, df_diet, by = "SEQN")


# ---------------------------------------------------------------------
# 5. Variable Prep & Survey Weighting
# ---------------------------------------------------------------------

# Mapping self-rated dental health to 1-5 scale (1 = Excellent, 5 = Poor)
# Scale 2 year exam weights by dividing by 4, as 4 cycles
master_data_model <- master_data %>%
  mutate(
    RIDAGEYR = as.numeric(RIDAGEYR),
    RIAGENDR = as.factor(RIAGENDR),
    scaled_wt = as.numeric(WTMEC2YR) / 4,
    
    dental_health_rating = case_when(
      OHQ845 == "Excellent" ~ 1,
      OHQ845 == "Very good" ~ 2,
      OHQ845 == "Good" ~ 3,
      OHQ845 == "Fair" ~ 4,
      OHQ845 == "Poor" ~ 5,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(scaled_wt), !is.na(SDMVPSU), !is.na(SDMVSTRA),
    !is.na(dental_health_rating), !is.na(cumulative_ACB_score),
    !is.na(RIDAGEYR), !is.na(RIAGENDR)
  )

# Setting up NHANES survey design object using PSU, Strata, and scaled weights
nhanes_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~scaled_wt, nest = TRUE, data = master_data_model
)


# ---------------------------------------------------------------------
# 6. Baseline Survey Regressions (GLMs)
# ---------------------------------------------------------------------

# Testing association between ACB score & self rated oral health
model_oral_impact <- svyglm(
  dental_health_rating ~ cumulative_ACB_score + RIDAGEYR + RIAGENDR, 
  design = nhanes_design, family = gaussian()
)

# Testing association between ACB score & daily caloric intake
master_data_model$caloric_intake <- parse_number(master_data_model$DR1TKCAL)
master_data_model_nutr <- subset(master_data_model, !is.na(caloric_intake))

nhanes_design_nutr <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~scaled_wt, nest = TRUE, data = master_data_model_nutr
)

model_nutrition <- svyglm(
  caloric_intake ~ cumulative_ACB_score + RIDAGEYR + RIAGENDR,
  design = nhanes_design_nutr, family = gaussian()
)


# ---------------------------------------------------------------------
# 7. Machine Learning: XGBoost with Survey Weights
# ---------------------------------------------------------------------

# Cleaning up column names for modeling matrix
ml_data_clean <- master_data_model_nutr %>%
  mutate(
    acb_score      = cumulative_ACB_score,
    age            = RIDAGEYR,
    gender_male    = ifelse(RIAGENDR == "Male" | RIAGENDR == 1, 1, 0),
    oral_health    = dental_health_rating,
    caloric_intake = caloric_intake,
    scaled_wt      = scaled_wt
  ) %>%
  select(acb_score, age, gender_male, oral_health, caloric_intake, scaled_wt) %>%
  drop_na()

# Test Split (80/20)
set.seed(42)
train_idx <- createDataPartition(ml_data_clean$caloric_intake, p = 0.8, list = FALSE)

train_df <- ml_data_clean[train_idx, ]
test_df  <- ml_data_clean[-train_idx, ]

train_x <- model.matrix(caloric_intake ~ acb_score + age + gender_male + oral_health - 1, data = train_df)
train_y <- train_df$caloric_intake
train_w <- train_df$scaled_wt

test_x  <- model.matrix(caloric_intake ~ acb_score + age + gender_male + oral_health - 1, data = test_df)
test_y  <- test_df$caloric_intake
test_w  <- test_df$scaled_wt

# Passing survey weights directly into DMatrix object
dtrain_w <- xgb.DMatrix(data = train_x, label = train_y, weight = train_w)
dtest_w  <- xgb.DMatrix(data = test_x, label = test_y, weight = test_w)

params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  max_depth = 4,
  eta = 0.05,
  subsample = 0.8
)

set.seed(42)
xgb_clean_model <- xgb.train(
  params = params, data = dtrain_w, nrounds = 150, 
  evals = list(train = dtrain_w, test = dtest_w), early_stopping_rounds = 15, verbose = 0
)

preds <- predict(xgb_clean_model, dtest_w)


# ---------------------------------------------------------------------
# 8. SHAP Plots & Saving Figures
# ---------------------------------------------------------------------

shp_clean <- shapviz(xgb_clean_model, X_pred = test_x, X = as.data.frame(test_x))
colnames(shp_clean$S) <- c("Cumulative ACB Score", "Age (Years)", "Gender (Male)", "Oral Health Status")
colnames(shp_clean$X) <- c("Cumulative ACB Score", "Age (Years)", "Gender (Male)", "Oral Health Status")

# Figure 1: Variable Importance Bar Plot
sv_importance(shp_clean, kind = "bar", fill = "steelblue") +
  ggtitle("Global Feature Importance: Drivers of Caloric Intake") +
  labs(x = "Mean Absolute SHAP Value", y = "Feature") +
  theme_minimal()
ggsave("Figure1_Global_Feature_Hierarchy.png", width = 8, height = 6, dpi = 300)

# Figure 2: ACB vs Oral Health Interaction
sv_dependence(shp_clean, v = "Cumulative ACB Score", color_var = "Oral Health Status", alpha = 0.5, size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "line", color = "black", linewidth = 1, group = 1) +
  ggtitle("ACB Impact Stratified by Self-Rated Oral Health") +
  labs(x = "Cumulative ACB Score", y = "SHAP Value", color = "Oral Health") +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_minimal()
ggsave("Figure2_Cascade_Interaction.png", width = 8, height = 6, dpi = 300)

# Figure 3: Tipping Point Plot
sv_dependence(shp_clean, v = "Cumulative ACB Score", alpha = 0.3, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "line", color = "darkblue", linewidth = 1.3, group = 1) +
  stat_summary(fun = mean, geom = "point", color = "darkblue", size = 2.5) +
  ggtitle("Caloric Drop Threshold by ACB Score") +
  labs(x = "Cumulative ACB Score", y = "SHAP Value") +
  theme_minimal()
ggsave("Figure3_Tipping_Point.png", width = 8, height = 6, dpi = 300)


# ---------------------------------------------------------------------
# 9. Model Metrics & Population Risk Estimates
# ---------------------------------------------------------------------

# Calculating survey weighted evaluation metrics on test set
weighted_mean_y <- sum(test_y * test_w) / sum(test_w)
weighted_r2 <- 1 - (sum(test_w * (test_y - preds)^2) / sum(test_w * (test_y - weighted_mean_y)^2))
weighted_rmse <- sqrt(sum(test_w * (test_y - preds)^2) / sum(test_w))
weighted_mae <- sum(test_w * abs(test_y - preds)) / sum(test_w)

cat("R2:", round(weighted_r2, 4), "| RMSE:", round(weighted_rmse, 2), "| MAE:", round(weighted_mae, 2), "\n")

# Calculating population percentage  (High ACB + Poor Oral Health)
risk_exposure_summary <- master_data_model_nutr %>%
  summarize(
    total_weighted_pop = sum(scaled_wt),
    pct_high_acb = sum(scaled_wt[cumulative_ACB_score >= 2]) / sum(scaled_wt) * 100,
    pct_poor_oral = sum(scaled_wt[dental_health_rating >= 4]) / sum(scaled_wt) * 100,
    pct_double_hit = sum(scaled_wt[cumulative_ACB_score >= 2 & dental_health_rating >= 4]) / sum(scaled_wt) * 100
  )
print(risk_exposure_summary)

# Exporting & Saving
write.csv(risk_exposure_summary, "Population_Risk_Exposure_Summary.csv", row.names = FALSE)
saveRDS(xgb_clean_model, "saved_xgb_model.rds")