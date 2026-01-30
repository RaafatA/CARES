"""
Machine learning optimization of wheat–berseem intercropping
under different water qualities and irrigation levels
------------------------------------------------------------

Algorithms:
- Random Forest
- XGBoost

Key features:
- Season-stratified bootstrapping
- Scenario-based optimization
- SHAP interpretability

Author: [Your Name]
"""

# =========================
# 1. IMPORT LIBRARIES
# =========================

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.utils import resample

from xgboost import XGBRegressor

import shap
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


# =========================
# 2. LOAD DATA
# =========================
# Your dataframe must already be cleaned and in long format
# Each row = plot × season

df = pd.read_csv("field_experiment_data.csv")

# -------------------------
# REQUIRED COLUMNS
# -------------------------
# Season
# CroppingSystem   -> Mono, Mixed
# WaterQuality     -> Fresh, Fish, Inorganic, Mixture
# FieldCapacity    -> 100, 50
# Soil_*           -> soil variables
# Weather_*        -> weather summaries
# AppliedWater
# TARGET VARIABLE  -> e.g. SystemYield / LER / WaterProductivity


# =========================
# 3. DEFINE TARGET & FEATURES
# =========================

TARGET = "SystemYield"  # change to LER or WaterProductivity

categorical_features = [
    "CroppingSystem",
    "WaterQuality",
    "FieldCapacity"
]

numeric_features = [
    c for c in df.columns
    if c.startswith("Soil_") or c.startswith("Weather_")
] + ["AppliedWater"]

X = df[categorical_features + numeric_features]
y = df[TARGET]
seasons = df["Season"]


# =========================
# 4. PREPROCESSING PIPELINE
# =========================

preprocessor = ColumnTransformer(
    transformers=[
        ("cat", OneHotEncoder(drop="first"), categorical_features),
        ("num", "passthrough", numeric_features)
    ]
)

# =========================
# 5. MODELS
# =========================

rf_model = RandomForestRegressor(
    n_estimators=600,
    min_samples_leaf=3,
    random_state=42,
    n_jobs=-1
)

xgb_model = XGBRegressor(
    n_estimators=500,
    learning_rate=0.05,
    max_depth=4,
    subsample=0.8,
    colsample_bytree=0.8,
    objective="reg:squarederror",
    random_state=42
)

rf_pipeline = Pipeline([
    ("prep", preprocessor),
    ("model", rf_model)
])

xgb_pipeline = Pipeline([
    ("prep", preprocessor),
    ("model", xgb_model)
])


# =========================
# 6. SEASON-STRATIFIED BOOTSTRAP FUNCTION
# =========================

def season_bootstrap(data, season_col="Season"):
    """
    Resample plots WITHIN each season to preserve
    seasonal structure.
    """
    boot = []
    for s in data[season_col].unique():
        season_data = data[data[season_col] == s]
        boot.append(
            resample(season_data, replace=True, n_samples=len(season_data))
        )
    return pd.concat(boot)


# =========================
# 7. SCENARIO GRID (OPTIMIZATION SPACE)
# =========================

scenario_grid = (
    df[categorical_features]
    .drop_duplicates()
    .reset_index(drop=True)
)

# Attach representative (mean) soil & weather conditions
for col in numeric_features:
    scenario_grid[col] = df[col].mean()


# =========================
# 8. BOOTSTRAP LOOP
# =========================

B = 1000  # number of bootstrap iterations

rf_preds = []
xgb_preds = []

for b in range(B):

    # ---- Bootstrap sampling ----
    boot_df = season_bootstrap(df)

    Xb = boot_df[categorical_features + numeric_features]
    yb = boot_df[TARGET]

    # ---- Fit models ----
    rf_pipeline.fit(Xb, yb)
    xgb_pipeline.fit(Xb, yb)

    # ---- Predict scenarios ----
    rf_preds.append(rf_pipeline.predict(scenario_grid))
    xgb_preds.append(xgb_pipeline.predict(scenario_grid))

# Convert to arrays
rf_preds = np.array(rf_preds)
xgb_preds = np.array(xgb_preds)


# =========================
# 9. DECISION TABLE
# =========================

decision_table = scenario_grid.copy()

# ---- Random Forest summaries ----
decision_table["RF_Mean"] = rf_preds.mean(axis=0)
decision_table["RF_Lower95"] = np.percentile(rf_preds, 2.5, axis=0)
decision_table["RF_Upper95"] = np.percentile(rf_preds, 97.5, axis=0)

# ---- XGBoost summaries ----
decision_table["XGB_Mean"] = xgb_preds.mean(axis=0)
decision_table["XGB_Lower95"] = np.percentile(xgb_preds, 2.5, axis=0)
decision_table["XGB_Upper95"] = np.percentile(xgb_preds, 97.5, axis=0)

# ---- Average rank across models ----
decision_table["MeanRank"] = (
    decision_table[["RF_Mean", "XGB_Mean"]]
    .mean(axis=1)
    .rank(ascending=False)
)

decision_table.sort_values("MeanRank", inplace=True)

decision_table.to_csv("optimization_decision_table.csv", index=False)


# =========================
# 10. SHAP INTERPRETATION
# =========================

# ---- Random Forest SHAP ----
rf_pipeline.fit(X, y)

X_transformed = rf_pipeline.named_steps["prep"].transform(X)
rf_model_fitted = rf_pipeline.named_steps["model"]

rf_explainer = shap.TreeExplainer(rf_model_fitted)
rf_shap = rf_explainer.shap_values(X_transformed)

feature_names = rf_pipeline.named_steps["prep"].get_feature_names_out()

shap.summary_plot(
    rf_shap,
    X_transformed,
    feature_names=feature_names,
    show=False
)
plt.title("SHAP Summary – Random Forest")
plt.tight_layout()
plt.show()


# ---- XGBoost SHAP ----
xgb_pipeline.fit(X, y)

X_transformed = xgb_pipeline.named_steps["prep"].transform(X)
xgb_model_fitted = xgb_pipeline.named_steps["model"]

xgb_explainer = shap.TreeExplainer(xgb_model_fitted)
xgb_shap = xgb_explainer.shap_values(X_transformed)

shap.summary_plot(
    xgb_shap,
    X_transformed,
    feature_names=feature_names,
    show=False
)
plt.title("SHAP Summary – XGBoost")
plt.tight_layout()
plt.show()


# =========================
# 11. FINAL MESSAGE
# =========================

print("Analysis complete.")
print("Decision table saved as 'optimization_decision_table.csv'")
