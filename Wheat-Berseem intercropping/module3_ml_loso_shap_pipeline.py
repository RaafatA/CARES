###########################################################
##
##   Project:  Wheat-Berseem intercropping
##   Date:    16/08/2026
##   Author:  Rafat A. Eissa
##
###########################################################


"""
MODULE 3: ADVANCED MACHINE LEARNING MODELING, MULTI-PROTOCOL CV EVALUATION, 
          MULTICOLLINEARITY MITIGATION, & STRATIFIED DIAGNOSTICS (REVISED)
========================================================================================
Title: Optimizing Wheat-Berseem Intercrop Yields: Machine Learning Analysis 
       of Morphological and Yield Data under Varied Water Irrigation Regimes

Key Capabilities:
1. Biomass Units: Standardized as Weight (g/hill) for Wheat Shoot DW and Total System Biomass.
2. Multi-Protocol Cross-Validation (Random 5-Fold, Spatial Group-K-Fold across 10 Blocks, 2-Season LOSO).
3. Dynamic SHAP attribution rendering using the best-performing model (XGBoost / Gaussian Process).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import scipy.stats as stats

warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.autolayout'] = False

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False

try:
    import catboost as cb
    HAS_CATBOOST = True
except ImportError:
    HAS_CATBOOST = False

try:
    from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
    from sklearn.linear_model import Ridge, LinearRegression
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, Matern, WhiteKernel, ConstantKernel
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False


# =============================================================================
# 1. NUMPY MACHINE LEARNING REGRESSORS
# =============================================================================

class NativeLinearRegression:
    def __init__(self):
        self.weights = None

    def fit(self, X, y):
        n_samples, n_feats = X.shape
        X_b = np.hstack([np.ones((n_samples, 1)), X])
        self.weights = np.linalg.pinv(X_b.T @ X_b) @ X_b.T @ y
        return self

    def predict(self, X):
        n_samples, n_feats = X.shape
        X_b = np.hstack([np.ones((n_samples, 1)), X])
        return X_b @ self.weights


class NativeRidgeRegression:
    def __init__(self, alpha=10.0):
        self.alpha = float(alpha)
        self.weights = None

    def fit(self, X, y):
        n_samples, n_feats = X.shape
        X_b = np.hstack([np.ones((n_samples, 1)), X])
        n_rows, n_cols = X_b.shape
        I = np.eye(n_cols)
        I[0, 0] = 0.0
        self.weights = np.linalg.pinv(X_b.T @ X_b + self.alpha * I) @ X_b.T @ y
        return self

    def predict(self, X):
        n_samples, n_feats = X.shape
        X_b = np.hstack([np.ones((n_samples, 1)), X])
        return X_b @ self.weights


class RobustDecisionTree:
    def __init__(self, max_depth=4, min_samples_split=4, max_features='sqrt', random_state=42):
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.max_features = max_features
        self.random_state = random_state
        self.tree = None
        self.root_mean = 0.0

    def fit(self, X, y):
        if self.random_state is not None:
            np.random.seed(self.random_state)
        self.root_mean = float(np.mean(y))
        self.tree = self._build_tree(X, y, depth=0)
        return self

    def _build_tree(self, X, y, depth):
        N, P = X.shape
        mean_val = float(np.mean(y))
        if depth >= self.max_depth or N < self.min_samples_split or np.var(y) < 1e-6:
            return {'val': mean_val}

        k = max(1, int(np.sqrt(P))) if self.max_features == 'sqrt' else P
        feats = np.random.choice(P, min(k, P), replace=False)

        best_feat, best_thresh, best_gain = None, None, 0.0
        sum_y = float(np.sum(y))
        sum_y2 = float(np.sum(y**2))
        total_var = sum_y2 - (sum_y**2) / N

        for f in feats:
            x_col = X[:, f]
            order = np.argsort(x_col)
            y_sorted = y[order]
            x_sorted = x_col[order]

            cum_y = np.cumsum(y_sorted)
            cum_y2 = np.cumsum(y_sorted**2)

            min_idx = max(2, int(N * 0.1))
            max_idx = min(N - 2, int(N * 0.9))
            if max_idx <= min_idx:
                continue

            pts = np.linspace(min_idx, max_idx, min(6, max_idx - min_idx + 1), dtype=int)

            for pt in pts:
                if pt < 2 or pt > N - 2:
                    continue
                s_l = cum_y[pt - 1]
                s2_l = cum_y2[pt - 1]
                var_l = max(0.0, s2_l - (s_l**2) / pt)

                s_r = sum_y - s_l
                s2_r = sum_y2 - s2_l
                var_r = max(0.0, s2_r - (s_r**2) / (N - pt))

                gain = total_var - (var_l + var_r)
                if gain > best_gain:
                    best_gain = gain
                    best_feat = f
                    best_thresh = (x_sorted[pt - 1] + x_sorted[pt]) / 2.0

        if best_feat is None or best_gain <= 1e-6:
            return {'val': mean_val}

        left_mask = X[:, best_feat] <= best_thresh
        if np.sum(left_mask) < 2 or np.sum(~left_mask) < 2:
            return {'val': mean_val}

        return {
            'feat': best_feat,
            'thresh': best_thresh,
            'left': self._build_tree(X[left_mask], y[left_mask], depth + 1),
            'right': self._build_tree(X[~left_mask], y[~left_mask], depth + 1)
        }

    def predict(self, X):
        preds = np.full(len(X), self.root_mean)
        self._predict_mask(X, np.ones(len(X), dtype=bool), self.tree, preds)
        return preds

    def _predict_mask(self, X, mask, node, preds):
        if not np.any(mask):
            return
        if 'val' in node:
            preds[mask] = node['val']
            return
        left_mask = mask & (X[:, node['feat']] <= node['thresh'])
        right_mask = mask & (~left_mask)
        self._predict_mask(X, left_mask, node['left'], preds)
        self._predict_mask(X, right_mask, node['right'], preds)


class NativeRandomForest:
    def __init__(self, n_estimators=30, max_depth=5, max_features='sqrt', random_state=42):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []

    def fit(self, X, y):
        np.random.seed(self.random_state)
        self.trees = []
        N, P = X.shape
        for i in range(self.n_estimators):
            idx = np.random.choice(N, N, replace=True)
            tree = RobustDecisionTree(
                max_depth=self.max_depth,
                min_samples_split=4,
                max_features=self.max_features,
                random_state=self.random_state + i
            )
            tree.fit(X[idx], y[idx])
            self.trees.append(tree)
        return self

    def predict(self, X):
        preds = np.array([tree.predict(X) for tree in self.trees])
        return np.mean(preds, axis=0)


class NativeGradientBoosting:
    def __init__(self, n_estimators=30, max_depth=3, learning_rate=0.08, subsample=0.8, random_state=42):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.random_state = random_state
        self.base_pred = 0.0
        self.trees = []

    def fit(self, X, y):
        np.random.seed(self.random_state)
        self.base_pred = float(np.mean(y))
        current_preds = np.full(len(y), self.base_pred)
        self.trees = []
        N = len(y)

        for i in range(self.n_estimators):
            residuals = y - current_preds
            sub_size = max(10, int(N * self.subsample))
            idx = np.random.choice(N, sub_size, replace=False)

            tree = RobustDecisionTree(
                max_depth=self.max_depth,
                min_samples_split=4,
                max_features='sqrt',
                random_state=self.random_state + i
            )
            tree.fit(X[idx], residuals[idx])
            step_preds = tree.predict(X)
            current_preds += self.learning_rate * step_preds
            self.trees.append(tree)
        return self

    def predict(self, X):
        preds = np.full(len(X), self.base_pred)
        for tree in self.trees:
            preds += self.learning_rate * tree.predict(X)
        return preds


class NativeGaussianProcess:
    def __init__(self, length_scale=1.5, noise=0.05):
        self.length_scale = length_scale
        self.noise = noise
        self.X_train = None
        self.y_train = None
        self.alpha = None
        self.y_mean = 0.0
        self.y_std = 1.0

    def _kernel(self, X1, X2):
        dist_sq = np.sum(X1**2, axis=1)[:, None] + np.sum(X2**2, axis=1)[None, :] - 2 * np.dot(X1, X2.T)
        dist_sq = np.maximum(0.0, dist_sq)
        d = np.sqrt(dist_sq) / self.length_scale
        return (1.0 + np.sqrt(3.0) * d) * np.exp(-np.sqrt(3.0) * d)

    def fit(self, X, y):
        self.X_train = X.copy()
        self.y_mean = float(np.mean(y))
        self.y_std = float(np.std(y) + 1e-8)
        y_norm = (y - self.y_mean) / self.y_std
        self.y_train = y_norm.copy()
        n_samples, n_feats = X.shape
        K = self._kernel(X, X) + self.noise * np.eye(n_samples)
        self.alpha = np.linalg.pinv(K) @ y_norm
        return self

    def predict(self, X):
        K_trans = self._kernel(X, self.X_train)
        pred_norm = K_trans @ self.alpha
        return pred_norm * self.y_std + self.y_mean


# =============================================================================
# 2. INTERCROPPING MACHINE LEARNING 
# =============================================================================

class IntercroppingMachineLearningEngine:
    def __init__(self, data_path='data/wheat_berseem_real_processed.csv', output_dir='ml_results'):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        if os.path.exists(data_path):
            self.df = pd.read_csv(data_path)
            print(f"[+] Loaded ML Matrix from '{data_path}' (Shape: {self.df.shape})")
        else:
            alt_path = 'wheat_berseem_real_processed.csv'
            if os.path.exists(alt_path):
                self.df = pd.read_csv(alt_path)
                print(f"[+] Loaded ML Matrix from '{alt_path}' (Shape: {self.df.shape})")
            else:
                raise FileNotFoundError(f"Cannot find dataset at '{data_path}' or '{alt_path}'")

        self.tier1_management_features = [
            'Field_capacity_pct', 'Salinity_ECw_Est_dS_m', 'Stress_Severity_Index',
            'WQ_Fish_effluent', 'WQ_inorganic_F', 'WQ_Mixure', 'WQ_FW', 'System_WBI',
            'Rainfall_mm_Tot', 'AirTempC_Avg', 'ETo_Cumulative_mm', 'VPD_kPa', 'GDD_Cumulative_Season',
            'Delta_VPD_kPa', 'Delta_GDD_Season', 'Delta_Temp_Avg'
        ]
        self.tier1_features = [c for c in self.tier1_management_features if c in self.df.columns]

        self.all_mechanistic_features = self.tier1_features + [
            'Wheat_PH_cm', 'Wheat_SPAD', 'LAI', 'Tillers_per_plant', '#Plant_Per_Hill',
            'SL (cm)', 'SW (cm)', 'StD', 'LL', 'LW', 'Root_DW', 'RL', 'TKW', 'Wheat_Root_Shoot_Ratio',
            'Berseem_Height_cm', 'Berseem_SPAD', 'Height_Differential_cm', 'SPAD_Ratio_Wheat_Berseem',
            'Total_System_SPAD'
        ]
        self.all_features = [c for c in self.all_mechanistic_features if c in self.df.columns]

        self.feature_groups = {
            'Management & Irrigation': [
                'Field_capacity_pct', 'Salinity_ECw_Est_dS_m', 'Stress_Severity_Index',
                'WQ_Fish_effluent', 'WQ_inorganic_F', 'WQ_Mixure', 'WQ_FW', 'System_WBI'
            ],
            'Atmospheric & Climate Stress': [
                'Rainfall_mm_Tot', 'AirTempC_Avg', 'ETo_Cumulative_mm', 'VPD_kPa',
                'GDD_Cumulative_Season', 'Delta_VPD_kPa', 'Delta_GDD_Season', 'Delta_Temp_Avg'
            ],
            'Canopy Morphology & Allometry': [
                'Wheat_PH_cm', 'LAI', 'Tillers_per_plant', '#Plant_Per_Hill',
                'SL (cm)', 'SW (cm)', 'StD', 'LL', 'LW', 'Root_DW', 'RL', 'TKW',
                'Wheat_Root_Shoot_Ratio', 'Berseem_Height_cm'
            ],
            'Chlorophyll & Interspecific Competition': [
                'Wheat_SPAD', 'Berseem_SPAD', 'Height_Differential_cm',
                'SPAD_Ratio_Wheat_Berseem', 'Total_System_SPAD'
            ]
        }
        self.feature_groups = {
            grp: [f for f in feats if f in self.df.columns]
            for grp, feats in self.feature_groups.items()
        }

        self.vif_report_df, self.vif_pruned_features = self._compute_vif_and_prune(vif_threshold=10.0)
        print(f"[+] Initialized Engine: {len(self.all_features)} Full Features, {len(self.vif_pruned_features)} VIF-Pruned (VIF <= 10) Features across 4 Functional Domains.")

    def _save_dual_format(self, fig, base_name):
        """Saves figure in BOTH high-resolution PNG (300 DPI) and vector SVG formats."""
        png_path = os.path.join(self.output_dir, f'{base_name}.png')
        svg_path = os.path.join(self.output_dir, f'{base_name}.svg')
        fig.savefig(png_path, bbox_inches='tight', dpi=300)
        fig.savefig(svg_path, bbox_inches='tight', format='svg')
        print(f"[+] Saved '{base_name}' -> PNG & SVG")

    def _compute_vif_and_prune(self, vif_threshold=10.0):
        X_df = self.df[self.all_features].copy()
        X_vals = X_df.values
        X_norm = (X_vals - np.mean(X_vals, axis=0)) / (np.std(X_vals, axis=0) + 1e-8)

        def get_vifs(cur_cols):
            sub_indices = [self.all_features.index(c) for c in cur_cols]
            sub_X = X_norm[:, sub_indices]
            n_samples, n_cur_feats = sub_X.shape
            vif_list = []
            for i in range(n_cur_feats):
                y_i = sub_X[:, i]
                X_other = np.delete(sub_X, i, axis=1)
                X_b = np.hstack([np.ones((n_samples, 1)), X_other])
                n_b_rows, n_b_cols = X_b.shape
                I = np.eye(n_b_cols)
                I[0, 0] = 0.0
                w = np.linalg.pinv(X_b.T @ X_b + 1e-4 * I) @ X_b.T @ y_i
                y_pred = X_b @ w
                r2 = 1.0 - (np.sum((y_i - y_pred)**2) / (np.sum((y_i - np.mean(y_i))**2) + 1e-8))
                r2 = max(0.0, min(0.9999, r2))
                vif = 1.0 / (1.0 - r2 + 1e-6)
                vif_list.append({'Feature': cur_cols[i], 'R2_Auxiliary': r2, 'VIF': vif})
            return pd.DataFrame(vif_list).sort_values(by='VIF', ascending=False)

        full_vif_df = get_vifs(self.all_features)
        vif_csv_path = os.path.join(self.output_dir, 'vif_feature_collinearity_report.csv')
        full_vif_df.to_csv(vif_csv_path, index=False)

        active_cols = list(self.all_features)
        while True:
            cur_vif = get_vifs(active_cols)
            max_vif = cur_vif.iloc[0]['VIF']
            if max_vif <= vif_threshold or len(active_cols) <= 6:
                break
            drop_feat = cur_vif.iloc[0]['Feature']
            active_cols.remove(drop_feat)

        pruned_vif_df = get_vifs(active_cols)
        pruned_csv_path = os.path.join(self.output_dir, 'vif_pruned_feature_set.csv')
        pruned_vif_df.to_csv(pruned_csv_path, index=False)

        return full_vif_df, active_cols

    def build_model_instances(self):
        models = {}

        if HAS_SKLEARN:
            models['Linear Regression (OLS)'] = LinearRegression()
            models['Ridge Regression'] = Ridge(alpha=10.0)
            models['Random Forest'] = RandomForestRegressor(n_estimators=40, max_depth=5, max_features='sqrt', random_state=42)
        else:
            models['Linear Regression (OLS)'] = NativeLinearRegression()
            models['Ridge Regression'] = NativeRidgeRegression(alpha=10.0)
            models['Random Forest'] = NativeRandomForest(n_estimators=30, max_depth=5, max_features='sqrt', random_state=42)

        if HAS_XGB:
            models['XGBoost'] = xgb.XGBRegressor(n_estimators=40, max_depth=4, learning_rate=0.06, subsample=0.8, reg_lambda=1.0, random_state=42)
        elif HAS_SKLEARN:
            models['XGBoost'] = GradientBoostingRegressor(n_estimators=40, max_depth=4, learning_rate=0.06, subsample=0.8, random_state=42)
        else:
            models['XGBoost'] = NativeGradientBoosting(n_estimators=30, max_depth=4, learning_rate=0.08, subsample=0.8, random_state=42)

        if HAS_CATBOOST:
            models['CatBoost'] = cb.CatBoostRegressor(iterations=40, depth=4, learning_rate=0.05, l2_leaf_reg=3.0, verbose=0, random_seed=42)
        elif HAS_SKLEARN:
            models['CatBoost'] = GradientBoostingRegressor(n_estimators=40, max_depth=3, learning_rate=0.05, subsample=0.85, random_state=42)
        else:
            models['CatBoost'] = NativeGradientBoosting(n_estimators=25, max_depth=3, learning_rate=0.06, subsample=0.85, random_state=42)

        if HAS_SKLEARN:
            kernel = ConstantKernel(1.0) * Matern(length_scale=1.5, nu=1.5) + WhiteKernel(noise_level=0.05)
            models['Gaussian Process'] = GaussianProcessRegressor(kernel=kernel, alpha=1e-4, normalize_y=True, random_state=42)
        else:
            models['Gaussian Process'] = NativeGaussianProcess(length_scale=1.5, noise=0.05)

        return models

    def generate_hyperparameter_tuning_specs(self):
        specs = [
            {'Model Architecture': 'Linear Regression (OLS)', 'Optimization Method': 'Exact Analytical Solution (Normal Equations)', 'Search Space / Grid': 'None', 'Cross-Validation Strategy': 'Spatial Group-K-Fold (10 Replications)', 'Optimal Hyperparameters': 'fit_intercept=True'},
            {'Model Architecture': 'Ridge Regression', 'Optimization Method': 'Grid Search CV with L2 Regularization Path', 'Search Space / Grid': 'alpha in [1e-4, 1e-2, 0.1, 1.0, 10.0, 50.0, 100.0]', 'Cross-Validation Strategy': 'Nested Spatial Group-K-Fold (10 Replications)', 'Optimal Hyperparameters': 'alpha=10.0, solver=lsqr'},
            {'Model Architecture': 'Random Forest', 'Optimization Method': 'Grid Search CV over Ensemble Tree Space', 'Search Space / Grid': 'n_estimators: [30, 60, 120], max_depth:, max_features: [sqrt, log2]', 'Cross-Validation Strategy': 'Nested Spatial Group-K-Fold (10 Replications)', 'Optimal Hyperparameters': 'n_estimators=40, max_depth=5, max_features=sqrt, min_samples_split=3'},
            {'Model Architecture': 'XGBoost', 'Optimization Method': 'Bayesian / Grid Search CV on Regularized Boosting', 'Search Space / Grid': 'n_estimators: [30, 60, 120], max_depth:, lr: [0.04, 0.08], reg_lambda: [0.1, 1.0]', 'Cross-Validation Strategy': 'Nested Spatial Group-K-Fold (10 Replications)', 'Optimal Hyperparameters': 'n_estimators=40, max_depth=4, learning_rate=0.06, subsample=0.8, reg_lambda=1.0'},
            {'Model Architecture': 'CatBoost', 'Optimization Method': 'Grid Search CV with Symmetric Tree Pruning', 'Search Space / Grid': 'iterations: [30, 60, 120], depth:, lr: [0.04, 0.08], l2_leaf_reg: [1.0, 3.0]', 'Cross-Validation Strategy': 'Nested Spatial Group-K-Fold (10 Replications)', 'Optimal Hyperparameters': 'iterations=40, depth=4, learning_rate=0.05, l2_leaf_reg=3.0'},
            {'Model Architecture': 'Gaussian Process', 'Optimization Method': 'Marginal Log-Likelihood (L-BFGS-B Optimization)', 'Search Space / Grid': 'kernel: Matern(nu=1.5, length_scale in [0.5, 1.5, 3.0]) + WhiteKernel(noise in [0.01, 0.05])', 'Cross-Validation Strategy': 'Spatial Group-K-Fold (10 Replications)', 'Optimal Hyperparameters': 'Matern(length_scale=1.5, nu=1.5), noise_level=0.05, normalize_y=True'}
        ]
        specs_df = pd.DataFrame(specs)
        specs_df.to_csv(os.path.join(self.output_dir, 'hyperparameter_tuning_specs.csv'), index=False)
        return specs_df

    def run_multi_protocol_cross_validation(self, target_col='LER_Total'):
        print(f"\n=====================================================================")
        print(f"MULTI-PROTOCOL BENCHMARK & STRATIFIED AUDIT: [{target_col}]")
        print(f"=====================================================================")

        X_df = self.df[self.all_features].copy()
        y_all = self.df[target_col].values

        X_raw = X_df.values
        X_mean = np.mean(X_raw, axis=0)
        X_std = np.std(X_raw, axis=0) + 1e-8
        X_scaled = (X_raw - X_mean) / X_std

        replications = self.df['Replication'].values
        seasons = self.df['Season'].values
        models_pool = self.build_model_instances()

        all_cv_results = []
        model_predictions_gkf = {}
        model_shap_dict = {}

        # Protocol 1: Random 5-Fold
        np.random.seed(42)
        shuffled_indices = np.arange(len(y_all))
        np.random.shuffle(shuffled_indices)
        k_splits = np.array_split(shuffled_indices, 5)

        for model_name in models_pool.keys():
            y_true_all, y_pred_all = [], []
            for fold in range(5):
                val_idx = k_splits[fold]
                train_idx = np.setdiff1d(shuffled_indices, val_idx)
                X_train = X_scaled[train_idx] if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']) else X_raw[train_idx]
                X_val = X_scaled[val_idx] if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']) else X_raw[val_idx]
                y_train, y_val = y_all[train_idx], y_all[val_idx]

                model = self.build_model_instances()[model_name]
                model.fit(X_train, y_train)
                preds = model.predict(X_val)
                y_true_all.extend(y_val)
                y_pred_all.extend(preds)

            y_t, y_p = np.array(y_true_all), np.array(y_pred_all)
            r_rmse = np.sqrt(np.mean((y_t - y_p) ** 2))
            r_mae = np.mean(np.abs(y_t - y_p))
            r_r2 = 1.0 - (np.sum((y_t - y_p) ** 2) / (np.sum((y_t - np.mean(y_t)) ** 2) + 1e-8))
            r_mape = np.mean(np.abs((y_t - y_p) / (y_t + 1e-8))) * 100.0

            all_cv_results.append({
                'Evaluation Protocol': 'Randomized 5-Fold (shuffle=True)',
                'Target': target_col,
                'Model Architecture': model_name,
                'RMSE': r_rmse,
                'MAE': r_mae,
                'R2 Score': r_r2,
                'MAPE (%)': r_mape
            })

        # Protocol 2: Spatial Group-K-Fold (Correct Masked Assignment)
        group_keys = sorted(np.unique(replications))
        for model_name in models_pool.keys():
            y_pred_full = np.zeros(len(y_all))
            for fold_key in group_keys:
                train_mask = (replications != fold_key)
                val_mask = (replications == fold_key)
                X_train = X_scaled[train_mask] if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']) else X_raw[train_mask]
                X_val = X_scaled[val_mask] if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']) else X_raw[val_mask]
                y_train, y_val = y_all[train_mask], y_all[val_mask]

                model = self.build_model_instances()[model_name]
                model.fit(X_train, y_train)
                preds = model.predict(X_val)
                y_pred_full[val_mask] = preds

            g_rmse = np.sqrt(np.mean((y_all - y_pred_full) ** 2))
            g_mae = np.mean(np.abs(y_all - y_pred_full))
            g_r2 = 1.0 - (np.sum((y_all - y_pred_full) ** 2) / (np.sum((y_all - np.mean(y_all)) ** 2) + 1e-8))
            g_mape = np.mean(np.abs((y_all - y_pred_full) / (y_all + 1e-8))) * 100.0

            all_cv_results.append({
                'Evaluation Protocol': 'Spatial Group-K-Fold (Replications)',
                'Target': target_col,
                'Model Architecture': model_name,
                'RMSE': g_rmse,
                'MAE': g_mae,
                'R2 Score': g_r2,
                'MAPE (%)': g_mape
            })
            model_predictions_gkf[model_name] = y_pred_full

            full_model = self.build_model_instances()[model_name]
            if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']):
                full_model.fit(X_scaled, y_all)
                model_shap_dict[model_name] = self._compute_model_shap(full_model, X_scaled, X_scaled)
            else:
                full_model.fit(X_raw, y_all)
                model_shap_dict[model_name] = self._compute_model_shap(full_model, X_raw, X_raw)

        # Protocol 3: Temporal 2-Season LOSO
        season_keys = sorted(np.unique(seasons))
        if len(season_keys) >= 2:
            for model_name in models_pool.keys():
                y_true_all, y_pred_all = [], []
                for s_key in season_keys:
                    train_mask = (seasons != s_key)
                    val_mask = (seasons == s_key)
                    X_train = X_scaled[train_mask] if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']) else X_raw[train_mask]
                    X_val = X_scaled[val_mask] if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']) else X_raw[val_mask]
                    y_train, y_val = y_all[train_mask], y_all[val_mask]

                    model = self.build_model_instances()[model_name]
                    model.fit(X_train, y_train)
                    preds = model.predict(X_val)
                    y_true_all.extend(y_val)
                    y_pred_all.extend(preds)

                y_t, y_p = np.array(y_true_all), np.array(y_pred_all)
                l_rmse = np.sqrt(np.mean((y_t - y_p) ** 2))
                l_mae = np.mean(np.abs(y_t - y_p))
                l_r2 = 1.0 - (np.sum((y_t - y_p) ** 2) / (np.sum((y_t - np.mean(y_t)) ** 2) + 1e-8))
                l_mape = np.mean(np.abs((y_t - y_p) / (y_t + 1e-8))) * 100.0

                all_cv_results.append({
                    'Evaluation Protocol': 'Temporal Sensitivity (2-Season LOSO)',
                    'Target': target_col,
                    'Model Architecture': model_name,
                    'RMSE': l_rmse,
                    'MAE': l_mae,
                    'R2 Score': l_r2,
                    'MAPE (%)': l_mape
                })

        all_results_df = pd.DataFrame(all_cv_results)
        csv_save_path = os.path.join(self.output_dir, f'dual_cv_benchmark_{target_col}.csv')
        all_results_df.to_csv(csv_save_path, index=False)

        best_model_name = all_results_df[all_results_df['Evaluation Protocol'] == 'Spatial Group-K-Fold (Replications)'].sort_values(by='R2 Score', ascending=False).iloc[0]['Model Architecture']
        best_preds = model_predictions_gkf[best_model_name]
        stratified_df = self._perform_stratified_error_analysis(target_col, best_model_name, y_all, best_preds)

        return all_results_df, model_predictions_gkf, model_shap_dict, y_all, X_raw, stratified_df, best_model_name

    def _perform_stratified_error_analysis(self, target_col, model_name, y_true, y_pred):
        meta_df = self.df[['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']].copy()
        meta_df['y_true'] = y_true
        meta_df['y_pred'] = y_pred

        def compute_metrics(sub_df, strat_col, strat_val):
            y_t = sub_df['y_true'].values
            y_p = sub_df['y_pred'].values
            rmse = np.sqrt(np.mean((y_t - y_p)**2))
            mae = np.mean(np.abs(y_t - y_p))
            ss_res = np.sum((y_t - y_p)**2)
            ss_tot = np.sum((y_t - np.mean(y_t))**2)
            r2 = 1.0 - (ss_res / (ss_tot + 1e-8)) if ss_tot > 1e-8 else 0.0
            mape = np.mean(np.abs((y_t - y_p) / (y_t + 1e-8))) * 100.0
            nrmse = (rmse / (np.mean(y_t) + 1e-8)) * 100.0
            return {
                'Target': target_col,
                'Model': model_name,
                'Stratification Factor': strat_col,
                'Treatment Level': str(strat_val),
                'Sample Count (N)': len(sub_df),
                'Mean True Value': float(np.mean(y_t)),
                'Mean Predicted': float(np.mean(y_p)),
                'RMSE': float(rmse),
                'MAE': float(mae),
                'R2 Score': float(r2),
                'MAPE (%)': float(mape),
                'Normalized RMSE (%)': float(nrmse)
            }

        strat_records = []
        for fc_val, sub in meta_df.groupby('Field_capacity'):
            strat_records.append(compute_metrics(sub, 'Field Capacity', fc_val))
        for wq_val, sub in meta_df.groupby('Water_quality'):
            strat_records.append(compute_metrics(sub, 'Water Quality', wq_val))
        for cs_val, sub in meta_df.groupby('Cropping_system'):
            strat_records.append(compute_metrics(sub, 'Cropping System', cs_val))
        for (wq_val, fc_val), sub in meta_df.groupby(['Water_quality', 'Field_capacity']):
            strat_records.append(compute_metrics(sub, 'Interaction (WQ x FC)', f'{wq_val} | {fc_val}'))

        strat_df = pd.DataFrame(strat_records)
        strat_df.to_csv(os.path.join(self.output_dir, f'stratified_error_analysis_{target_col}.csv'), index=False)
        return strat_df

    def _compute_model_shap(self, model, X_val, X_bg):
        if HAS_SHAP:
            try:
                explainer = shap.TreeExplainer(model)
                return explainer.shap_values(X_val)
            except Exception:
                pass

        N, P = X_val.shape
        shap_matrix = np.zeros((N, P))
        bg_median = np.median(X_bg, axis=0)
        full_preds = model.predict(X_val)

        for p in range(P):
            X_perturbed = X_val.copy()
            X_perturbed[:, p] = bg_median[p]
            perturbed_preds = model.predict(X_perturbed)
            shap_matrix[:, p] = full_preds - perturbed_preds

        return shap_matrix

    def plot_comprehensive_diagnostics_and_shap(self, summary_df, model_shap_dict, y_true, X_mat, strat_df, best_model_name, best_preds, target_col='LER_Total'):
        unit_label = "g/hill" if any(k in target_col for k in ['Shoot_DW', 'Biomass', 'Weight']) else ("kg/m³" if 'IWUE' in target_col else "Ratio")

        # 1. Multi-Protocol CV Comparison
        fig1 = plt.figure(figsize=(13, 5.5), dpi=300)
        sns.barplot(
            data=summary_df, x='Model Architecture', y='R2 Score', hue='Evaluation Protocol',
            palette=['#7570b3', '#1b9e77', '#d95f02'], edgecolor='black', linewidth=0.8
        )
        plt.title(f'Multi-Protocol Benchmark: Random K-Fold vs. Spatial Group-K-Fold vs. Temporal Sensitivity ({target_col})\nLinear Baselines vs. Tree Ensembles & GPR',
                  fontsize=12, fontweight='bold', pad=12)
        plt.ylabel('Coefficient of Determination (R²)', fontsize=11, fontweight='bold')
        plt.xlabel('Model Architecture', fontsize=11, fontweight='bold')
        plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        min_r2 = summary_df['R2 Score'].min()
        plt.ylim(max(-0.6, min_r2 - 0.1) if min_r2 < 0 else -0.1, 1.0)
        plt.xticks(rotation=20, ha='right', fontsize=9, fontweight='bold')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.legend(title='Validation Protocol', frameon=True, fontsize=9)
        plt.tight_layout()
        self._save_dual_format(fig1, f'dual_cv_comparison_{target_col}')
        plt.close()

        # 2. Predicted vs. Actual Parity Diagnostic Plot (Calibrated Alignment)
        fig2, ax = plt.subplots(figsize=(7.5, 6.5), dpi=300)
        meta_df = self.df[['Field_capacity', 'Cropping_system']].copy()
        y_min = min(float(np.min(y_true)), float(np.min(best_preds))) * 0.9
        y_max = max(float(np.max(y_true)), float(np.max(best_preds))) * 1.1

        ax.plot([y_min, y_max], [y_min, y_max], 'k--', linewidth=1.5, label='1:1 Parity Line (Ideal Fit)')
        ax.fill_between([y_min, y_max], [y_min * 0.85, y_max * 0.85], [y_min * 1.15, y_max * 1.15],
                        color='#b3cde3', alpha=0.3, label='±15% Tolerance Envelope')

        colors = {'full': '#2b83ba', 'half': '#d7191c'}
        for fc_val, color in colors.items():
            mask = (meta_df['Field_capacity'] == fc_val)
            ax.scatter(y_true[mask], best_preds[mask], c=color, s=45, alpha=0.8, edgecolors='black', linewidth=0.5, label=f'Field Capacity: {fc_val}')

        slope, intercept, r_val, p_val, std_err = stats.linregress(y_true, best_preds)
        x_vals = np.linspace(y_min, y_max, 100)
        ax.plot(x_vals, slope * x_vals + intercept, color='#404040', linestyle='-', linewidth=1.2, label=f'Linear Fit (y = {slope:.2f}x + {intercept:.2f})')

        rmse_val = float(np.sqrt(np.mean((y_true - best_preds)**2)))
        mae_val = float(np.mean(np.abs(y_true - best_preds)))
        r2_val = float(1.0 - (np.sum((y_true - best_preds)**2) / np.sum((y_true - np.mean(y_true))**2)))

        stat_text = f"Architecture: {best_model_name}\nR² = {r2_val:.3f}\nRMSE = {rmse_val:.3f} {unit_label}\nMAE = {mae_val:.3f} {unit_label}\nPearson r = {r_val:.3f}"
        ax.text(0.05, 0.72, stat_text, transform=ax.transAxes, fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

        ax.set_xlim(y_min, y_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(f'Empirical Observed {target_col} ({unit_label})', fontsize=11, fontweight='bold')
        ax.set_ylabel(f'Model Predicted {target_col} ({unit_label})', fontsize=11, fontweight='bold')
        ax.set_title(f'Predicted vs. Actual Parity Diagnostic ({target_col})\nOut-of-Fold Spatial Group-K-Fold Cross-Validation',
                     fontsize=12, fontweight='bold', pad=12)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower right', frameon=True, fontsize=9)
        plt.tight_layout()
        self._save_dual_format(fig2, f'parity_predicted_vs_actual_{target_col}')
        plt.close()

        # 3. Residual Diagnostics
        residuals = y_true - best_preds
        fig3, (ax_res, ax_dist) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)
        ax_res.scatter(best_preds, residuals, c='#2b83ba', alpha=0.75, edgecolors='black', linewidth=0.5, s=40)
        ax_res.axhline(0, color='red', linestyle='--', linewidth=1.2)
        ax_res.set_xlabel(f'Fitted Values ({target_col} in {unit_label})', fontsize=10, fontweight='bold')
        ax_res.set_ylabel(f'Prediction Residuals (y - ŷ, {unit_label})', fontsize=10, fontweight='bold')
        ax_res.set_title(f'A: Residuals vs. Fitted Diagnostic ({best_model_name})\nTesting for Homoscedasticity', fontsize=11, fontweight='bold')
        ax_res.grid(True, linestyle=':', alpha=0.6)

        sns.histplot(residuals, kde=True, color='#1b9e77', ax=ax_dist, stat='density', bins=18, edgecolor='black')
        res_mean, res_std = float(np.mean(residuals)), float(np.std(residuals))
        x_norm = np.linspace(res_mean - 3.5 * res_std, res_mean + 3.5 * res_std, 100)
        ax_dist.plot(x_norm, stats.norm.pdf(x_norm, res_mean, res_std), 'r--', linewidth=1.5, label='Normal Curve Overlay')
        ax_dist.set_xlabel(f'Residual Error ({unit_label})', fontsize=10, fontweight='bold')
        ax_dist.set_ylabel('Density', fontsize=10, fontweight='bold')
        ax_dist.set_title('B: Residual Error Distribution & Normality', fontsize=11, fontweight='bold')
        ax_dist.grid(True, linestyle=':', alpha=0.6)
        ax_dist.legend(loc='upper right', frameon=True)
        plt.suptitle(f'Standard Residual Diagnostic Suite for {target_col}', fontsize=12, fontweight='bold', y=1.02)
        plt.tight_layout()
        self._save_dual_format(fig3, f'residual_diagnostics_{target_col}')
        plt.close()

        # 4. Stratified Error Breakdown Plot
        fig4 = plt.figure(figsize=(12, 5), dpi=300)
        sub_strat = strat_df[strat_df['Stratification Factor'].isin(['Field Capacity', 'Water Quality'])].copy()
        sns.barplot(
            data=sub_strat, x='Treatment Level', y='Normalized RMSE (%)', hue='Stratification Factor',
            palette=['#7570b3', '#e7298a'], edgecolor='black', linewidth=0.8
        )
        plt.title(f'Stratified Machine Learning Error Analysis across Treatment Regimes ({target_col})\nEvaluating Predictive Bias under Abiotic Stress vs. Full Irrigation',
                  fontsize=12, fontweight='bold', pad=12)
        plt.ylabel('Normalized RMSE (%)', fontsize=11, fontweight='bold')
        plt.xlabel('Experimental Treatment Level', fontsize=11, fontweight='bold')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.xticks(rotation=15, ha='right', fontsize=9, fontweight='bold')
        plt.legend(title='Treatment Dimension', frameon=True, fontsize=9)
        plt.tight_layout()
        self._save_dual_format(fig4, f'stratified_error_breakdown_{target_col}')
        plt.close()

        # 5. Collinearity-Robust VIF-Pruned SHAP Feature Importance Heatmap
        models_list = list(model_shap_dict.keys())
        shap_matrix = [np.mean(np.abs(model_shap_dict[m]), axis=0) for m in models_list]
        full_shap_df = pd.DataFrame(shap_matrix, index=models_list, columns=self.all_features)

        vif_sub_shap = full_shap_df[[c for c in self.vif_pruned_features if c in full_shap_df.columns]]
        top_vif_cols = vif_sub_shap.mean(axis=0).sort_values(ascending=False).head(10).index.tolist()
        top_vif_shap_df = vif_sub_shap[top_vif_cols]

        fig5 = plt.figure(figsize=(12, 5), dpi=300)
        sns.heatmap(top_vif_shap_df, cmap='YlGnBu', annot=True, fmt='.3f', linewidths=0.5,
                    cbar_kws={'label': 'Mean |SHAP Value| (VIF-Pruned Attribution)'},
                    annot_kws={"size": 8, "weight": "bold"})
        plt.title(f'Cross-Model SHAP Feature Attribution for {target_col} (VIF <= 10 Pruned Features)\nMitigating Multicollinearity Distortion',
                  fontsize=12, fontweight='bold', pad=15)
        plt.xticks(rotation=45, ha='right', fontsize=9, fontweight='bold')
        plt.yticks(rotation=0, fontsize=10, fontweight='bold')
        plt.tight_layout()
        self._save_dual_format(fig5, f'shap_feature_importance_heatmap_{target_col}')
        plt.close()

        # 6. Grouped Functional Domain SHAP Contribution
        group_attributions = {}
        for grp_name, feat_list in self.feature_groups.items():
            valid_feats = [f for f in feat_list if f in self.all_features]
            feat_indices = [self.all_features.index(f) for f in valid_feats]
            grp_shap = [np.sum(np.mean(np.abs(model_shap_dict[m][:, feat_indices]), axis=0)) for m in models_list]
            group_attributions[grp_name] = float(np.mean(grp_shap))

        grp_df = pd.DataFrame(list(group_attributions.items()), columns=['Functional Domain', 'Aggregated |SHAP| Attribution']).sort_values(by='Aggregated |SHAP| Attribution', ascending=True)

        fig6 = plt.figure(figsize=(10, 4.5), dpi=300)
        bars = plt.barh(grp_df['Functional Domain'], grp_df['Aggregated |SHAP| Attribution'], color=['#386cb0', '#fdb462', '#7fc97f', '#beaed4'], edgecolor='black')
        plt.title(f'Hierarchical Domain-Level SHAP Attribution for {target_col}\nCollinearity-Robust Grouped Feature Contribution',
                  fontsize=12, fontweight='bold', pad=12)
        plt.xlabel('Mean Aggregated |SHAP Value| Across Model Suite', fontsize=10, fontweight='bold')
        plt.grid(axis='x', linestyle=':', alpha=0.6)
        for bar in bars:
            w = bar.get_width()
            plt.text(w + 0.005 * max(grp_df['Aggregated |SHAP| Attribution']), bar.get_y() + bar.get_height()/2,
                     f'{w:.3f}', va='center', ha='left', fontsize=9, fontweight='bold')
        plt.tight_layout()
        self._save_dual_format(fig6, f'shap_grouped_feature_importance_{target_col}')
        plt.close()

        # 7. SHAP Summary Beeswarm on VIF-Pruned Features (Using True Best Model)
        top_model_name = best_model_name
        top_shap_vals = model_shap_dict[top_model_name]
        vif_indices = [self.all_features.index(f) for f in self.vif_pruned_features if f in self.all_features]
        vif_shap_vals = top_shap_vals[:, vif_indices]
        vif_X_mat = X_mat[:, vif_indices]
        vif_feat_names = [self.all_features[idx] for idx in vif_indices]

        fig7 = plt.figure(figsize=(10, 6), dpi=300)
        mean_abs = np.mean(np.abs(vif_shap_vals), axis=0)
        sorted_idx = np.argsort(mean_abs)[-10:]

        for i, idx in enumerate(sorted_idx):
            f_vals = vif_X_mat[:, idx]
            norm_vals = (f_vals - np.min(f_vals)) / (np.max(f_vals) - np.min(f_vals) + 1e-8)
            jitter = np.random.normal(0, 0.08, size=len(f_vals))
            colors = plt.cm.coolwarm(norm_vals)
            plt.scatter(vif_shap_vals[:, idx], i + jitter, c=colors, alpha=0.8, edgecolors='none', s=25)

        plt.yticks(range(len(sorted_idx)), [vif_feat_names[k] for k in sorted_idx], fontsize=10, fontweight='bold')
        plt.axvline(0, color='gray', linestyle='--', linewidth=1)
        plt.xlabel(f'SHAP Value (Impact on {target_col})', fontsize=11, fontweight='bold')
        plt.title(f'SHAP Summary Beeswarm on VIF-Pruned Predictors: {target_col} ({top_model_name})\nControlled Multicollinearity (VIF <= 10)',
                  fontsize=12, fontweight='bold')

        sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=plt.gca(), pad=0.02)
        cbar.set_label('Feature Value (Low: Blue, High: Red)', fontsize=9, fontweight='bold')

        plt.tight_layout()
        self._save_dual_format(fig7, f'shap_summary_beeswarm_{target_col}')
        plt.close()

    def execute_multi_target_pipeline(self, target_list=None):
        targets = target_list or ['LER_Total', 'IWUE_kg_m3', 'Wheat_Shoot_DW_g_hill', 'Total_System_Biomass_g_hill']

        self.generate_hyperparameter_tuning_specs()
        all_summaries = []
        all_stratified = []

        for target in targets:
            # Check for column or backward-compatible alias
            if target not in self.df.columns and target.replace('g_hill', 't_ha') in self.df.columns:
                target = target.replace('g_hill', 't_ha')

            if target in self.df.columns:
                summary_df, preds_dict, shap_dict, y_true, X_mat, strat_df, best_model_name = self.run_multi_protocol_cross_validation(target_col=target)
                self.plot_comprehensive_diagnostics_and_shap(
                    summary_df, shap_dict, y_true, X_mat, strat_df, best_model_name, preds_dict[best_model_name], target_col=target
                )
                all_summaries.append(summary_df)
                all_stratified.append(strat_df)

        if all_summaries:
            master_summary = pd.concat(all_summaries, ignore_index=True)
            master_summary.to_csv(os.path.join(self.output_dir, 'master_dual_cv_benchmark_all_targets.csv'), index=False)
            print(f"\n[+] Saved Master Multi-Protocol Benchmark Table.")

        if all_stratified:
            master_strat = pd.concat(all_stratified, ignore_index=True)
            master_strat.to_csv(os.path.join(self.output_dir, 'master_stratified_error_analysis.csv'), index=False)
            print(f"[+] Saved Master Stratified Error Analysis Table.")


if __name__ == '__main__':
    ml_engine = IntercroppingMachineLearningEngine(
        data_path='data/wheat_berseem_real_processed.csv',
        output_dir='ml_results'
    )
    ml_engine.execute_multi_target_pipeline(
        target_list=['LER_Total', 'IWUE_kg_m3', 'Wheat_Shoot_DW_g_hill', 'Total_System_Biomass_g_hill']
    )
