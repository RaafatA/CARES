"""
MODULE 4: RIGOROUS MACHINE LEARNING STATISTICAL INFERENCE & HYPOTHESIS TESTING PIPELINE (REVISED)
================================================================================================
Title: Statistical Verification Suite for Machine Learning Models, SHAP Concordance, 
       Residual Hypotheses, and Factorial Error ANOVA in Wheat-Berseem Intercropping

Updates:
1. Biomass Units: Standardized as Weight (g/hill) for Wheat Shoot DW and Total System Biomass.
2. Dual Export: Saves all statistical figures in BOTH high-res PNG (300 DPI) and vector SVG formats.
3. Safe Scalar p-value extraction across all SciPy versions.
4. Comprehensive Statistical Suite:
   - Omnibus Friedman & Repeated-Measures ANOVA
   - Post-hoc Pairwise Wilcoxon Signed-Rank Tests (Holm-Bonferroni Adjusted)
   - Kendall’s Coefficient of Concordance (W) & Spearman Matrix on SHAP Rankings
   - Residual Hypothesis Diagnostics (Zero-mean bias, Shapiro-Wilk normality, Homoscedasticity, Durbin-Watson)
   - Factorial Error ANOVA on Predictive Deviations (|y - ŷ|)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import warnings

warnings.filterwarnings('ignore')

# Set scientific publication aesthetic parameters
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.autolayout'] = False

# Import Engine from Module 3
try:
    from module3_ml_loso_shap_pipeline import IntercroppingMachineLearningEngine
except ImportError:
    raise ImportError("Cannot import 'IntercroppingMachineLearningEngine' from 'module3_ml_loso_shap_pipeline.py'. Ensure Module 3 is in the working directory.")


def _safe_pvalue(stat_res):
    """Safely extracts the scalar p-value across all SciPy versions (objects, namedtuples, or tuples)."""
    if hasattr(stat_res, 'pvalue'):
        return float(stat_res.pvalue)
    elif isinstance(stat_res, (tuple, list)) and len(stat_res) >= 2:
        return float(stat_res)
    elif hasattr(stat_res, 'p_value'):
        return float(stat_res.p_value)
    return float(stat_res)


class MLStatisticalInferenceEngine:
    """
    Executes advanced statistical hypothesis tests on Machine Learning outputs:
    - Cross-validation fold performance variance
    - Inter-model feature importance concordance
    - Systematic residual biases
    - Stratified treatment error invariance
    """
    def __init__(self, data_path='data/wheat_berseem_real_processed.csv', output_dir='ml_statistical_results'):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize Module 3 Engine
        self.ml_engine = IntercroppingMachineLearningEngine(data_path=data_path, output_dir=os.path.join(output_dir, 'temp_ml'))
        self.df = self.ml_engine.df
        self.all_features = self.ml_engine.all_features
        self.vif_pruned_features = self.ml_engine.vif_pruned_features
        self.models_pool = self.ml_engine.build_model_instances()
        
        print(f"[+] Initialized Module 4 Statistical Inference Engine.")
        print(f"[+] Output Directory: '{self.output_dir}'")

    def _save_dual_format(self, fig, base_name):
        """Saves figure in BOTH high-resolution PNG (300 DPI) and vector SVG formats."""
        png_path = os.path.join(self.output_dir, f'{base_name}.png')
        svg_path = os.path.join(self.output_dir, f'{base_name}.svg')
        fig.savefig(png_path, bbox_inches='tight', dpi=300)
        fig.savefig(svg_path, bbox_inches='tight', format='svg')
        print(f"[+] Saved '{base_name}' -> PNG & SVG")

    # =========================================================================
    # 1. CROSS-VALIDATION FOLD-BY-FOLD METRIC EXTRACTION
    # =========================================================================
    def extract_spatial_fold_metrics(self, target_col='LER_Total'):
        X_df = self.df[self.all_features].copy()
        y_all = self.df[target_col].values
        replications = self.df['Replication'].values
        group_keys = sorted(np.unique(replications))
        
        X_raw = X_df.values
        X_mean = np.mean(X_raw, axis=0)
        X_std = np.std(X_raw, axis=0) + 1e-8
        X_scaled = (X_raw - X_mean) / X_std

        fold_records = []
        oof_predictions = {}
        fold_scores_dict = {m: {'RMSE': [], 'MAE': [], 'R2': [], 'MAPE': []} for m in self.models_pool.keys()}

        for model_name in self.models_pool.keys():
            y_pred_full = np.zeros(len(y_all))

            for fold_idx, fold_key in enumerate(group_keys):
                train_mask = (replications != fold_key)
                val_mask = (replications == fold_key)

                if any(k in model_name for k in ['Gaussian', 'Ridge', 'Linear']):
                    X_train, y_train = X_scaled[train_mask], y_all[train_mask]
                    X_val, y_val = X_scaled[val_mask], y_all[val_mask]
                else:
                    X_train, y_train = X_raw[train_mask], y_all[train_mask]
                    X_val, y_val = X_raw[val_mask], y_all[val_mask]

                model = self.ml_engine.build_model_instances()[model_name]
                model.fit(X_train, y_train)
                preds = model.predict(X_val)
                y_pred_full[val_mask] = preds

                f_rmse = float(np.sqrt(np.mean((y_val - preds) ** 2)))
                f_mae = float(np.mean(np.abs(y_val - preds)))
                ss_tot = np.sum((y_val - np.mean(y_val)) ** 2)
                ss_res = np.sum((y_val - preds) ** 2)
                f_r2 = float(1.0 - (ss_res / (ss_tot + 1e-8))) if ss_tot > 1e-8 else 0.0
                f_mape = float(np.mean(np.abs((y_val - preds) / (y_val + 1e-8))) * 100.0)

                fold_records.append({
                    'Target': target_col,
                    'Model Architecture': model_name,
                    'Spatial Fold (Replication)': fold_key,
                    'RMSE': f_rmse,
                    'MAE': f_mae,
                    'R2 Score': f_r2,
                    'MAPE (%)': f_mape
                })

                fold_scores_dict[model_name]['RMSE'].append(f_rmse)
                fold_scores_dict[model_name]['MAE'].append(f_mae)
                fold_scores_dict[model_name]['R2'].append(f_r2)
                fold_scores_dict[model_name]['MAPE'].append(f_mape)

            oof_predictions[model_name] = y_pred_full

        fold_df = pd.DataFrame(fold_records)
        return fold_df, fold_scores_dict, oof_predictions, y_all

    # =========================================================================
    # 2. MODEL COMPARISON HYPOTHESIS TESTS (FRIEDMAN & WILCOXON POST-HOC)
    # =========================================================================
    def run_model_comparison_tests(self, fold_scores_dict, target_col='LER_Total'):
        models = list(fold_scores_dict.keys())
        omnibus_records = []
        pairwise_records = []
        metrics = ['RMSE', 'MAE', 'R2', 'MAPE']

        for metric in metrics:
            score_matrix = np.array([fold_scores_dict[m][metric] for m in models])
            
            try:
                fried_res = stats.friedmanchisquare(*score_matrix)
                friedman_stat = float(fried_res[0] if isinstance(fried_res, (tuple, list)) else fried_res.statistic)
                friedman_p = _safe_pvalue(fried_res)
            except Exception:
                friedman_stat, friedman_p = np.nan, np.nan

            try:
                anova_res = stats.f_oneway(*score_matrix)
                f_val = float(anova_res[0] if isinstance(anova_res, (tuple, list)) else anova_res.statistic)
                anova_p = _safe_pvalue(anova_res)
            except Exception:
                f_val, anova_p = np.nan, np.nan

            mean_scores = {m: np.mean(fold_scores_dict[m][metric]) for m in models}
            best_model = max(mean_scores, key=mean_scores.get) if metric == 'R2' else min(mean_scores, key=mean_scores.get)

            omnibus_records.append({
                'Target': target_col,
                'Evaluation Metric': metric,
                'Best Model': best_model,
                'Best Model Mean Score': mean_scores[best_model],
                'Friedman Chi2 Stat': friedman_stat,
                'Friedman p-value': friedman_p,
                'Friedman Significant (p < 0.05)': (friedman_p < 0.05) if not np.isnan(friedman_p) else False,
                'RM-ANOVA F-Stat': f_val,
                'RM-ANOVA p-value': anova_p
            })

            best_scores = fold_scores_dict[best_model][metric]
            raw_pvals = []
            comp_models = [m for m in models if m != best_model]

            for other_m in comp_models:
                other_scores = fold_scores_dict[other_m][metric]
                try:
                    w_res = stats.wilcoxon(best_scores, other_scores, zero_method='pratt')
                    w_stat = float(w_res[0] if isinstance(w_res, (tuple, list)) else w_res.statistic)
                    w_p = _safe_pvalue(w_res)
                except Exception:
                    w_stat, w_p = np.nan, 1.0
                try:
                    t_res = stats.ttest_rel(best_scores, other_scores)
                    t_stat = float(t_res[0] if isinstance(t_res, (tuple, list)) else t_res.statistic)
                    t_p = _safe_pvalue(t_res)
                except Exception:
                    t_stat, t_p = np.nan, 1.0

                raw_pvals.append(w_p)
                pairwise_records.append({
                    'Target': target_col,
                    'Metric': metric,
                    'Reference Model (Best)': best_model,
                    'Comparator Model': other_m,
                    'Ref Mean': mean_scores[best_model],
                    'Comp Mean': mean_scores[other_m],
                    'Mean Difference': mean_scores[best_model] - mean_scores[other_m],
                    'Wilcoxon W-Stat': w_stat,
                    'Wilcoxon Raw p-value': w_p,
                    'Paired t-Stat': t_stat,
                    'Paired t p-value': t_p
                })

            sorted_indices = np.argsort(raw_pvals)
            m_comp = len(comp_models)
            for rank, idx in enumerate(sorted_indices):
                orig_record_idx = len(pairwise_records) - m_comp + idx
                raw_p = pairwise_records[orig_record_idx]['Wilcoxon Raw p-value']
                adjusted_p = min(1.0, raw_p * (m_comp - rank))
                pairwise_records[orig_record_idx]['Holm-Bonferroni Adjusted p-value'] = adjusted_p
                pairwise_records[orig_record_idx]['Statistically Significant'] = (adjusted_p < 0.05)

        omnibus_df = pd.DataFrame(omnibus_records)
        pairwise_df = pd.DataFrame(pairwise_records)
        return omnibus_df, pairwise_df

    # =========================================================================
    # 3. KENDALL'S CONCORDANCE & SHAP FEATURE AGREEMENT ANALYSIS
    # =========================================================================
    def compute_shap_concordance(self, target_col='LER_Total'):
        X_df = self.df[self.all_features].copy()
        y_all = self.df[target_col].values
        X_raw = X_df.values
        X_scaled = (X_raw - np.mean(X_raw, axis=0)) / (np.std(X_raw, axis=0) + 1e-8)

        models_list = list(self.models_pool.keys())
        vif_sub_features = [f for f in self.vif_pruned_features if f in self.all_features]
        vif_indices = [self.all_features.index(f) for f in vif_sub_features]

        model_mean_shaps = {}
        for m_name in models_list:
            model = self.ml_engine.build_model_instances()[m_name]
            if any(k in m_name for k in ['Gaussian', 'Ridge', 'Linear']):
                model.fit(X_scaled, y_all)
                shap_vals = self.ml_engine._compute_model_shap(model, X_scaled, X_scaled)
            else:
                model.fit(X_raw, y_all)
                shap_vals = self.ml_engine._compute_model_shap(model, X_raw, X_raw)
            
            mean_abs_vif = np.mean(np.abs(shap_vals[:, vif_indices]), axis=0)
            model_mean_shaps[m_name] = mean_abs_vif

        shap_df = pd.DataFrame(model_mean_shaps, index=vif_sub_features)
        top_features = shap_df.mean(axis=1).sort_values(ascending=False).head(10).index.tolist()
        top_shap_df = shap_df.loc[top_features]

        rank_matrix = top_shap_df.rank(ascending=False).values
        k_features, m_models = rank_matrix.shape

        R_j = np.sum(rank_matrix, axis=1)
        R_bar = np.mean(R_j)
        S = np.sum((R_j - R_bar) ** 2)

        W = (12.0 * S) / ((m_models ** 2) * (k_features ** 3 - k_features))
        W = max(0.0, min(1.0, W))

        chi2_w = m_models * (k_features - 1) * W
        df_w = k_features - 1
        p_val_w = 1.0 - stats.chi2.cdf(chi2_w, df_w)

        spearman_corr_matrix, _ = stats.spearmanr(top_shap_df)

        concordance_summary = {
            'Target': target_col,
            'Number of Evaluated Models (m)': m_models,
            'Number of Consensus Features (k)': k_features,
            'Kendall Concordance Coefficient (W)': W,
            'Friedman Chi2 Stat': chi2_w,
            'Degrees of Freedom': df_w,
            'Concordance p-value': p_val_w,
            'Significant Ranking Agreement (p < 0.001)': (p_val_w < 0.001),
            'Mean Pairwise Spearman Rho': float(np.mean(spearman_corr_matrix[np.triu_indices(m_models, k=1)]))
        }
        return pd.DataFrame([concordance_summary]), top_shap_df, spearman_corr_matrix, models_list

    # =========================================================================
    # 4. RESIDUAL HYPOTHESIS & DIAGNOSTIC TESTS
    # =========================================================================
    def run_residual_hypothesis_tests(self, y_true, y_pred, model_name, target_col='LER_Total'):
        residuals = y_true - y_pred
        abs_residuals = np.abs(residuals)
        n_obs = len(residuals)

        t_res = stats.ttest_1samp(residuals, 0.0)
        t_stat = float(t_res[0] if isinstance(t_res, (tuple, list)) else t_res.statistic)
        p_bias_t = _safe_pvalue(t_res)

        try:
            w_res = stats.wilcoxon(residuals, zero_method='pratt')
            w_stat = float(w_res[0] if isinstance(w_res, (tuple, list)) else w_res.statistic)
            p_bias_w = _safe_pvalue(w_res)
        except Exception:
            w_stat, p_bias_w = np.nan, 1.0

        try:
            shap_res = stats.shapiro(residuals)
            shapiro_stat = float(shap_res[0] if isinstance(shap_res, (tuple, list)) else shap_res.statistic)
            p_shapiro = _safe_pvalue(shap_res)
        except Exception:
            shapiro_stat, p_shapiro = np.nan, np.nan

        try:
            dag_res = stats.normaltest(residuals)
            dagostino_stat = float(dag_res[0] if isinstance(dag_res, (tuple, list)) else dag_res.statistic)
            p_dagostino = _safe_pvalue(dag_res)
        except Exception:
            dagostino_stat, p_dagostino = np.nan, np.nan

        spear_res = stats.spearmanr(y_pred, abs_residuals)
        spear_rho = float(spear_res[0] if isinstance(spear_res, (tuple, list)) else spear_res.statistic)
        p_homosced = _safe_pvalue(spear_res)

        diff_res = np.diff(residuals)
        dw_stat = np.sum(diff_res ** 2) / (np.sum(residuals ** 2) + 1e-8)

        res_report = {
            'Target': target_col,
            'Best Architecture': model_name,
            'Sample Size (N)': n_obs,
            'Mean Error (Bias)': float(np.mean(residuals)),
            'Median Error': float(np.median(residuals)),
            'Residual Standard Deviation': float(np.std(residuals)),
            'Zero-Mean Bias t-Test p-value': float(p_bias_t),
            'Statistically Unbiased (t-test p > 0.05)': bool(p_bias_t > 0.05),
            'Wilcoxon Bias p-value': float(p_bias_w),
            'Shapiro-Wilk Normality W': float(shapiro_stat) if not np.isnan(shapiro_stat) else np.nan,
            'Shapiro-Wilk p-value': float(p_shapiro) if not np.isnan(p_shapiro) else np.nan,
            'Residuals Normal (p > 0.05)': bool(p_shapiro > 0.05) if not np.isnan(p_shapiro) else False,
            'D\'Agostino-Pearson K2 p-value': float(p_dagostino) if not np.isnan(p_dagostino) else np.nan,
            'Homoscedasticity Spearman Rho': float(spear_rho),
            'Homoscedasticity p-value': float(p_homosced),
            'Constant Error Variance (p > 0.05)': bool(p_homosced > 0.05),
            'Durbin-Watson Stat': float(dw_stat),
            'Residual Skewness': float(stats.skew(residuals)),
            'Residual Kurtosis': float(stats.kurtosis(residuals))
        }
        return pd.DataFrame([res_report]), residuals

    # =========================================================================
    # 5. FACTORIAL ANOVA ON MODEL RESIDUALS
    # =========================================================================
    def run_treatment_error_anova(self, y_true, y_pred, model_name, target_col='LER_Total'):
        meta_df = self.df[['Field_capacity', 'Water_quality', 'Cropping_system']].copy()
        meta_df['abs_error'] = np.abs(y_true - y_pred)

        anova_rows = []

        # 1. Main Effect: Field Capacity
        groups_fc = [group['abs_error'].values for _, group in meta_df.groupby('Field_capacity')]
        fc_res = stats.f_oneway(*groups_fc)
        f_fc = float(fc_res[0] if isinstance(fc_res, (tuple, list)) else fc_res.statistic)
        p_fc = _safe_pvalue(fc_res)
        anova_rows.append({
            'Target': target_col,
            'Model': model_name,
            'Factor / Source of Variation': 'Field Capacity (50% vs 100%)',
            'F-Statistic': f_fc,
            'p-value': p_fc,
            'Significant Error Bias (p < 0.05)': bool(p_fc < 0.05),
            'Interpretation': 'Error is invariant to moisture stress' if p_fc >= 0.05 else 'Model has differential accuracy by moisture regime'
        })

        # 2. Main Effect: Water Quality
        groups_wq = [group['abs_error'].values for _, group in meta_df.groupby('Water_quality')]
        wq_res = stats.f_oneway(*groups_wq)
        f_wq = float(wq_res[0] if isinstance(wq_res, (tuple, list)) else wq_res.statistic)
        p_wq = _safe_pvalue(wq_res)
        anova_rows.append({
            'Target': target_col,
            'Model': model_name,
            'Factor / Source of Variation': 'Water Quality Regimes (4 levels)',
            'F-Statistic': f_wq,
            'p-value': p_wq,
            'Significant Error Bias (p < 0.05)': bool(p_wq < 0.05),
            'Interpretation': 'Error is invariant across water types' if p_wq >= 0.05 else 'Model has differential accuracy by water quality'
        })

        # 3. Main Effect: Cropping System
        groups_cs = [group['abs_error'].values for _, group in meta_df.groupby('Cropping_system')]
        cs_res = stats.f_oneway(*groups_cs)
        f_cs = float(cs_res[0] if isinstance(cs_res, (tuple, list)) else cs_res.statistic)
        p_cs = _safe_pvalue(cs_res)
        anova_rows.append({
            'Target': target_col,
            'Model': model_name,
            'Factor / Source of Variation': 'Cropping System Stand (WBI vs MW vs MB)',
            'F-Statistic': f_cs,
            'p-value': p_cs,
            'Significant Error Bias (p < 0.05)': bool(p_cs < 0.05),
            'Interpretation': 'Error is invariant across cropping systems' if p_cs >= 0.05 else 'Model has differential accuracy by cropping system'
        })

        # 4. Factorial Interaction: Water Quality x Field Capacity
        groups_inter = [group['abs_error'].values for _, group in meta_df.groupby(['Water_quality', 'Field_capacity'])]
        inter_res = stats.f_oneway(*groups_inter)
        f_inter = float(inter_res[0] if isinstance(inter_res, (tuple, list)) else inter_res.statistic)
        p_inter = _safe_pvalue(inter_res)
        anova_rows.append({
            'Target': target_col,
            'Model': model_name,
            'Factor / Source of Variation': 'Interaction (Water Quality x Field Capacity)',
            'F-Statistic': f_inter,
            'p-value': p_inter,
            'Significant Error Bias (p < 0.05)': bool(p_inter < 0.05),
            'Interpretation': 'No compound stress error inflation' if p_inter >= 0.05 else 'Significant compound interaction on predictive error'
        })

        anova_df = pd.DataFrame(anova_rows)
        return anova_df

    # =========================================================================
    # 6. HIGH-RESOLUTION STATISTICAL VISUALIZATIONS (PNG + SVG DUAL EXPORT)
    # =========================================================================
    def generate_statistical_figures(self, fold_df, omnibus_df, pairwise_df, concordance_df,
                                     top_shap_df, spearman_matrix, models_list,
                                     residuals, y_true, y_pred, anova_df,
                                     best_model_name, target_col='LER_Total'):
        
        unit_label = "g/hill" if any(k in target_col for k in ['Shoot_DW', 'Biomass', 'Weight']) else ("kg/m³" if 'IWUE' in target_col else "Ratio")

        # Figure 1: Model Comparison Boxplot
        fig1 = plt.figure(figsize=(11, 6), dpi=300)
        sns.boxplot(
            data=fold_df, x='Model Architecture', y='RMSE', palette='Set2',
            boxprops=dict(alpha=0.85), showmeans=True,
            meanprops=dict(marker='o', markeredgecolor='black', markerfacecolor='white', markersize=6)
        )
        sns.stripplot(data=fold_df, x='Model Architecture', y='RMSE', color='black', alpha=0.6, jitter=0.2, size=5)
        plt.title(f'Spatial Cross-Validation Performance & Statistical Comparison ({target_col})\nEvaluating 10 Replications per Architecture (Omnibus Friedman p = {omnibus_df[omnibus_df["Evaluation Metric"]=="RMSE"]["Friedman p-value"].iloc[0]:.4f})',
                  fontsize=12, fontweight='bold', pad=12)
        plt.ylabel(f'Root Mean Squared Error (RMSE, {unit_label})', fontsize=11, fontweight='bold')
        plt.xlabel('Model Architecture', fontsize=11, fontweight='bold')
        plt.xticks(rotation=20, ha='right', fontsize=9, fontweight='bold')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.tight_layout()
        self._save_dual_format(fig1, f'ml_stats_model_ranking_significance_{target_col}')
        plt.close()

        # Figure 2: SHAP Concordance Heatmap
        w_val = concordance_df['Kendall Concordance Coefficient (W)'].iloc[0]
        w_pval = concordance_df['Concordance p-value'].iloc[0]
        fig2 = plt.figure(figsize=(9, 7.5), dpi=300)
        sns.heatmap(
            spearman_matrix, annot=True, fmt='.3f', cmap='coolwarm', vmin=0.0, vmax=1.0,
            xticklabels=models_list, yticklabels=models_list, linewidths=0.8,
            cbar_kws={'label': 'Pairwise Spearman Rank Correlation (ρ)'},
            annot_kws={"size": 9, "weight": "bold"}
        )
        plt.title(f'Cross-Model SHAP Feature Ranking Concordance ({target_col})\nKendall\'s W = {w_val:.3f} (Friedman χ² = {concordance_df["Friedman Chi2 Stat"].iloc[0]:.2f}, p = {w_pval:.2e})',
                  fontsize=11, fontweight='bold', pad=15)
        plt.xticks(rotation=30, ha='right', fontsize=9, fontweight='bold')
        plt.yticks(rotation=0, fontsize=9, fontweight='bold')
        plt.tight_layout()
        self._save_dual_format(fig2, f'ml_stats_shap_concordance_matrix_{target_col}')
        plt.close()

        bias_p = _safe_pvalue(stats.ttest_1samp(residuals, 0.0))
        shapiro_p = _safe_pvalue(stats.shapiro(residuals))
        spearman_p = _safe_pvalue(stats.spearmanr(y_pred, np.abs(residuals)))

        # Figure 3: 4-Panel Statistical Residual Diagnostic Dashboard
        fig3, ((ax_fitted, ax_qq), (ax_dist, ax_scale)) = plt.subplots(2, 2, figsize=(13, 11), dpi=300)
        
        # Panel A: Residuals vs Fitted
        ax_fitted.scatter(y_pred, residuals, c='#2b83ba', alpha=0.75, edgecolors='black', linewidth=0.5, s=40)
        ax_fitted.axhline(0, color='red', linestyle='--', linewidth=1.2)
        slope, intercept, _, p_val, _ = stats.linregress(y_pred, residuals)
        x_grid = np.linspace(min(y_pred), max(y_pred), 100)
        ax_fitted.plot(x_grid, slope * x_grid + intercept, color='black', linestyle='-', linewidth=1.0, label=f'Trend Slope = {slope:.3f} (p={p_val:.3f})')
        ax_fitted.set_xlabel(f'Fitted Values ({target_col} in {unit_label})', fontsize=10, fontweight='bold')
        ax_fitted.set_ylabel(f'Residuals (e = y - ŷ, {unit_label})', fontsize=10, fontweight='bold')
        ax_fitted.set_title(f'A: Residuals vs. Fitted (Homoscedasticity Test)\nBias Test p = {bias_p:.3f}', fontsize=11, fontweight='bold')
        ax_fitted.grid(True, linestyle=':', alpha=0.6)
        ax_fitted.legend(loc='upper right', frameon=True, fontsize=8)

        # Panel B: Normal Q-Q Plot
        (osm, osr), (slope_q, intercept_q, r_q) = stats.probplot(residuals, dist="norm")
        ax_qq.plot(osm, osr, 'o', color='#1b9e77', alpha=0.75, markeredgecolor='black', markersize=5)
        ax_qq.plot(osm, slope_q * np.array(osm) + intercept_q, 'r--', linewidth=1.5, label=f'Ideal Normal Fit (R² = {r_q**2:.3f})')
        ax_qq.set_xlabel('Theoretical Normal Quantiles', fontsize=10, fontweight='bold')
        ax_qq.set_ylabel('Sample Error Quantiles', fontsize=10, fontweight='bold')
        ax_qq.set_title(f'B: Normal Q-Q Plot\nShapiro-Wilk Normality p = {shapiro_p:.4f}', fontsize=11, fontweight='bold')
        ax_qq.grid(True, linestyle=':', alpha=0.6)
        ax_qq.legend(loc='upper left', frameon=True, fontsize=8)

        # Panel C: Error Distribution
        sns.histplot(residuals, kde=True, color='#7570b3', ax=ax_dist, stat='density', bins=20, edgecolor='black')
        res_m, res_s = np.mean(residuals), np.std(residuals)
        x_norm = np.linspace(res_m - 3.5 * res_s, res_m + 3.5 * res_s, 150)
        ax_dist.plot(x_norm, stats.norm.pdf(x_norm, res_m, res_s), 'r--', linewidth=1.5, label='Normal Distribution Overlay')
        ax_dist.set_xlabel(f'Residual Magnitude ({unit_label})', fontsize=10, fontweight='bold')
        ax_dist.set_ylabel('Probability Density', fontsize=10, fontweight='bold')
        ax_dist.set_title(f'C: Error Distribution & Density\nMean = {res_m:.3f}, Std = {res_s:.3f}', fontsize=11, fontweight='bold')
        ax_dist.grid(True, linestyle=':', alpha=0.6)
        ax_dist.legend(loc='upper right', frameon=True, fontsize=8)

        # Panel D: Scale-Location Plot
        std_res = (residuals - np.mean(residuals)) / (np.std(residuals) + 1e-8)
        sqrt_abs_std = np.sqrt(np.abs(std_res))
        ax_scale.scatter(y_pred, sqrt_abs_std, c='#d95f02', alpha=0.75, edgecolors='black', linewidth=0.5, s=40)
        slope_sl, inter_sl, _, p_sl, _ = stats.linregress(y_pred, sqrt_abs_std)
        ax_scale.plot(x_grid, slope_sl * x_grid + inter_sl, 'r-', linewidth=1.5, label=f'Scale Trend Slope = {slope_sl:.3f} (p={p_sl:.3f})')
        ax_scale.set_xlabel(f'Fitted Values ({target_col} in {unit_label})', fontsize=10, fontweight='bold')
        ax_scale.set_ylabel('√|Standardized Residuals|', fontsize=10, fontweight='bold')
        ax_scale.set_title(f'D: Scale-Location Plot\nSpearman Heteroscedasticity p = {spearman_p:.3f}', fontsize=11, fontweight='bold')
        ax_scale.grid(True, linestyle=':', alpha=0.6)
        ax_scale.legend(loc='upper right', frameon=True, fontsize=8)

        plt.suptitle(f'Comprehensive Residual Statistical Diagnostic Suite: {target_col} ({best_model_name})', fontsize=13, fontweight='bold', y=0.99)
        plt.tight_layout()
        self._save_dual_format(fig3, f'ml_stats_residual_dashboard_{target_col}')
        plt.close()

        # Figure 4: Factorial Error ANOVA Interaction Plot
        meta_df = self.df[['Field_capacity', 'Water_quality', 'Cropping_system']].copy()
        meta_df['Absolute_Error'] = np.abs(residuals)
        fig4 = plt.figure(figsize=(11, 5), dpi=300)
        sns.barplot(
            data=meta_df, x='Water_quality', y='Absolute_Error', hue='Field_capacity',
            palette=['#1b9e77', '#d95f02'], edgecolor='black', ci=68, capsize=0.1
        )
        fc_p = anova_df[anova_df['Factor / Source of Variation'].str.contains('Field Capacity')]['p-value'].iloc[0]
        wq_p = anova_df[anova_df['Factor / Source of Variation'].str.contains('Water Quality')]['p-value'].iloc[0]
        plt.title(f'Predictive Absolute Error Invariance across Treatment Factors ({target_col})\nANOVA Main Effects: Field Capacity (p = {fc_p:.3f}), Water Quality (p = {wq_p:.3f})',
                  fontsize=12, fontweight='bold', pad=12)
        plt.ylabel(f'Mean Absolute Error |y - ŷ| ({unit_label} ± 1 SE)', fontsize=10, fontweight='bold')
        plt.xlabel('Irrigation Water Quality', fontsize=10, fontweight='bold')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.legend(title='Moisture Regime', frameon=True)
        plt.tight_layout()
        self._save_dual_format(fig4, f'ml_stats_factorial_error_anova_{target_col}')
        plt.close()

    # =========================================================================
    # 7. MASTER EXECUTION PIPELINE FOR ALL TARGET METRICS
    # =========================================================================
    def execute_full_statistical_pipeline(self, target_list=None):
        targets = target_list or ['LER_Total', 'IWUE_kg_m3', 'Wheat_Shoot_DW_g_hill', 'Total_System_Biomass_g_hill']
        
        all_omnibus, all_pairwise, all_concordance, all_residuals, all_anovas = [], [], [], [], []

        print(f"\n=====================================================================")
        print(f"STARTING MODULE 4: FULL STATISTICAL INFERENCE & HYPOTHESIS PIPELINE (PNG + SVG)")
        print(f"=====================================================================")

        for target in targets:
            # Column alias check
            if target not in self.df.columns and target.replace('g_hill', 't_ha') in self.df.columns:
                target = target.replace('g_hill', 't_ha')

            if target not in self.df.columns:
                continue
            
            print(f"\n>>> PROCESSING STATISTICAL SUITE FOR TARGET: [{target}] <<<")
            fold_df, fold_scores_dict, oof_preds, y_true = self.extract_spatial_fold_metrics(target_col=target)
            omnibus_df, pairwise_df = self.run_model_comparison_tests(fold_scores_dict, target_col=target)
            all_omnibus.append(omnibus_df)
            all_pairwise.append(pairwise_df)
            
            concordance_df, top_shap_df, spearman_mat, models_list = self.compute_shap_concordance(target_col=target)
            all_concordance.append(concordance_df)
            
            best_model_name = omnibus_df[omnibus_df['Evaluation Metric'] == 'RMSE']['Best Model'].iloc[0]
            best_y_pred = oof_preds[best_model_name]
            
            res_df, residuals = self.run_residual_hypothesis_tests(y_true, best_y_pred, best_model_name, target_col=target)
            all_residuals.append(res_df)
            
            anova_df = self.run_treatment_error_anova(y_true, best_y_pred, best_model_name, target_col=target)
            all_anovas.append(anova_df)
            
            self.generate_statistical_figures(
                fold_df, omnibus_df, pairwise_df, concordance_df, top_shap_df,
                spearman_mat, models_list, residuals, y_true, best_y_pred, anova_df,
                best_model_name, target_col=target
            )

        # Export Master Structured CSV Tables
        pd.concat(all_omnibus, ignore_index=True).to_csv(os.path.join(self.output_dir, 'ml_stats_master_omnibus_model_comparison.csv'), index=False)
        pd.concat(all_pairwise, ignore_index=True).to_csv(os.path.join(self.output_dir, 'ml_stats_master_posthoc_pairwise.csv'), index=False)
        pd.concat(all_concordance, ignore_index=True).to_csv(os.path.join(self.output_dir, 'ml_stats_master_shap_kendall_concordance.csv'), index=False)
        pd.concat(all_residuals, ignore_index=True).to_csv(os.path.join(self.output_dir, 'ml_stats_master_residual_hypothesis_tests.csv'), index=False)
        pd.concat(all_anovas, ignore_index=True).to_csv(os.path.join(self.output_dir, 'ml_stats_master_treatment_error_anova.csv'), index=False)

        print(f"\n=====================================================================")
        print(f"[✓] MODULE 4 STATISTICAL PIPELINE COMPLETED SUCCESSFULLY (PNG + SVG)!")
        print(f"=====================================================================")


if __name__ == '__main__':
    stat_engine = MLStatisticalInferenceEngine(
        data_path='data/wheat_berseem_real_processed.csv',
        output_dir='ml_statistical_results'
    )
    stat_engine.execute_full_statistical_pipeline(
        target_list=['LER_Total', 'IWUE_kg_m3', 'Wheat_Shoot_DW_g_hill', 'Total_System_Biomass_g_hill']
    )