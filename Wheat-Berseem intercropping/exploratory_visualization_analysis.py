"""
MODULE 2: EXPLORATORY DATA ANALYSIS (EDA) & PUBLICATION-QUALITY VISUALIZATION
=============================================================================
Title: Optimizing Wheat-Berseem Intercrop Yields: Machine Learning Analysis 
       of Morphological and Yield Data under Varied Water Irrigation Regimes

Generates 8 high-resolution (300 DPI) publication figures:
- Fig 1: Agro-Meteorological & Atmospheric Stress Dynamics across 2 Seasons
- Fig 2: Multi-Factorial Yield & Biomass Performance under Abiotic Stress Regimes
- Fig 3: Intercropping System Efficiency & Interspecific Competition Isoclines (LER, ATER, CR, Aggressivity)
- Fig 4: Irrigation Water Use Efficiency (IWUE) & System Productivity Index (SPI) vs Stress Severity (SSI)
- Fig 5: Berseem Clover Multi-Cutting Regrowth & Temporal Resilience Dynamics
- Fig 6: Comprehensive Clustered Correlation Heatmap across Morphological, Weather, & Yield Metrics
- Fig 7: Canopy Allometry, Root-to-Shoot Allocation, & Interspecific Chlorophyll Dynamics
- Fig 8: Principal Component Analysis (PCA) Multivariate Treatment Clustering & Phenotypic Loadings
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set scientific publication style
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.color'] = '#e0e0e0'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.6


class IntercroppingExploratoryVisualizer:
    def __init__(self, processed_data_path, daily_weather_path=None, output_dir='figures'):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 1. Ingest Processed Matrix
        resolved_matrix_path = self._resolve_path([
            processed_data_path,
            'data/wheat_berseem_real_processed.csv',
            'wheat_berseem_real_processed.csv'
        ])
        
        if resolved_matrix_path is None:
            raise FileNotFoundError(f"Could not locate processed data at '{processed_data_path}' or local folders.")
            
        print(f"[+] Ingesting processed experimental matrix from: '{resolved_matrix_path}'")
        self.df = pd.read_csv(resolved_matrix_path)
        
        # 2. Ingest Daily Weather Time-Series with Auto-Path Resolution
        resolved_weather_path = self._resolve_path([
            daily_weather_path,
            'weather_daily_data.csv',
            'data/weather_daily_data.csv',
            'weather_data.csv',
            'data/weather_data.csv'
        ])
        
        self.df_weather_daily = None
        if resolved_weather_path:
            print(f"[+] Ingesting daily weather time-series from: '{resolved_weather_path}'")
            raw_w = pd.read_csv(resolved_weather_path)
            raw_w = raw_w.dropna(subset=['Date'])
            raw_w = raw_w[~raw_w['Date'].astype(str).str.startswith((',', 'mbar', 'Average', 'Total'))]
            raw_w['Date'] = pd.to_datetime(raw_w['Date'], errors='coerce')
            raw_w = raw_w.dropna(subset=['Date']).sort_values('Date').reset_index(drop=True)
            for c in raw_w.columns:
                if c != 'Date': 
                    raw_w[c] = pd.to_numeric(raw_w[c], errors='coerce')
            
            raw_w['Season'] = [1 if (d.year == 2021 or (d.year == 2022 and d.month <= 4)) else 2 for d in raw_w['Date']]
            
            # Daily Vapor Pressure Deficit (VPD) and Growing Degree Days (GDD)
            es_tmax = 0.61078 * np.exp((17.27 * raw_w['Air_Temp_C_Max']) / (raw_w['Air_Temp_C_Max'] + 237.3))
            es_tmin = 0.61078 * np.exp((17.27 * raw_w['Air_Temp_C_Min']) / (raw_w['Air_Temp_C_Min'] + 237.3))
            es = (es_tmax + es_tmin) / 2.0
            ea = (es_tmin * (raw_w['RH_Max'] / 100.0) + es_tmax * (raw_w['RH_Min'] / 100.0)) / 2.0
            raw_w['VPD_Daily_kPa'] = np.maximum(0.0, es - ea)
            raw_w['GDD_Daily'] = np.maximum(0.0, ((raw_w['Air_Temp_C_Max'] + raw_w['Air_Temp_C_Min']) / 2.0) - 5.0)
            self.df_weather_daily = raw_w
        else:
            print("[!] Daily weather CSV not found on disk. Falling back to seasonal weather metrics in processed matrix.")

        # Color palettes
        self.wq_palette = {
            'FW': '#2b83ba',
            'inorganic-F': '#4daf4a',
            'Mixure': '#fdae61',
            'Fish_effluent': '#d7191c'
        }
        self.sys_palette = {
            'MW': '#7570b3',
            'WBI': '#1b9e77'
        }

    def _resolve_path(self, candidates):
        """Finds the first existing candidate file path."""
        for p in candidates:
            if p and os.path.exists(p):
                return p
        return None

    # =========================================================================
    # FIGURE 1: AGRO-METEOROLOGICAL & ATMOSPHERIC DYNAMICS
    # =========================================================================
    def plot_fig1_weather_dynamics(self):
        """Figure 1: Multi-Seasonal Agro-Meteorological Dynamics (VPD, GDD, Temperature, Rainfall)."""
        if self.df_weather_daily is not None:
            # Mode A: High-resolution daily time-series
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), dpi=300, sharex=False)
            
            for s, color, label in [(1, '#1f77b4', 'Season 1 (2021-2022)'), (2, '#d62728', 'Season 2 (2022-2023)')]:
                s_data = self.df_weather_daily[self.df_weather_daily['Season'] == s].reset_index(drop=True)
                days = np.arange(len(s_data))
                
                # Panel A: Mean, Max, Min Temperature
                ax1.plot(days, s_data['AirTempC_Avg'], color=color, label=f"{label} Mean Temp", linewidth=1.8)
                ax1.fill_between(days, s_data['Air_Temp_C_Min'], s_data['Air_Temp_C_Max'], color=color, alpha=0.15)
                
                # Panel B: Vapor Pressure Deficit (VPD)
                ax2.plot(days, s_data['VPD_Daily_kPa'], color=color, label=f"{label} VPD", linewidth=1.8)
                
                # Panel C: Cumulative GDD
                ax3.plot(days, np.cumsum(s_data['GDD_Daily']), color=color, label=f"{label} Cumulative GDD", linewidth=2.0)

            ax1.set_ylabel('Air Temperature (°C)', fontsize=11, fontweight='bold')
            ax1.set_title('(A) Daily Thermal Dynamics & Min-Max Ranges across Growing Seasons', fontsize=12, fontweight='bold', loc='left')
            ax1.legend(loc='upper right', frameon=True)
            ax1.grid(True)
            
            ax2.set_ylabel('VPD (kPa)', fontsize=11, fontweight='bold')
            ax2.set_title('(B) Atmospheric Evaporative Demand (Vapor Pressure Deficit, VPD)', fontsize=12, fontweight='bold', loc='left')
            ax2.legend(loc='upper left', frameon=True)
            ax2.grid(True)
            
            ax3.set_xlabel('Days After Sowing (DAS - November to April)', fontsize=11, fontweight='bold')
            ax3.set_ylabel('Thermal Time (°C·days)', fontsize=11, fontweight='bold')
            ax3.set_title('(C) Accumulated Growing Degree Days (GDD Base 5.0 °C)', fontsize=12, fontweight='bold', loc='left')
            ax3.legend(loc='upper left', frameon=True)
            ax3.grid(True)

        else:
            # Mode B: Fallback using seasonal aggregated weather columns in processed dataset
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 9), dpi=300)
            
            # Cumulative Rainfall
            sns.barplot(data=self.df, x='Season', y='Rainfall_mm_Tot', ax=ax1, palette=['#1f77b4', '#d62728'], edgecolor='black')
            ax1.set_title('(A) Total Cumulative Precipitation (mm)', fontsize=11, fontweight='bold', loc='left')
            ax1.set_ylabel('Rainfall (mm)', fontsize=10, fontweight='bold')
            ax1.grid(axis='y')

            # Mean Air Temperature
            sns.barplot(data=self.df, x='Season', y='AirTempC_Avg', ax=ax2, palette=['#1f77b4', '#d62728'], edgecolor='black')
            ax2.set_title('(B) Mean Growing Season Temperature (°C)', fontsize=11, fontweight='bold', loc='left')
            ax2.set_ylabel('Air Temperature (°C)', fontsize=10, fontweight='bold')
            ax2.grid(axis='y')

            # Seasonal VPD
            sns.barplot(data=self.df, x='Season', y='VPD_kPa', ax=ax3, palette=['#1f77b4', '#d62728'], edgecolor='black')
            ax3.set_title('(C) Mean Vapor Pressure Deficit (kPa)', fontsize=11, fontweight='bold', loc='left')
            ax3.set_ylabel('VPD (kPa)', fontsize=10, fontweight='bold')
            ax3.grid(axis='y')

            # Cumulative GDD
            sns.barplot(data=self.df, x='Season', y='GDD_Cumulative_Season', ax=ax4, palette=['#1f77b4', '#d62728'], edgecolor='black')
            ax4.set_title('(D) Cumulative Thermal Time (GDD °C·days)', fontsize=11, fontweight='bold', loc='left')
            ax4.set_ylabel('Thermal Time (°C·days)', fontsize=10, fontweight='bold')
            ax4.grid(axis='y')

        plt.suptitle('Figure 1: Multi-Seasonal Agro-Meteorological Profiles and Thermal Accumulation', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig1_AgroMeteorological_Dynamics.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 1: '{save_path}'")

    # =========================================================================
    # FIGURE 2: MULTI-FACTORIAL YIELD & BIOMASS PERFORMANCE
    # =========================================================================
    def plot_fig2_biomass_yield_performance(self):
        """Figure 2: Faceted Bar/Box Plots of Wheat Biomass, Berseem Forage, and Total System Biomass."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
        
        # 1. Wheat Shoot Dry Weight by System and Field Capacity
        sns.barplot(data=self.df, x='Water_quality', y='Wheat_Shoot_DW_t_ha', hue='Field_capacity', 
                    palette=['#2ca02c', '#d62728'], ax=ax1, errorbar='se', capsize=0.1, edgecolor='black')
        ax1.set_title('(A) Wheat Shoot Dry Weight by Irrigation Regime', fontsize=12, fontweight='bold', loc='left')
        ax1.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Wheat Shoot Dry Weight (t/ha)', fontsize=11, fontweight='bold')
        ax1.grid(axis='y')

        # 2. Berseem Total Dry Weight (Intercropped WBI plots)
        wbi_df = self.df[self.df['Cropping_system'] == 'WBI']
        sns.barplot(data=wbi_df, x='Water_quality', y='Berseem_Dry_Weight_t_ha', hue='Field_capacity',
                    palette=['#1f77b4', '#ff7f0e'], ax=ax2, errorbar='se', capsize=0.1, edgecolor='black')
        ax2.set_title('(B) Berseem Cumulative Forage Dry Matter (WBI)', fontsize=12, fontweight='bold', loc='left')
        ax2.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Berseem Dry Weight (t/ha)', fontsize=11, fontweight='bold')
        ax2.grid(axis='y')

        # 3. Total System Biomass (Wheat + Berseem) by Cropping System
        sns.boxplot(data=self.df, x='Cropping_system', y='Total_System_Biomass_t_ha', hue='Field_capacity',
                    palette=['#7fc97f', '#beaed4'], ax=ax3)
        ax3.set_title('(C) Total System Biomass Productivity by Cropping Pattern', fontsize=12, fontweight='bold', loc='left')
        ax3.set_xlabel('Cropping System (MW vs WBI)', fontsize=11, fontweight='bold')
        ax3.set_ylabel('Total System Biomass (t/ha)', fontsize=11, fontweight='bold')
        ax3.grid(axis='y')

        # 4. Seasonal Yield Comparison (Season 1 vs Season 2)
        sns.barplot(data=self.df, x='Water_quality', y='Total_System_Biomass_t_ha', hue='Season',
                    palette=['#386cb0', '#f0027f'], ax=ax4, errorbar='se', capsize=0.1, edgecolor='black')
        ax4.set_title('(D) Total System Biomass Stability Across 2 Growing Seasons', fontsize=12, fontweight='bold', loc='left')
        ax4.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax4.set_ylabel('Total Biomass (t/ha)', fontsize=11, fontweight='bold')
        ax4.grid(axis='y')

        plt.suptitle('Figure 2: Crop Growth and Aboveground Biomass Responses to Water Quality & Deficit Regimes', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig2_Biomass_Yield_Performance.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 2: '{save_path}'")

    # =========================================================================
    # FIGURE 3: INTERCROPPING COMPETITION ISOCLINES & EFFICIENCY
    # =========================================================================
    def plot_fig3_intercropping_indices(self):
        """Figure 3: pLER Isocline Diagram, LER_Total & ATER, and Competitive Ratio / Aggressivity."""
        wbi_df = self.df[self.df['Cropping_system'] == 'WBI'].copy()
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
        
        # 1. Bivariate LER Isocline (pLER_Wheat vs pLER_Berseem)
        for wq, color in self.wq_palette.items():
            sub = wbi_df[wbi_df['Water_quality'] == wq]
            ax1.scatter(sub['pLER_Berseem'], sub['pLER_Wheat'], color=color, label=wq, s=60, edgecolors='black', alpha=0.85)

        # Draw LER = 1.0 Isocline
        x_iso = np.linspace(0, 1.2, 100)
        ax1.plot(x_iso, 1.0 - x_iso, 'k--', linewidth=1.5, label='LER = 1.0 (Break-even)')
        ax1.plot(x_iso, 1.2 - x_iso, 'g:', linewidth=1.5, label='LER = 1.2 (+20% Gain)')
        ax1.set_xlim(0, 1.2)
        ax1.set_ylim(0, 2.5)
        ax1.set_xlabel('Partial LER Berseem (pLER_B)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Partial LER Wheat (pLER_W)', fontsize=11, fontweight='bold')
        ax1.set_title('(A) Land Equivalent Ratio (LER) Bivariate Phase Isocline', fontsize=12, fontweight='bold', loc='left')
        ax1.legend(loc='upper right', frameon=True, fontsize=9)
        ax1.grid(True)

        # 2. Total LER vs Area-Time Equivalent Ratio (ATER)
        ler_ater_df = wbi_df.melt(id_vars=['Water_quality', 'Field_capacity'], value_vars=['LER_Total', 'ATER'], 
                                  var_name='Index_Type', value_name='Index_Value')
        sns.boxplot(data=ler_ater_df, x='Water_quality', y='Index_Value', hue='Index_Type',
                    palette=['#1b9e77', '#d95f02'], ax=ax2)
        ax2.axhline(1.0, color='red', linestyle='--', linewidth=1.5)
        ax2.set_title('(B) Land Equivalent Ratio (LER) vs Area-Time Equivalent Ratio (ATER)', fontsize=12, fontweight='bold', loc='left')
        ax2.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax2.set_ylabel('System Ratio Value', fontsize=11, fontweight='bold')
        ax2.grid(axis='y')

        # 3. Competitive Ratio for Wheat (CR_Wheat)
        sns.barplot(data=wbi_df, x='Water_quality', y='CR_Wheat', hue='Field_capacity',
                    palette=['#7570b3', '#e7298a'], ax=ax3, errorbar='se', capsize=0.1, edgecolor='black')
        ax3.axhline(1.0, color='black', linestyle=':', linewidth=1.2, label='Neutral Competition (CR=1)')
        ax3.set_title('(C) Interspecific Competitive Ratio (CR_Wheat)', fontsize=12, fontweight='bold', loc='left')
        ax3.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax3.set_ylabel('CR_Wheat (Weighted Ratio)', fontsize=11, fontweight='bold')
        ax3.grid(axis='y')
        ax3.legend(loc='upper right', fontsize=9)

        # 4. Aggressivity of Wheat (Aggressivity_Wheat)
        sns.barplot(data=wbi_df, x='Water_quality', y='Aggressivity_Wheat', hue='Field_capacity',
                    palette=['#e6ab02', '#a6761d'], ax=ax4, errorbar='se', capsize=0.1, edgecolor='black')
        ax4.axhline(0.0, color='red', linestyle='--', linewidth=1.2, label='Zero Aggressivity (A=0)')
        ax4.set_title('(D) Aggressivity Index of Wheat over Berseem (A_W)', fontsize=12, fontweight='bold', loc='left')
        ax4.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax4.set_ylabel('Aggressivity Value', fontsize=11, fontweight='bold')
        ax4.grid(axis='y')
        ax4.legend(loc='upper right', fontsize=9)

        plt.suptitle('Figure 3: System Land-Use Efficiency and Interspecific Species Competition Dynamics', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig3_Intercropping_Competition_Indices.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 3: '{save_path}'")

    # =========================================================================
    # FIGURE 4: WATER PRODUCTIVITY & RESOURCE EFFICIENCY (IWUE)
    # =========================================================================
    def plot_fig4_water_use_efficiency(self):
        """Figure 4: Irrigation Water Use Efficiency (IWUE) and System Productivity Index (SPI)."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

        # 1. IWUE across Cropping System & Field Capacity
        sns.barplot(data=self.df, x='Water_quality', y='IWUE_kg_m3', hue='Cropping_system',
                    palette=self.sys_palette, ax=ax1, errorbar='se', capsize=0.1, edgecolor='black')
        ax1.set_title('(A) Irrigation Water Use Efficiency (IWUE) by Cropping System', fontsize=12, fontweight='bold', loc='left')
        ax1.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax1.set_ylabel('IWUE (kg dry matter / m³ water)', fontsize=11, fontweight='bold')
        ax1.grid(axis='y')

        # 2. IWUE under Full vs Deficit Irrigation
        sns.boxplot(data=self.df, x='Field_capacity', y='IWUE_kg_m3', hue='Cropping_system',
                    palette=self.sys_palette, ax=ax2)
        ax2.set_title('(B) IWUE under Full (100% FC) vs Deficit (50% FC) Irrigation', fontsize=12, fontweight='bold', loc='left')
        ax2.set_xlabel('Field Capacity Regime', fontsize=11, fontweight='bold')
        ax2.set_ylabel('IWUE (kg/m³)', fontsize=11, fontweight='bold')
        ax2.grid(axis='y')

        # 3. System Productivity Index (SPI) across Treatments
        sns.barplot(data=self.df, x='Water_quality', y='SPI', hue='Field_capacity',
                    palette=['#41b6c4', '#253494'], ax=ax3, errorbar='se', capsize=0.1, edgecolor='black')
        ax3.set_title('(C) System Productivity Index (SPI Standardized Yield)', fontsize=12, fontweight='bold', loc='left')
        ax3.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax3.set_ylabel('SPI (Standardized t/ha)', fontsize=11, fontweight='bold')
        ax3.grid(axis='y')

        # 4. IWUE vs Stress Severity Index (SSI) Regression
        sns.scatterplot(data=self.df, x='Stress_Severity_Index', y='IWUE_kg_m3', hue='Cropping_system',
                        style='Field_capacity', s=75, palette=self.sys_palette, ax=ax4, edgecolor='black')
        sns.regplot(data=self.df, x='Stress_Severity_Index', y='IWUE_kg_m3', scatter=False, ax=ax4, color='gray', line_kws={'linestyle': '--'})
        ax4.set_title('(D) IWUE Trajectory across Compounded Stress Severity (SSI)', fontsize=12, fontweight='bold', loc='left')
        ax4.set_xlabel('Stress Severity Index (ECw / %FC)', fontsize=11, fontweight='bold')
        ax4.set_ylabel('IWUE (kg/m³)', fontsize=11, fontweight='bold')
        ax4.grid(True)

        plt.suptitle('Figure 4: Water Use Productivity and Resource Transformation Efficiency', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig4_Water_Use_Efficiency.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 4: '{save_path}'")

    # =========================================================================
    # FIGURE 5: BERSEEM MULTI-CUTTING REGROWTH & RESILIENCE
    # =========================================================================
    def plot_fig5_berseem_cutting_dynamics(self):
        """Figure 5: Multi-Cut Forage Regrowth Trajectory, SPAD, and Cutting Recovery Ratios."""
        cut_cols = [c for c in self.df.columns if 'Berseem_' in c and '_DW_t_ha' in c and 'Cut' in c]
        if not cut_cols:
            return

        wbi_df = self.df[self.df['Cropping_system'] == 'WBI'].copy()
        melted_cuts = wbi_df.melt(id_vars=['Season', 'Replication', 'Field_capacity', 'Water_quality'],
                                  value_vars=cut_cols, var_name='Cut_Name', value_name='Cut_DW_t_ha')
        melted_cuts['Cut_Label'] = melted_cuts['Cut_Name'].apply(lambda x: x.replace('Berseem_', '').replace('_DW_t_ha', ''))

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

        # 1. Regrowth Yield Trajectory across Cuttings by Water Quality
        sns.lineplot(data=melted_cuts, x='Cut_Label', y='Cut_DW_t_ha', hue='Water_quality', style='Water_quality',
                     markers=True, dashes=False, palette=self.wq_palette, ax=ax1, errorbar='se', linewidth=2.2, markersize=8)
        ax1.set_title('(A) Forage Regrowth Trajectory Across Successive Cuttings', fontsize=12, fontweight='bold', loc='left')
        ax1.set_xlabel('Harvest Cutting Number', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Forage Dry Weight (t/ha / cut)', fontsize=11, fontweight='bold')
        ax1.grid(True)

        # 2. Regrowth by Field Capacity (Full vs Deficit)
        sns.lineplot(data=melted_cuts, x='Cut_Label', y='Cut_DW_t_ha', hue='Field_capacity', style='Field_capacity',
                     markers=True, dashes=False, palette=['#31a354', '#e6550d'], ax=ax2, errorbar='se', linewidth=2.2, markersize=8)
        ax2.set_title('(B) Forage Regrowth Dynamics under 100% vs 50% Field Capacity', fontsize=12, fontweight='bold', loc='left')
        ax2.set_xlabel('Harvest Cutting Number', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Forage Dry Weight (t/ha / cut)', fontsize=11, fontweight='bold')
        ax2.grid(True)

        # 3. Cutting Recovery Ratio (Cut_3 / Cut_1)
        if 'Berseem_Cutting_Recovery_Ratio' in wbi_df.columns:
            sns.barplot(data=wbi_df, x='Water_quality', y='Berseem_Cutting_Recovery_Ratio', hue='Field_capacity',
                        palette=['#74c476', '#fd8d3c'], ax=ax3, errorbar='se', capsize=0.1, edgecolor='black')
            ax3.axhline(1.0, color='red', linestyle='--', linewidth=1.2, label='Equal Regrowth (RR=1.0)')
            ax3.set_title('(C) Forage Cutting Recovery Ratio (Final Cut / Initial Cut)', fontsize=12, fontweight='bold', loc='left')
            ax3.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
            ax3.set_ylabel('Cutting Recovery Ratio', fontsize=11, fontweight='bold')
            ax3.grid(axis='y')
            ax3.legend(loc='upper right', fontsize=9)

        # 4. Cutting Yield Stability Coefficient of Variation (CV)
        if 'Berseem_Cut_Stability_CV' in wbi_df.columns:
            sns.boxplot(data=wbi_df, x='Water_quality', y='Berseem_Cut_Stability_CV', hue='Field_capacity',
                        palette=['#9ecae1', '#fc9272'], ax=ax4)
            ax4.set_title('(D) Temporal Forage Production Stability (Cutting CV)', fontsize=12, fontweight='bold', loc='left')
            ax4.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
            ax4.set_ylabel('Cutting Stability CV (Lower = More Stable)', fontsize=11, fontweight='bold')
            ax4.grid(axis='y')

        plt.suptitle('Figure 5: Berseem Clover Multi-Harvest Regrowth Dynamics and Seasonal Yield Stability', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig5_Berseem_Cutting_Dynamics.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 5: '{save_path}'")

    # =========================================================================
    # FIGURE 6: CLUSTERED CORRELATION HEATMAP
    # =========================================================================
    def plot_fig6_correlation_heatmap(self):
        """Figure 6: High-Resolution Clustered Pearson Correlation Heatmap across Traits and Metrics."""
        selected_traits = [
            'Wheat_Shoot_DW_t_ha', 'Berseem_Dry_Weight_t_ha', 'Total_System_Biomass_t_ha',
            'pLER_Wheat', 'pLER_Berseem', 'LER_Total', 'ATER', 'CR_Wheat', 'IWUE_kg_m3', 'SPI',
            'Wheat_PH_cm', 'Wheat_SPAD', 'Wheat_LAI', 'Tillers_per_plant', 'Root_DW', 'TKW',
            'Berseem_Height_cm', 'Berseem_SPAD', 'Height_Differential_cm', 'SPAD_Ratio_Wheat_Berseem',
            'Stress_Severity_Index', 'VPD_kPa', 'GDD_Cumulative_Season'
        ]
        available_traits = [c for c in selected_traits if c in self.df.columns]
        
        corr_matrix = self.df[available_traits].corr(method='pearson')
        
        plt.figure(figsize=(15, 12), dpi=300)
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
        
        sns.heatmap(corr_matrix, mask=mask, cmap='coolwarm', vmin=-1.0, vmax=1.0, center=0,
                    square=True, linewidths=0.5, cbar_kws={"shrink": 0.8, "label": "Pearson Correlation Coefficient (r)"},
                    annot=True, fmt='.2f', annot_kws={"size": 7, "weight": "bold"})
        
        plt.title('Figure 6: Clustered Pearson Correlation Matrix linking Morphological Traits, Weather Indices, & System Efficiencies', 
                  fontsize=13, fontweight='bold', pad=20)
        plt.xticks(rotation=45, ha='right', fontsize=9, fontweight='bold')
        plt.yticks(rotation=0, fontsize=9, fontweight='bold')
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig6_Correlation_Heatmap.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 6: '{save_path}'")

    # =========================================================================
    # FIGURE 7: ALLOMETRIC TRAIT ALLOCATION & CANOPY PLASTICITY
    # =========================================================================
    def plot_fig7_allometric_plasticity(self):
        """Figure 7: Canopy Height Differential, Root-to-Shoot Ratios, and Chlorophyll Dynamics."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

        # 1. Canopy Height Differential (Wheat - Berseem) vs LER_Total
        wbi_df = self.df[self.df['Cropping_system'] == 'WBI']
        sns.scatterplot(data=wbi_df, x='Height_Differential_cm', y='LER_Total', hue='Water_quality', 
                        style='Field_capacity', s=80, palette=self.wq_palette, ax=ax1, edgecolor='black')
        sns.regplot(data=wbi_df, x='Height_Differential_cm', y='LER_Total', scatter=False, ax=ax1, color='black', line_kws={'linestyle': '--'})
        ax1.set_title('(A) Canopy Shading Differential (ΔH) vs Land Equivalent Ratio (LER)', fontsize=12, fontweight='bold', loc='left')
        ax1.set_xlabel('Canopy Height Differential ΔH (Wheat - Berseem, cm)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('LER Total', fontsize=11, fontweight='bold')
        ax1.grid(True)

        # 2. Wheat Root-to-Shoot Ratio (RSR) across Water Quality & Deficit
        sns.boxplot(data=self.df, x='Water_quality', y='Wheat_Root_Shoot_Ratio', hue='Field_capacity',
                    palette=['#a1d99b', '#fcae91'], ax=ax2)
        ax2.set_title('(B) Belowground Biomass Partitioning (Wheat Root:Shoot Ratio)', fontsize=12, fontweight='bold', loc='left')
        ax2.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Wheat Root-to-Shoot Ratio (RSR)', fontsize=11, fontweight='bold')
        ax2.grid(axis='y')

        # 3. SPAD Ratio (Wheat SPAD / Berseem SPAD) across Treatments
        sns.barplot(data=wbi_df, x='Water_quality', y='SPAD_Ratio_Wheat_Berseem', hue='Field_capacity',
                    palette=['#9ecae1', '#fc9272'], ax=ax3, errorbar='se', capsize=0.1, edgecolor='black')
        ax3.axhline(1.0, color='red', linestyle='--', linewidth=1.2, label='Equal Greenness (Ratio=1.0)')
        ax3.set_title('(C) Interspecific Chlorophyll SPAD Ratio (Wheat / Berseem)', fontsize=12, fontweight='bold', loc='left')
        ax3.set_xlabel('Irrigation Water Quality', fontsize=11, fontweight='bold')
        ax3.set_ylabel('SPAD Ratio (Wheat : Berseem)', fontsize=11, fontweight='bold')
        ax3.grid(axis='y')
        ax3.legend(loc='upper right', fontsize=9)

        # 4. Specific Leaf Area (SLA) vs Stress Severity Index (SSI)
        if 'Specific_Leaf_Area_cm2_g' in self.df.columns:
            sns.scatterplot(data=self.df, x='Stress_Severity_Index', y='Specific_Leaf_Area_cm2_g', hue='Cropping_system',
                            style='Field_capacity', s=80, palette=self.sys_palette, ax=ax4, edgecolor='black')
            ax4.set_title('(D) Specific Leaf Area (SLA) Plasticity under Compounded Stress', fontsize=12, fontweight='bold', loc='left')
            ax4.set_xlabel('Stress Severity Index (SSI)', fontsize=11, fontweight='bold')
            ax4.set_ylabel('Specific Leaf Area (cm²/g)', fontsize=11, fontweight='bold')
            ax4.grid(True)

        plt.suptitle('Figure 7: Canopy Architecture, Root Partitioning, and Photosynthetic Plasticity under Abiotic Stress', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig7_Allometric_Plasticity.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 7: '{save_path}'")

    # =========================================================================
    # FIGURE 8: PRINCIPAL COMPONENT ANALYSIS (PCA) MULTIVARIATE ORDINATION
    # =========================================================================
    def plot_fig8_pca_multivariate_analysis(self):
        """Figure 8: PCA Score Biplot with Treatment Clustering and Trait Loading Vectors."""
        pca_features = [
            'Wheat_Shoot_DW_t_ha', 'Total_System_Biomass_t_ha', 'LER_Total', 'IWUE_kg_m3',
            'Wheat_PH_cm', 'Wheat_SPAD', 'Wheat_LAI', 'Tillers_per_plant', 'Root_DW', 'TKW',
            'Stress_Severity_Index'
        ]
        available_cols = [c for c in pca_features if c in self.df.columns]
        
        X = self.df[available_cols].values
        X_std = (X - np.mean(X, axis=0)) / (np.std(X, axis=0) + 1e-8)
        
        # SVD PCA
        U, S, Vt = np.linalg.svd(X_std, full_matrices=False)
        scores = X_std @ Vt.T
        var_exp = (S ** 2) / np.sum(S ** 2) * 100.0
        
        # Explicit 1D coordinate vectors to avoid dimension mismatch
        dim_pc1 = 0
        dim_pc2 = 1
        
        pc1_scores = scores[:, dim_pc1]
        pc2_scores = scores[:, dim_pc2]
        pc1_variance = var_exp[dim_pc1]
        pc2_variance = var_exp[dim_pc2]
        vt_row_pc1 = Vt[dim_pc1, :]
        vt_row_pc2 = Vt[dim_pc2, :]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5), dpi=300)
        
        # 1. PCA Score Biplot by Cropping System
        for sys_code, color in self.sys_palette.items():
            mask = (self.df['Cropping_system'] == sys_code).values
            ax1.scatter(pc1_scores[mask], pc2_scores[mask], c=color, label=f"System: {sys_code}", s=60, edgecolors='black', alpha=0.8)
            
        # Draw feature loading vectors
        for i, feat in enumerate(available_cols):
            ax1.arrow(0, 0, vt_row_pc1[i] * 3.5, vt_row_pc2[i] * 3.5, color='#333333', alpha=0.7, head_width=0.08, linewidth=1.2)
            ax1.text(vt_row_pc1[i] * 3.8, vt_row_pc2[i] * 3.8, feat.replace('_t_ha', '').replace('_cm', '').replace('_kg_m3', ''), 
                     fontsize=8, fontweight='bold', color='#111111')

        ax1.set_xlabel(f"Principal Component 1 ({pc1_variance:.1f}% Variance)", fontsize=11, fontweight='bold')
        ax1.set_ylabel(f"Principal Component 2 ({pc2_variance:.1f}% Variance)", fontsize=11, fontweight='bold')
        ax1.set_title('(A) PCA Score Biplot & Trait Loading Vectors', fontsize=12, fontweight='bold', loc='left')
        ax1.legend(loc='upper right', frameon=True)
        ax1.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        ax1.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax1.grid(True)
        ax1.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax1.grid(True)

        # 2. PCA Scores Colored by Water Quality
        for wq, color in self.wq_palette.items():
            mask = (self.df['Water_quality'] == wq).values
            ax2.scatter(pc1_scores[mask], pc2_scores[mask], c=color, label=f"WQ: {wq}", s=60, edgecolors='black', alpha=0.85)

        ax2.set_xlabel(f"Principal Component 1 ({pc1_variance:.1f}% Variance)", fontsize=11, fontweight='bold')
        ax2.set_ylabel(f"Principal Component 2 ({pc2_variance:.1f}% Variance)", fontsize=11, fontweight='bold')
        ax2.set_title('(B) Multivariate Treatment Ordination by Water Quality', fontsize=12, fontweight='bold', loc='left')
        ax2.legend(loc='upper right', frameon=True)
        ax2.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        ax2.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax2.grid(True)

        plt.suptitle('Figure 8: Principal Component Analysis (PCA) Ordination and Multivariate Trait Loadings', fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, 'Fig8_PCA_Multivariate_Ordination.png')
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
        print(f"[+] Saved Fig 8: '{save_path}'")

    # =========================================================================
    # EXECUTION PIPELINE
    # =========================================================================
    def generate_all_publication_figures(self):

        self.plot_fig1_weather_dynamics()
        self.plot_fig2_biomass_yield_performance()
        self.plot_fig3_intercropping_indices()
        self.plot_fig4_water_use_efficiency()
        self.plot_fig5_berseem_cutting_dynamics()
        self.plot_fig6_correlation_heatmap()
        self.plot_fig7_allometric_plasticity()
        self.plot_fig8_pca_multivariate_analysis()
        print(f"\n[+] All 8 Publication Figures Successfully Generated in: '{os.path.abspath(self.output_dir)}/'")


if __name__ == '__main__':
    visualizer = IntercroppingExploratoryVisualizer(
        processed_data_path='data/wheat_berseem_real_processed.csv',
        daily_weather_path='weather_daily_data.csv',
        output_dir='figures'
    )
    visualizer.generate_all_publication_figures()
