"""
MODULE 1: PRODUCTION DATA INGESTION, ROBUST CLEANING, & FEATURE ENGINEERING
===========================================================================
Engineered to process real-world experimental field data:
- Wheat Trait Dataset (Handles missing values, typos, extra trailing columns, 18 traits)
- Berseem Multi-Cutting Dataset (Handles long-format Cut_No, missing values, cuts 1..3)
- Weather Daily Dataset (Handles unit header rows, daily VPD & GDD calculation, seasonal aggregation)
- Handles Single-Season or Multi-Season datasets dynamically
"""

import os
import re
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings('ignore')


class WheatBerseemFeatureEngine:
    """
    Robust Data Ingestion, Cleaning, Missing Data Imputation, and Feature Engineering Engine.
    """
    def __init__(self, wheat_duration_days=150, berseem_duration_days=180, gdd_base_temp=5.0):
        self.t_w = wheat_duration_days
        self.t_b = berseem_duration_days
        self.T_total = max(wheat_duration_days, berseem_duration_days)
        self.gdd_base_temp = gdd_base_temp
        
        # Salinity ECw approximation map (dS/m) for stress severity indexing
        self.wq_ec_map = {
            'FW': 0.8,
            'inorganic-F': 1.5,
            'Mixure': 3.5,
            'Fish_effluent': 6.8
        }

    # =========================================================================
    # 1. ROBUST WEATHER DATA INGESTION & CODING
    # =========================================================================
    def process_weather_data(self, weather_source):
        """
        Parses daily weather time-series, filters header metadata, assigns seasons,
        calculates daily VPD (Tetens) and GDD, and aggregates seasonal climate metrics.
        """
        print("[+] Processing Micro-Meteorological & Weather Observations...")
        if isinstance(weather_source, pd.DataFrame):
            df_raw = weather_source.copy()
        elif isinstance(weather_source, str):
            if weather_source.endswith('.csv'):
                df_raw = pd.read_csv(weather_source)
            else:
                df_raw = pd.read_excel(weather_source)
        else:
            raise TypeError("Weather source must be a file path or DataFrame.")

        # Filter out metadata / unit description rows
        df_raw = df_raw.dropna(subset=['Date'])
        df_raw = df_raw[~df_raw['Date'].astype(str).str.startswith((',', 'mbar', 'Average', 'Total'))]

        # Parse Date
        df_raw['Date'] = pd.to_datetime(df_raw['Date'], errors='coerce')
        df_raw = df_raw.dropna(subset=['Date']).sort_values('Date').reset_index(drop=True)

        # Convert numeric columns
        num_cols = [c for c in df_raw.columns if c != 'Date']
        for c in num_cols:
            df_raw[c] = pd.to_numeric(df_raw[c], errors='coerce')

        # Detect Season based on dates (Nov 2021 - Apr 2022 = Season 1; Nov 2022 - Apr 2023 = Season 2)
        def assign_season(dt):
            if dt >= pd.Timestamp('2021-11-01') and dt <= pd.Timestamp('2022-04-30'):
                return 1
            elif dt >= pd.Timestamp('2022-11-01') and dt <= pd.Timestamp('2023-04-30'):
                return 2
            else:
                return 1 if dt.year == 2021 or (dt.year == 2022 and dt.month <= 4) else 2

        df_raw['Season'] = df_raw['Date'].apply(assign_season)

        # Calculate daily Vapor Pressure Deficit (VPD in kPa) via Tetens formulation
        t_max = df_raw['Air_Temp_C_Max']
        t_min = df_raw['Air_Temp_C_Min']
        rh_max = df_raw['RH_Max']
        rh_min = df_raw['RH_Min']

        es_tmax = 0.61078 * np.exp((17.27 * t_max) / (t_max + 237.3))
        es_tmin = 0.61078 * np.exp((17.27 * t_min) / (t_min + 237.3))
        es = (es_tmax + es_tmin) / 2.0
        ea = (es_tmin * (rh_max / 100.0) + es_tmax * (rh_min / 100.0)) / 2.0
        df_raw['VPD_Daily_kPa'] = np.maximum(0.0, es - ea)

        # Calculate daily Growing Degree Days (GDD in deg C - days) with base temp = 5.0 C
        df_raw['GDD_Daily'] = np.maximum(0.0, ((t_max + t_min) / 2.0) - self.gdd_base_temp)

        # Seasonal Weather Named Aggregations (Guarantees flat 1D column headers)
        weath_named_aggs = {
            'Barometric Pressure_mbar_Avg': ('Barometric Pressure_mbar_Avg', 'mean'),
            'Rainfall_mm_Tot': ('Rainfall_mm_Tot', 'sum'),
            'AirTempC_Avg': ('AirTempC_Avg', 'mean'),
            'Air_Temp_C_Max': ('Air_Temp_C_Max', 'mean'),
            'Air_Temp_C_Min': ('Air_Temp_C_Min', 'mean'),
            'RH_Max': ('RH_Max', 'mean'),
            'RH_Min': ('RH_Min', 'mean'),
            'SolarW_Avg': ('SolarW_Avg', 'mean'),
            'SolarMJ_Tot': ('SolarMJ_Tot', 'sum'),
            'WindSpeed_ms_S_WVT': ('WindSpeed_ms_S_WVT', 'mean'),
            'ETo_Daily': ('ETo_Daily', 'mean'),
            'ETo_Cumulative_mm': ('ETo_Daily', 'sum'),
            'VPD_kPa': ('VPD_Daily_kPa', 'mean'),
            'VPD_Max_kPa': ('VPD_Daily_kPa', 'max'),
            'GDD_Daily': ('GDD_Daily', 'mean'),
            'GDD_Cumulative_Season': ('GDD_Daily', 'sum')
        }
        seasonal_weather = df_raw.groupby('Season').agg(**weath_named_aggs).reset_index()

        print(f"[+] Weather Data Aggregated across {len(seasonal_weather)} Seasons.")
        return seasonal_weather

    # =========================================================================
    # 2. ROBUST WHEAT DATA INGESTION & IMPUTATION
    # =========================================================================
    def load_and_clean_wheat_data(self, wheat_source):
        """
        Loads wheat data robustly, cleans outliers, standardizes treatment names,
        and imputes missing values using treatment-level medians.
        """
        print("[+] Ingesting and Cleaning Wheat Dataset...")
        if isinstance(wheat_source, pd.DataFrame):
            df = wheat_source.copy()
        elif isinstance(wheat_source, str):
            if wheat_source.endswith('.csv'):
                with open(wheat_source, 'r') as f:
                    lines = f.readlines()
                header_cols = [c.strip() for c in lines[0].split(',') if c.strip()]
                n_cols = len(header_cols)
                parsed_rows = []
                for line in lines[1:]:
                    if not line.strip(): continue
                    parts = [p.strip() for p in line.strip().split(',')]
                    row_data = parts[:n_cols] if len(parts) >= n_cols else parts + [''] * (n_cols - len(parts))
                    parsed_rows.append(row_data)
                df = pd.DataFrame(parsed_rows, columns=header_cols)
            else:
                df = pd.read_excel(wheat_source)
        else:
            raise TypeError("Wheat source must be a file path or DataFrame.")

        # Clean factor columns
        df['Season'] = pd.to_numeric(df['Season'], errors='coerce').fillna(2).astype(int)

        df['Cropping_system'] = df['Cropping_system'].astype(str).str.strip()
        df['Cropping_system'] = df['Cropping_system'].replace({'WM': 'MW', 'Wheat_monocrop': 'MW', 'Wheat-clover_intercrop': 'WBI'})

        df['Field_capacity'] = df['Field_capacity'].astype(str).str.strip().str.lower()
        df['Field_capacity'] = df['Field_capacity'].replace({'100': 'full', '50': 'half', 'full': 'full', 'half': 'half'})

        df['Water_quality'] = df['Water_quality'].astype(str).str.strip().replace({
            'Fish-effluent': 'Fish_effluent',
            'Inorganic_F': 'inorganic-F',
            'Inorganic-F': 'inorganic-F',
            'Mixture': 'Mixure',
            'Mixure': 'Mixure',
            'FW': 'FW',
            'FW ': 'FW'
        })

        # Generate Replication / Hill ID if missing
        if 'Replication' not in df.columns:
            df['Replication'] = df.groupby(['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']).cumcount() + 1
        else:
            df['Replication'] = pd.to_numeric(df['Replication'], errors='coerce')
            df['Replication'] = df['Replication'].fillna(df.groupby(['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']).cumcount() + 1).astype(int)

        # Standardize feature names
        df.rename(columns={
            'PH': 'Wheat_PH_cm', 'SPAD': 'Wheat_SPAD', 'Shoot_DW': 'Wheat_Shoot_DW_t_ha'
        }, inplace=True)

        # Outlier Detection and Correction
        df['Wheat_PH_cm'] = pd.to_numeric(df['Wheat_PH_cm'], errors='coerce')
        df.loc[df['Wheat_PH_cm'] > 200, 'Wheat_PH_cm'] = np.nan

        df['#Plant_Per_Hill'] = pd.to_numeric(df['#Plant_Per_Hill'], errors='coerce')
        df.loc[df['#Plant_Per_Hill'] > 50, '#Plant_Per_Hill'] = np.nan

        # Convert numerical traits and impute missing data using Treatment Medians
        num_w_traits = [
            '#Plant_Per_Hill', '#Tillers_Per_Hill', 'Tillers_per_plant', 'SL (cm)', 'SW (cm)',
            'Wheat_SPAD', 'Wheat_PH_cm', 'StD', 'LL', 'LW', 'LA', 'LAI', 
            'Wheat_Shoot_DW_t_ha', 'Root_DW', 'RL', 'RW', 'HKW', 'TKW'
        ]
        treatment_keys = ['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']
        
        for col in num_w_traits:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                # 1st tier: Treatment-level median imputation
                df[col] = df.groupby(treatment_keys)[col].transform(lambda s: s.fillna(s.median()))
                # 2nd tier: Factor-level median imputation
                df[col] = df.groupby(['Cropping_system', 'Field_capacity'])[col].transform(lambda s: s.fillna(s.median()))
                # 3rd tier: Global column median
                df[col] = df[col].fillna(df[col].median())

        print(f"[+] Cleaned Wheat Dataset: {df.shape} plots across Seasons: {sorted(df['Season'].unique())}")
        return df

    # =========================================================================
    # 3. ROBUST BERSEEM DATA INGESTION & MULTI-CUT RESHAPING
    # =========================================================================
    def load_and_clean_berseem_data(self, berseem_source):
        """
        Ingests berseem data, handles repeated measures / Cut_No, cleans typos,
        imputes missing values with cut-specific treatment medians, and pivots to wide format.
        """
        print("[+] Ingesting and Processing Berseem Dataset...")
        if isinstance(berseem_source, pd.DataFrame):
            df = berseem_source.copy()
        elif isinstance(berseem_source, str):
            if berseem_source.endswith('.csv'):
                df = pd.read_csv(berseem_source)
            else:
                df = pd.read_excel(berseem_source)
        else:
            raise TypeError("Berseem source must be a file path or DataFrame.")

        # Harmonize column names
        rename_map = {
            'Cut No': 'Cut_No', 'Cut_no': 'Cut_No',
            'Fertilizer': 'Water_quality', 'Irrigation': 'Field_capacity',
            'Cultivation system': 'Cropping_system', 'Cultivation_system': 'Cropping_system',
            'Cropping_System': 'Cropping_system', 'Rep': 'Replication', 'Block': 'Replication',
            'Height': 'Berseem_Height_cm', 'SPAD': 'Berseem_SPAD', 'Dry_weight': 'Berseem_Dry_Weight_t_ha'
        }
        df.rename(columns=rename_map, inplace=True)

        df['Season'] = pd.to_numeric(df['Season'], errors='coerce').fillna(2).astype(int)

        df['Cropping_system'] = df['Cropping_system'].astype(str).str.strip().replace({'WM': 'MW'})
        df['Field_capacity'] = df['Field_capacity'].astype(str).str.strip().str.lower()
        df['Field_capacity'] = df['Field_capacity'].replace({'100': 'full', '50': 'half', 'full': 'full', 'half': 'half'})

        df['Water_quality'] = df['Water_quality'].astype(str).str.strip().replace({
            'Fish-effluent': 'Fish_effluent', 'Inorganic_F': 'inorganic-F', 'Inorganic-F': 'inorganic-F',
            'Mixture': 'Mixure', 'Mixure': 'Mixure', 'FW': 'FW', 'FW ': 'FW'
        })

        if 'Replication' in df.columns:
            df['Replication'] = pd.to_numeric(df['Replication'], errors='coerce').fillna(1).astype(int)
        else:
            df['Replication'] = df.groupby(['Season', 'Cut_No', 'Cropping_system', 'Field_capacity', 'Water_quality']).cumcount() + 1

        # Impute missing values for Height, SPAD, and Dry_weight at the Treatment x Cut level
        num_b_traits = ['Berseem_Height_cm', 'Berseem_SPAD', 'Berseem_Dry_Weight_t_ha']
        impute_keys = ['Season', 'Cut_No', 'Cropping_system', 'Field_capacity', 'Water_quality'] if 'Cut_No' in df.columns else ['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']

        for col in num_b_traits:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                df[col] = df.groupby(impute_keys)[col].transform(lambda s: s.fillna(s.median()))
                df[col] = df.groupby(['Cropping_system', 'Field_capacity'])[col].transform(lambda s: s.fillna(s.median()))
                df[col] = df[col].fillna(0.0)

        # Check if Cut_No exists -> Long-to-Wide Repeated Measures Aggregation
        plot_keys = ['Season', 'Replication', 'Cropping_system', 'Field_capacity', 'Water_quality']
        
        if 'Cut_No' in df.columns:
            cuts = sorted(df['Cut_No'].unique())
            print(f"[+] Found {len(cuts)} Berseem Cuttings: {cuts}. Aggregating Multi-Cut Matrix...")

            # 1. Direct, Flat Summary Aggregations (Completely avoids MultiIndex tuple issues)
            summary_dw = df.groupby(plot_keys)['Berseem_Dry_Weight_t_ha'].agg(['sum', 'mean', 'max', 'std']).reset_index()
            summary_dw.rename(columns={
                'sum': 'Berseem_Dry_Weight_t_ha',
                'mean': 'Berseem_Mean_Cut_Dry_Weight_t_ha',
                'max': 'Berseem_Peak_Cut_Dry_Weight_t_ha',
                'std': 'Berseem_Cut_Dry_Weight_Std'
            }, inplace=True)

            summary_ht = df.groupby(plot_keys)['Berseem_Height_cm'].agg(['mean', 'max']).reset_index()
            summary_ht.rename(columns={
                'mean': 'Berseem_Height_cm',
                'max': 'Berseem_Peak_Height_cm'
            }, inplace=True)

            summary_spad = df.groupby(plot_keys)['Berseem_SPAD'].agg(['mean', 'min']).reset_index()
            summary_spad.rename(columns={
                'mean': 'Berseem_SPAD',
                'min': 'Berseem_Min_SPAD'
            }, inplace=True)

            summary = summary_dw.merge(summary_ht, on=plot_keys).merge(summary_spad, on=plot_keys)

            # 2. Pivot Individual Cuts to Wide Columns
            pivot_dw = df.pivot(index=plot_keys, columns='Cut_No', values='Berseem_Dry_Weight_t_ha').reset_index()
            pivot_dw.columns = [f"Berseem_{col}_DW_t_ha" if col not in plot_keys else col for col in pivot_dw.columns]

            pivot_ht = df.pivot(index=plot_keys, columns='Cut_No', values='Berseem_Height_cm').reset_index()
            pivot_ht.columns = [f"Berseem_{col}_Height_cm" if col not in plot_keys else col for col in pivot_ht.columns]

            pivot_spad = df.pivot(index=plot_keys, columns='Cut_No', values='Berseem_SPAD').reset_index()
            pivot_spad.columns = [f"Berseem_{col}_SPAD" if col not in plot_keys else col for col in pivot_spad.columns]

            # Merge Pivots with Summary
            pivoted_berseem = summary.merge(pivot_dw, on=plot_keys).merge(pivot_ht, on=plot_keys).merge(pivot_spad, on=plot_keys)

            # 3. Dynamic Cut Regrowth & Resilience Indices
            dw_cut_cols = [c for c in pivot_dw.columns if 'Berseem_' in c and '_DW_t_ha' in c]
            if len(dw_cut_cols) >= 2:
                first_cut, last_cut = dw_cut_cols[0], dw_cut_cols[-1]
                pivoted_berseem['Berseem_Cutting_Recovery_Ratio'] = np.where(
                    pivoted_berseem[first_cut] > 0,
                    pivoted_berseem[last_cut] / pivoted_berseem[first_cut], 0.0
                )
                pivoted_berseem['Berseem_Cut_Stability_CV'] = np.where(
                    pivoted_berseem['Berseem_Mean_Cut_Dry_Weight_t_ha'] > 0,
                    pivoted_berseem['Berseem_Cut_Dry_Weight_Std'].fillna(0) / pivoted_berseem['Berseem_Mean_Cut_Dry_Weight_t_ha'], 0.0
                )

            return pivoted_berseem
        else:
            return df

    # =========================================================================
    # 4. DATASET INTEGRATION & COMPETITION METRICS
    # =========================================================================
    def integrate_and_calculate_all_metrics(self, df_w, df_b, df_weather=None):
        """
        Integrates Wheat, Berseem, and Weather datasets; computes all agronomic,
        competition, water productivity, atmospheric, and canopy interaction indices.
        """
        join_keys = ['Season', 'Replication', 'Cropping_system', 'Field_capacity', 'Water_quality']
        print(f"[+] Performing Relational Join on Primary Keys: {join_keys}...")
        merged = pd.merge(df_w, df_b, on=join_keys, how='outer', suffixes=('_wheat', '_berseem'))

        # Fill NaNs for Monoculture Wheat (MW has 0 Berseem traits)
        berseem_cols = [c for c in merged.columns if 'Berseem_' in c]
        for col in berseem_cols:
            merged[col] = merged[col].fillna(0.0)

        # Fill NaNs for Wheat features
        num_w_traits = [
            '#Plant_Per_Hill', '#Tillers_Per_Hill', 'Tillers_per_plant', 'SL (cm)', 'SW (cm)',
            'Wheat_SPAD', 'Wheat_PH_cm', 'StD', 'LL', 'LW', 'LA', 'LAI', 
            'Wheat_Shoot_DW_t_ha', 'Root_DW', 'RL', 'RW', 'HKW', 'TKW'
        ]
        for col in num_w_traits:
            if col in merged.columns:
                merged[col] = merged[col].fillna(0.0)

        # Merge Weather Data
        if df_weather is not None:
            if 'Season' in df_weather.columns:
                print("[+] Merging Weather Observations by Season...")
                merged = pd.merge(merged, df_weather, on='Season', how='left')

        # ---------------------------------------------------------------------
        # Monoculture Baseline Matching (Supports 1 Season or Multi-Season)
        # ---------------------------------------------------------------------
        match_keys = ['Season', 'Field_capacity', 'Water_quality']
        sole_w = merged[merged['Cropping_system'] == 'MW']

        if len(sole_w) > 0 and 'Wheat_Shoot_DW_t_ha' in merged.columns:
            ref_w = sole_w.groupby(match_keys)['Wheat_Shoot_DW_t_ha'].mean().reset_index().rename(
                columns={'Wheat_Shoot_DW_t_ha': 'MW_Baseline_Yield_Ref'}
            )
            merged = merged.merge(ref_w, on=match_keys, how='left')
            merged['pLER_Wheat'] = np.where(merged['MW_Baseline_Yield_Ref'] > 0, merged['Wheat_Shoot_DW_t_ha'] / merged['MW_Baseline_Yield_Ref'], 0.0)
        else:
            merged['MW_Baseline_Yield_Ref'] = merged['Wheat_Shoot_DW_t_ha'].max()
            merged['pLER_Wheat'] = merged['Wheat_Shoot_DW_t_ha'] / (merged['MW_Baseline_Yield_Ref'] + 1e-8)

        # Baseline Berseem Reference
        merged['MB_Baseline_Yield_Ref'] = np.where(
            merged['Cropping_system'] == 'WBI',
            merged['Berseem_Dry_Weight_t_ha'] * 2.5,
            merged['Berseem_Dry_Weight_t_ha']
        )
        merged['pLER_Berseem'] = np.where(
            merged['MB_Baseline_Yield_Ref'] > 0,
            merged['Berseem_Dry_Weight_t_ha'] / merged['MB_Baseline_Yield_Ref'], 0.0
        )

        # Total Land Equivalent Ratio (LER) & Area-Time Equivalent Ratio (ATER)
        merged['LER_Total'] = merged['pLER_Wheat'] + merged['pLER_Berseem']
        merged['ATER'] = ((merged['pLER_Wheat'] * self.t_w) + (merged['pLER_Berseem'] * self.t_b)) / self.T_total

        # Sowing Proportions: Wheat = 67% (0.67), Berseem = 33% (0.33)
        z_w, z_b = 0.67, 0.33

        # Competitive Ratio (CR_Wheat) & Aggressivity (Aggressivity_Wheat)
        merged['CR_Wheat'] = np.where(
            (merged['pLER_Berseem'] > 0) & (merged['Cropping_system'] == 'WBI'),
            (merged['pLER_Wheat'] / merged['pLER_Berseem']) * (z_b / z_w), 0.0
        )
        merged['Aggressivity_Wheat'] = np.where(
            merged['Cropping_system'] == 'WBI',
            (merged['pLER_Wheat'] / z_w) - (merged['pLER_Berseem'] / z_b), 0.0
        )

        # System Productivity Index (SPI)
        mean_mw = merged['MW_Baseline_Yield_Ref'].mean()
        mean_mb = merged['MB_Baseline_Yield_Ref'].replace(0, np.nan).mean()
        merged['SPI'] = merged['Wheat_Shoot_DW_t_ha'] + (merged['Berseem_Dry_Weight_t_ha'] * (mean_mw / (mean_mb + 1e-8)))

        # Relative Crowding Coefficient (RCC / K)
        merged['K_Wheat'] = np.where(
            (merged['Cropping_system'] == 'WBI') & (merged['MW_Baseline_Yield_Ref'] > merged['Wheat_Shoot_DW_t_ha']),
            (merged['Wheat_Shoot_DW_t_ha'] * z_b) / ((merged['MW_Baseline_Yield_Ref'] - merged['Wheat_Shoot_DW_t_ha']) * z_w + 1e-8), 1.0
        )
        merged['K_Berseem'] = np.where(
            (merged['Cropping_system'] == 'WBI') & (merged['MB_Baseline_Yield_Ref'] > merged['Berseem_Dry_Weight_t_ha']),
            (merged['Berseem_Dry_Weight_t_ha'] * z_w) / ((merged['MB_Baseline_Yield_Ref'] - merged['Berseem_Dry_Weight_t_ha']) * z_b + 1e-8), 1.0
        )
        merged['RCC_Total'] = merged['K_Wheat'] * merged['K_Berseem']

        # Monetary Advantage Index (MAI)
        p_w, p_b = 250.0, 180.0
        combined_val = (merged['Wheat_Shoot_DW_t_ha'] * p_w) + (merged['Berseem_Dry_Weight_t_ha'] * p_b)
        merged['MAI'] = np.where(
            merged['LER_Total'] > 0,
            combined_val * ((merged['LER_Total'] - 1.0) / merged['LER_Total']), 0.0
        )

        # Water Productivity (IWUE & CWP)
        merged['Total_System_Biomass_t_ha'] = merged['Wheat_Shoot_DW_t_ha'] + merged['Berseem_Dry_Weight_t_ha']
        if 'Applied_Water_m3_ha' not in merged.columns:
            merged['Applied_Water_m3_ha'] = np.where(merged['Field_capacity'] == 'half', 2750.0, 5500.0)

        merged['IWUE_kg_m3'] = np.where(
            merged['Applied_Water_m3_ha'] > 0,
            (merged['Total_System_Biomass_t_ha'] * 1000.0) / merged['Applied_Water_m3_ha'], 0.0
        )

        if 'ETo_Daily' in merged.columns:
            seasonal_eto_m3 = merged['ETo_Daily'] * 150.0 * 10.0
            merged['CWP_kg_m3'] = (merged['Total_System_Biomass_t_ha'] * 1000.0) / (seasonal_eto_m3 + 1e-8)

        # Compounded Stress Severity Index (SSI = ECw / FC_ratio)
        merged['Field_capacity_pct'] = np.where(merged['Field_capacity'] == 'half', 50.0, 100.0)
        merged['Salinity_ECw_Est_dS_m'] = merged['Water_quality'].map(self.wq_ec_map).fillna(1.0)
        merged['Stress_Severity_Index'] = merged['Salinity_ECw_Est_dS_m'] / (merged['Field_capacity_pct'] / 100.0)

        # Treatment One-Hot Encodings
        merged['WQ_Fish_effluent'] = np.where(merged['Water_quality'] == 'Fish_effluent', 1, 0)
        merged['WQ_inorganic_F'] = np.where(merged['Water_quality'] == 'inorganic-F', 1, 0)
        merged['WQ_Mixure'] = np.where(merged['Water_quality'] == 'Mixure', 1, 0)
        merged['WQ_FW'] = np.where(merged['Water_quality'] == 'FW', 1, 0)
        merged['System_WBI'] = np.where(merged['Cropping_system'] == 'WBI', 1, 0)

        # Canopy & Allometric Interactions
        if 'Wheat_PH_cm' in merged.columns and 'Berseem_Height_cm' in merged.columns:
            merged['Height_Differential_cm'] = merged['Wheat_PH_cm'] - merged['Berseem_Height_cm']

        if 'Wheat_SPAD' in merged.columns and 'Berseem_SPAD' in merged.columns:
            merged['SPAD_Ratio_Wheat_Berseem'] = np.where(merged['Berseem_SPAD'] > 0, merged['Wheat_SPAD'] / merged['Berseem_SPAD'], 0.0)
            merged['Total_System_SPAD'] = merged['Wheat_SPAD'] + merged['Berseem_SPAD']

        if 'Root_DW' in merged.columns and 'Wheat_Shoot_DW_t_ha' in merged.columns:
            merged['Wheat_Root_Shoot_Ratio'] = np.where(merged['Wheat_Shoot_DW_t_ha'] > 0, merged['Root_DW'] / (merged['Wheat_Shoot_DW_t_ha'] + 1e-8), 0.0)

        if 'SL (cm)' in merged.columns and 'SW (cm)' in merged.columns:
            merged['Spike_Shape_Factor'] = np.where(merged['SW (cm)'] > 0, merged['SL (cm)'] / merged['SW (cm)'], 0.0)

        if 'LL' in merged.columns and 'LW' in merged.columns:
            merged['Leaf_Aspect_Ratio'] = np.where(merged['LW'] > 0, merged['LL'] / merged['LW'], 0.0)

        if 'LA' in merged.columns and 'Wheat_Shoot_DW_t_ha' in merged.columns:
            merged['Specific_Leaf_Area_cm2_g'] = np.where(merged['Wheat_Shoot_DW_t_ha'] > 0, (merged['LA'] * 100.0) / (merged['Wheat_Shoot_DW_t_ha'] + 1e-8), 0.0)

        print(f"[+] Full Processing Complete! Matrix Shape: {merged.shape} with {len(merged.columns)} Features.")
        return merged

    # =========================================================================
    # 5. END-TO-END PIPELINE EXECUTION
    # =========================================================================
    def run_pipeline(self, wheat_source, berseem_source, weather_source=None, save_csv_path=None):
        """
        Executes end-to-end loading, cleaning, imputation, weather coding, and metric calculation.
        """
        print("=====================================================================")
        print("MODULE 1: ROBUST DATA LOADING, CLEANING & METRICS ENGINE")
        print("=====================================================================")
        
        # 1. Process Weather Data
        df_weath_summary = self.process_weather_data(weather_source) if weather_source is not None else None

        # 2. Ingest & Clean Wheat Data
        df_w_clean = self.load_and_clean_wheat_data(wheat_source)

        # 3. Ingest & Clean Berseem Data (Handles repeated cuts)
        df_b_clean = self.load_and_clean_berseem_data(berseem_source)

        # 4. Integrate & Calculate All Metrics
        final_df = self.integrate_and_calculate_all_metrics(df_w_clean, df_b_clean, df_weath_summary)

        # 5. Export to CSV
        if save_csv_path:
            os.makedirs(os.path.dirname(os.path.abspath(save_csv_path)), exist_ok=True)
            final_df.to_csv(save_csv_path, index=False)
            print(f"[+] Exported Clean Machine Learning Dataset to: '{save_csv_path}'")

        return final_df
