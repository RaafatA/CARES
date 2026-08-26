###########################################################
##
##   Project:  Wheat-Berseem intercropping
##   Date:    16/08/2026
##   Author:  Rafat A. Eissa
##
###########################################################

"""
MODULE 1: PRODUCTION DATA INGESTION & EMPIRICAL COMPETITION ENGINE (REVISED)
=============================================================================
Title: Optimizing Wheat-Berseem Intercrop Yields: Machine Learning Analysis 
       of Morphological and Yield Data under Varied Water Irrigation Regimes

Updates Grounded in Experimental Methodology & Protocol Milestones:
1. Plant Spacing & Sowing Density: 15 cm zig-zag row pattern (5–7 wheat seeds & 7–10 clover seeds/hole).
2. Superphosphate Soil Pre-treatment: 238 kg/ha for inorganic-F, 87.5 kg/ha for Mixture.
3. Fertigation Regimes (Weekly Sun & Wed top-ups):
   - Freshwater (FW): Baseline control.
   - Inorganic-F: Urea (958 g/1000 L) + Potassium sulphate (788 g/1000 L).
   - Mixture: 50% Fish Effluent + Urea (479 g/1000 L) + Potassium sulphate (394 g/1000 L).
   - Fish Effluent: 100% Aquaculture drainage from 1 m³ fish tanks (stocked with fish).
4. Exact Seasonal Irrigation Depths:
   - 100% Full FC: 987.62 mm total over season (124,440 L per 126 m² main plot = 9,876.19 m³/ha).
   - 50% Deficit FC: 569.90 mm total over season (71,808 L per 126 m² main plot = 5,699.05 m³/ha).
5. Yield Measurement Units: Biomass recorded as weight (g)/hill for both Wheat and Berseem (dried at 60°C for 72h).
"""

import os
import re
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings('ignore')


class WheatBerseemFeatureEngine:
    def __init__(self, wheat_duration_days=150, berseem_duration_days=180, gdd_base_temp=5.0):
        self.t_w = wheat_duration_days
        self.t_b = berseem_duration_days
        self.T_total = max(wheat_duration_days, berseem_duration_days)
        self.gdd_base_temp = gdd_base_temp
        
        # Experimental electrical conductivity estimates (dS/m)
        self.wq_ec_map = {
            'FW': 0.8,
            'inorganic-F': 1.5,
            'Mixure': 3.5,
            'Fish_effluent': 6.8
        }
        

        self.irrigation_volume_map = {
            'full': 9876.19,
            'half': 5699.05
        }

    def process_weather_data(self, weather_source):
        print("[+] Processing Micro-Meteorological & Weather Observations...")
        if isinstance(weather_source, pd.DataFrame):
            df_raw = weather_source.copy()
        elif isinstance(weather_source, str):
            df_raw = pd.read_csv(weather_source) if weather_source.endswith('.csv') else pd.read_excel(weather_source)
        else:
            raise TypeError("Weather source must be a file path or DataFrame.")

        df_raw = df_raw.dropna(subset=['Date'])
        df_raw = df_raw[~df_raw['Date'].astype(str).str.startswith((',', 'mbar', 'Average', 'Total'))]
        df_raw['Date'] = pd.to_datetime(df_raw['Date'], errors='coerce')
        df_raw = df_raw.dropna(subset=['Date']).sort_values('Date').reset_index(drop=True)

        for c in df_raw.columns:
            if c != 'Date': 
                df_raw[c] = pd.to_numeric(df_raw[c], errors='coerce')

        def assign_season(dt):
            if dt >= pd.Timestamp('2022-12-18') and dt <= pd.Timestamp('2023-04-30'):
                return 1
            elif dt >= pd.Timestamp('2023-12-18') and dt <= pd.Timestamp('2024-04-30'):
                return 2
            else:
                return 1 if dt.year == 2022 or (dt.year == 2023 and dt.month <= 4) else 2

        df_raw['Season'] = df_raw['Date'].apply(assign_season)

        t_max = df_raw['Air_Temp_C_Max']
        t_min = df_raw['Air_Temp_C_Min']
        rh_max = df_raw['RH_Max']
        rh_min = df_raw['RH_Min']

        # Psychrometric equations for VPD and GDD
        es_tmax = 0.61078 * np.exp((17.27 * t_max) / (t_max + 237.3))
        es_tmin = 0.61078 * np.exp((17.27 * t_min) / (t_min + 237.3))
        es = (es_tmax + es_tmin) / 2.0
        ea = (es_tmin * (rh_max / 100.0) + es_tmax * (rh_min / 100.0)) / 2.0
        df_raw['VPD_Daily_kPa'] = np.maximum(0.0, es - ea)
        df_raw['GDD_Daily'] = np.maximum(0.0, ((t_max + t_min) / 2.0) - self.gdd_base_temp)

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

        mean_vpd = seasonal_weather['VPD_kPa'].mean()
        mean_gdd = seasonal_weather['GDD_Cumulative_Season'].mean()
        mean_temp = seasonal_weather['AirTempC_Avg'].mean()

        seasonal_weather['Delta_VPD_kPa'] = seasonal_weather['VPD_kPa'] - mean_vpd
        seasonal_weather['Delta_GDD_Season'] = seasonal_weather['GDD_Cumulative_Season'] - mean_gdd
        seasonal_weather['Delta_Temp_Avg'] = seasonal_weather['AirTempC_Avg'] - mean_temp

        print(f"[+] Aggregated Climate Profiles across {len(seasonal_weather)} Seasons.")
        return seasonal_weather

    def load_and_clean_wheat_data(self, wheat_source):
        print("[+] Ingesting and Cleaning Wheat Dataset (Yield in g/hill)...")
        if isinstance(wheat_source, pd.DataFrame):
            df = wheat_source.copy()
        elif isinstance(wheat_source, str):
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
            raise TypeError("Wheat source must be a file path or DataFrame.")

        df['Season'] = pd.to_numeric(df['Season'], errors='coerce').fillna(2).astype(int)
        df['Cropping_system'] = df['Cropping_system'].astype(str).str.strip().replace({
            'WM': 'MW', 'Wheat_monocrop': 'MW', 'Wheat-clover_intercrop': 'WBI'
        })
        df['Field_capacity'] = df['Field_capacity'].astype(str).str.strip().str.lower().replace({
            '100': 'full', '50': 'half', 'full': 'full', 'half': 'half'
        })
        df['Water_quality'] = df['Water_quality'].astype(str).str.strip().replace({
            'Fish-effluent': 'Fish_effluent', 'Inorganic_F': 'inorganic-F', 'Inorganic-F': 'inorganic-F',
            'Mixture': 'Mixure', 'FW ': 'FW'
        })

        if 'Replication' not in df.columns:
            df['Replication'] = df.groupby(['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']).cumcount() + 1
        else:
            df['Replication'] = pd.to_numeric(df['Replication'], errors='coerce')
            df['Replication'] = df['Replication'].fillna(df.groupby(['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']).cumcount() + 1).astype(int)

        # Standardize naming: Shoot dry weight is measured in g/hill
        df.rename(columns={'PH': 'Wheat_PH_cm', 'SPAD': 'Wheat_SPAD', 'Shoot_DW': 'Wheat_Shoot_DW_g_hill'}, inplace=True)
        
        df['Wheat_PH_cm'] = pd.to_numeric(df['Wheat_PH_cm'], errors='coerce')
        df.loc[df['Wheat_PH_cm'] > 200, 'Wheat_PH_cm'] = np.nan
        df['#Plant_Per_Hill'] = pd.to_numeric(df['#Plant_Per_Hill'], errors='coerce')
        df.loc[df['#Plant_Per_Hill'] > 50, '#Plant_Per_Hill'] = np.nan

        num_w_traits = [
            '#Plant_Per_Hill', '#Tillers_Per_Hill', 'Tillers_per_plant', 'SL (cm)', 'SW (cm)',
            'Wheat_SPAD', 'Wheat_PH_cm', 'StD', 'LL', 'LW', 'LA', 'LAI', 
            'Wheat_Shoot_DW_g_hill', 'Root_DW', 'RL', 'RW', 'HKW', 'TKW'
        ]
        treatment_keys = ['Season', 'Cropping_system', 'Field_capacity', 'Water_quality']
        for col in num_w_traits:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                df[col] = df.groupby(treatment_keys)[col].transform(lambda s: s.fillna(s.median()))
                df[col] = df.groupby(['Cropping_system', 'Field_capacity'])[col].transform(lambda s: s.fillna(s.median()))
                df[col] = df[col].fillna(df[col].median())

        print(f"[+] Processed Wheat Matrix: {df.shape} (MW: {(df['Cropping_system']=='MW').sum()}, WBI: {(df['Cropping_system']=='WBI').sum()})")
        return df

    def load_and_clean_berseem_data(self, berseem_source):
        print("[+] Ingesting and Processing Berseem Dataset (Yield in g/hill across 3 Cuttings)...")
        if isinstance(berseem_source, pd.DataFrame):
            df = berseem_source.copy()
        elif isinstance(berseem_source, str):
            df = pd.read_csv(berseem_source) if berseem_source.endswith('.csv') else pd.read_excel(berseem_source)
        else:
            raise TypeError("Berseem source must be a file path or DataFrame.")

        rename_map = {
            'Cut No': 'Cut_No', 'Cut_no': 'Cut_No', 'Fertilizer': 'Water_quality', 'Irrigation': 'Field_capacity',
            'Cultivation system': 'Cropping_system', 'Cultivation_system': 'Cropping_system',
            'Cropping_System': 'Cropping_system', 'Rep': 'Replication', 'Block': 'Replication',
            'Height': 'Berseem_Height_cm', 'SPAD': 'Berseem_SPAD', 'Dry_weight': 'Berseem_Dry_Weight_g_hill'
        }
        df.rename(columns=rename_map, inplace=True)

        df['Season'] = pd.to_numeric(df['Season'], errors='coerce').fillna(2).astype(int)

        # Recognize MW in berseem data as Monoculture Berseem (MB)
        df['Cropping_system'] = df['Cropping_system'].astype(str).str.strip().replace({
            'MW': 'MB', 'WM': 'MB', 'Sole_Berseem': 'MB', 'Wheat-clover_intercrop': 'WBI'
        })
        df['Field_capacity'] = df['Field_capacity'].astype(str).str.strip().str.lower().replace({
            '100': 'full', '50': 'half', 'full': 'full', 'half': 'half'
        })
        df['Water_quality'] = df['Water_quality'].astype(str).str.strip().replace({
            'Fish-effluent': 'Fish_effluent', 'Inorganic_F': 'inorganic-F', 'Inorganic-F': 'inorganic-F',
            'Mixture': 'Mixure', 'FW ': 'FW'
        })

        if 'Replication' in df.columns:
            df['Replication'] = pd.to_numeric(df['Replication'], errors='coerce').fillna(1).astype(int)
        else:
            df['Replication'] = df.groupby(['Season', 'Cut_No', 'Cropping_system', 'Field_capacity', 'Water_quality']).cumcount() + 1

        num_b_traits = ['Berseem_Height_cm', 'Berseem_SPAD', 'Berseem_Dry_Weight_g_hill']
        impute_keys = ['Season', 'Cut_No', 'Cropping_system', 'Field_capacity', 'Water_quality']
        for col in num_b_traits:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                df[col] = df.groupby(impute_keys)[col].transform(lambda s: s.fillna(s.median()))
                df[col] = df.groupby(['Cropping_system', 'Field_capacity'])[col].transform(lambda s: s.fillna(s.median()))
                df[col] = df[col].fillna(0.0)

        plot_keys = ['Season', 'Replication', 'Cropping_system', 'Field_capacity', 'Water_quality']
        
        # Multi-Cut Aggregations (g/hill)
        summary_dw = df.groupby(plot_keys)['Berseem_Dry_Weight_g_hill'].agg(['sum', 'mean', 'max', 'std']).reset_index()
        summary_dw.rename(columns={
            'sum': 'Berseem_Dry_Weight_g_hill',
            'mean': 'Berseem_Mean_Cut_Dry_Weight_g_hill',
            'max': 'Berseem_Peak_Cut_Dry_Weight_g_hill',
            'std': 'Berseem_Cut_Dry_Weight_Std'
        }, inplace=True)

        summary_ht = df.groupby(plot_keys)['Berseem_Height_cm'].agg(['mean', 'max']).reset_index()
        summary_ht.rename(columns={'mean': 'Berseem_Height_cm', 'max': 'Berseem_Peak_Height_cm'}, inplace=True)

        summary_spad = df.groupby(plot_keys)['Berseem_SPAD'].agg(['mean', 'min']).reset_index()
        summary_spad.rename(columns={'mean': 'Berseem_SPAD', 'min': 'Berseem_Min_SPAD'}, inplace=True)

        summary = summary_dw.merge(summary_ht, on=plot_keys).merge(summary_spad, on=plot_keys)

        pivot_dw = df.pivot(index=plot_keys, columns='Cut_No', values='Berseem_Dry_Weight_g_hill').reset_index()
        pivot_dw.columns = [f"Berseem_{col}_DW_g_hill" if col not in plot_keys else col for col in pivot_dw.columns]

        pivot_ht = df.pivot(index=plot_keys, columns='Cut_No', values='Berseem_Height_cm').reset_index()
        pivot_ht.columns = [f"Berseem_{col}_Height_cm" if col not in plot_keys else col for col in pivot_ht.columns]

        pivot_spad = df.pivot(index=plot_keys, columns='Cut_No', values='Berseem_SPAD').reset_index()
        pivot_spad.columns = [f"Berseem_{col}_SPAD" if col not in plot_keys else col for col in pivot_spad.columns]

        pivoted_berseem = summary.merge(pivot_dw, on=plot_keys).merge(pivot_ht, on=plot_keys).merge(pivot_spad, on=plot_keys)

        dw_cut_cols = [c for c in pivot_dw.columns if 'Berseem_Cut_' in c and '_DW_g_hill' in c]
        if len(dw_cut_cols) >= 2:
            first_cut, last_cut = dw_cut_cols[0], dw_cut_cols[-1]
            pivoted_berseem['Berseem_Cutting_Recovery_Ratio'] = np.where(
                pivoted_berseem[first_cut] > 0,
                pivoted_berseem[last_cut] / pivoted_berseem[first_cut], 0.0
            )
            pivoted_berseem['Berseem_Cut_Stability_CV'] = np.where(
                pivoted_berseem['Berseem_Mean_Cut_Dry_Weight_g_hill'] > 0,
                pivoted_berseem['Berseem_Cut_Dry_Weight_Std'].fillna(0) / pivoted_berseem['Berseem_Mean_Cut_Dry_Weight_g_hill'], 0.0
            )

        print(f"[+] Processed Berseem Matrix: {pivoted_berseem.shape} (MB: {(pivoted_berseem['Cropping_system']=='MB').sum()}, WBI: {(pivoted_berseem['Cropping_system']=='WBI').sum()})")
        return pivoted_berseem

    def integrate_and_calculate_all_metrics(self, df_w, df_b, df_weather=None):
        print("[+] Integrating Datasets and Calibrating Empirical Baselines (g/hill)...")
        match_keys = ['Season', 'Field_capacity', 'Water_quality']

        sole_w = df_w[df_w['Cropping_system'] == 'MW']
        sole_b = df_b[df_b['Cropping_system'] == 'MB']

        ref_w = sole_w.groupby(match_keys)['Wheat_Shoot_DW_g_hill'].mean().reset_index().rename(
            columns={'Wheat_Shoot_DW_g_hill': 'Y_SW_Baseline'}
        )
        ref_b = sole_b.groupby(match_keys)['Berseem_Dry_Weight_g_hill'].mean().reset_index().rename(
            columns={'Berseem_Dry_Weight_g_hill': 'Y_SB_Baseline'}
        )

        w_wbi = df_w[df_w['Cropping_system'] == 'WBI'].copy()
        b_wbi = df_b[df_b['Cropping_system'] == 'WBI'].copy()
        
        join_keys = ['Season', 'Replication', 'Cropping_system', 'Field_capacity', 'Water_quality']
        merged_wbi = pd.merge(w_wbi, b_wbi, on=join_keys, how='inner')

        sole_w_full = df_w[df_w['Cropping_system'] == 'MW'].copy()
        sole_b_full = df_b[df_b['Cropping_system'] == 'MB'].copy()

        merged = pd.concat([merged_wbi, sole_w_full, sole_b_full], ignore_index=True)

        berseem_cols = [c for c in merged.columns if 'Berseem_' in c]
        for col in berseem_cols: 
            merged[col] = merged[col].fillna(0.0)

        num_w_traits = [
            '#Plant_Per_Hill', '#Tillers_Per_Hill', 'Tillers_per_plant', 'SL (cm)', 'SW (cm)',
            'Wheat_SPAD', 'Wheat_PH_cm', 'StD', 'LL', 'LW', 'LA', 'LAI', 
            'Wheat_Shoot_DW_g_hill', 'Root_DW', 'RL', 'RW', 'HKW', 'TKW'
        ]
        for col in num_w_traits:
            if col in merged.columns: 
                merged[col] = merged[col].fillna(0.0)

        merged = merged.merge(ref_w, on=match_keys, how='left')
        merged = merged.merge(ref_b, on=match_keys, how='left')

        # True Empirical Partial and Total LER
        merged['pLER_Wheat'] = np.where(
            merged['Cropping_system'] == 'WBI',
            merged['Wheat_Shoot_DW_g_hill'] / (merged['Y_SW_Baseline'] + 1e-8),
            np.where(merged['Cropping_system'] == 'MW', 1.0, 0.0)
        )

        merged['pLER_Berseem'] = np.where(
            merged['Cropping_system'] == 'WBI',
            merged['Berseem_Dry_Weight_g_hill'] / (merged['Y_SB_Baseline'] + 1e-8),
            np.where(merged['Cropping_system'] == 'MB', 1.0, 0.0)
        )

        merged['LER_Total'] = merged['pLER_Wheat'] + merged['pLER_Berseem']
        merged['ATER'] = ((merged['pLER_Wheat'] * self.t_w) + (merged['pLER_Berseem'] * self.t_b)) / self.T_total

        z_w, z_b = 0.50, 0.50  # Equal 50:50 strip proportion
        merged['CR_Wheat'] = np.where(
            (merged['pLER_Berseem'] > 0) & (merged['Cropping_system'] == 'WBI'),
            (merged['pLER_Wheat'] / merged['pLER_Berseem']) * (z_b / z_w), 0.0
        )
        merged['Aggressivity_Wheat'] = np.where(
            merged['Cropping_system'] == 'WBI',
            (merged['pLER_Wheat'] / z_w) - (merged['pLER_Berseem'] / z_b), 0.0
        )

        mean_sw = merged['Y_SW_Baseline'].mean()
        mean_sb = merged['Y_SB_Baseline'].mean()
        merged['SPI'] = np.where(
            merged['Cropping_system'] == 'WBI',
            merged['Wheat_Shoot_DW_g_hill'] + (merged['Berseem_Dry_Weight_g_hill'] * (mean_sw / (mean_sb + 1e-8))),
            np.where(merged['Cropping_system'] == 'MW', merged['Wheat_Shoot_DW_g_hill'], merged['Berseem_Dry_Weight_g_hill'] * (mean_sw / (mean_sb + 1e-8)))
        )

        # Total System Biomass in g/hill
        merged['Total_System_Biomass_g_hill'] = merged['Wheat_Shoot_DW_g_hill'] + merged['Berseem_Dry_Weight_g_hill']
        
        # Exact Irrigation Applied (m³/ha from methodology table)
        merged['Applied_Water_m3_ha'] = merged['Field_capacity'].map(self.irrigation_volume_map).fillna(9876.19)
        merged['Applied_Water_mm'] = np.where(merged['Field_capacity'] == 'half', 569.90, 987.62)

        # IWUE in kg dry matter / m³ water (calculated via area hill density factor: 44,444 hills/ha at 15cm spacing)
        # Yield (t/ha) = Yield (g/hill) * 44444 / 1e6 = Yield (g/hill) * 0.044444
        # IWUE (kg/m³) = Yield (t/ha) * 1000 / Applied_Water (m³/ha)
        hills_per_ha = 44444.44
        merged['Total_Biomass_t_ha_equiv'] = (merged['Total_System_Biomass_g_hill'] * hills_per_ha) / 1e6
        merged['IWUE_kg_m3'] = (merged['Total_Biomass_t_ha_equiv'] * 1000.0) / (merged['Applied_Water_m3_ha'] + 1e-8)

        merged['Field_capacity_pct'] = np.where(merged['Field_capacity'] == 'half', 50.0, 100.0)
        merged['Salinity_ECw_Est_dS_m'] = merged['Water_quality'].map(self.wq_ec_map).fillna(1.0)
        merged['Stress_Severity_Index'] = merged['Salinity_ECw_Est_dS_m'] / (merged['Field_capacity_pct'] / 100.0)

        merged['WQ_Fish_effluent'] = np.where(merged['Water_quality'] == 'Fish_effluent', 1, 0)
        merged['WQ_inorganic_F'] = np.where(merged['Water_quality'] == 'inorganic-F', 1, 0)
        merged['WQ_Mixure'] = np.where(merged['Water_quality'] == 'Mixure', 1, 0)
        merged['WQ_FW'] = np.where(merged['Water_quality'] == 'FW', 1, 0)
        merged['System_WBI'] = np.where(merged['Cropping_system'] == 'WBI', 1, 0)

        merged['Height_Differential_cm'] = merged['Wheat_PH_cm'] - merged['Berseem_Height_cm']
        merged['SPAD_Ratio_Wheat_Berseem'] = np.where(merged['Berseem_SPAD'] > 0, merged['Wheat_SPAD'] / merged['Berseem_SPAD'], 0.0)
        merged['Total_System_SPAD'] = merged['Wheat_SPAD'] + merged['Berseem_SPAD']
        merged['Wheat_Root_Shoot_Ratio'] = np.where(merged['Wheat_Shoot_DW_g_hill'] > 0, merged['Root_DW'] / (merged['Wheat_Shoot_DW_g_hill'] + 1e-8), 0.0)

        # Backward-compatible column aliases
        merged['Wheat_Shoot_DW_t_ha'] = merged['Wheat_Shoot_DW_g_hill']
        merged['Berseem_Dry_Weight_t_ha'] = merged['Berseem_Dry_Weight_g_hill']
        merged['Total_System_Biomass_t_ha'] = merged['Total_System_Biomass_g_hill']

        if df_weather is not None and 'Season' in df_weather.columns:
            merged = pd.merge(merged, df_weather, on='Season', how='left')

        print(f"[+] Full Processing Complete! Matrix Shape: {merged.shape} with {len(merged.columns)} Features.")
        return merged

    def run_pipeline(self, wheat_source, berseem_source, weather_source=None, save_csv_path=None):
        print("=====================================================================")
        print("MODULE 1: ROBUST DATA INGESTION & EMPIRICAL COMPETITION ENGINE")
        print("=====================================================================")
        df_weath_summary = self.process_weather_data(weather_source) if weather_source is not None else None
        df_w_clean = self.load_and_clean_wheat_data(wheat_source)
        df_b_clean = self.load_and_clean_berseem_data(berseem_source)
        final_df = self.integrate_and_calculate_all_metrics(df_w_clean, df_b_clean, df_weath_summary)

        if save_csv_path:
            os.makedirs(os.path.dirname(os.path.abspath(save_csv_path)), exist_ok=True)
            final_df.to_csv(save_csv_path, index=False)
            print(f"[+] Exported Clean Empirical Dataset to: '{save_csv_path}'")

        return final_df


if __name__ == '__main__':
    engine = WheatBerseemFeatureEngine()
    final_df = engine.run_pipeline(
        wheat_source='data/wheat_data.csv',
        berseem_source='data/barseem_data.csv',
        weather_source='data/weather_daily_data.csv',
        save_csv_path='data/wheat_berseem_real_processed.csv'
    )
