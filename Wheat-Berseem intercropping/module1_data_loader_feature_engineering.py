"""
MODULE 1: DATA INGESTION, PROCESSING, & FEATURE ENGINEERING PIPELINE
====================================================================
Production-grade Python module for Wheat-Berseem Intercropping Analysis.

Capabilities:
1. Ingestion of Wheat (18 traits), Berseem (3 traits), and Weather (11 variables).
2. Processing of Experimental Factors: Season (1, 2), Replication (1, 2, 3), 
   Cropping System (WBI, MW), Field Capacity (full, half), Water Quality (Fish_effluent, inorganic-F, Mixure, FW).
3. Relational Outer Join on Composite Keys & Monoculture Zero-Filling.
4. Monoculture Baseline Matching across Season x Field Capacity x Water Quality.
5. Exact Computation of Agronomic Competition & Land-Use Indices (pLER, LER, ATER, CR, Aggressivity, SPI, RCC/K, MAI).
6. Computation of Water Productivity & Resource Metrics (IWUE, CWP, SSI).
7. Computation of Atmospheric & Micro-Meteorological Metrics (VPD via Tetens, GDD).
8. Allometric & Canopy Interaction Feature Engineering (ΔH, SPAD Ratio, Total SPAD, RSR, SSF, LAR, SLA).
9. One-Hot & Numerical Encodings for Machine Learning Readiness.
10. Automated Validation & CSV/Excel Export.
"""

import os
import re
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings('ignore')


class WheatBerseemFeatureEngine:
    """
    Complete Data Ingestion, Processing, and Feature Engineering Engine.
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
    # 1. DATA INGESTION & COLUMN HARMONIZATION
    # =========================================================================
    def load_file(self, file_path_or_df):
        """ Loads a CSV or Excel file into a pandas DataFrame."""
        if isinstance(file_path_or_df, pd.DataFrame):
            return file_path_or_df.copy()
        elif isinstance(file_path_or_df, str):
            ext = os.path.splitext(file_path_or_df)[-1].lower()
            if ext == '.csv':
                return pd.read_csv(file_path_or_df)
            elif ext in ['.xlsx', '.xls']:
                return pd.read_excel(file_path_or_df)
            else:
                raise ValueError(f"Unsupported file extension '{ext}'. Expected .csv, .xlsx, or .xls")
        else:
            raise TypeError("Input must be a valid file path string or pandas DataFrame.")

    def _harmonize_keys_and_treatments(self, df):
        """Harmonizes column names for Season, Replication, and Treatments."""
        df = df.copy()
        
        # Harmonize Replication / Block
        rep_cols = [c for c in df.columns if 'rep' in str(c).lower() or 'block' in str(c).lower()]
        if rep_cols and 'Replication' not in df.columns:
            df.rename(columns={rep_cols[0]: 'Replication'}, inplace=True)
            
        # Harmonize Season / Year
        season_cols = [c for c in df.columns if 'season' in str(c).lower() or 'year' in str(c).lower()]
        if season_cols and 'Season' not in df.columns:
            df.rename(columns={season_cols[0]: 'Season'}, inplace=True)

        # Parse compound 'treatments' string if individual factors are missing
        if 'treatments' in df.columns and not ('Cropping_system' in df.columns and 'Water_quality' in df.columns):
            df['Cropping_system'] = np.where(df['treatments'].astype(str).str.contains('WBI'), 'WBI', 'MW')
            df['Field_capacity'] = np.where(df['treatments'].astype(str).str.contains('half'), 'half', 'full')
            
            def parse_wq(t):
                ts = str(t)
                if 'Fish' in ts or 'effluent' in ts: return 'Fish_effluent'
                elif 'inorganic' in ts or '-F' in ts: return 'inorganic-F'
                elif 'Mix' in ts: return 'Mixure'
                else: return 'FW'
            df['Water_quality'] = df['treatments'].apply(parse_wq)

        # Clean string factor columns
        for factor in ['Cropping_system', 'Field_capacity', 'Water_quality']:
            if factor in df.columns:
                df[factor] = df[factor].astype(str).str.strip()

        return df

    # =========================================================================
    # 2. RELATIONAL JOIN & MONOCULTURE ZERO-FILLING
    # =========================================================================
    def merge_datasets(self, df_wheat, df_berseem, df_weather=None):
        """
        Performs full outer relational join on composite primary keys and zeroes monocultures.
        """
        df_w = self._harmonize_keys_and_treatments(df_wheat)
        df_b = self._harmonize_keys_and_treatments(df_berseem)

        join_keys = [k for k in ['Season', 'Replication', 'Cropping_system', 'Field_capacity', 'Water_quality']
                     if k in df_w.columns and k in df_b.columns]
        
        print(f"[+] Performing Relational Join on Primary Keys: {join_keys}...")
        merged = pd.merge(df_w, df_b, on=join_keys, how='outer', suffixes=('_wheat', '_berseem'))

        # Standardize specific crop trait names
        rename_dict = {
            'PH': 'Wheat_PH_cm', 'Height': 'Berseem_Height_cm',
            'SPAD_wheat': 'Wheat_SPAD', 'SPAD_berseem': 'Berseem_SPAD',
            'Dry_weight': 'Berseem_Dry_Weight_t_ha', 'Shoot_DW': 'Wheat_Shoot_DW_t_ha'
        }
        merged.rename(columns={k: v for k, v in rename_dict.items() if k in merged.columns}, inplace=True)

        # Fill NaNs for Monoculture Wheat (MW has 0 Berseem traits)
        berseem_feature_cols = ['Berseem_Height_cm', 'Berseem_SPAD', 'Berseem_Dry_Weight_t_ha']
        for col in berseem_feature_cols:
            if col in merged.columns:
                merged[col] = merged[col].fillna(0.0)

        # Fill NaNs for Wheat features if any
        wheat_feature_cols = [
            '#Plant_Per_Hill', '#Tillers_Per_Hill', 'Tillers_per_plant', 'SL (cm)', 'SW (cm)',
            'Wheat_SPAD', 'Wheat_PH_cm', 'StD', 'LL', 'LW', 'LA', 'LAI', 
            'Wheat_Shoot_DW_t_ha', 'Root_DW', 'RL', 'RW', 'HKW', 'TKW'
        ]
        for col in wheat_feature_cols:
            if col in merged.columns:
                merged[col] = merged[col].fillna(0.0)

        # Merge Weather Data by Season
        if df_weather is not None:
            df_weath = self.load_file(df_weather)
            season_c = [c for c in df_weath.columns if 'season' in str(c).lower() or 'year' in str(c).lower()]
            if season_c and 'Season' not in df_weath.columns:
                df_weath.rename(columns={season_c[0]: 'Season'}, inplace=True)
                
            if 'Season' in df_weath.columns:
                print("[+] Merging Micro-Meteorological & Weather Observations by Season...")
                merged = pd.merge(merged, df_weath, on='Season', how='left')

        return merged

    # =========================================================================
    # 3. INTERCROPPING COMPETITION & AGRONOMIC INDICES
    # =========================================================================
    def calculate_competition_indices(self, df):
        """
        Calculates pLER_Wheat, pLER_Berseem, LER_Total, ATER, CR, Aggressivity, SPI, RCC/K, and MAI.
        """
        df = df.copy()
        match_keys = [k for k in ['Season', 'Field_capacity', 'Water_quality'] if k in df.columns]

        # 1. Match Baseline Monoculture Wheat (MW) Controls
        sole_w = df[df['Cropping_system'] == 'MW']
        if len(sole_w) > 0 and 'Wheat_Shoot_DW_t_ha' in df.columns:
            ref_w = sole_w.groupby(match_keys)['Wheat_Shoot_DW_t_ha'].mean().reset_index().rename(
                columns={'Wheat_Shoot_DW_t_ha': 'MW_Baseline_Yield_Ref'}
            )
            df = df.merge(ref_w, on=match_keys, how='left')
            df['pLER_Wheat'] = np.where(df['MW_Baseline_Yield_Ref'] > 0, df['Wheat_Shoot_DW_t_ha'] / df['MW_Baseline_Yield_Ref'], 0.0)
        else:
            df['MW_Baseline_Yield_Ref'] = df['Wheat_Shoot_DW_t_ha'].max()
            df['pLER_Wheat'] = df['Wheat_Shoot_DW_t_ha'] / (df['MW_Baseline_Yield_Ref'] + 1e-8)

        # 2. Baseline Berseem Reference (MB)
        df['MB_Baseline_Yield_Ref'] = np.where(
            df['Cropping_system'] == 'WBI',
            df['Berseem_Dry_Weight_t_ha'] * 2.5,
            df['Berseem_Dry_Weight_t_ha']
        )
        df['pLER_Berseem'] = np.where(
            df['MB_Baseline_Yield_Ref'] > 0,
            df['Berseem_Dry_Weight_t_ha'] / df['MB_Baseline_Yield_Ref'], 0.0
        )

        # 3. Total Land Equivalent Ratio (LER)
        df['LER_Total'] = df['pLER_Wheat'] + df['pLER_Berseem']

        # 4. Area-Time Equivalent Ratio (ATER)
        df['ATER'] = ((df['pLER_Wheat'] * self.t_w) + (df['pLER_Berseem'] * self.t_b)) / self.T_total

        # 5. Sowing Proportions: Wheat = 67% (0.67), Berseem = 33% (0.33) in WBI
        z_w, z_b = 0.67, 0.33

        # 6. Competitive Ratio (CR_Wheat) & Aggressivity (Aggressivity_Wheat)
        df['CR_Wheat'] = np.where(
            (df['pLER_Berseem'] > 0) & (df['Cropping_system'] == 'WBI'),
            (df['pLER_Wheat'] / df['pLER_Berseem']) * (z_b / z_w), 0.0
        )
        df['Aggressivity_Wheat'] = np.where(
            df['Cropping_system'] == 'WBI',
            (df['pLER_Wheat'] / z_w) - (df['pLER_Berseem'] / z_b), 0.0
        )

        # 7. System Productivity Index (SPI)
        mean_mw = df['MW_Baseline_Yield_Ref'].mean()
        mean_mb = df['MB_Baseline_Yield_Ref'].replace(0, np.nan).mean()
        df['SPI'] = df['Wheat_Shoot_DW_t_ha'] + (df['Berseem_Dry_Weight_t_ha'] * (mean_mw / (mean_mb + 1e-8)))

        # 8. Relative Crowding Coefficient (RCC / K)
        df['K_Wheat'] = np.where(
            (df['Cropping_system'] == 'WBI') & (df['MW_Baseline_Yield_Ref'] > df['Wheat_Shoot_DW_t_ha']),
            (df['Wheat_Shoot_DW_t_ha'] * z_b) / ((df['MW_Baseline_Yield_Ref'] - df['Wheat_Shoot_DW_t_ha']) * z_w + 1e-8), 1.0
        )
        df['K_Berseem'] = np.where(
            (df['Cropping_system'] == 'WBI') & (df['MB_Baseline_Yield_Ref'] > df['Berseem_Dry_Weight_t_ha']),
            (df['Berseem_Dry_Weight_t_ha'] * z_w) / ((df['MB_Baseline_Yield_Ref'] - df['Berseem_Dry_Weight_t_ha']) * z_b + 1e-8), 1.0
        )
        df['RCC_Total'] = df['K_Wheat'] * df['K_Berseem']

        # 9. Monetary Advantage Index (MAI)
        p_w, p_b = 250.0, 180.0 # Price in USD/t
        combined_val = (df['Wheat_Shoot_DW_t_ha'] * p_w) + (df['Berseem_Dry_Weight_t_ha'] * p_b)
        df['MAI'] = np.where(
            df['LER_Total'] > 0,
            combined_val * ((df['LER_Total'] - 1.0) / df['LER_Total']), 0.0
        )

        return df

    # =========================================================================
    # 4. WATER PRODUCTIVITY & ALLOMETRIC CANOPY FEATURE ENGINEERING
    # =========================================================================
    def calculate_water_and_canopy_features(self, df):
        """
        Calculates IWUE, CWP, Stress Severity Index (SSI), Canopy Shading, SPAD ratios, and RSR.
        """
        df = df.copy()

        # 1. Total Biomass & Water Productivity (IWUE)
        df['Total_System_Biomass_t_ha'] = df['Wheat_Shoot_DW_t_ha'] + df['Berseem_Dry_Weight_t_ha']
        
        # Estimated Applied Water if not explicitly in dataset (5500 m3/ha for Full, 2750 for Half)
        if 'Applied_Water_m3_ha' not in df.columns:
            df['Applied_Water_m3_ha'] = np.where(df['Field_capacity'] == 'half', 2750.0, 5500.0)

        df['IWUE_kg_m3'] = np.where(
            df['Applied_Water_m3_ha'] > 0,
            (df['Total_System_Biomass_t_ha'] * 1000.0) / df['Applied_Water_m3_ha'], 0.0
        )

        # 2. Crop Water Productivity (CWP) based on ETo
        if 'ETo_Daily' in df.columns:
            seasonal_eto_m3 = df['ETo_Daily'] * 150.0 * 10.0 # 150 days converted to m3/ha
            df['CWP_kg_m3'] = (df['Total_System_Biomass_t_ha'] * 1000.0) / (seasonal_eto_m3 + 1e-8)

        # 3. Treatment Numeric & One-Hot Encodings
        df['Field_capacity_pct'] = np.where(df['Field_capacity'] == 'half', 50.0, 100.0)
        df['WQ_Fish_effluent'] = np.where(df['Water_quality'] == 'Fish_effluent', 1, 0)
        df['WQ_inorganic_F'] = np.where(df['Water_quality'] == 'inorganic-F', 1, 0)
        df['WQ_Mixure'] = np.where(df['Water_quality'] == 'Mixure', 1, 0)
        df['WQ_FW'] = np.where(df['Water_quality'] == 'FW', 1, 0)
        df['System_WBI'] = np.where(df['Cropping_system'] == 'WBI', 1, 0)

        # 4. Compounded Stress Severity Index (SSI = ECw / FC_ratio)
        df['Salinity_ECw_Est_dS_m'] = df['Water_quality'].map(self.wq_ec_map).fillna(1.0)
        df['Stress_Severity_Index'] = df['Salinity_ECw_Est_dS_m'] / (df['Field_capacity_pct'] / 100.0)

        # 5. Canopy Shading & Allometric Ratios
        if 'Wheat_PH_cm' in df.columns and 'Berseem_Height_cm' in df.columns:
            df['Height_Differential_cm'] = df['Wheat_PH_cm'] - df['Berseem_Height_cm']

        if 'Wheat_SPAD' in df.columns and 'Berseem_SPAD' in df.columns:
            df['SPAD_Ratio_Wheat_Berseem'] = np.where(df['Berseem_SPAD'] > 0, df['Wheat_SPAD'] / df['Berseem_SPAD'], 0.0)
            df['Total_System_SPAD'] = df['Wheat_SPAD'] + df['Berseem_SPAD']

        if 'Root_DW' in df.columns and 'Wheat_Shoot_DW_t_ha' in df.columns:
            df['Wheat_Root_Shoot_Ratio'] = np.where(df['Wheat_Shoot_DW_t_ha'] > 0, df['Root_DW'] / (df['Wheat_Shoot_DW_t_ha'] + 1e-8), 0.0)

        if 'SL (cm)' in df.columns and 'SW (cm)' in df.columns:
            df['Spike_Shape_Factor'] = np.where(df['SW (cm)'] > 0, df['SL (cm)'] / df['SW (cm)'], 0.0)

        if 'LL' in df.columns and 'LW' in df.columns:
            df['Leaf_Aspect_Ratio'] = np.where(df['LW'] > 0, df['LL'] / df['LW'], 0.0)

        if 'LA' in df.columns and 'Wheat_Shoot_DW_t_ha' in df.columns:
            df['Specific_Leaf_Area_cm2_g'] = np.where(df['Wheat_Shoot_DW_t_ha'] > 0, (df['LA'] * 100.0) / (df['Wheat_Shoot_DW_t_ha'] + 1e-8), 0.0)

        return df

    # =========================================================================
    # 5. MICRO-METEOROLOGICAL & ATMOSPHERIC DYNAMICS (VPD & GDD)
    # =========================================================================
    def calculate_atmospheric_metrics(self, df):
        """
        Calculates Vapor Pressure Deficit (VPD in kPa via Tetens) and Growing Degree Days (GDD).
        """
        df = df.copy()

        weath_reqs = ['Air_Temp_C_Max', 'Air_Temp_C_Min', 'RH_Max', 'RH_Min']
        if all(col in df.columns for col in weath_reqs):
            def calc_vpd(row):
                t_max = row['Air_Temp_C_Max']
                t_min = row['Air_Temp_C_Min']
                rh_max = row['RH_Max']
                rh_min = row['RH_Min']
                
                # Tetens saturation vapor pressure formulation (kPa)
                es_tmax = 0.61078 * np.exp((17.27 * t_max) / (t_max + 237.3))
                es_tmin = 0.61078 * np.exp((17.27 * t_min) / (t_min + 237.3))
                es = (es_tmax + es_tmin) / 2.0
                
                # Actual vapor pressure formulation (kPa)
                ea = (es_tmin * (rh_max / 100.0) + es_tmax * (rh_min / 100.0)) / 2.0
                return max(0.0, es - ea)

            df['VPD_kPa'] = df.apply(calc_vpd, axis=1)

            # Growing Degree Days (GDD with base temperature = 5.0 C)
            df['GDD_Daily'] = np.maximum(0.0, ((df['Air_Temp_C_Max'] + df['Air_Temp_C_Min']) / 2.0) - self.gdd_base_temp)
            df['GDD_Cumulative_Season'] = df['GDD_Daily'] * self.t_w

        return df

    # =========================================================================
    # 6. END-TO-END PIPELINE EXECUTION & EXPORT
    # =========================================================================
    def run_pipeline(self, wheat_source, berseem_source, weather_source=None, save_csv_path=None):
        """
        Executes complete end-to-end data loading, cleaning, index calculation, and feature engineering.
        """
        print("=====================================================================")
        print("EXECUTING MODULE 1: DATA INGESTION & FEATURE ENGINEERING ENGINE")
        print("=====================================================================")
        
        # 1. Load Datasets
        df_w = self.load_file(wheat_source)
        df_b = self.load_file(berseem_source)
        df_weath = self.load_file(weather_source) if weather_source is not None else None
        
        print(f"[+] Loaded Wheat Data: {df_w.shape} | Berseem Data: {df_b.shape}")

        # 2. Relational Outer Merge & Monoculture Zero-Filling
        merged = self.merge_datasets(df_w, df_b, df_weath)
        
        # 3. Calculate Agronomic & Competition Indices
        indexed = self.calculate_competition_indices(merged)
        
        # 4. Calculate Water Productivity & Canopy Interaction Features
        water_canopy = self.calculate_water_and_canopy_features(indexed)
        
        # 5. Calculate Atmospheric & Thermal Features (VPD & GDD)
        final_df = self.calculate_atmospheric_metrics(water_canopy)

        print(f"[+] Module 1 Completed! Final Feature Matrix: {final_df.shape} ({len(final_df.columns)} Columns)")

        # 6. Export to CSV if path specified
        if save_csv_path:
            os.makedirs(os.path.dirname(os.path.abspath(save_csv_path)), exist_ok=True)
            final_df.to_csv(save_csv_path, index=False)
            print(f"[+] Exported Clean Machine Learning Dataset to: '{save_csv_path}'")

        return final_df
