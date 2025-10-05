#!/usr/bin/env python3
"""
Data Overview Visualizations
Create comprehensive visualizations of the temperature sensor data to illustrate the analysis task.
"""

import pandas as pd
import numpy as np
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.signal import correlate
import warnings
warnings.filterwarnings('ignore')

# Sensor name mapping dictionary
SENSOR_NAME_MAPPING = {
    'ams': '3D printer closet',
    'bed': 'Bedroom',
    'couch': 'Living room',
    'desk': 'Office',
    'freezer': 'Freezer',
    'fridge': 'Fridge',
    'it': 'IT closet',
    'storage': 'NAS',
    'weather': 'Outdoors',
    'wine_fridge': 'Wine fridge'
}

def load_temperature_data():
    """Load and process temperature sensor data."""
    print("Loading temperature sensor data...")

    # Load data
    df = pd.read_csv('Data/history.csv')

    # Filter for temperature sensors only
    temp_sensors = [
        'sensor.ams_temperature',
        'sensor.bed_temperature',
        'sensor.couch_temperature',
        'sensor.desk_temperature',
        'sensor.freezer_temperature',
        'sensor.fridge_temperature',
        'sensor.it_temperature',
        'sensor.storage_temperature',
        'sensor.weather_temperature_sensor',
        'sensor.wine_fridge_temperature'
    ]

    df = df[df['entity_id'].isin(temp_sensors)]

    # Convert timestamp and state
    df['timestamp'] = pd.to_datetime(df['last_changed'])
    # Convert from UTC to California time (Pacific Time)
    df['timestamp'] = df['timestamp'].dt.tz_convert('US/Pacific')
    df['temperature'] = pd.to_numeric(df['state'], errors='coerce')

    # Create sensor short names
    df['sensor_short'] = df['entity_id'].str.replace('sensor.', '').str.replace('_temperature', '').str.replace('_sensor', '')

    # Apply sensor name mapping
    df['sensor'] = df['sensor_short'].map(SENSOR_NAME_MAPPING)

    # Remove outliers (beyond 3 standard deviations)
    df = df.groupby('sensor').apply(
        lambda x: x[np.abs(x['temperature'] - x['temperature'].mean()) <= 3 * x['temperature'].std()]
    ).reset_index(drop=True)

    print(f"Loaded {len(df)} temperature readings from {df['sensor'].nunique()} sensors")
    print(f"Date range: {df['timestamp'].min()} to {df['timestamp'].max()}")

    return df

def create_sensor_overview_plot(df):
    """Create overview plot showing all temperature sensors over time."""
    print("Creating sensor overview time series plot...")

    # Use all raw data points (no sampling/averaging)
    df_sample = df.copy()

    # Create facet grouping for sensors
    def assign_facet_group(sensor):
        if sensor == 'NAS':
            return 'NAS'
        elif sensor in ['Fridge', 'Freezer', 'Wine fridge']:
            return 'Refrigeration'
        else:
            return 'Environmental Sensors'

    df_sample['facet_group'] = df_sample['sensor'].apply(assign_facet_group)

    # Define facet order (top to bottom)
    facet_order = ['Environmental Sensors', 'Refrigeration', 'NAS']
    df_sample['facet_group'] = pd.Categorical(df_sample['facet_group'], categories=facet_order, ordered=True)

    # Define first week of September for highlighting (Pacific time)
    sept_start = pd.Timestamp('2025-09-01', tz='US/Pacific')
    sept_week_end = pd.Timestamp('2025-09-07 23:59:59', tz='US/Pacific')

    # Create plot with vertical facets and free Y scales
    p = (ggplot(df_sample, aes(x='timestamp', y='temperature', color='sensor')) +
         geom_rect(aes(xmin=sept_start, xmax=sept_week_end, ymin=-float('inf'), ymax=float('inf')),
                   fill='lightgrey', alpha=0.3, inherit_aes=False) +
         geom_line(alpha=0.6, size=0.5) +
         facet_wrap('~facet_group', scales='free_y', ncol=1) +
         labs(
             title='Temperature Sensor Data Overview',
             subtitle='All raw temperature readings with faceted groupings (independent Y-axis ranges)',
             x='Date',
             y='Temperature (°C)',
             color='Sensor'
         ) +
         theme_nxn() +
         theme(
             figure_size=(8, 6),
             legend_position='right',
             legend_title=element_text(size=12, weight='bold'),
             legend_text=element_text(size=10),
             axis_text_x=element_text(rotation=45, hjust=1),
             plot_title=element_text(size=16, weight='bold'),
             plot_subtitle=element_text(size=12),
             strip_text=element_text(size=12, weight='bold'),
             panel_spacing=0.05
         ) +
         scale_color_manual(values=get_nxn_palette()[:10]) +
         scale_x_datetime(date_breaks='1 month', date_labels='%b %Y') +
         guides(color=guide_legend(title='Sensor', ncol=1))
    )

    p.save('sensor_overview_timeseries.png', dpi=300, bbox_inches='tight')
    return p

def create_september_week_plot(df):
    """Create detailed plot showing first week of September."""
    print("Creating September week detailed plot...")

    # Filter data for first week of September (Pacific time)
    sept_start = pd.Timestamp('2025-09-01', tz='US/Pacific')
    sept_week_end = pd.Timestamp('2025-09-07 23:59:59', tz='US/Pacific')

    df_sept = df[(df['timestamp'] >= sept_start) & (df['timestamp'] <= sept_week_end)].copy()

    if df_sept.empty:
        print("No data found for first week of September")
        return None

    print(f"September week data: {len(df_sept)} readings from {df_sept['sensor'].nunique()} sensors")

    # Create facet grouping for sensors
    def assign_facet_group(sensor):
        if sensor == 'NAS':
            return 'NAS'
        elif sensor in ['Fridge', 'Freezer', 'Wine fridge']:
            return 'Refrigeration'
        else:
            return 'Environmental Sensors'

    df_sept['facet_group'] = df_sept['sensor'].apply(assign_facet_group)

    # Define facet order (top to bottom)
    facet_order = ['Environmental Sensors', 'Refrigeration', 'NAS']
    df_sept['facet_group'] = pd.Categorical(df_sept['facet_group'], categories=facet_order, ordered=True)

    # Create plot with vertical facets and free Y scales
    p = (ggplot(df_sept, aes(x='timestamp', y='temperature', color='sensor')) +
         geom_line(alpha=0.8, size=0.8) +
         facet_wrap('~facet_group', scales='free_y', ncol=1) +
         labs(
             title='Temperature Sensor Data: First Week of September 2025',
             subtitle='High-resolution view of all raw temperature readings (Sep 1-7, 2025)',
             x='Date and Time',
             y='Temperature (°C)',
             color='Sensor'
         ) +
         theme_nxn() +
         theme(
             figure_size=(8, 6),
             legend_position='right',
             legend_title=element_text(size=12, weight='bold'),
             legend_text=element_text(size=10),
             axis_text_x=element_text(rotation=45, hjust=1),
             plot_title=element_text(size=16, weight='bold'),
             plot_subtitle=element_text(size=12),
             strip_text=element_text(size=12, weight='bold'),
             panel_spacing=0.05
         ) +
         scale_color_manual(values=get_nxn_palette()[:10]) +
         scale_x_datetime(date_breaks='1 day', date_labels='%b %d') +
         guides(color=guide_legend(title='Sensor', ncol=1))
    )

    p.save('sensor_overview_september_week.png', dpi=300, bbox_inches='tight')
    return p

def create_correlation_heatmap(df):
    """Create correlation heatmap between sensors."""
    print("Creating sensor correlation heatmap...")

    # Pivot data to get sensors as columns
    pivot_df = df.pivot_table(
        index='timestamp',
        columns='sensor',
        values='temperature',
        aggfunc='mean'
    )

    # Calculate correlation matrix
    corr_matrix = pivot_df.corr()

    # Create heatmap using matplotlib/seaborn for better control
    plt.figure(figsize=(6, 5))

    # Create custom colormap using nxn palette
    nxn_colors = get_nxn_palette()
    cmap = sns.blend_palette([nxn_colors[1], 'white', nxn_colors[0]], as_cmap=True)

    # Create heatmap
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))  # Mask upper triangle
    sns.heatmap(
        corr_matrix,
        mask=mask,
        annot=True,
        fmt='.2f',
        cmap=cmap,
        center=0,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8}
    )

    plt.title('Temperature Sensor Correlation Matrix', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Sensor', fontweight='bold')
    plt.ylabel('Sensor', fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig('sensor_correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_sensor_statistics_summary(df):
    """Create summary statistics visualization."""
    print("Creating sensor statistics summary...")

    # Calculate statistics by sensor
    stats_df = df.groupby('sensor').agg({
        'temperature': ['mean', 'std', 'min', 'max', 'count']
    }).round(2)

    stats_df.columns = ['Mean', 'Std Dev', 'Min', 'Max', 'Count']
    stats_df = stats_df.reset_index()

    # Melt for plotting
    stats_melted = stats_df.melt(
        id_vars=['sensor'],
        value_vars=['Mean', 'Std Dev', 'Min', 'Max'],
        var_name='statistic',
        value_name='value'
    )

    # Create faceted plot
    p = (ggplot(stats_melted, aes(x='sensor', y='value', fill='sensor')) +
         geom_col(alpha=0.8) +
         facet_wrap('~statistic', scales='free_y', ncol=2) +
         labs(
             title='Temperature Sensor Statistics Summary',
             subtitle='Mean, standard deviation, min and max temperatures by sensor',
             x='Sensor',
             y='Temperature (°C)'
         ) +
         theme_nxn() +
         theme(
             figure_size=(8, 5),
             legend_position='none',
             strip_text=element_text(size=11, weight='bold'),
             axis_text_x=element_text(rotation=45, hjust=1)
         ) +
         scale_fill_manual(values=get_nxn_palette()[:10])
    )

    p.save('sensor_statistics_summary.png', dpi=300, bbox_inches='tight')

    # Also save the stats table
    stats_df.to_csv('sensor_statistics_table.csv', index=False)
    print("Saved sensor statistics to sensor_statistics_table.csv")

    return p

def create_cross_correlation_analysis(df):
    """Create cross-correlation analysis between key sensor pairs."""
    print("Creating cross-correlation analysis...")

    # Pivot data
    pivot_df = df.pivot_table(
        index='timestamp',
        columns='sensor',
        values='temperature',
        aggfunc='mean'
    ).fillna(method='ffill').fillna(method='bfill')

    # Select interesting sensor pairs for cross-correlation
    sensor_pairs = [
        ('Outdoors', '3D printer closet'),
        ('Outdoors', 'Living room'),
        ('Outdoors', 'Office'),
        ('3D printer closet', 'Living room'),
        ('Living room', 'Office'),
        ('Fridge', 'Freezer')
    ]

    cross_corr_results = []

    for sensor1, sensor2 in sensor_pairs:
        if sensor1 in pivot_df.columns and sensor2 in pivot_df.columns:
            # Get clean data
            s1 = pivot_df[sensor1].dropna()
            s2 = pivot_df[sensor2].dropna()

            # Align timestamps
            common_idx = s1.index.intersection(s2.index)
            if len(common_idx) > 100:  # Minimum data points
                s1_aligned = s1.loc[common_idx]
                s2_aligned = s2.loc[common_idx]

                # Calculate cross-correlation
                correlation = correlate(s1_aligned, s2_aligned, mode='full')
                lags = np.arange(-len(s2_aligned) + 1, len(s1_aligned))

                # Find peak correlation and lag
                max_corr_idx = np.argmax(np.abs(correlation))
                max_corr = correlation[max_corr_idx]
                lag_at_max = lags[max_corr_idx]

                # Also calculate simple Pearson correlation
                pearson_corr, _ = pearsonr(s1_aligned, s2_aligned)

                cross_corr_results.append({
                    'sensor_pair': f"{sensor1} → {sensor2}",
                    'max_cross_correlation': max_corr,
                    'lag_hours': lag_at_max,
                    'pearson_correlation': pearson_corr,
                    'data_points': len(common_idx)
                })

    # Create results dataframe
    cross_corr_df = pd.DataFrame(cross_corr_results)

    if not cross_corr_df.empty:
        # Create visualization
        p = (ggplot(cross_corr_df, aes(x='sensor_pair', y='pearson_correlation', fill='sensor_pair')) +
             geom_col(alpha=0.8) +
             geom_text(aes(label='pearson_correlation'),
                      format_string='{:.2f}',
                      position=position_nudge(y=0.02)) +
             labs(
                 title='Cross-Correlation Analysis Between Sensor Pairs',
                 subtitle='Pearson correlation coefficients for selected sensor combinations',
                 x='Sensor Pair',
                 y='Pearson Correlation'
             ) +
             theme_nxn() +
             theme(
                 figure_size=(8, 4),
                 legend_position='none',
                 axis_text_x=element_text(rotation=45, hjust=1)
             ) +
             scale_fill_manual(values=get_nxn_palette()[:len(cross_corr_df)])
        )

        p.save('cross_correlation_analysis.png', dpi=300, bbox_inches='tight')

        # Save results
        cross_corr_df.to_csv('cross_correlation_results.csv', index=False)
        print("Saved cross-correlation results to cross_correlation_results.csv")

    return cross_corr_df

def create_sensor_location_diagram():
    """Create a conceptual sensor location diagram."""
    print("Creating sensor location/topology diagram...")

    # Define sensor locations and types using proper names
    sensor_info = {
        'Outdoors': {'type': 'External', 'location': 'Outside', 'x': 0, 'y': 5},
        '3D printer closet': {'type': 'Room', 'location': 'AMS Room', 'x': 2, 'y': 3},
        'Living room': {'type': 'Living', 'location': 'Living Room', 'x': 4, 'y': 3},
        'Office': {'type': 'Office', 'location': 'Desk Area', 'x': 6, 'y': 3},
        'Bedroom': {'type': 'Bedroom', 'location': 'Bedroom', 'x': 8, 'y': 3},
        'IT closet': {'type': 'Tech', 'location': 'IT Room', 'x': 2, 'y': 1},
        'NAS': {'type': 'Storage', 'location': 'Storage', 'x': 4, 'y': 1},
        'Fridge': {'type': 'Kitchen', 'location': 'Refrigerator', 'x': 6, 'y': 1},
        'Freezer': {'type': 'Kitchen', 'location': 'Freezer', 'x': 7, 'y': 1},
        'Wine fridge': {'type': 'Kitchen', 'location': 'Wine Fridge', 'x': 8, 'y': 1}
    }

    # Create DataFrame
    sensor_df = pd.DataFrame.from_dict(sensor_info, orient='index').reset_index()
    sensor_df.columns = ['sensor', 'type', 'location', 'x', 'y']

    # Create plot
    p = (ggplot(sensor_df, aes(x='x', y='y', color='type', size='type')) +
         geom_point(alpha=0.8) +
         geom_text(aes(label='sensor'),
                  position=position_nudge(y=0.2),
                  size=8,
                  fontweight='bold') +
         labs(
             title='Temperature Sensor Network Topology',
             subtitle='Conceptual layout showing sensor locations and types',
             x='Relative X Position',
             y='Relative Y Position',
             color='Sensor Type',
             size='Sensor Type'
         ) +
         theme_nxn() +
         theme(
             figure_size=(8, 4),
             panel_grid_major=element_blank(),
             panel_grid_minor=element_blank(),
             axis_ticks=element_blank()
         ) +
         scale_color_manual(values=get_nxn_palette()[:len(sensor_df['type'].unique())]) +
         scale_size_manual(values=[4] * len(sensor_df['type'].unique())) +
         xlim(-1, 9) +
         ylim(0, 6)
    )

    p.save('sensor_topology_diagram.png', dpi=300, bbox_inches='tight')

    # Save sensor info
    sensor_df.to_csv('sensor_location_info.csv', index=False)

    return p

def main():
    """Main function to create all visualizations."""
    print("=== Temperature Sensor Data Visualization Suite ===")

    # Load data
    df = load_temperature_data()

    # Create visualizations
    print("\n1. Creating overview time series...")
    overview_plot = create_sensor_overview_plot(df)

    print("\n1b. Creating September week detailed plot...")
    september_plot = create_september_week_plot(df)

    print("\n2. Creating correlation heatmap...")
    create_correlation_heatmap(df)

    print("\n3. Creating statistics summary...")
    stats_plot = create_sensor_statistics_summary(df)

    print("\n4. Creating cross-correlation analysis...")
    cross_corr_results = create_cross_correlation_analysis(df)

    print("\n5. Creating sensor topology diagram...")
    topology_plot = create_sensor_location_diagram()

    print("\n=== Visualization Complete ===")
    print("Generated files:")
    print("- sensor_overview_timeseries.png (with September week highlighted)")
    print("- sensor_overview_september_week.png (detailed September 1-7 view)")
    print("- sensor_correlation_heatmap.png")
    print("- sensor_statistics_summary.png")
    print("- cross_correlation_analysis.png")
    print("- sensor_topology_diagram.png")
    print("- sensor_statistics_table.csv")
    print("- cross_correlation_results.csv")
    print("- sensor_location_info.csv")

if __name__ == "__main__":
    main()