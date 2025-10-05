#!/usr/bin/env python3
"""
Clean Temperature Sensor Granger Causality Analysis
Streamlined analysis focused only on temperature sensor data with proper preprocessing.
"""

import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests, adfuller
from statsmodels.tsa.vector_ar.var_model import VAR
import networkx as nx
import matplotlib.pyplot as plt
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette
import seaborn as sns
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
    """Load and preprocess temperature sensor data."""
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
    df['temperature'] = pd.to_numeric(df['state'], errors='coerce')

    # Create sensor short names
    df['sensor'] = df['entity_id'].str.replace('sensor.', '').str.replace('_temperature', '').str.replace('_sensor', '')

    # Remove outliers (beyond 3 standard deviations)
    df = df.groupby('sensor').apply(
        lambda x: x[np.abs(x['temperature'] - x['temperature'].mean()) <= 3 * x['temperature'].std()]
    ).reset_index(drop=True)

    print(f"Loaded {len(df)} temperature readings from {df['sensor'].nunique()} sensors")
    print(f"Date range: {df['timestamp'].min()} to {df['timestamp'].max()}")

    return df

def prepare_time_series_data(df, resample_freq='1H'):
    """Prepare time series data for Granger causality analysis."""
    print(f"Preparing time series data with {resample_freq} frequency...")

    # Pivot to get sensors as columns
    pivot_df = df.pivot_table(
        index='timestamp',
        columns='sensor',
        values='temperature',
        aggfunc='mean'
    )

    # Resample to regular intervals
    pivot_df = pivot_df.resample(resample_freq).mean()

    # Forward fill then backward fill missing values
    pivot_df = pivot_df.fillna(method='ffill').fillna(method='bfill')

    # Remove any remaining NaN values
    pivot_df = pivot_df.dropna()

    print(f"Prepared data shape: {pivot_df.shape}")
    print(f"Date range: {pivot_df.index.min()} to {pivot_df.index.max()}")
    print(f"Sensors: {list(pivot_df.columns)}")

    return pivot_df

def test_stationarity(data, sensor_name):
    """Test stationarity using Augmented Dickey-Fuller test."""
    adf_result = adfuller(data.dropna())
    return {
        'sensor': sensor_name,
        'adf_statistic': adf_result[0],
        'p_value': adf_result[1],
        'critical_values': adf_result[4],
        'is_stationary': adf_result[1] < 0.05
    }

def make_stationary(data):
    """Make time series stationary using first differences."""
    return data.diff().dropna()

def analyze_stationarity(pivot_df):
    """Analyze stationarity of all sensors."""
    print("Testing stationarity of temperature time series...")

    stationarity_results = []
    for sensor in pivot_df.columns:
        result = test_stationarity(pivot_df[sensor], sensor)
        stationarity_results.append(result)

    stationarity_df = pd.DataFrame(stationarity_results)

    print("\nStationarity Test Results:")
    print(f"Stationary sensors: {stationarity_df['is_stationary'].sum()}/{len(stationarity_df)}")
    print(stationarity_df[['sensor', 'adf_statistic', 'p_value', 'is_stationary']])

    return stationarity_df

def perform_granger_causality_tests(data, max_lag=5):
    """Perform pairwise Granger causality tests."""
    print(f"Performing Granger causality tests with max lag {max_lag}...")

    sensors = list(data.columns)
    n_sensors = len(sensors)
    results = []

    for i, cause_sensor in enumerate(sensors):
        for j, effect_sensor in enumerate(sensors):
            if i != j:  # Don't test sensor against itself
                print(f"Testing {cause_sensor} → {effect_sensor}")

                # Create bivariate time series
                test_data = data[[effect_sensor, cause_sensor]].dropna()

                if len(test_data) > max_lag * 2:  # Ensure enough data
                    try:
                        # Perform Granger causality test
                        gc_result = grangercausalitytests(test_data, max_lag, verbose=False)

                        # Extract results for each lag
                        for lag in range(1, max_lag + 1):
                            if lag in gc_result:
                                lag_results = gc_result[lag][0]

                                # Get F-statistic and p-value from ssr_ftest
                                f_stat = lag_results['ssr_ftest'][0]
                                p_value = lag_results['ssr_ftest'][1]

                                results.append({
                                    'cause_sensor': cause_sensor,
                                    'effect_sensor': effect_sensor,
                                    'lag': lag,
                                    'f_statistic': f_stat,
                                    'p_value': p_value,
                                    'significant': p_value < 0.05,
                                    'data_points': len(test_data)
                                })
                    except Exception as e:
                        print(f"Error testing {cause_sensor} → {effect_sensor}: {e}")
                        continue

    results_df = pd.DataFrame(results)

    if not results_df.empty:
        # Apply professional sensor names
        results_df['cause_sensor'] = results_df['cause_sensor'].map(SENSOR_NAME_MAPPING)
        results_df['effect_sensor'] = results_df['effect_sensor'].map(SENSOR_NAME_MAPPING)

        print(f"\nGranger Causality Results:")
        print(f"Total tests: {len(results_df)}")
        print(f"Significant relationships: {results_df['significant'].sum()}")
        print(f"Significance rate: {results_df['significant'].mean():.1%}")

    return results_df

def select_optimal_lag_results(results_df):
    """Select optimal lag for each sensor pair based on lowest p-value."""
    if results_df.empty:
        return results_df

    # Group by sensor pair and select the lag with minimum p-value
    optimal_results = results_df.loc[
        results_df.groupby(['cause_sensor', 'effect_sensor'])['p_value'].idxmin()
    ].copy()

    print(f"\nOptimal lag results: {len(optimal_results)} sensor pairs")
    print(f"Significant relationships: {optimal_results['significant'].sum()}")

    return optimal_results

def create_causality_network_plot(results_df, title_suffix=""):
    """Create network visualization of causal relationships with consistent color scheme."""
    print("Creating causality network visualization...")

    if results_df.empty:
        print("No results to plot")
        return None

    # Filter for significant relationships with F > 400
    sig_results = results_df[(results_df['significant']) & (results_df['f_statistic'] > 400)].copy()

    if sig_results.empty:
        print("No significant relationships found with F > 400")
        return None

    print(f"Displaying {len(sig_results)} relationships with F > 400 (filtered from {results_df['significant'].sum()} total significant)")

    # Create network graph
    G = nx.DiGraph()

    # Add nodes (all sensors)
    all_sensors = sorted(set(results_df['cause_sensor'].unique()) | set(results_df['effect_sensor'].unique()))
    G.add_nodes_from(all_sensors)

    # Add edges for significant relationships
    for _, row in sig_results.iterrows():
        G.add_edge(
            row['cause_sensor'],
            row['effect_sensor'],
            weight=row['f_statistic'],
            p_value=row['p_value'],
            lag=row['lag']
        )

    # Create color mapping matching time series plots (alphabetical order with palette cycling)
    palette = get_nxn_palette()
    sensor_color_map = {sensor: palette[i % len(palette)] for i, sensor in enumerate(all_sensors)}
    node_colors = [sensor_color_map[sensor] for sensor in G.nodes()]

    # Create matplotlib plot with web-optimized sizing
    fig, ax = plt.subplots(figsize=(10, 8))

    # Import numpy for positioning
    import numpy as np

    # Position nodes using hierarchical layout based on network role
    # Calculate in-degree and out-degree for each node
    in_degree = dict(G.in_degree())
    out_degree = dict(G.out_degree())

    # Assign layers based on network role (source vs sink)
    # Set as node attributes for multipartite_layout
    for node in G.nodes():
        # Classify: more outputs = left (layer 0), more inputs = right (layer 2)
        if out_degree[node] > in_degree[node]:
            G.nodes[node]['layer'] = 0  # Source layer
        elif in_degree[node] > out_degree[node]:
            G.nodes[node]['layer'] = 2  # Sink layer
        else:
            G.nodes[node]['layer'] = 1  # Intermediate layer

    # Use multipartite layout with node attributes
    pos = nx.multipartite_layout(G, subset_key='layer', scale=3, align='horizontal')

    # Draw edges FIRST (before nodes) with manual arrow positioning
    from matplotlib.patches import FancyArrowPatch

    edges = G.edges(data=True)
    edge_weights = [edge[2]['weight'] for edge in edges]

    if edge_weights:
        # Normalize edge weights for visualization with larger minimum
        min_weight = min(edge_weights)
        max_weight = max(edge_weights)
        normalized_weights = {(u, v): (w - min_weight) / (max_weight - min_weight + 1e-6) * 4 + 1.5
                             for (u, v, d), w in zip(edges, edge_weights)}

        # Draw each edge manually with proper margins
        # Calculate node radius in data coordinates from node_size
        # node_size is in points^2, need to convert to data coordinates
        # Get axis limits to calculate scale
        ax_width = 6  # Approximate width in data coordinates (from multipartite scale=3)
        fig_width_inches = 10  # Updated figure width
        points_per_data_unit = (fig_width_inches * 72) / ax_width  # 72 points per inch
        node_size_points = 2000  # Our node_size parameter
        node_radius_points = np.sqrt(node_size_points / np.pi)  # Radius from area
        node_radius = node_radius_points / points_per_data_unit

        for (u, v, d) in edges:
            # Get positions
            x1, y1 = pos[u]
            x2, y2 = pos[v]

            # Calculate angle and shorten edge to avoid node overlap
            dx = x2 - x1
            dy = y2 - y1
            dist = np.sqrt(dx**2 + dy**2)

            # Shorten both ends
            margin = node_radius * 1.1  # 10% larger than node radius
            x1_new = x1 + (dx / dist) * margin
            y1_new = y1 + (dy / dist) * margin
            x2_new = x2 - (dx / dist) * margin
            y2_new = y2 - (dy / dist) * margin

            # Create curved arrow with triangular arrowhead
            arrow = FancyArrowPatch(
                (x1_new, y1_new), (x2_new, y2_new),
                arrowstyle='-|>',  # Triangular arrowhead
                connectionstyle='arc3,rad=0.1',
                mutation_scale=25,
                linewidth=normalized_weights[(u, v)],
                alpha=1.0,
                color='#666666',
                zorder=1
            )
            ax.add_patch(arrow)

    # Draw nodes AFTER edges (on top) with consistent colors
    nodes = nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=2000,
        alpha=0.8,
        ax=ax
    )
    nodes.set_zorder(2)

    # Draw labels with larger font and white background for readability
    labels = nx.draw_networkx_labels(
        G, pos,
        font_size=11,
        font_weight='bold',
        ax=ax,
        bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=2)
    )

    # Add title
    ax.set_title(f'Temperature Sensor Granger Causality Network{title_suffix}',
                 fontsize=16, fontweight='bold', pad=20)
    ax.axis('off')

    plt.tight_layout()

    # Save plot
    filename = f'temperature_granger_network{title_suffix.lower().replace(" ", "_")}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Network statistics
    print(f"\nNetwork Statistics:")
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")
    print(f"Network density: {nx.density(G):.1%}")

    return G

def create_results_heatmap(results_df, title_suffix=""):
    """Create heatmap of Granger causality results."""
    print("Creating results heatmap...")

    if results_df.empty:
        print("No results to plot")
        return None

    # Create pivot table for F-statistics
    f_stat_matrix = results_df.pivot_table(
        index='cause_sensor',
        columns='effect_sensor',
        values='f_statistic',
        fill_value=0
    )

    # Create significance mask
    sig_matrix = results_df.pivot_table(
        index='cause_sensor',
        columns='effect_sensor',
        values='significant',
        fill_value=False
    )

    # Create plot
    plt.figure(figsize=(12, 10))

    # Create colormap
    nxn_colors = get_nxn_palette()
    cmap = sns.blend_palette([nxn_colors[1], 'white', nxn_colors[0]], as_cmap=True)

    # Use log scale for better visualization of large F-statistic range
    # But keep original values for annotations
    log_f_stat_matrix = np.log10(f_stat_matrix + 1)  # +1 to handle zeros

    # Create heatmap with log scale
    sns.heatmap(
        log_f_stat_matrix,
        annot=f_stat_matrix.applymap(lambda x: f'{x:.0f}' if x < 100 else f'{x:.0e}'),
        fmt='',
        cmap=cmap,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8, "label": "log₁₀(F-Statistic + 1)"}
    )

    # Add significance markers
    for i, cause in enumerate(f_stat_matrix.index):
        for j, effect in enumerate(f_stat_matrix.columns):
            if sig_matrix.loc[cause, effect]:
                plt.text(j + 0.5, i + 0.8, '*',
                        fontsize=20, fontweight='bold',
                        ha='center', va='center', color='red')

    plt.title(f'Granger Causality F-Statistics Heatmap{title_suffix}',
              fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Effect Sensor (Y)', fontweight='bold')
    plt.ylabel('Cause Sensor (X)', fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Add note about significance
    plt.figtext(0.02, 0.02, '* indicates p < 0.05', fontsize=10, style='italic')

    plt.tight_layout()

    # Save plot
    filename = f'temperature_granger_heatmap{title_suffix.lower().replace(" ", "_")}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def generate_summary_report(stationarity_df, granger_results, filename="temperature_granger_report.txt"):
    """Generate comprehensive summary report."""
    print("Generating summary report...")

    with open(filename, 'w') as f:
        f.write("Temperature Sensor Granger Causality Analysis Report\n")
        f.write("=" * 55 + "\n\n")

        f.write("Data Overview:\n")
        f.write(f"- Number of sensors: {len(stationarity_df)}\n")
        f.write(f"- Sensors: {', '.join(stationarity_df['sensor'].tolist())}\n\n")

        f.write("Stationarity Analysis:\n")
        f.write(f"- Stationary sensors: {stationarity_df['is_stationary'].sum()}/{len(stationarity_df)}\n")
        non_stationary = stationarity_df[~stationarity_df['is_stationary']]['sensor'].tolist()
        if non_stationary:
            f.write(f"- Non-stationary sensors: {', '.join(non_stationary)}\n")
        f.write("\n")

        if not granger_results.empty:
            f.write("Granger Causality Results:\n")
            f.write(f"- Total sensor pairs tested: {len(granger_results)}\n")
            f.write(f"- Significant relationships: {granger_results['significant'].sum()}\n")
            f.write(f"- Significance rate: {granger_results['significant'].mean():.1%}\n\n")

            # Top relationships by F-statistic
            top_results = granger_results.nlargest(10, 'f_statistic')
            f.write("Top 10 Relationships by F-Statistic:\n")
            f.write("-" * 45 + "\n")
            for _, row in top_results.iterrows():
                sig_marker = "*" if row['significant'] else ""
                f.write(f"{row['cause_sensor']} → {row['effect_sensor']}: "
                       f"F={row['f_statistic']:.2f}, p={row['p_value']:.4f}, "
                       f"lag={row['lag']}{sig_marker}\n")
            f.write("\n")

            # Significant relationships
            sig_results = granger_results[granger_results['significant']]
            if not sig_results.empty:
                f.write("Significant Causal Relationships (p < 0.05):\n")
                f.write("-" * 40 + "\n")
                for _, row in sig_results.iterrows():
                    f.write(f"{row['cause_sensor']} → {row['effect_sensor']}: "
                           f"F={row['f_statistic']:.2f}, p={row['p_value']:.4f}, "
                           f"lag={row['lag']} hours\n")
        else:
            f.write("No Granger causality results available.\n")

def main():
    """Main analysis function."""
    print("=== Temperature Sensor Granger Causality Analysis ===")

    # Load and prepare data
    df = load_temperature_data()
    pivot_df = prepare_time_series_data(df, resample_freq='1H')

    # Analyze stationarity
    stationarity_df = analyze_stationarity(pivot_df)

    # Apply first differences to make series stationary
    print("\nApplying first differences to ensure stationarity...")
    diff_df = pivot_df.diff().dropna()

    # Verify stationarity after differencing
    print("Testing stationarity after differencing...")
    diff_stationarity = []
    for sensor in diff_df.columns:
        result = test_stationarity(diff_df[sensor], sensor)
        diff_stationarity.append(result)
    diff_stationarity_df = pd.DataFrame(diff_stationarity)
    print(f"Stationary after differencing: {diff_stationarity_df['is_stationary'].sum()}/{len(diff_stationarity_df)}")

    # Perform Granger causality analysis on differenced data
    granger_results = perform_granger_causality_tests(diff_df, max_lag=5)

    if not granger_results.empty:
        # Select optimal lag results
        optimal_results = select_optimal_lag_results(granger_results)

        # Create visualizations
        network = create_causality_network_plot(optimal_results, " (Stationary Data)")
        create_results_heatmap(optimal_results, " (Stationary Data)")

        # Save results
        optimal_results.to_csv('temperature_granger_results.csv', index=False)
        print("Saved results to temperature_granger_results.csv")

    # Generate report
    final_results = optimal_results if not granger_results.empty else pd.DataFrame()
    generate_summary_report(diff_stationarity_df, final_results)

    print("\n=== Analysis Complete ===")
    print("Generated files:")
    print("- temperature_granger_network_stationary_data.png")
    print("- temperature_granger_heatmap_stationary_data.png")
    print("- temperature_granger_results.csv")
    print("- temperature_granger_report.txt")

if __name__ == "__main__":
    main()