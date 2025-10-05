#!/usr/bin/env python3
"""
Granger Causality Illustration
Visualize how Granger causality works by comparing predictions from AR-only vs. AR+cause models.
"""

import pandas as pd
import numpy as np
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette
from sklearn.linear_model import LinearRegression
from statsmodels.tsa.stattools import grangercausalitytests
import warnings
warnings.filterwarnings('ignore')

# Sensor name mapping
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

def load_september_week_data():
    """Load and prepare September week data."""
    print("Loading September week data...")

    # Load data
    df = pd.read_csv('Data/history.csv')

    # Filter for temperature sensors
    temp_sensors = [
        'sensor.couch_temperature',  # Living room
        'sensor.desk_temperature'     # Office
    ]

    df = df[df['entity_id'].isin(temp_sensors)]

    # Convert timestamp and state
    df['timestamp'] = pd.to_datetime(df['last_changed'])
    df['timestamp'] = df['timestamp'].dt.tz_convert('US/Pacific')
    df['temperature'] = pd.to_numeric(df['state'], errors='coerce')

    # Create sensor names
    df['sensor_short'] = df['entity_id'].str.replace('sensor.', '').str.replace('_temperature', '')
    df['sensor'] = df['sensor_short'].map(SENSOR_NAME_MAPPING)

    # Filter for September week
    sept_start = pd.Timestamp('2025-09-01', tz='US/Pacific')
    sept_week_end = pd.Timestamp('2025-09-07 23:59:59', tz='US/Pacific')
    df = df[(df['timestamp'] >= sept_start) & (df['timestamp'] <= sept_week_end)]

    # Pivot to get sensors as columns
    pivot_df = df.pivot_table(
        index='timestamp',
        columns='sensor',
        values='temperature',
        aggfunc='mean'
    )

    # Resample to hourly
    pivot_df = pivot_df.resample('1H').mean()
    pivot_df = pivot_df.fillna(method='ffill').fillna(method='bfill')

    print(f"Loaded data shape: {pivot_df.shape}")
    print(f"Date range: {pivot_df.index.min()} to {pivot_df.index.max()}")

    return pivot_df

def create_lag_features(data, target_col, predictor_cols, max_lag=3):
    """Create lagged features for time series prediction."""
    df_features = pd.DataFrame(index=data.index)

    # Add target variable
    df_features['target'] = data[target_col]

    # Add lagged features for each predictor
    for col in predictor_cols:
        for lag in range(1, max_lag + 1):
            df_features[f'{col}_lag{lag}'] = data[col].shift(lag)

    # Drop rows with NaN (from lagging)
    df_features = df_features.dropna()

    return df_features

def fit_models(data, target_col='Office', max_lag=3):
    """Fit baseline and Granger models."""
    print(f"\nFitting models with max_lag={max_lag}...")

    # Baseline model: Office ~ Office lags only
    baseline_features = create_lag_features(
        data,
        target_col=target_col,
        predictor_cols=[target_col],
        max_lag=max_lag
    )

    # Granger model: Office ~ Office lags + Living room lags
    granger_features = create_lag_features(
        data,
        target_col=target_col,
        predictor_cols=[target_col, 'Living room'],
        max_lag=max_lag
    )

    # Fit baseline model
    baseline_X = baseline_features[[col for col in baseline_features.columns if col != 'target']]
    baseline_y = baseline_features['target']
    baseline_model = LinearRegression()
    baseline_model.fit(baseline_X, baseline_y)
    baseline_pred = baseline_model.predict(baseline_X)
    baseline_r2 = baseline_model.score(baseline_X, baseline_y)

    # Fit Granger model
    granger_X = granger_features[[col for col in granger_features.columns if col != 'target']]
    granger_y = granger_features['target']
    granger_model = LinearRegression()
    granger_model.fit(granger_X, granger_y)
    granger_pred = granger_model.predict(granger_X)
    granger_r2 = granger_model.score(granger_X, granger_y)

    # Calculate prediction errors
    baseline_error = np.sqrt(np.mean((baseline_y - baseline_pred) ** 2))
    granger_error = np.sqrt(np.mean((granger_y - granger_pred) ** 2))

    print(f"Baseline model R²: {baseline_r2:.4f}, RMSE: {baseline_error:.4f}°C")
    print(f"Granger model R²: {granger_r2:.4f}, RMSE: {granger_error:.4f}°C")
    print(f"Error reduction: {(baseline_error - granger_error):.4f}°C ({(1 - granger_error/baseline_error)*100:.1f}%)")

    # Run formal Granger causality test
    print("\nRunning formal Granger causality test...")
    test_data = data[['Office', 'Living room']].copy()
    # Need to make stationary for Granger test
    test_data_diff = test_data.diff().dropna()
    gc_result = grangercausalitytests(test_data_diff, max_lag, verbose=False)

    # Get F-statistic for optimal lag
    f_stats = []
    for lag in range(1, max_lag + 1):
        f_stat = gc_result[lag][0]['ssr_ftest'][0]
        p_value = gc_result[lag][0]['ssr_ftest'][1]
        f_stats.append((lag, f_stat, p_value))

    best_lag = min(f_stats, key=lambda x: x[2])
    print(f"Best lag: {best_lag[0]}, F-statistic: {best_lag[1]:.2f}, p-value: {best_lag[2]:.4e}")

    # Create results dataframe
    results = {
        'baseline': {
            'features': baseline_features,
            'model': baseline_model,
            'predictions': baseline_pred,
            'r2': baseline_r2,
            'rmse': baseline_error
        },
        'granger': {
            'features': granger_features,
            'model': granger_model,
            'predictions': granger_pred,
            'r2': granger_r2,
            'rmse': granger_error
        },
        'granger_test': {
            'f_statistic': best_lag[1],
            'p_value': best_lag[2],
            'optimal_lag': best_lag[0]
        }
    }

    return results

def create_visualization(data, results, max_lag=3):
    """Create two-facet visualization showing predictions and lag connections."""
    print("\nCreating visualization...")

    # Filter data to just two days (Sept 2-3)
    sept_2_start = pd.Timestamp('2025-09-02 00:00:00', tz='US/Pacific')
    sept_3_end = pd.Timestamp('2025-09-03 23:59:59', tz='US/Pacific')
    data = data[(data.index >= sept_2_start) & (data.index <= sept_3_end)]

    # Choose a specific time point to illustrate (Sept 3, 2pm)
    example_time = pd.Timestamp('2025-09-03 14:00:00', tz='US/Pacific')

    # Prepare data for plotting - filter predictions to 2-day window too
    plot_data = []

    # Add actual Office temperature
    for idx in data.index:
        plot_data.append({
            'timestamp': idx,
            'temperature': data.loc[idx, 'Office'],
            'series': 'Actual Office',
            'model': 'Baseline'
        })
        plot_data.append({
            'timestamp': idx,
            'temperature': data.loc[idx, 'Office'],
            'series': 'Actual Office',
            'model': 'Granger (with Living room)'
        })

    # Add Living room temperature (for reference)
    for idx in data.index:
        plot_data.append({
            'timestamp': idx,
            'temperature': data.loc[idx, 'Living room'],
            'series': 'Living room (reference)',
            'model': 'Baseline'
        })
        plot_data.append({
            'timestamp': idx,
            'temperature': data.loc[idx, 'Living room'],
            'series': 'Living room (reference)',
            'model': 'Granger (with Living room)'
        })

    # Add baseline predictions - filter to 2-day window
    baseline_features = results['baseline']['features']
    baseline_pred = results['baseline']['predictions']
    for idx, pred in zip(baseline_features.index, baseline_pred):
        if sept_2_start <= idx <= sept_3_end:
            plot_data.append({
                'timestamp': idx,
                'temperature': pred,
                'series': 'Predicted Office',
                'model': 'Baseline'
            })

    # Add Granger predictions - filter to 2-day window
    granger_features = results['granger']['features']
    granger_pred = results['granger']['predictions']
    for idx, pred in zip(granger_features.index, granger_pred):
        if sept_2_start <= idx <= sept_3_end:
            plot_data.append({
                'timestamp': idx,
                'temperature': pred,
                'series': 'Predicted Office',
                'model': 'Granger (with Living room)'
            })

    plot_df = pd.DataFrame(plot_data)

    # Create lag connection data for the example time point with right-angle paths
    lag_connections = []

    if example_time in baseline_features.index:
        # Baseline connections (Office lags only)
        baseline_pred_value = baseline_pred[baseline_features.index.get_loc(example_time)]
        for lag in range(1, max_lag + 1):
            lag_time = example_time - pd.Timedelta(hours=lag)
            if lag_time in data.index:
                y_start = data.loc[lag_time, 'Office']
                # Create path: up -> right -> down
                # Determine intermediate y-level (slightly above max to avoid overlap)
                y_intermediate = data['Office'].max() + (lag * 0.3)

                # Vertical segment up
                lag_connections.append({
                    'x': lag_time,
                    'xend': lag_time,
                    'y': y_start,
                    'yend': y_intermediate,
                    'model': 'Baseline',
                    'lag_type': f'Office(t-{lag})',
                    'segment': 'vertical_up'
                })
                # Horizontal segment
                lag_connections.append({
                    'x': lag_time,
                    'xend': example_time,
                    'y': y_intermediate,
                    'yend': y_intermediate,
                    'model': 'Baseline',
                    'lag_type': f'Office(t-{lag})',
                    'segment': 'horizontal'
                })
                # Vertical segment down
                lag_connections.append({
                    'x': example_time,
                    'xend': example_time,
                    'y': y_intermediate,
                    'yend': baseline_pred_value,
                    'model': 'Baseline',
                    'lag_type': f'Office(t-{lag})',
                    'segment': 'vertical_down'
                })

        # Granger connections (Office + Living room lags)
        granger_pred_value = granger_pred[granger_features.index.get_loc(example_time)]
        for lag in range(1, max_lag + 1):
            lag_time = example_time - pd.Timedelta(hours=lag)
            if lag_time in data.index:
                # Office lag path
                y_start_office = data.loc[lag_time, 'Office']
                y_intermediate_office = data['Office'].max() + (lag * 0.3)

                lag_connections.append({
                    'x': lag_time,
                    'xend': lag_time,
                    'y': y_start_office,
                    'yend': y_intermediate_office,
                    'model': 'Granger (with Living room)',
                    'lag_type': f'Office(t-{lag})',
                    'segment': 'vertical_up'
                })
                lag_connections.append({
                    'x': lag_time,
                    'xend': example_time,
                    'y': y_intermediate_office,
                    'yend': y_intermediate_office,
                    'model': 'Granger (with Living room)',
                    'lag_type': f'Office(t-{lag})',
                    'segment': 'horizontal'
                })
                lag_connections.append({
                    'x': example_time,
                    'xend': example_time,
                    'y': y_intermediate_office,
                    'yend': granger_pred_value,
                    'model': 'Granger (with Living room)',
                    'lag_type': f'Office(t-{lag})',
                    'segment': 'vertical_down'
                })

                # Living room lag path
                y_start_living = data.loc[lag_time, 'Living room']
                y_intermediate_living = data['Living room'].max() + (lag * 0.3)

                lag_connections.append({
                    'x': lag_time,
                    'xend': lag_time,
                    'y': y_start_living,
                    'yend': y_intermediate_living,
                    'model': 'Granger (with Living room)',
                    'lag_type': f'Living room(t-{lag})',
                    'segment': 'vertical_up'
                })
                lag_connections.append({
                    'x': lag_time,
                    'xend': example_time,
                    'y': y_intermediate_living,
                    'yend': y_intermediate_living,
                    'model': 'Granger (with Living room)',
                    'lag_type': f'Living room(t-{lag})',
                    'segment': 'horizontal'
                })
                lag_connections.append({
                    'x': example_time,
                    'xend': example_time,
                    'y': y_intermediate_living,
                    'yend': granger_pred_value,
                    'model': 'Granger (with Living room)',
                    'lag_type': f'Living room(t-{lag})',
                    'segment': 'vertical_down'
                })

    lag_df = pd.DataFrame(lag_connections)

    # Create color palette matching sensor_overview_september_week.png
    # Colors are assigned to sensors in alphabetical order (plotnine default behavior)
    # Palette has 7 colors, so it cycles for 10 sensors
    palette = get_nxn_palette()

    # All sensor names in alphabetical order
    all_sensors_alphabetical = [
        '3D printer closet', 'Bedroom', 'Freezer', 'Fridge',
        'IT closet', 'Living room', 'NAS', 'Office', 'Outdoors', 'Wine fridge'
    ]

    # Create color mapping: sensor name -> palette color by alphabetical position (with cycling)
    sensor_color_map = {sensor: palette[i % len(palette)] for i, sensor in enumerate(all_sensors_alphabetical)}

    # Simplify series names and separate color vs linetype
    plot_df['sensor'] = plot_df['series'].apply(lambda x: 'Office' if 'Office' in x else 'Living room')
    plot_df['kind'] = plot_df['series'].apply(lambda x: 'Prediction' if 'Predicted' in x else 'Observed')

    # Extract colors for our two sensors
    colors = {
        'Office': sensor_color_map['Office'],
        'Living room': sensor_color_map['Living room']
    }

    p = (ggplot(plot_df, aes(x='timestamp', y='temperature', color='sensor', linetype='kind')) +
         geom_line(size=1, alpha=0.8) +
         facet_wrap('~model', ncol=1, scales='free_y') +
         scale_color_manual(values=colors, name='Sensor') +
         scale_linetype_manual(values={'Observed': 'solid', 'Prediction': 'dashed'}, name='Kind') +
         labs(
             title='Granger Causality: Comparing Prediction Models',
             subtitle=' ',
             x='Date and Time',
             y='Temperature (°C)',
             color='Series'
         ) +
         scale_x_datetime(date_breaks='6 hours', date_labels='%b %d\n%H:%M') +
         theme_nxn() +
         theme(
             figure_size=(8, 6),
             legend_position='right',
             strip_text=element_text(size=12, weight='bold'),
             axis_text_x=element_text(rotation=45, hjust=1)
         )
    )

    # Add vertical line at example time point with annotate (vline doesn't work well with datetime)
    if example_time in baseline_features.index:
        # Get y range for the vertical line
        y_min = plot_df['temperature'].min()
        y_max = plot_df['temperature'].max()

        # Create vertical line data
        vline_data = pd.DataFrame({
            'x': [example_time, example_time],
            'y': [y_min, y_max],
            'model': ['Baseline', 'Granger (with Living room)']
        })

        p = p + geom_line(
            data=vline_data,
            mapping=aes(x='x', y='y'),
            linetype='dashed',
            color='darkgray',
            size=1.5,
            alpha=0.6,
            inherit_aes=False
        )

    # Add lag connections if they exist
    if not lag_df.empty:
        # Separate Office and Living room lags for better coloring
        lag_df['source'] = lag_df['lag_type'].apply(
            lambda x: 'Office lags' if 'Office' in x else 'Living room lags'
        )

        # Create different colors for lag sources (but don't add to legend)
        # Match arrow colors to their corresponding sensor colors
        lag_colors = {
            'Office lags': sensor_color_map['Office'],
            'Living room lags': sensor_color_map['Living room']
        }

        # Add arrow only on the final downward segment
        lag_df_no_arrow = lag_df[lag_df['segment'] != 'vertical_down']
        lag_df_arrow = lag_df[lag_df['segment'] == 'vertical_down']

        # Separate Office and Living room arrows to control plotting order
        # Plot Office arrows first, then Living room on top
        office_no_arrow = lag_df_no_arrow[lag_df_no_arrow['source'] == 'Office lags']
        living_no_arrow = lag_df_no_arrow[lag_df_no_arrow['source'] == 'Living room lags']
        office_arrow = lag_df_arrow[lag_df_arrow['source'] == 'Office lags']
        living_arrow = lag_df_arrow[lag_df_arrow['source'] == 'Living room lags']

        # Office segments without arrows (first layer)
        if not office_no_arrow.empty:
            p = p + geom_segment(
                data=office_no_arrow,
                mapping=aes(x='x', y='y', xend='xend', yend='yend'),
                color=lag_colors['Office lags'],
                size=1.0,
                alpha=1.0,
                inherit_aes=False,
                show_legend=False
            )

        # Office final segment with arrow
        if not office_arrow.empty:
            p = p + geom_segment(
                data=office_arrow,
                mapping=aes(x='x', y='y', xend='xend', yend='yend'),
                color=lag_colors['Office lags'],
                arrow=arrow(length=0.12, type='closed'),
                size=1.0,
                alpha=1.0,
                inherit_aes=False,
                show_legend=False
            )

        # Living room segments without arrows (on top of Office)
        if not living_no_arrow.empty:
            p = p + geom_segment(
                data=living_no_arrow,
                mapping=aes(x='x', y='y', xend='xend', yend='yend'),
                color=lag_colors['Living room lags'],
                size=1.0,
                alpha=1.0,
                inherit_aes=False,
                show_legend=False
            )

        # Living room final segment with arrow (topmost layer)
        if not living_arrow.empty:
            p = p + geom_segment(
                data=living_arrow,
                mapping=aes(x='x', y='y', xend='xend', yend='yend'),
                color=lag_colors['Living room lags'],
                arrow=arrow(length=0.12, type='closed'),
                size=1.0,
                alpha=1.0,
                inherit_aes=False,
                show_legend=False
            )

    # Statistics will be saved to markdown file instead of plotting
    baseline_rmse = results['baseline']['rmse']
    granger_rmse = results['granger']['rmse']
    f_stat = results['granger_test']['f_statistic']
    p_val = results['granger_test']['p_value']

    # Save statistics to markdown file
    stats_text = f"""# Granger Causality Analysis Results

## Living room → Office Relationship

**Time period**: September 2-3, 2025
**Prediction target**: Office temperature at Sept 3, 2pm
**Maximum lag**: {max_lag} hours

## Model Performance

### Baseline Model (Autoregressive only)
- **Predictors**: Office(t-1, t-2, t-3)
- **RMSE**: {baseline_rmse:.4f}°C
- **R²**: {results['baseline']['r2']:.4f}

### Granger Model (with Living room)
- **Predictors**: Office(t-1, t-2, t-3) + Living room(t-1, t-2, t-3)
- **RMSE**: {granger_rmse:.4f}°C
- **R²**: {results['granger']['r2']:.4f}

### Improvement
- **Error reduction**: {(baseline_rmse - granger_rmse):.4f}°C ({(1 - granger_rmse/baseline_rmse)*100:.1f}%)
- **F-statistic**: {f_stat:.2f}
- **p-value**: {p_val:.4e}

## Interpretation

The Granger causality test shows that adding Living room lagged values (t-1, t-2, t-3)
significantly improves Office temperature prediction beyond what Office's own history provides.
The F-statistic of {f_stat:.2f} with p < 0.0001 indicates this improvement is highly statistically
significant.

However, as documented in the project analysis, this statistical relationship should not be
interpreted as true causality due to:
- Shared external drivers (weather, HVAC)
- Common building thermal mass effects
- High network density in observational data
"""

    with open('granger_causality_illustration_stats.md', 'w') as f:
        f.write(stats_text)

    # Save plot
    p.save('granger_causality_illustration.png', dpi=300, bbox_inches='tight')
    print("Saved visualization to granger_causality_illustration.png")

    return p

def main():
    """Main execution function."""
    print("=== Granger Causality Illustration ===")

    # Load data
    data = load_september_week_data()

    # Fit models
    max_lag = 3
    results = fit_models(data, target_col='Office', max_lag=max_lag)

    # Create visualization
    create_visualization(data, results, max_lag=max_lag)

    print("\n=== Complete ===")

if __name__ == "__main__":
    main()
