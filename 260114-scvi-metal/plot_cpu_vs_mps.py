"""
Bar plot comparing CPU vs MPS training times.
"""
import pandas as pd
from plotnine import (
    ggplot, aes, geom_bar,
    labs, scale_fill_manual, coord_flip, scale_y_continuous,
    theme, element_text
)
from theme_nxn import theme_nxn, get_nxn_palette

OUTPUT_DIR = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/260113 - d - SCVI joint embedding"

# Load results
results_df = pd.read_csv(f"{OUTPUT_DIR}/cpu_vs_mps_results.csv")

# Get total time per accelerator
timing_df = results_df.groupby('accelerator')['total_training_time'].first().reset_index()
timing_df.columns = ['Accelerator', 'Time (s)']
timing_df['Accelerator'] = timing_df['Accelerator'].str.upper()

# Calculate speedup for label
cpu_time = timing_df[timing_df['Accelerator'] == 'CPU']['Time (s)'].values[0]
mps_time = timing_df[timing_df['Accelerator'] == 'MPS']['Time (s)'].values[0]
speedup = cpu_time / mps_time

palette = get_nxn_palette()

def format_minutes_seconds(x):
    """Format seconds as M:SS."""
    result = []
    for val in x:
        minutes = int(val // 60)
        seconds = int(val % 60)
        result.append(f'{minutes}:{seconds:02d}')
    return result

p = (
    ggplot(timing_df, aes(x='Accelerator', y='Time (s)', fill='Accelerator'))
    + geom_bar(stat='identity', width=0.6)
    + scale_fill_manual(values={'CPU': palette[0], 'MPS': palette[1]})
    + scale_y_continuous(labels=format_minutes_seconds, expand=(0, 0, 0.05, 0))
    + coord_flip()
    + labs(
        x='',
        y='Training time (min:sec)',
        title='SCVI training',
        subtitle='5 epochs, 90,852 cells'
    )
    + theme_nxn()
    + theme(
        figure_size=(8, 3),
        legend_position='none'
    )
)

output_path = f"{OUTPUT_DIR}/cpu_vs_mps.png"
p.save(output_path, dpi=150)
print(f"Plot saved to: {output_path}")
print(f"\nCPU: {cpu_time:.1f}s")
print(f"MPS: {mps_time:.1f}s")
print(f"Speedup: {speedup:.2f}x")
