"""
Plot training curves from scVI batch size benchmark with sqrt LR scaling.
"""
import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_point,
    labs, scale_x_continuous, scale_color_manual, scale_linetype_manual,
    theme, element_text
)
from theme_nxn import theme_nxn, get_nxn_palette
import matplotlib.pyplot as plt

OUTPUT_DIR = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/260113 - d - SCVI on MPS"

# Load results
results_df = pd.read_csv(f"{OUTPUT_DIR}/benchmark_sqrt_lr_scaling_results.csv")

# Reshape for plotting: train and validation as separate rows
train_data = results_df[['batch_size', 'epoch', 'cells_processed', 'train_loss']].copy()
train_data['loss_type'] = 'Train'
train_data = train_data.rename(columns={'train_loss': 'loss'})

val_data = results_df[['batch_size', 'epoch', 'cells_processed', 'validation_loss']].copy()
val_data['loss_type'] = 'Validation'
val_data = val_data.rename(columns={'validation_loss': 'loss'})

plot_data = pd.concat([train_data, val_data], ignore_index=True)
batch_order = ['128', '256', '512', '1024', '2048', '4096']
plot_data['batch_size'] = pd.Categorical(plot_data['batch_size'].astype(str), categories=batch_order, ordered=True)

# Get palette colors - need 6 for all batch sizes
palette = get_nxn_palette()
# Extend palette if needed
if len(palette) < 6:
    cmap = plt.cm.viridis
    palette = [cmap(i / 5) for i in range(6)]
    palette = ['#%02x%02x%02x' % (int(c[0]*255), int(c[1]*255), int(c[2]*255)) for c in palette]

color_map = {
    '128': palette[0],
    '256': palette[1],
    '512': palette[2],
    '1024': palette[3],
    '2048': palette[4],
    '4096': palette[5],
}

# Create plot
p = (
    ggplot(plot_data, aes(x='cells_processed', y='loss', color='batch_size', linetype='loss_type'))
    + geom_line(size=1)
    + geom_point(size=2.5)
    + scale_color_manual(values=color_map, name='Batch size')
    + scale_linetype_manual(values={'Train': 'solid', 'Validation': 'dashed'}, name='Loss type')
    + scale_x_continuous(
        labels=lambda x: [f'{int(v/1000)}k' for v in x]
    )
    + labs(
        x='Cells processed',
        y='Loss',
        title='SCVI training on MPS (LR = 0.001 × √(batch_size / 128))',
        subtitle='5 epochs, 90,852 cells'
    )
    + theme_nxn()
    + theme(
        figure_size=(10, 6),
        legend_position='right',
        plot_title=element_text(size=18),
        plot_subtitle=element_text(size=14),
        axis_title=element_text(size=14),
        axis_text=element_text(size=12),
        legend_title=element_text(size=13),
        legend_text=element_text(size=12)
    )
)

# Save plot
output_path = f"{OUTPUT_DIR}/training_curves_sqrt_lr_scaling.png"
p.save(output_path, dpi=150)
print(f"Plot saved to: {output_path}")
