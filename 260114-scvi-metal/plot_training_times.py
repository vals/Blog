import pandas as pd
from plotnine import ggplot, aes, geom_bar, labs, scale_y_continuous, theme, element_text
from theme_nxn import theme_nxn, get_nxn_palette

def format_time(seconds):
    """Format seconds as 'Xm Ys'"""
    m, s = divmod(int(seconds), 60)
    if m > 0:
        return f'{m}m {s}s'
    return f'{s}s'

# Load data and get unique batch size / training time pairs
df = pd.read_csv('benchmark_lr_scaling_results.csv')
df_unique = df.groupby('batch_size')['total_training_time'].first().reset_index()

# Create bar plot
p = (
    ggplot(df_unique, aes(x='factor(batch_size)', y='total_training_time'))
    + geom_bar(stat='identity', fill=get_nxn_palette()[0])
    + scale_y_continuous(expand=(0, 0, 0.05, 0), labels=lambda x: [format_time(v) for v in x])
    + labs(
        x='Batch Size',
        y='Training Time',
        title='SCVI training on MPS',
        subtitle='5 epochs, 90,852 cells'
    )
    + theme_nxn()
    + theme(
        figure_size=(8, 5),
        plot_title=element_text(size=14),
        plot_subtitle=element_text(size=11),
        axis_title=element_text(size=12),
        axis_text=element_text(size=11)
    )
)

p.save('mps_training_times.png', dpi=150)
print('Saved to mps_training_times.png')
