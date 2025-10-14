#!/bin/bash
# Create animation from volcano plots

# Check if ffmpeg is available
if ! command -v ffmpeg &> /dev/null; then
    echo "Error: ffmpeg not found. Install with: brew install ffmpeg"
    exit 1
fi

# Check if plots directory exists
if [ ! -d "plots" ] || [ ! "$(ls -A plots/volcano_epoch_*.png 2>/dev/null)" ]; then
    echo "Error: No plots found in 'plots/' directory"
    echo ""
    echo "Please run the analysis first:"
    echo "  1. python train_and_save_de.py    # Train model and save DE results"
    echo "  2. python create_plots.py         # Generate plots with fixed axes"
    echo "  3. ./create_animation.sh          # Create animation"
    exit 1
fi

# Count frames
n_frames=$(ls plots/volcano_epoch_*.png 2>/dev/null | wc -l)
echo "Creating animation from $n_frames volcano plots..."

# Create concat file with variable durations
concat_file="plots/concat_list.txt"
echo "Generating concat file with variable frame durations..."

# Clear/create concat file
> "$concat_file"

# Get sorted list of frames (just filenames)
frames=($(ls plots/volcano_epoch_*.png | sort -V | xargs -n 1 basename))

# Add frames with durations
last_idx=$((${#frames[@]} - 1))
for i in "${!frames[@]}"; do
    # First 3 frames: 1.5 seconds each
    # Next 5 frames: 0.5 seconds each
    # Remaining frames: 0.04 seconds (25 fps)
    if [ $i -lt 3 ]; then
        duration=1.0
    elif [ $i -lt 8 ]; then
        duration=0.5
    else
        duration=0.04
    fi

    echo "file '${frames[$i]}'" >> "$concat_file"
    echo "duration $duration" >> "$concat_file"
done

# Add last frame again without duration (required by ffmpeg)
echo "file '${frames[$last_idx]}'" >> "$concat_file"

# Create animation using concat demuxer (run from plots directory)
cd plots
ffmpeg -y -f concat -safe 0 -i "concat_list.txt" \
    -c:v libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
    de_evolution.mp4

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ MP4 animation created successfully!"
    echo "  Output: plots/de_evolution.mp4"
    echo "  Frames: $n_frames"
    echo "  First 3 frames: 1.5s each"
    echo "  Next 5 frames: 0.5s each"
    echo "  Remaining frames: 0.04s each (25 fps)"

    # Generate GIF with optimized palette (800px width)
    echo ""
    echo "Creating GIF version (800px wide)..."
    ffmpeg -y -i de_evolution.mp4 -vf "fps=25,scale=800:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=256[p];[s1][p]paletteuse=dither=bayer:bayer_scale=5" de_evolution.gif

    if [ $? -eq 0 ]; then
        echo "✓ GIF created successfully!"
        echo "  Output: plots/de_evolution.gif"

        # Show file sizes
        mp4_size=$(du -h de_evolution.mp4 | cut -f1)
        gif_size=$(du -h de_evolution.gif | cut -f1)
        echo ""
        echo "File sizes:"
        echo "  MP4: $mp4_size"
        echo "  GIF: $gif_size"
    else
        echo "✗ Error creating GIF (MP4 still available)"
    fi

    rm "concat_list.txt"
    cd ..
else
    cd ..
    echo ""
    echo "✗ Error creating animation"
    exit 1
fi
