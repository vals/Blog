#!/usr/bin/env python3
"""
Temperature Causality Analysis Runner
====================================

Simple script to run both basic and advanced causality analyses.

Usage:
    python run_analysis.py [--basic-only] [--advanced-only]

Author: Claude Code
Date: 2025-09-24
"""

import argparse
import sys
import os
from pathlib import Path

def install_requirements():
    """Check and install required packages."""
    try:
        import subprocess

        print("Checking requirements...")
        result = subprocess.run([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"],
                              capture_output=True, text=True)

        if result.returncode != 0:
            print("Warning: Some packages may not have installed correctly")
            print(result.stderr)
        else:
            print("Requirements check complete")

    except Exception as e:
        print(f"Could not install requirements: {e}")
        print("Please install manually: pip install -r requirements.txt")

def run_basic_analysis():
    """Run basic causality analysis."""
    print("\n" + "="*60)
    print("RUNNING BASIC CAUSALITY ANALYSIS")
    print("="*60)

    try:
        # Import and run basic analysis
        import importlib.util
        spec = importlib.util.spec_from_file_location("causal_analysis", "causal_temperature_analysis.py")
        causal_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(causal_module)

        # Run main function
        causal_module.main()

        print("\n✅ Basic analysis completed successfully!")

    except Exception as e:
        print(f"\n❌ Basic analysis failed: {e}")
        import traceback
        traceback.print_exc()

def run_advanced_analysis():
    """Run advanced causality analysis."""
    print("\n" + "="*60)
    print("RUNNING ADVANCED CAUSALITY ANALYSIS")
    print("="*60)

    try:
        # Import and run advanced analysis
        import importlib.util
        spec = importlib.util.spec_from_file_location("advanced_analysis", "advanced_causality_analysis.py")
        advanced_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(advanced_module)

        # Run main function
        advanced_module.main()

        print("\n✅ Advanced analysis completed successfully!")

    except Exception as e:
        print(f"\n❌ Advanced analysis failed: {e}")
        print("This is likely due to missing optional dependencies (tigramite, pyinform)")
        import traceback
        traceback.print_exc()

def main():
    """Main runner function."""
    parser = argparse.ArgumentParser(description="Run temperature causality analysis")
    parser.add_argument("--basic-only", action="store_true",
                       help="Run only basic analysis (Granger causality)")
    parser.add_argument("--advanced-only", action="store_true",
                       help="Run only advanced analysis (TE, PCMCI)")
    parser.add_argument("--skip-install", action="store_true",
                       help="Skip automatic package installation")

    args = parser.parse_args()

    print("Temperature Sensor Causality Analysis Runner")
    print("=" * 45)

    # Check if data file exists
    data_file = Path("Data/history.csv")
    if not data_file.exists():
        print(f"❌ Data file not found: {data_file}")
        print("Please ensure Data/history.csv exists in the current directory")
        return 1

    print(f"✅ Found data file: {data_file}")
    print(f"   File size: {data_file.stat().st_size / (1024*1024):.1f} MB")

    # Install requirements unless skipped
    if not args.skip_install:
        install_requirements()

    # Determine which analyses to run
    run_basic = not args.advanced_only
    run_advanced = not args.basic_only

    # Run analyses
    if run_basic:
        run_basic_analysis()

    if run_advanced:
        run_advanced_analysis()

    # Summary
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)

    generated_files = []

    # Check for generated files
    possible_files = [
        "temperature_time_series.png",
        "causal_network.png",
        "temperature_causality_report.txt",
        "pcmci_causal_network.png",
        "causality_heatmaps.png",
        "causality_methods_comparison.csv"
    ]

    for file in possible_files:
        if Path(file).exists():
            generated_files.append(file)

    if generated_files:
        print("Generated files:")
        for file in generated_files:
            file_size = Path(file).stat().st_size
            print(f"  ✅ {file} ({file_size:,} bytes)")
    else:
        print("⚠️  No output files found")

    print(f"\nTo explore results interactively, consider using:")
    print(f"  jupyter notebook temperature_causality_exploration.ipynb")

    return 0

if __name__ == "__main__":
    sys.exit(main())