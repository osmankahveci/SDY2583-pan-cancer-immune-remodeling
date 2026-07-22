#!/usr/bin/env python3
"""Run convergence, immunotype, PCA, and clustering analyses."""
import argparse
import os
from final_analysis.convergence_immunotype import run as run_convergence
from final_analysis.pca_clustering import run as run_pca

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=os.getenv("SDY2583_ALL10_MATRIX", ""))
    parser.add_argument("--output-dir", default=os.getenv("SDY2583_FINAL_RESULTS_DIR", "results/local-final"))
    args = parser.parse_args()
    run_convergence(args.matrix, args.output_dir)
    run_pca(args.matrix, args.output_dir)
