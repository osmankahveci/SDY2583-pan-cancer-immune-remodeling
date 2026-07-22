#!/usr/bin/env python3
"""Run convergence, immunotype, PCA, and clustering analyses."""
from final_analysis.convergence_immunotype import main as convergence_main
from final_analysis.pca_clustering import main as pca_main

if __name__ == "__main__":
    # The release orchestrator calls the modules separately so each receives its
    # own command-line arguments. This wrapper is retained as a documented
    # entry point and intentionally delegates through the shell orchestrator.
    raise SystemExit(
        "Use scripts/40_run_final_integrated_release.sh, or run "
        "python -m final_analysis.convergence_immunotype and "
        "python -m final_analysis.pca_clustering separately."
    )
