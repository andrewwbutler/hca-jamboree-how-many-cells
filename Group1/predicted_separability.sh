#!/bin/bash

# Example: bash predicted_separability.sh human_pancreas_activated_quiescent_stellate.loom activated_stellate 0.03 quiescent_stellate 0.02 predictions.tsv

Rscript --vanilla prediction_metrics.R $1 $2 $3 $4 $5

python predict-separability.py features-predicted-separability.csv $6


# Trial runs
# Rscript --vanilla prediction_metrics.R human_pancreas_activated_quiescent_stellate.loom activated_stellate 0.03 quiescent_stellate 0.02 predictions.tsv
# Rscript --vanilla prediction_metrics.R human_pancreas_gamma_delta.loom gamma 0.03 delta 0.02 predictions.tsv
