# Create Interaction Heatmaps

**Dependencies:** python 3.8, numpy 1.19.2, pandas 1.1.3, seaborn 0.10.1, matplotlib 3.2.2

Creates two heatmaps of SHAP interaction values. To create the heatmap of averages, the SHAP interaction values for each TF pair are first averaged across all positive samples for a given model, then these averages were averaged again over all models. To create the heatmap of maximums, the SHAP interaction values at each positive sample are averaged across all models and the maximum (across all samples) is taken for each TF pair.

Example: For a model with 200 TFs (features) and 1000 DE genes (2000 samples after accounting for negative set), we have `200 * 200 * 2000` SHAP interaction values organized into 200-by-200 matrices for each of the 2000 samples. If we have 100 such models, we can consolidate all of these SHAP interaction values into a 4D array with the dimensions `(100, 2000, 200, 200)`. 

- All negative samples are discarded, leaving an array with the dimensions `(100, 1000, 200, 200)`.
- Heatmap of averages: Take the average across all samples to get an array with the dimensions `(100, 200, 200)`. Take the average across all models to get an array with the dimensions `(200, 200)`.
- Heatmap of maximums: Take the average across all models to get an array with the dimensions `(1000, 200, 200)`. Take the maximum across all samples to get an array with the dimensions `(200, 200)`.