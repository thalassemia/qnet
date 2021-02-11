# Interaction Heatmaps

Splits task of creating heatmap into two steps:
1. [interactMatrix.py](interactMatrix.py): Compile and save two separate matrices of average and maximum interaction values for each set of 100 models as described [here](../info.md#interaction-heatmaps) (time-consuming as it involves sequentially loading 100 saved python shelves for each set).
1. [interactMap.py](interactMap.py): Load saved interaction matrices and create heatmaps (much faster than first step).
    - All heatmaps are clustered along both their rows and columns using Euclidean distance and Ward linkage, before the resulting dendrograms are arranged using dendsort.