# Epilepsy-Brain-Graph-Analysis
This project focuses on identifying statistically significant changes between 18 EEG channels in both the time and frequency domains through univariate and bivariate feature extraction. The dataset used is the [CHB-MIT](https://physionet.org/content/chbmit/1.0.0/), which contains epileptic EEG recordings from 22 subjects. The preprocessing pipeline includes the following steps:
1. High-pass filtering at 0.5 Hz
2. Re-referencing to a common average
3. Cleanline filtering at 60 Hz
4. Independent Component Analysis (ICA) to remove noise and artifacts
Graph learning techniques are then employed to construct adjacency matrices, and graphical features are analyzed to uncover trends and changes across different seizure phases and frequency bands.
