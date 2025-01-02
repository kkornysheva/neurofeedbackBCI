# **cqEEG**

Analysis of the (already preprocessed) competitive queuing EEG dataset, including TFRs, permutation tests, and linear mixed-effects models.
The goal was to investigate the relationship between behavior (RT, movement time) and power features to find possible targets for a EEG-neurofeedback BCI. 

## Folder Structure 
```
│projectdir          <- Project's main folder
│
├── data             <- Output from scripts (for example figures, table for further GLM analysis, ...)
│   ├── eeg...       <- each of these are run with different parameters (I deleted most to clean up space)
│
├── matlab           <- matlab scripts and functions
│   ├── functions    <- functions 
│   ├── figures      <- figures
│
├── README.md        <- Top-level README
```