# fc-fMRI_processing_1stLevel_analysis
Scripts to run processing and 1st-level functional connectivity analyses in CONN/Matlab

"dg_fcmri_rs_wholeGroup.m" runs default MNI processing available in CONN, and 1st-level ROI-to-ROI and Seed-to-Voxel analyses on RESTING-STATE DATASETS. Ultimately, it allows to open CONN's GUI to visualize results and perform 2nd-level analysis. The script runs these processes with a WHOLE group of patients and controls. Without 2nd-level covariates (groups), the script allows a Quality Check and explores FC tendencies in the sample measured.

"dg_fcmri_wl_wholeGroup.m" runs the same processing and 1st-level analysis but calculates FC within phonological or visual conditions from a replication of the language paradigm of Ebner et al. (2011) [https://doi.org/10.1016/j.neuroimage.2011.06.048]. As a result, it is possible to run with the GUI a 2nd-level analysis to compare FC between patients and controls during phonological processing or visual pattern discrimination. The script includes a mask to calculate FC only in ROIs relevant for phonological awareness and reading, downloaded from Neurosynth: (https://neurosynth.org/analyses/terms/phonological/) AND (https://neurosynth.org/analyses/terms/reading/).

"dg_fcmri_wl_wholeGroup.m" runs the same processing and 1st-level analysis but calculates FC within encoding or retrieval conditions from a replication of the verbal working memory paradigm of Siffredi et al. (2017) [https://doi.org/10.1371/journal.pone.0179959]. As a result, it is possible to run with the GUI a 2nd-level analysis to compare FC between patients and controls during encoding or retrieval. The script includes a mask to calculate FC only in ROIs relevant for verbal working memory, downloaded from Neurosynth:(https://neurosynth.org/analyses/terms/verbal%20working/).


