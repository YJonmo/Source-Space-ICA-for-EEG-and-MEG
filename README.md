# Source-Space-ICA-for-EEG
This Matlab code is the implementation of the Source-space ICA (Jonmohamadi et al. (2014, NeuroImage)) for EEG soure reconstruction. 
The article is available on:
http://nzbri.org/resources/publications/130/Jonmohamadi_NeuroImage_2014.pdf

Dependencies:

fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/

How to use:

1) Download the fieldtrip toolbox and add its path to your MATLAB (folders and subfolders) 
2) Download the contents of this repository and put them in a folder
3) Run the code SourceSpaceICA_Demo.m 

The code 'SourceSpaceICA_Demo.m' simulates 4 EEG sources and add white noise to them. Then applies source-space ICA to find the sources and at the end it plots their brain space maps with two different approaches. You need to provide which component to be plotted. By default it plots component 1 (Current_comp = 1). 

A few hints:
To use this code for your study you need to provide the EEG, the head model, and electrode positions yourself and use the codes. 
You may use the prebuilt forward model shown here, but it is important to have the subject specific electrode postions (elec.elecpos) from your own set up. 
You need to make sure your EEG data has the same structure as simulated EEG shown here otherwise it may not work (i.e., Data.trial, Data.time, Data.fsample, Data.label).
When calling the source-space ICA, you need to know that it requires substantial amount of RAM from your computer.
To reduce the required RAM, you can reduce the resolution for the scanning grid, or reduce the length of the data duration by or reduce the sample rate of the date using which can be done using cfg.grid.resolution before calling ft_prepare_leadfield function, cfg.NoTrials and cfg.ReSampleFs before calling Source_Space_ICA_Beta.

For any other question please contact Yaqub Jonmohamadi at y.jonmo@auckland.ac.nz
