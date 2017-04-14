# Source-Space-ICA-for-EEG
This Matlab code is the implementation of the Source-space ICA (Jonmohamadi et al. (2014 and 2016)) for EEG and MEG soure reconstruction. 
The articles are available at:

http://nzbri.org/resources/publications/130/Jonmohamadi_NeuroImage_2014.pdf

https://pdfs.semanticscholar.org/f306/06490cfa7e3404ef5288ca24652bcdca7dfd.pdf

Dependencies:

fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/
*Above codes are tested on fieltrip version 20170410.

How to use:

1) Download the fieldtrip toolbox and add its path to your MATLAB (folders and subfolders) 
2) Download the contents of this repository and put them in a folder
3) For EEG, run the code SourceSpaceICA_Demo_EEG.m
4) For MEG, run the code SourceSpaceICA_Demo_MEG.m

The code 'SourceSpaceICA_Demo_EEG.m' and 'SourceSpaceICA_Demo_MEG.m' simulates 4 sources and add white noise to them. Then applies source-space ICA to find the sources and finally it plots their brain space maps with two different approaches. You need to provide which component to be plotted. By default it plots component 1 (Current_comp = 1). 

A few hints:
To use this code for your study you need to provide the EEG (or MEG), the head model, and sensor information (elec for EEG or grad for MEG) yourself. 
You may use the prebuilt forward model shown here, but it is important to have the subject specific sensor file from your own set up. 
You need to make sure your EEG oor MEG data has the same structure as simulated data shown here otherwise it may not work (i.e., Data.trial, Data.time, Data.fsample, Data.label).
When calling the source-space ICA, you need to know that it requires substantial amount of RAM from your computer.
To reduce the required RAM, you can reduce the resolution of the scanning grid, or reduce the length of the data or reduce the sampling rate of the date, which can be done using cfg.grid.resolution before calling ft_prepare_leadfield function, cfg.NoTrials and cfg.ReSampleFs before calling Source_Space_ICA_Beta.

For any other question please contact Yaqub Jonmohamadi at y.jonmo@auckland.ac.nz
