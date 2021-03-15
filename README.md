# Source-Space ICA for EEG and MEG
Source-space ICA in brief:

Source-space ICA is an approach for brain source localization and time-course reconstruction via application of singular value decomposition and ICA (independent component analysis) post minimum-variance beamforming of sensor-space EEG or MEG. Compared with the beamforming and power mapping approach, source-space ICA is far more robust to artifacts as it separates the activity of different brain sources and artifacts and provides a unique tomographic map for each. That is, the artifacts have their own maps and cannot interfere with the maps of other sources. Source-space ICA is also superior to sensor-space ICA on accuracy and spatial resolution for localization of sources. Each component identified by source-space ICA has its own tomographic map which shows the extent to which each voxel has contributed to that component.
The articles are available at:

http://nzbri.org/resources/publications/130/Jonmohamadi_NeuroImage_2014.pdf

https://pdfs.semanticscholar.org/f306/06490cfa7e3404ef5288ca24652bcdca7dfd.pdf

About this repository:

This repository contains Matlab codes and sample sensor templates for implementation of the Source-space ICA (Jonmohamadi et al. (2014 and 2016)) for EEG and MEG source reconstruction. 

Dependencies:

fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/

*The codes are tested on fieldtrip version 20170410. So this version or newer ones are recommended. 

How to use (easy!):

1) Download the fieldtrip toolbox and add its path to your MATLAB (folders and subfolders) 
2) Download the contents of this repository and put them in a folder and add its path to MATLAB 
3) For EEG, run the code SourceSpaceICA_Demo_EEG.m
4) For MEG, run the code SourceSpaceICA_Demo_MEG.m

The code 'SourceSpaceICA_Demo_EEG.m' and 'SourceSpaceICA_Demo_MEG.m' simulates 4 sources and add white noise to them. Then applies source-space ICA to find the sources and finally it plots their brain space maps with two different approaches. You need to provide which component to be plotted. By default it plots component 1 (Current_comp = 1). 


A few hints:

To use this code for your study you need to provide the EEG (or MEG), the head model, and sensor information (elec for EEG or grad for MEG) yourself. 
You may use the prebuilt head model shown here, but it is less than ideal.
You need to make sure your EEG or MEG data has the same structure as simulated data shown here otherwise it may not work, i.e., the object 'Data' which is used in these codes should have 'Data.trial', 'Data.time', 'Data.fsample', 'Data.label'). In other words, you need to replace the 'Data' in these codes with your own 'Data'.
When calling the source-space ICA, you need to know that it requires substantial amount of RAM from your computer.
To reduce the required RAM, you can reduce the resolution of the scanning grid, or reduce the length of the data or reduce the sampling rate of the date, which can be done using cfg.grid.resolution before calling ft_prepare_leadfield function, cfg.NoTrials and cfg.ReSampleFs before calling Source_Space_ICA_Beta.

For any other question please contact Yaqub Jonmohamadi at y.jonmo@gmail.com
