% This code demonstrates the application of source-space ICA for EEG source
% reconstruction
% The reference article is Jonmohamadi et al. (2014, NeuroImage)

%% Using the prebuild BEM on T1 provided by the fieldtrip
% FieldTripPath = uigetdir;
% addpath(genpath(FieldTripPath))   % When GUI popped up, highlight the fieldtrip folder and press open
load(['/standard_bem.mat']); %template boundary element model
load(['/standard_mri.mat']); %template mri

elec = [] ; 
elec.elecpos = [-29     0    30   -55   -34    36    56   -70   -64   -50   -27     0    30    52    68    73   -84   -81   -77   -60   -34     0    35    62    80    82    84   -84   -80   -65   -36     0    38  67    83    85   -86   -85   -80   -64   -36     0    38    67    83    86    86   -72   -67   -53   -29     0    32    56    68    73   -55   -37     0    37    56   -29     0    30
     84    88    85    69    77    78    70    42    48    53    57    59    58    54    50    44    15    14    19    23    26    27    26    24    20    15    14   -16   -14   -12   -10    -9   -10   -11   -13   -15   -47   -46   -47   -47   -47   -47   -47   -47   -46   -46   -47   -73   -76   -79   -81   -81   -80   -79   -76   -73   -98  -101  -102  -101   -98  -112  -115  -112
     -7    -2    -7   -11    21    22   -11   -11    17    42    60    66    60    41    16   -12   -50   -11    24    56    80    89    79    56    24   -11   -51    -9    29    64    90   100    88 64    29    -9   -46    -7    31    66    91    99    91    66    31    -7   -46    -2    28    56    75    83    77    57    28    -3     3    37    51    36     3     9    15     9]' ;

%Chans = 1:size(elec.elecpos,1);
for I = 1:size(elec.elecpos,1)
    elec.label{I,1} = strcat('Chan',num2str(I));
end

%Generate sourcemodel and leadfields
cfg = [];
cfg.headmodel = vol;        % Came from standard_bem.mat
cfg.elec = elec;
cfg.grid.resolution = 8 ;   % use a 3-D grid with a 8mm resolution
cfg.grid.unit       = 'mm';
cfg.channel = elec.label; %Only generate leadfields for good channels else will bug out later
grid = ft_prepare_leadfield(cfg);


%% Create a few dipoles with custom timecourses and locations
Fsample = 100 ;
SourceMagnitude = 1000;

cfg      = [];
cfg.vol  = vol;               % see the calculated volume conduction
cfg.elec = elec;              % elec file which contains the label and the position of the electrodes
cfg.dip.pos = [0 0.5 30];
cfg.dip.mom = [1 0 0]';       % note, it should be transposed
cfg.fsample = Fsample;            % Hz
time = (0:Fsample-1)/Fsample;           % manually create a time axis
signal1 = sin(8*time*2*pi)*SourceMagnitude;   % manually create a signal
cfg.dip.signal = {signal1, signal1};  % two trials
raw1 = ft_dipolesimulation(cfg);
figure
plot(signal1)
title('One trial time-coure for source 1') 


cfg.dip.pos = [-31 -60.5 -10];
cfg.dip.mom = [0.5 0.2 0.5]';       
cfg.fsample = Fsample;            % Hz
signal2 = (sin(4*time*2*pi)*SourceMagnitude);   % manually create a signal
cfg.dip.signal = {signal2, signal2, signal2, signal2};  % three trials
raw2 = ft_dipolesimulation(cfg);
figure
plot(signal2)
title('One trial time-coure for source 2') 

cfg.dip.pos = [-29 40 11.1];
cfg.dip.mom = [0 0 1]';       % note, it should be transposed
cfg.fsample = Fsample;            % Hz
signal3 = (sin(16*time*2*pi)*SourceMagnitude)./((time+1).^(10));   % manually create a signal
cfg.dip.signal = {signal3, signal3, signal3, signal3, signal3};  % three trials
raw3 = ft_dipolesimulation(cfg);
figure
plot(signal3)
title('One trial time-coure for source 3') 

cfg.dip.pos = [ 9 -40 22];
cfg.dip.mom = [0.5 1 0.3]';       % note, it should be transposed
cfg.fsample = Fsample;            % Hz
signal4 = (sin(2*time*2*pi)*SourceMagnitude)./((time+1).^(6));   % manually create a signal
cfg.dip.signal = {signal4, signal4, signal4, signal4, signal4, signal4};  % three trials
raw4 = ft_dipolesimulation(cfg);
figure
plot(signal4)
title('One trial time-coure for source 4') 

Sources = [raw1, raw2, raw3, raw4];
%% Superimposing the simulated time courses and adding the white noise
SuperImposed = raw1;
No_Trials = 22
BackGroundNoise = randn(size(Sources(1).trial{1},1),size(Sources(1).trial{1},2)*No_Trials)./(SourceMagnitude/15);

SourceDelays = [0, 2, 6, 1]; % These are the delays that each source has, i.e., source number 3 starts at the 6th trial

for Trial_Index = 1:No_Trials
    SuperImposed.trial{Trial_Index} = BackGroundNoise(:, 1 + (Trial_Index - 1)*Fsample : Trial_Index*Fsample) ; 
    SuperImposed.time{Trial_Index} = time ; 
    
end

for Source_Index = 1:length(SourceDelays)
    for Trial_Index = SourceDelays(Source_Index) + 1:size(Sources(Source_Index).trial,2) + SourceDelays(Source_Index)
        SuperImposed.trial{Trial_Index} = SuperImposed.trial{Trial_Index} + Sources(Source_Index).trial{Trial_Index-SourceDelays(Source_Index)}  +Sources(Source_Index).trial{Trial_Index-SourceDelays(Source_Index)};     
    end
end
%% Plottig the simulated data 
cfg = [];
Data = SuperImposed;
cfg.viewmode = 'vertical';
cfg.colorgroups = ones([1 length(Data.label)]);
cfg.channelcolormap = [0 0 0];  %Above two lines plot in blue (events will be black)
cfg.ploteventlabels = 'colorvalue';
cfg.ylim = [-0.1 .1]; %sets yscale
cfg.blocksize = 10;
cfg.continuous              = 'yes';
ft_databrowser(cfg,Data);

%% Up to this point was all providing simulations. In reality you need to provide the sensor data, the head model, and electrode positions yourself and use the codes below for localization of the sources
%% You may use the forward model shown here, but it is important to have the electrode positions (elec.elecpos) from your own set ups.
%% You need to make sure your real data has the same format as 'SuperImposed' otherwise it may not work.
%% Calling the source-space ICA
% When calling the source-space ICA, you need to know that it requires substantial amount of RAM from your computer.
% To reduce the required RAM, you can reduce the resolution for the
% scanning grid, or reduce the length of the data duration by cfg.NoTrials
% and also reducing the sample rate of the date using the cfg.ReSampleFs
cfg = [] ; 
cfg.NoTrials = 18 ; 
cfg.vol = vol;
cfg.elec = elec;
% cfg.ReSampleFs = 100;  
cfg.grid = grid;
cfg.NumComp = 20;                      % This reduces the number of principal components (here to 20). 
                                            % You may not need this, but it is recommend for small suration EEG
SensorData = SuperImposed
[SourceSpaceStuff] = Source_Space_ICA_Beta(cfg, SensorData);
              
No_Vox = size(SourceSpaceStuff.SpatialICs_Maps,1) ;

%% Plot the temporal ICs 
cfg = [];
cfg.viewmode = 'vertical';
cfg.colorgroups = ones([1 length(SourceSpaceStuff.TemporalICs.label)]);
cfg.channelcolormap = [0 0 0];  %Above two lines plot in blue (events will be black)
cfg.ploteventlabels = 'colorvalue';
cfg.ylim = [-0.2 .2]; %sets yscale
cfg.blocksize = 10 ;
cfg.continuous              = 'yes';
ft_databrowser(cfg,SourceSpaceStuff.TemporalICs);

%% Function for 3D plotting of the components. Here you need to provide which component you are interested to look at by providing number for 'Current_comp' 
positions = grid.pos(grid.inside,:);
Current_comp =  1 ;  
Map = SourceSpaceStuff.SpatialICs_Maps(:,Current_comp) ; 
FigHandle = figure('Position', [1000, 500, 550, 170]);
plot(Map/max(Map))
ylim([0 1.05])
xlim([0 No_Vox+10])
set(gcf,'Color',[1 1 1])
set(gca,'Color',[1 1 1])
xlabel('Voxel')
ylabel('Normalized intensitycl')

figure
scatter3(positions(:,1),positions(:,2),positions(:,3),(Map/max(Map)*50),Map,'filled')
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'Color',[0.8 0.8 0.8])
set(gcf,'Color',[0.8 0.8 0.8])
colorbar;

%% Arrow plot for showing the direction of the sources 
hold on 
quiver3(positions(:,1),positions(:,2),positions(:,3),[SourceSpaceStuff.SpatialICs(1:3:No_Vox*3,Current_comp)] , [SourceSpaceStuff.SpatialICs(2:3:No_Vox*3,Current_comp)] , [SourceSpaceStuff.SpatialICs(3:3:No_Vox*3,Current_comp)],2.0)

%% Map on the MRI image.  This is a more convensional way to observe the source maps
% Current_comp = 1 ; 

Inside_count = 0;
source = grid;
for Vox_Index = 1 : size(grid.inside,1)
    if  grid.inside(Vox_Index) == 1
        Inside_count = Inside_count + 1;
        source.avg.Map(Vox_Index) = SourceSpaceStuff.SpatialICs_Maps(Inside_count,Current_comp);
    else
        source.avg.Map(Vox_Index) = 0;
    end
end
    
cfg              = [];
cfg.parameter    = 'Map';
cfg.interpmethod = 'nearest';
source_diff_int  = ft_sourceinterpolate(cfg, source, mri);

%And plot the results
source_diff_int.coordsys = 'mni';
cfg = [];
cfg.funparameter    = 'Map';
cfg.funcolorlim   = 'maxabs';
%cfg.method = 'slice';
cfg.method = 'ortho';
ft_sourceplot(cfg,source_diff_int);
