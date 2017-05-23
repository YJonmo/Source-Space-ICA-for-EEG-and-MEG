% This code demonstrates the application of source-space ICA for EEG source
% reconstruction
% The reference article is Jonmohamadi and Jones (2016, J Neural Eng)
%% Create a single shell volume conductor
% load T1 MRI
load(['/standard_mri.mat']); %template mri 
cfg = [];
cfg.write      = 'no';
[segmentedmri] = ft_volumesegment(cfg, mri);
% mri = ft_volumereslice([], mri);
cfg = [];
cfg.method = 'singleshell';
%load SegmentedMRI_T1
headmodel = ft_prepare_headmodel(cfg, segmentedmri);

%% Generate sourcemodel and leadfields
% The file dataFIC contains information on sensor positions and direction
% (from fieldtrip)
load dataFIC
grad = dataFIC.grad ; 

cfg                 = [];
cfg.grad            = grad;
cfg.headmodel       = headmodel;
%cfg.reducerank      = 2;
cfg.channel         = {'MEG'};  % Channels 1 to 151 are the MEG channels
cfg.grid.resolution = 8;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'mm';
[grid] = ft_prepare_leadfield(cfg);

figure
scatter3(grad.pnt(1:151,1),grad.pnt(1:151,2),grad.pnt(1:151,3),20)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('positions of the simulated coils')

%% Creating a few dipoles with custom timecourses and locations
Fsample = 100 ;
SourceMagnitude = 200000000000;

cfg = [];
cfg.vol = headmodel;
cfg.channel         = {'MEG'};
cfg.grad = dataFIC.grad;              % grad 
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
No_Trials = 22 ;
BackGroundNoise = randn(size(Sources(1).trial{1},1),size(Sources(1).trial{1},2)*No_Trials)./10;

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
cfg.ylim = [-0.5 0.5]; %sets yscale
cfg.blocksize = 10;
cfg.continuous              = 'yes';
ft_databrowser(cfg,Data);


%% Calling the source-space ICA
cfg = [] ; 
cfg.NoTrials = 18 ; 
cfg.vol = headmodel;
cfg.grad = dataFIC.grad;
% cfg.ReSampleFs = 100;  
cfg.grid = grid;
cfg.ReduceRankBy = 100;                      % This reduces the number of principal components (here by 20). 
                                            % You may not need this, but it is recommend for short duration data
SensorData = SuperImposed
[SourceSpaceStuff] = Source_Space_ICA_Beta(cfg, Data);
              
No_Vox = size(SourceSpaceStuff.SpacialICs_Maps,1) ;


%% Plot the temporal ICs 
cfg = [];
cfg.viewmode = 'vertical';
cfg.colorgroups = ones([1 length(SourceSpaceStuff.TemporalICs.label)]);
cfg.channelcolormap = [0 0 0];  %Above two lines plot in blue (events will be black)
cfg.ploteventlabels = 'colorvalue';
cfg.ylim = [-0.5 0.5] ; %sets yscale
cfg.blocksize = 10 ;
cfg.continuous              = 'yes';
ft_databrowser(cfg,SourceSpaceStuff.TemporalICs);

%% Function for 3D plotting of the components. Here you need to provide which component you are interested to look at by providing number for 'Current_comp' 

Current_comp =  1 ;

positions = grid.pos(grid.inside,:);  
Map = SourceSpaceStuff.SpacialICs_Maps(:,Current_comp) ; 
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
quiver3(positions(:,1),positions(:,2),positions(:,3),[SourceSpaceStuff.SpacialICs(1:3:No_Vox*3,Current_comp)] , [SourceSpaceStuff.SpacialICs(2:3:No_Vox*3,Current_comp)] , [SourceSpaceStuff.SpacialICs(3:3:No_Vox*3,Current_comp)],2.0)


%% Map on the MRI image.  This is a more convensional way to observe the source maps
% Current_comp = 1 ; 

Inside_count = 0;
source = grid;
for Vox_Index = 1 : size(grid.inside,1)
    if  grid.inside(Vox_Index) == 1
        Inside_count = Inside_count + 1;
        source.avg.Map(Vox_Index) = SourceSpaceStuff.SpacialICs_Maps(Inside_count,Current_comp);
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
