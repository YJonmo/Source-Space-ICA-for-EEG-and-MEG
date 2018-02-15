% This code demonstrates the application of source-space ICA for EEG source
% reconstruction
% The reference article is Jonmohamadi et al. (2014, NeuroImage)

%% Using the prebuild BEM on T1 provided by the fieldtrip
% FieldTripPath = uigetdir;
% addpath(genpath(FieldTripPath))   % When GUI popped up, highlight the fieldtrip folder and press open
load(['/standard_bem.mat']); %template boundary element model
load(['/standard_mri.mat']); %template mri

load('Subject1_faces_scramb')
%load('Subject2_faces_scramb')

%Generate sourcemodel and leadfields
cfg = [];
cfg.headmodel = vol;        % Came from standard_bem.mat
cfg.elec = data.elec;
cfg.grid.resolution = 8 ;   % use a 3-D grid with a 8mm resolution
cfg.grid.unit       = 'mm';
cfg.channel = elec.label; %Only generate leadfields for good channels else will bug out later
grid = ft_prepare_leadfield(cfg);

%% Calling the source-space ICA
% When calling the source-space ICA, you need to know that it requires substantial amount of RAM from your computer.
% To reduce the required RAM, you can reduce the resolution for the
% scanning grid, or reduce the length of the data duration by cfg.NoTrials
% and also reducing the sample rate of the date using the cfg.ReSampleFs
cfg = [] ; 
cfg.NoTrials = 1 ; 
cfg.vol = vol;
cfg.elec = data.elec;
% cfg.ReSampleFs = 100;  
cfg.grid = grid;
cfg.ReduceRankBy = 57 ;                      % This reduces the number of principal components (here by 20). 
                                            % You may not need this, but it is recommend for small suration EEG
SensorData = data
[SourceSpaceStuff] = Source_Space_ICA_Beta(cfg, SensorData);
              
No_Vox = size(SourceSpaceStuff.SpacialICs_Maps,1) ;

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
Current_comp =  6 ;  
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
