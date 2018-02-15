function [SourceSpace] = Source_Space_ICA_Beta(cfg, dataToRun)


% Sample using syntax:
% cfg = [] ; 
% cfg.NoTrials = 100 ;
% cfg.vol = vol;
% cfg.elec = elec;      % for EEG
% cfg.grad = grad;      % for MEG    
% cfg.grid = grid;      
% cfg.ReSampleFs = 100;
% cfg.NumComp = 20;
% SensorData = data_simulation_Trialed
% [SourceSpaceStuff] = Source_Space_ICA(cfg, SensorData);


% This function performs beamforming on the input data and then performs SVD to separate the Spatial and temporal subspaces, then
% performs ICA on the temporal subspace and finally estimates the Spatial
% maps of the temporal ICs. 
% The necessary parts of the cfg are the 
% cfg.vol which contains volume conduction and cfg.elec for EEG or cfg.grad for MEG 
% for the location of the sensor and their labels. 

% The optional inputs are  
% cfg.ReSampleFs (almost necessary to reduce the size of the required
% memory) when the sampling rate is high (e.g., 500)
% cfg.NoTrials which indicates how many trials of the given data shall be processed. If not give, all the trials will be processed
% cfg.NumComp which shows how many of the components shall be considered
% by SVD. If not given nothing will be removed. 

% The outputs of this function are:
% SourceSpace.data which is the source-space projected data via the vector
% beamfomer
% SourceSpace.TemporalPCs:      the temporal principal components ;
% SourceSpace.TemporalICs:      the temporal independent components ;
% SourceSpace.SpatialPCs:       the Spatial principal components  (x,y,z format);
% SourceSpace.SpatialICs:       the Spatial independent components (x, y, z format) ;
% SourceSpace.SpatialPCs_Maps:  the principal Spatial maps for the purpose plotting;
% SourceSpace.SpatialICs_Maps:  the independent Spatial maps for the purpose of plotting;
% SourceSpace.MixingMatrix:     the mixing matrix found by the ICA on the temporal data ;
% SourceSpace.filters:          the beamformer coefficients ;




%%
if  isfield(cfg, 'NoTrials') 
    trials = 1:cfg.NoTrials;
else
    trials = 1:length(dataToRun.trial)
end

%% Data cannot be too long for the source-space ICA (or it can contain all the trials but down-sampled)
% This part may reduce the data length as provided by the user.
data_shortened = dataToRun;
data_shortened = rmfield(data_shortened, 'trial')
data_shortened = rmfield(data_shortened, 'time')
for Trial_Index = trials
    data_shortened.trial{Trial_Index} = dataToRun.trial{Trial_Index};
    data_shortened.time{Trial_Index} = dataToRun.time{Trial_Index}; 
end

%% Down sample the data if user provided
if  isfield(cfg, 'ReSampleFs') 
    cfg2 = [];
    cfg2.resamplefs = cfg.ReSampleFs;
    data_shortened = ft_resampledata(cfg2, data_shortened)
end

%% Vector LCMV beamformer 
cfg2 = [];
cfg2.covariance = 'yes';
cfg2.covariancewindow = 'all';
timelock = ft_timelockanalysis(cfg2, data_shortened);

% EEG or MEG 
cfg2= [];
cfg2.headmodel = cfg.vol;
if  (isfield(cfg, 'elec'))      % EEG
    cfg2.elec = cfg.elec;
elseif(isfield(cfg, 'grad'))    % MEG
    cfg2.grad = cfg.grad;
else
    error ('please provide sensor file (elec for EEG or grad for MEG')
end
cfg2.grid = cfg.grid;
cfg2.method = 'lcmv';
cfg2.lcmv.projectnoise='yes'; %needed for neural activity index
%cfg2.lcmv.fixedori = 'yes';  %Project onto largest variance orientation
cfg2.lcmv.keepfilter = 'yes'; %Keep the beamformer weights
cfg2.lcmv.lambda = '5%'; %Regularise a little
source = ft_sourceanalysis(cfg2, timelock);

filters = cell2mat(source.avg.filter(source.inside == 1)); 

%% Beamformer's weight normalization
No_Vox = size(filters,1)/3;
for Vox_Index = 1:No_Vox*3
    filters(Vox_Index,:) = filters(Vox_Index,:)/norm(filters(Vox_Index,:));
end
    
%% Converting to the continuous data
%Source_space = data_shortened;

data_continuous = data_shortened;
data_continuous = rmfield(data_continuous, 'trial');
data_continuous = rmfield(data_continuous, 'time');
template_time = (0: (data_shortened.time{1}(2) - data_shortened.time{1}(1)) : 1/data_shortened.fsample*size(data_shortened.time{1},2)) ; 
data_continuous.trial{1}  = zeros(size(data_shortened.trial{1},1),  size(data_shortened.trial{1},2) * length(data_shortened.trial));
data_continuous.time{1}  = zeros(1,                             size(data_shortened.time{1},2)  * length(data_shortened.time));
for Trial_Index = 1:length(data_shortened.trial)
    data_continuous.trial{1}(:,(Trial_Index-1)*size(data_shortened.trial{Trial_Index},2)+ 1: Trial_Index*size(data_shortened.trial{Trial_Index},2)) = data_shortened.trial{Trial_Index}(:,:); 
    data_continuous.time{1}(1,(Trial_Index-1)*size(data_shortened.time{Trial_Index},2)+ 1: Trial_Index*  size(data_shortened.time{Trial_Index},2)) = template_time(1:end-1) + (Trial_Index-1)*template_time(end) ;
end

%Source_space = trial2continuous(data_shortened);
Source_space = data_continuous ; 
Source_space.continuous{1} = zeros(size(filters,1), size(data_shortened.trial{1},2)*length(data_shortened.trial));
Source_space.trial = data_shortened.trial;

for Trial_Index = 1:length(data_shortened.trial)
    ft_progress(Trial_Index/length(data_shortened.trial), 'Processing trial %d from %d', Trial_Index, length(data_shortened.trial));
    Source_space.continuous{1}(:, (Trial_Index-1)*size(data_shortened.time{Trial_Index},2) + 1: Trial_Index*size(data_shortened.time{Trial_Index},2)) = (filters * data_shortened.trial{Trial_Index}); 
    Source_space.trial{Trial_Index} = filters * data_shortened.trial{Trial_Index};
end
ft_progress('close')

                    
%% Source-space ICA

rank_source = rank(Source_space.continuous{1})  

%[U Sig V] = svd(Source_space.continuous{1}) ; 
[V U Eigens] = pca(Source_space.continuous{1}) ; 
% It is recommended to reduce the rank when the data is short or when the
% band passed data is used as MATLAB sometimes over estimates the number
% the rank of the data matrix
if  isfield(cfg, 'NumComp')
    rank_source = cfg.NumComp ;
end

SpatialPCs  = U(:,1:rank_source) ;      % Spatial subspace
TemporalSubSpace  = V(:,1:rank_source) ;      % Temporal subspace  
%Sig_D= Sig(1:rank_source,1:rank_source) ; 
%TemporalPCs = Sig_D*TemporalSubSpace' ;   % Principal components (time-courses)
TemporalPCs = TemporalSubSpace' ;   % Principal components (time-courses)


data_TemporalSubSpace = Source_space;
data_TemporalSubSpace = rmfield(data_TemporalSubSpace, 'trial')
data_TemporalSubSpace.trial{1} = TemporalPCs;
data_TemporalSubSpace = rmfield(data_TemporalSubSpace, 'label')
%data_TemporalSubSpace.label = Source_space.label(1:rank_source);    % this is just to fool the fieldtrip to continue working

PCs = 1:rank_source;
for I = 1:rank_source
    data_TemporalSubSpace.label{I} = strcat('PC',num2str(I));
end

cfg2 = [];
[TemporalICs] = ft_componentanalysis(cfg2, data_TemporalSubSpace) ;          % Running ICA for the temporal subsapce

%G = U_D*comp.unmixing;
Mixing = pinv(TemporalICs.unmixing);                                   
SpatialICs = SpatialPCs*Mixing ;                                                % Calculating the spatial maps of the ICs

for I = 1:rank_source
    TemporalICs.label{I} = strcat('IC',num2str(I));
end

%% Producing the 3D maps from the Spatial ICs and PCs
% No_Vox = size(SpatialICs,1)/3
SpatialICs_Maps = zeros(No_Vox, rank_source);
SpatialPCs_Maps = zeros(No_Vox, rank_source);
for Curren_comp = 1:rank_source
    SpatialICs_Maps(:,Curren_comp) = sqrt(SpatialICs(1:3:No_Vox*3,Curren_comp).^2 + SpatialICs(2:3:No_Vox*3,Curren_comp).^2 + SpatialICs(3:3:No_Vox*3,Curren_comp).^2) ;
    SpatialPCs_Maps(:,Curren_comp) = sqrt(SpatialPCs(1:3:No_Vox*3,Curren_comp).^2 + SpatialPCs(2:3:No_Vox*3,Curren_comp).^2 + SpatialPCs(3:3:No_Vox*3,Curren_comp).^2) ;
end

%% Returning the processed data
SourceSpace.data = Source_space;
SourceSpace.TemporalPCs = data_TemporalSubSpace ;
SourceSpace.TemporalICs = TemporalICs ;
SourceSpace.SpatialPCs = SpatialPCs ;
SourceSpace.SpatialICs = SpatialICs ;
SourceSpace.SpatialPCs_Maps = SpatialPCs_Maps ;
SourceSpace.SpatialICs_Maps = SpatialICs_Maps ;
SourceSpace.MixingMatrix = Mixing ;
SourceSpace.filters = filters ;
