% Takes in real and simulated features and compares them, to classify epileptic brain state and find optimal
% ABG values (most similar set of simulated features)

%% HEADER

clear; clc; close all;
PATH='C:\Users\idallmer\Desktop\PhD\Projects\ABG\Study';
PATH_vars=[PATH,'\Scripts&Vars\2_onset patterns\public\vars'];
PATH_code=[PATH,'\Scripts&Vars\2_onset patterns\public\functions'];

addpath(PATH_code);
cd(PATH_vars);

%% load real data features
dataset='epilepsiae'; 
load([PATH_vars, '/FEATURES-',dataset],'FEATURES','pars');
%FEATURES is a matlab structure with dimensions (file_ind,channel_ind)
% with fields like this:
%              error: ''
%             fileID: 'Seiz121_5900-1'
%         featuresF0: [22×11 double]
%     featuresF0_raw: [22×11 double]

savestring=['RESULT-',dataset];

%% set parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pars.dataset=dataset;
pars.parsweepnum=2; 
% parsweep 2: abnd=[2 14], step=0.25 bbnd=[1 25],gbnd=[1 30], step=0.5; single iteration, fs model 512, fs features 512, dt scale of noise

%data selection
pars.cluster_datatype='sim';             %chose whether prototype dataset is 'sim' or 'real'
pars.type_dataset='model1_features';     %dataset for prototype generation and PCA projection (see Dallmer-Zerbe et al., 2023)

%feature preprocessing
pars.oldonly=false;                 %use Wendling 2005 features only
pars.pcafeat=true;                  %use PCA during feature preprocessing or not
pars.n_comp=4;                      %if yes, set number of PCA components

%prototype generation
pars.nclust=4;                      %chose number of clusters for prototype generation between 1 and 6
pars.labeling='sim';                %chose whether to use 'sim' or 'real' dataset for labeling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERATE PROTOTYPES for BRAIN STATE CLASSIFICATION
prototypes = generate_prototypes(pars, PATH_vars); %see Dallmer-Zerbe et al. 2023, GITHUB REPOSITORY FOR EPILEPSY BRAIN STATE ESTIMATION


%% LOAD SIM DATA - PARSWEEP
PARSWEEP=load([PATH_vars,'\parsweep' num2str(pars.parsweepnum),'.mat']);
fss=PARSWEEP.sim_pars.fs_sim;

%normalize with regards to interictal segments (just in case it hasn't been done in a previous step)
load([PATH_vars,'\',pars.type_dataset],'featuresFS_raw','typeFS')
PARSWEEP.features=PARSWEEP.features_raw;
PARSWEEP.features=((PARSWEEP.features)...
    -mean(featuresFS_raw(typeFS==1,:),1))...
    ./std(featuresFS_raw(typeFS==1,:),[],1);

%project normalized feature values into protoype-derived PCA space
if isfield(prototypes,'PCAcoeff')
    PARSWEEP.features_prePCA=PARSWEEP.features;
    PARSWEEP.features=PARSWEEP.features*prototypes.PCAcoeff;
end

%% SET UP parallel processing 
FEATURES = convert_to_single(FEATURES);
PARSWEEP.features = single(PARSWEEP.features);

delete(gcp('nocreate'));  % close any open pool
parpool('local', 2);      % set number of workers
parfevalOnAll(@maxNumCompThreads, 0, 1);

% Define one template struct outside parfor
template = struct( ...
    'fileID', '', ...
    'error', '', ...
    'featuresF0', [], ...
    'featuresF0_raw', [], ...
    'featuresF0_prepca', [], ...
    'minerr', [], ...
    'minerrtype', [], ...
    'minerrABG', [], ...
    'minerrparsABG', [], ...
    'minerrfeatures', [] ...
);

% Preallocate the final result
nFiles = size(FEATURES,1);
nChannels = size(FEATURES,2);

% set up constants
cPARSWEEP = parallel.pool.Constant(PARSWEEP);
cPrototypes = parallel.pool.Constant(prototypes);
cFEATURES = parallel.pool.Constant(FEATURES);


%% LOOP THROUGH REAL DATA and fit ABG
parfor file_ind = 1:nFiles
    prototypes_local = cPrototypes.Value;
    PARSWEEP_local = cPARSWEEP.Value;
    FEATURES_local = cFEATURES.Value;
    
    % Each worker gets its own independent copy
    tempResult = repmat(template, 1, nChannels);

    for channel_ind = 1:nChannels
        disp(['File ', num2str(file_ind),', Chan ', num2str(channel_ind),' of ',num2str(nChannels)]);

        % Fill in fileID
        tempResult(channel_ind).fileID = FEATURES_local(file_ind, channel_ind).fileID;

        % Skip empty data
        if isempty(FEATURES_local(file_ind, channel_ind).featuresF0)
            tempResult(channel_ind).error = 'no data';
            continue;
        end

        % Extract features
        featuresF0 = FEATURES_local(file_ind, channel_ind).featuresF0;
        featuresF0_raw = FEATURES_local(file_ind, channel_ind).featuresF0_raw;

        % PCA projection
        if isfield(prototypes_local, 'PCAcoeff')
            tempResult(channel_ind).featuresF0_prepca = featuresF0(:, prototypes_local.featureind, :);
            featuresF0 = featuresF0 * prototypes_local.PCAcoeff;
        end

        % Brain State Classification
        [minerrtype, minerr] = classify_type(featuresF0, prototypes_local.centroids);

        % ABG Optimization
        [minerrpars, minerrABG, minerrABGfeatures] = ...
            optimize_abg(featuresF0, PARSWEEP_local.features, PARSWEEP_local.pars);

        % Save results
        tempResult(channel_ind).error = '';
        tempResult(channel_ind).featuresF0 = featuresF0;
        tempResult(channel_ind).featuresF0_raw = featuresF0_raw;
        tempResult(channel_ind).minerr = minerr;
        tempResult(channel_ind).minerrtype = minerrtype; % THIS IS YOUR CLASSIFIED BRAIN STATE 
        tempResult(channel_ind).minerrABG = minerrABG;
        tempResult(channel_ind).minerrparsABG = minerrpars; % THIS IS YOUR FITTED ABG
        tempResult(channel_ind).minerrfeatures = minerrABGfeatures;
    end

    % Safe sliced assignment
    parsave([PATH_vars,'\tempResult\tempResult_',num2str(file_ind)], tempResult);

end

% Save results after parfor completes
RESULT(nFiles, nChannels) = template;
for file_ind = 1:nFiles
    load([PATH_vars,'\tempResult\tempResult_',num2str(file_ind)])
    RESULT(file_ind, :) = tempResult;
    clear tempResult
end
save(fullfile(PATH_vars, savestring), 'RESULT', 'pars', '-v7.3');



% --- helper function to save single file output ---
function parsave(fname, tempResult)
    save(fname, 'tempResult', '-v7.3');
end

% --- helper function to comvert Features to singles ---
function out = convert_to_single(in)
    out = in;
    for i = 1:numel(in)
        if isfield(in(i),'featuresF0')
            out(i).featuresF0 = single(in(i).featuresF0);
        end
        if isfield(in(i),'featuresF0_raw')
            out(i).featuresF0_raw = single(in(i).featuresF0_raw);
        end
    end
end