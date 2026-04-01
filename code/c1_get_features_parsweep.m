% Simulate WENDLING model and extract 11 signal features per [ABG] parameter set

%% SET PATHS
clear; clc; close all; 

PATH='...';
PATH_vars=[PATH,'\vars'];
PATH_code=[PATH,'\functions'];

addpath(PATH_code);
cd(PATH_vars);

%% SET PARAMETERS
parsweepnum=2; 

% sim pars
duration=5;             %simulation length in seconds
noise=[];               %0; %[] will generate random noise with mean 90, std 30
fs_sim=512;             %sampling rate for time series generation
fs_feat=512;            %sampling rate for features
dt_scale=true;          %rescale noise as in (Fietkiewicz & Loparo, 2016)

%store in structure
sim_pars.duration = duration; 
sim_pars.noise = noise; 
sim_pars.fs_sim = fs_sim; 
sim_pars.fs_feat = fs_feat; 
sim_pars.dt_scale = dt_scale;

% sweep bounds
sim_pars.abnd=[2 14]; sim_pars.astep=0.25;
sim_pars.bbnd=[1 25]; sim_pars.bstep=0.5;
sim_pars.gbnd=[1 30]; sim_pars.gstep=0.5;


%%%%%%%%%%%%%%%%%%
sim_pars.type_dataset='model1_features';     %dataset for prototype generation and PCA projection (see Dallmer-Zerbe et al., 2023)
load([PATH_vars,'\',sim_pars.type_dataset],'featuresFS_raw','typeFS')

%% PARSWEEP

% Create all combinations of (A,B,G) 
sim_pars.avals=sim_pars.abnd(1):sim_pars.astep:sim_pars.abnd(end);
sim_pars.bvals=sim_pars.bbnd(1):sim_pars.bstep:sim_pars.bbnd(end);
sim_pars.gvals=sim_pars.gbnd(1):sim_pars.gstep:sim_pars.gbnd(end);

[A,B,G] = ndgrid(sim_pars.avals, sim_pars.bvals, sim_pars.gvals);
pars = [A(:), B(:), G(:)];
N = size(pars,1);
clear A B G

%initialize vars
features_raw=zeros(N,size(featuresFS_raw,2));

% Start parallel pool if needed 
if isempty(gcp('nocreate'))
    parpool;  % or parpool(8) for 8 workers
end

tic
parfor parset = 1:N
    a = pars(parset,1);
    b = pars(parset,2); 
    g = pars(parset,3);
    
    %simulate data
    d=Wendl2005stimnoise(duration,a,b,g,0,noise,fs_sim,dt_scale,0); close
    
    %calculate features (best identified in errfuncdefinition.m)
    features_raw(parset,:)=calc_feat(d,fs_feat);
    
end
toc
clear parset N duration a b g noise fs_sim dt_scale


%% Z NORMALIZE FEATURES WITH REGARD TO INTERICTAL
features=(features_raw-mean(featuresFS_raw(typeFS==1,:),1))...
    ./std(featuresFS_raw(typeFS==1,:),[],1);

%% SAVE
save([PATH_vars,'\parsweep', num2str(parsweepnum)]);
%%%

