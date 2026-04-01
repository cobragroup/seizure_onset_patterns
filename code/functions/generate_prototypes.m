function [prototypes, X, real_type] = generate_prototypes(pars, PATHIN)

disp('GENERATING PROTOTYPES...')


%% PREPARATIONS

%load features of selected prototype dataset
cd(PATHIN);
load(pars.type_dataset);

if strcmp(pars.cluster_datatype,'real')
  
    %create X structure containing the prototype dataset
    X.features=featuresF0; %featuresF0 contain REAL data features
    X.data=DATA;
    X.datapars=datapars;
    X.fs=fs;
    real_type=typeF0';
    
elseif strcmp(pars.cluster_datatype,'sim')
   
    %create X structure containing the prototype dataset
    X.features=featuresFS;
    X.data=DATA;
    X.datapars=simpars;
    X.fs=fs;
    real_type=typeFS';
    
end


%% FEATURE PREPROCESSING BEFORE CLUSTERING

[X,N_features,featurenames,featureind,coeff]=...
    preprocess_features(X,pars);


%% FEATURE CLUSTERING

%cluster
clusters={'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5','cluster 6'};
prototypes.type=clusters(1:pars.nclust);
[cluster_type_pred, centroids] = kmeans(X.features,pars.nclust,'Replicates',100);

%test classification procedure using the centroids
cluster_type_test=classify_type(X.features,centroids);
disp('test type classification: summed classification error');
disp(sum(cluster_type_test-cluster_type_pred));


%% CENTROID LABELING 

%label helping dataset based on user-set parameters
if strcmp(pars.labeling,'sim')%labels according to overlap with wendling type
    load('model1_features');
    cluster_type_pred=classify_type(featuresFS(:,featureind)*coeff,centroids);
    type={'intICTAL-like','preONSET-like','ONSET-like','ICTAL-like'};
    cluster_type_real=[ones(1,100),repmat(2,1,100),...
        repmat(3,1,100),repmat(4,1,100)]';
elseif strcmp(pars.labeling,'real')%labels according to overlap with real type
    cluster_type_pred=classify_type(featuresF0(:,featureind)*coeff,centroids);
    cluster_type_real=typeF0';
    type={'intICTAL','preONSET','ONSET','ICTAL'};
end

%label centroids based on helping dataset
[centroids_labeled,unique_type_match,type_merged]=...
    type_labeling(centroids,cluster_type_real,cluster_type_pred);


%% COLLECT NEW PROTOTYPE STRUCTURE

prototypes.type=type;
prototypes.centroids=centroids_labeled;
prototypes.labeling=pars.labeling;
prototypes.type_n=sum(~isnan(centroids_labeled(:,1)));
prototypes.typematch_match=unique_type_match;
prototypes.typematch_merged=type_merged;
prototypes.N_features=N_features;
prototypes.featurenames=featurenames;
if pars.pcafeat==true
    if strcmp(pars.cluster_datatype,'sim')
        prototypes.featurenames_prePCA=simpars.featurenames(featureind);
        prototypes.normfeat=simpars.normalizefeatures;
    elseif strcmp(pars.cluster_datatype,'real')
        prototypes.featurenames_prePCA=datapars.featurenames(featureind);
        prototypes.normfeat=datapars.normalizefeatures;
    end
    prototypes.PCAcoeff=coeff;
end
prototypes.featureind=featureind;

disp('NEW PROTOTYPES');
disp(prototypes);
end