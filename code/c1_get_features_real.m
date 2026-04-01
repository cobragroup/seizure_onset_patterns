% Segments real EEG data provided into 5 second segments and calculates 11 signal features for each segment 

%% HEADER

clear; clc; close all;
PATH='...';
PATH_vars=[PATH,'\vars'];
PATH_code=[PATH,'\functions'];
PATH_data=[PATH,'\data'];

addpath(PATH_code);
cd(PATH_vars);

%% set PARS and savestring

pars.segmlength=5;
pars.feature_names={'b0power', 'b1power', 'b2power',...
        'b3power', 'b4power', 'alphdiff', 'spikeabs', 'sigmean', 'sigvar', 'autocorrel', 'linelen'};
pars.norm_segm=1:30/pars.segmlength; %[]; %segments that should be used for feature normalization (here: interictal segments, use '[]' for all)

savestring='FEATURES-epilepsiae';

%% SCRIPT

%get a list of available seizure files in your data folder
cd(PATH_data);
files=dir('Seiz*.mat');

%prepare FEATURE struct
FEATURES(length(files),1)=struct;

%loop through files
for file_ind=1:length(files) 
   
    %% load data file
    load([PATH_data,'/',files(file_ind).name],'data','info');   
%   seizure data from EPILEPSIAE database have been extracted into files called e.g. "Seiz133_5900.mat" 
%   containing "data(channel_ind,time_ind)" array and "info" structure of this format:
%               file_name: {'59000102_0028.data'}
%         seiz_start_time: 2006-07-11 23:30:42.429
%           seiz_end_time: 2006-07-11 23:31:56.990
%             seiz_length: 74.5605
%      seiz_onset_pattern: {'r'}
%     seiz_classification: {'UC'}
%              patient_id: '5900'
%                      fs: 1024
%        seiz_start_index: 674232
%          seiz_end_index: 750582
%              chan_names: {1×98 cell}
%          chan_names_soz: {1×17 cell}
%                     SOZ: [1×98 logical]
    
    %loop through channels
    for channel_ind=1:size(data,1)
                
        %check for data errors
        if isempty(data(channel_ind,:))
            FEATURES(file_ind,channel_ind).error='no data';
            continue;
        elseif size(data,2)<80*info.fs %60 seconds before seizure + 20 seconds
            FEATURES(file_ind,channel_ind).error='data too short';
            continue;
        else
            FEATURES(file_ind,channel_ind).error='';
        end
        
        fileID=[files(file_ind).name(1:end-4),'-',num2str(channel_ind)];
        FEATURES(file_ind,channel_ind).fileID=fileID;
       
        %get data per channel
        d=data(channel_ind,:);
        
        %% CENTER DATA
        d=d-mean(d);
        
        %% DATA SEGMENTATION and F0 FEATURE VECTORS FOR EACH SEGMENT (F0 is just some naming convention)
        
        datalength=length(d)/info.fs; %in seconds
        Nsegm=floor(datalength/(pars.segmlength));
        
        %loop through segments
        featuresF0=zeros(Nsegm,length(pars.feature_names));
        timeind=1;
        for segm=1:Nsegm 
            
            %select data by segment time window
            d_s=d(timeind:timeind+pars.segmlength*info.fs-1);
            
            %calculate features 
            F0=calc_feat(d_s,info.fs);
            featuresF0(segm,:)=F0;
            
            %prepare next iteration
            timeind=timeind+(pars.segmlength)*info.fs;
            clear d_s F0
        end
        clear timeind segm
        
        
        %% select & normalize features
        featuresF0_raw=featuresF0;
        if isempty(pars.norm_segm)
            featuresF0=(featuresF0-mean(featuresF0,'omitnan'))./std(featuresF0,'omitnan');
        else
            featuresF0=featuresF0(:,1:length(pars.feature_names));
            featuresF0=(featuresF0-mean(featuresF0(pars.norm_segm,:),'omitnan'))./std(featuresF0(pars.norm_segm,:),'omitnan');   
        end
        
        %% save result
        FEATURES(file_ind,channel_ind).featuresF0=featuresF0;
        FEATURES(file_ind,channel_ind).featuresF0_raw=featuresF0_raw;

        clear featuresF0 featuresF0_raw datalength Nsegm
    end
    
    %% save
    save([PATH_vars, '/',savestring],'FEATURES','pars',"-v7.3")
    
end
