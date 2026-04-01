% Get outcome vars

clear; close all; clc;
%% set PATHs
PATH='...';
PATHOUT=[PATH,'\vars'];
addpath([PATH,'\functions']);

%% load vars
load([PATH,'vars\RESULT-epilepsiae']); 
% RESULT is a Matlab structure with fields (file_ind,channel_ind) that
% with fields like this, containing the model fitting results:
%                fileID: 'Seiz121_5900-1'
%                 error: ''
%            featuresF0: [22×4 single]
%        featuresF0_raw: [22×11 single]
%     featuresF0_prepca: [22×11 single]
%                minerr: [22×1 double]
%            minerrtype: [22×1 double]
%             minerrABG: [22×1 double]
%         minerrparsABG: [22×3 double]
%        minerrfeatures: [22×4 single]

load([PATH,'\Scripts&Vars\2_onset patterns\vars\FILEINFO-epilepsiae'],'FILEINFO'); %file info, containing clinical metadata
% FILEINFO is a Matlab structure with fields (file_ind,1) that
% with fields like this:
%          patientID: {'5900'}
%             seizID: 121
%         n_channels: 98
%                SOZ: [1×98 logical]
%                 fs: 1024
%              SzInd: [61440 1.0403e+05]
%           SzLength: 41.5947
%     SzOnsetPattern: {'r'}
           
seiz_onset=60/5; %in segm


%% single seizure: RESULTS per channel
for file_ind=file_set
    
    soz=FILEINFO(file_ind).SOZ;    
    result_length=length(RESULT(file_ind,1).minerrtype);
    file_name=replace(RESULT(file_ind,1).fileID(1:end-2),'_','-');
    
    type=reshape([RESULT(file_ind,1:length(soz)).minerrtype],result_length,[])';    
    for channel_ind=1:length(soz)
        A(channel_ind,:)=RESULT(file_ind,channel_ind).minerrparsABG(:,1);
        B(channel_ind,:)=RESULT(file_ind,channel_ind).minerrparsABG(:,2);
        G(channel_ind,:)=RESULT(file_ind,channel_ind).minerrparsABG(:,3);
    end

    clear A B G type file_name soz result_length
end


%% AVERAGE ACROSS INTERVALS AND CHANNELS
for file_ind=1:size(FILEINFO,1)
    
    if or(isempty([RESULT(file_ind,:).minerrtype]),FILEINFO(file_ind).n_channels==0)
        mean_A_all(file_ind,:)=nan(4,1);
        mean_B_all(file_ind,:)=nan(4,1);
        mean_G_all(file_ind,:)=nan(4,1);
        mean_type_all(file_ind,:)=nan(4,1);
        
        mean_A_soz(file_ind,:)=nan(4,1);
        mean_B_soz(file_ind,:)=nan(4,1);
        mean_G_soz(file_ind,:)=nan(4,1);
        mean_type_soz(file_ind,:)=nan(4,1);
        
        mean_A_other(file_ind,:)=nan(4,1);
        mean_B_other(file_ind,:)=nan(4,1);
        mean_G_other(file_ind,:)=nan(4,1);
        mean_type_other(file_ind,:)=nan(4,1);
        continue
    end
    
    %get soz
    soz=FILEINFO(file_ind).SOZ;
    
    %extract results from time windows of different epileptic states ("types")
    TYPES=get_types(RESULT(file_ind,:),seiz_onset,seiz_onset+4);
    
    %% average over intervals
    mean_A= [mean(reshape([TYPES.a_intictal]',[],length(soz)),1);...
        mean(reshape([TYPES.a_preonset]',[],length(soz)),1);...
        mean(reshape([TYPES.a_onset]',[],length(soz)),1);...
        mean(reshape([TYPES.a_ictal]',[],length(soz)),1)];
    
    mean_B= [mean(reshape([TYPES.b_intictal]',[],length(soz)),1);...
        mean(reshape([TYPES.b_preonset]',[],length(soz)),1);...
        mean(reshape([TYPES.b_onset]',[],length(soz)),1);...
        mean(reshape([TYPES.b_ictal]',[],length(soz)),1)];
    
    mean_G= [mean(reshape([TYPES.g_intictal]',[],length(soz)),1);...
        mean(reshape([TYPES.g_preonset]',[],length(soz)),1);...
        mean(reshape([TYPES.g_onset]',[],length(soz)),1);...
        mean(reshape([TYPES.g_ictal]',[],length(soz)),1)];
    
    mean_type= [mean(reshape([TYPES.types_intictal]',[],length(soz)),1);...
        mean(reshape([TYPES.types_preonset]',[],length(soz)),1);...
        mean(reshape([TYPES.types_onset]',[],length(soz)),1);...
        mean(reshape([TYPES.types_ictal]',[],length(soz)),1)];
    
    A_all(file_ind)={mean_A};
    B_all(file_ind)={mean_B};
    G_all(file_ind)={mean_G};
    type_all(file_ind)={mean_type};    
    
    A_soz(file_ind)={mean_A(:,soz)};
    B_soz(file_ind)={mean_B(:,soz)};
    G_soz(file_ind)={mean_G(:,soz)};
    type_soz(file_ind)={mean_type(:,soz)};
    
    A_other(file_ind)={mean_A(:,~soz)};
    B_other(file_ind)={mean_B(:,~soz)};
    G_other(file_ind)={mean_G(:,~soz)};
    type_other(file_ind)={mean_type(:,~soz)};
    

    %% average over channels
    mean_A_all(file_ind,:)=mean(mean_A,2,'omitnan');
    mean_B_all(file_ind,:)=mean(mean_B,2,'omitnan');
    mean_G_all(file_ind,:)=mean(mean_G,2,'omitnan');
    mean_type_all(file_ind,:)=mean(mean_type,2,'omitnan');
    
    mean_A_soz(file_ind,:)=mean(mean_A(:,soz),2,'omitnan');
    mean_B_soz(file_ind,:)=mean(mean_B(:,soz),2,'omitnan');
    mean_G_soz(file_ind,:)=mean(mean_G(:,soz),2,'omitnan');
    mean_type_soz(file_ind,:)=mean(mean_type(:,soz),2,'omitnan');
    
    mean_A_other(file_ind,:)=mean(mean_A(:,~soz),2,'omitnan');
    mean_B_other(file_ind,:)=mean(mean_B(:,~soz),2,'omitnan');
    mean_G_other(file_ind,:)=mean(mean_G(:,~soz),2,'omitnan');
    mean_type_other(file_ind,:)=mean(mean_type(:,~soz),2,'omitnan');
    
    clear file_ind result_length soz TYPES mean_A mean_B mean_G mean_type
    
end


cd(PATHOUT)
save('OutcomeVars-epilepsiae',...
    'mean_A_all','mean_B_all','mean_G_all','mean_type_all',...
    'mean_A_soz','mean_A_other','mean_B_soz','mean_B_other',...
    'mean_G_soz','mean_G_other','mean_type_soz','mean_type_other',...       
    'A_all','B_all','G_all','type_all',...
    'A_soz','A_other','B_soz','B_other',...
    'G_soz','G_other','type_soz','type_other',...
    'FILEINFO','RESULT','-v7.3')

