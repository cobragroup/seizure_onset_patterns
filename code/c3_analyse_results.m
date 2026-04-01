% Main Analysis of Manuscript ....

clear; close all; clc
PATH='...';
addpath(genpath([PATH,'functions']));
load([PATH,'\functions\myCustomColormap.mat']);

PATHIN=[PATH,'vars'];
load([PATHIN,'\OutcomeVars-epilepsiae'])

% OutcomeVars-epilepsiae contains prepared outcome variables, see previous
% script
% save('OutcomeVars-epilepsiae',...
%     'mean_A_all','mean_B_all','mean_G_all','mean_type_all',...
%     'mean_A_soz','mean_A_other','mean_B_soz','mean_B_other',...
%     'mean_G_soz','mean_G_other','mean_type_soz','mean_type_other',...       
%     'A_all','B_all','G_all','type_all',...
%     'A_soz','A_other','B_soz','B_other',...
%     'G_soz','G_other','type_soz','type_other',...
%     'FILEINFO','RESULT','-v7.3')

%% Analysis 1: GRAND AVERAGE ACROSS SEIZURES PER PATIENT

[patient_names,~,I]=unique([FILEINFO.patientID]);
patient_names=patient_names(~contains(patient_names,'error')); %exclude error patients
for pat_ind=1:length(patient_names)
    
    x=[FILEINFO(I==pat_ind).SzOnsetPattern];
    pat_most_common_onset{pat_ind}=mode([x{:}]);
    n_seiz(pat_ind)=sum(I==pat_ind);
    
    mean_A_all_pat(pat_ind,:)=mean(mean_A_all(I==pat_ind,:),1, 'omitnan');
    mean_B_all_pat(pat_ind,:)=mean(mean_B_all(I==pat_ind,:),1, 'omitnan');
    mean_G_all_pat(pat_ind,:)=mean(mean_G_all(I==pat_ind,:),1, 'omitnan');
    
    mean_A_soz_pat(pat_ind,:)=mean(mean_A_soz(I==pat_ind,:),1, 'omitnan');
    mean_B_soz_pat(pat_ind,:)=mean(mean_B_soz(I==pat_ind,:),1, 'omitnan');
    mean_G_soz_pat(pat_ind,:)=mean(mean_G_soz(I==pat_ind,:),1, 'omitnan');
    
    mean_A_other_pat(pat_ind,:)=mean(mean_A_other(I==pat_ind,:),1, 'omitnan');
    mean_B_other_pat(pat_ind,:)=mean(mean_B_other(I==pat_ind,:),1, 'omitnan');
    mean_G_other_pat(pat_ind,:)=mean(mean_G_other(I==pat_ind,:),1, 'omitnan');
    
    std_A_all_pat(pat_ind,:)=std(mean_A_all(I==pat_ind,:),1, 'omitnan');
    std_B_all_pat(pat_ind,:)=std(mean_B_all(I==pat_ind,:),1, 'omitnan');
    std_G_all_pat(pat_ind,:)=std(mean_G_all(I==pat_ind,:),1, 'omitnan');
    
    std_A_soz_pat(pat_ind,:)=std(mean_A_soz(I==pat_ind,:),1, 'omitnan');
    std_B_soz_pat(pat_ind,:)=std(mean_B_soz(I==pat_ind,:),1, 'omitnan');
    std_G_soz_pat(pat_ind,:)=std(mean_G_soz(I==pat_ind,:),1, 'omitnan');
    std_type_soz_pat(pat_ind,:)=std(mean_type_soz(I==pat_ind,:),1, 'omitnan');
    
    std_A_other_pat(pat_ind,:)=std(mean_A_other(I==pat_ind,:),1, 'omitnan');
    std_B_other_pat(pat_ind,:)=std(mean_B_other(I==pat_ind,:),1, 'omitnan');
    std_G_other_pat(pat_ind,:)=std(mean_G_other(I==pat_ind,:),1, 'omitnan');
    
end 
clear x I


figure
subplot(3,1,1);
boxchart(mean_A_all_pat); ylabel('mean A','FontWeight','bold'); xticklabels({'intictal','preonset','onset','ictal'});
hold on; plot(1:4,mean(mean_A_all_pat,1),'Linewidth',2); legend({'per patient','average'},'Location','northwest')
subplot(3,1,2);
boxchart(mean_B_all_pat); ylabel('mean B','FontWeight','bold'); xticklabels({'intictal','preonset','onset','ictal'});
hold on; plot(1:4,mean(mean_B_all_pat,1),'Linewidth',2); 
subplot(3,1,3);
boxchart(mean_G_all_pat); ylabel('mean G','FontWeight','bold'); xticklabels({'intictal','preonset','onset','ictal'});
hold on; plot(1:4,mean(mean_G_all_pat,1),'Linewidth',2); 

[p_A,~,stats_A] = signrank(mean_A_all_pat(:,1),mean_A_all_pat(:,4),'tail','left')
[p_B,~,stats_B] = signrank(mean_B_all_pat(:,1),mean_B_all_pat(:,4),'tail','left')
[p_G,~,stats_G] = signrank(mean_G_all_pat(:,1),mean_G_all_pat(:,4),'tail','left')

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_A;p_B;p_G],0.05,'dep','yes');

disp('adjusted p vlaues:')
p_A_fdr = adj_p(1,:)
p_B_fdr = adj_p(2,:)
p_G_fdr = adj_p(3,:)


%% Analysis 2: AVERAGE ACROSS SEIZURES PER PATTERN
   
pattern_names={
'a: rhythmic alpha waves'                 
'b: rhythmic beta waves'
'l: low amp. fast activity'
'p: polyspikes'            
'r: repetitive spiking'
's: rhythmic sharp waves'	
't: rhythmic theta waves'		
};

onset_pattern=unique([FILEINFO.SzOnsetPattern]);

for pattern=1:length(onset_pattern)
    
    files_with_that_onset=find(strcmp([FILEINFO.SzOnsetPattern]',onset_pattern{pattern}));
    n_pattern(pattern)=length(files_with_that_onset)-sum(isnan(mean_A_other(files_with_that_onset,1)));
       
    mean_A_all_pattern(pattern,:)=mean(mean_A_all(files_with_that_onset,:),1,'omitnan');
    mean_B_all_pattern(pattern,:)=mean(mean_B_all(files_with_that_onset,:),1,'omitnan');
    mean_G_all_pattern(pattern,:)=mean(mean_G_all(files_with_that_onset,:),1,'omitnan');
    
    mean_A_soz_pattern(pattern,:)=mean(mean_A_soz(files_with_that_onset,:),1,'omitnan');
    mean_B_soz_pattern(pattern,:)=mean(mean_B_soz(files_with_that_onset,:),1,'omitnan');
    mean_G_soz_pattern(pattern,:)=mean(mean_G_soz(files_with_that_onset,:),1,'omitnan');
    
    mean_A_other_pattern(pattern,:)=mean(mean_A_other(files_with_that_onset,:),1,'omitnan');
    mean_B_other_pattern(pattern,:)=mean(mean_B_other(files_with_that_onset,:),1,'omitnan');
    mean_G_other_pattern(pattern,:)=mean(mean_G_other(files_with_that_onset,:),1,'omitnan');
   
    std_A_all_pattern(pattern,:)=std(mean_A_all(files_with_that_onset,:),1,'omitnan');
    std_B_all_pattern(pattern,:)=std(mean_B_all(files_with_that_onset,:),1,'omitnan');
    std_G_all_pattern(pattern,:)=std(mean_G_all(files_with_that_onset,:),1,'omitnan');
    
    std_A_soz_pattern(pattern,:)=std(mean_A_soz(files_with_that_onset,:),[],1,'omitnan');
    std_B_soz_pattern(pattern,:)=std(mean_B_soz(files_with_that_onset,:),[],1,'omitnan');
    std_G_soz_pattern(pattern,:)=std(mean_G_soz(files_with_that_onset,:),[],1,'omitnan');
    
    std_A_other_pattern(pattern,:)=std(mean_A_other(files_with_that_onset,:),[],1,'omitnan');
    std_B_other_pattern(pattern,:)=std(mean_B_other(files_with_that_onset,:),[],1,'omitnan');
    std_G_other_pattern(pattern,:)=std(mean_G_other(files_with_that_onset,:),[],1,'omitnan');
     
    
end
clear pattern

figure
subplot(3,1,1)
errorbar(mean_A_all_pattern',std_A_all_pattern'./sqrt(n_pattern),'linewidth',1.5); title('Seizure Onset Pattern'); 
ylabel('A (mean + SE)','FontWeight','bold'); xlim([0.5 4.5]); xticks(1:4); xticklabels({'intictal','preonset','onset','ictal'}); 
subplot(3,1,2)
errorbar(mean_B_all_pattern',std_B_all_pattern'./sqrt(n_pattern),'linewidth',1.5); 
ylabel('B (mean + SE)','FontWeight','bold'); xlim([0.5 4.5]); xticks(1:4); xticklabels({'intictal','preonset','onset','ictal'}); 
subplot(3,1,3)
errorbar(mean_G_all_pattern',std_G_all_pattern'./sqrt(n_pattern),'linewidth',1.5);
ylabel('G (mean + SE)','FontWeight','bold'); xlim([0.5 4.5]); xticks(1:4); xticklabels({'intictal','preonset','onset','ictal'}); 
legend(pattern_names,'Location','northwest');

%statistical analysis: classification performance, see later script 


%% Analysis 3: POST HOC "r" type vs "l" type seizures
pattern=[3,5];
n_pattern=2;
files1=find(strcmp([FILEINFO.SzOnsetPattern]',onset_pattern{pattern(1)})); %type "l"
files2=find(strcmp([FILEINFO.SzOnsetPattern]',onset_pattern{pattern(2)})); %type "r"
n_files1=sum(~isnan(nanmean(mean_A_soz(files1,:)-mean_A_other(files1,:),2)));
n_files2=sum(~isnan(nanmean(mean_A_soz(files2,:)-mean_A_other(files2,:),2)));

figure
subplot(3,2,1)
title('Excitation/inhibition changes in "l" and "r" type seizures'); 
b1=errorbar(mean_A_soz_pattern(pattern,:)',std_A_soz_pattern(pattern,:)'./sqrt(n_pattern),'linewidth',2); 
ylabel('A (mean + SE)','FontWeight','bold'); xlim([0.5 4.5]); xticks(1:4); xticklabels({'intictal','preonset','onset','ictal'}); y=ylim;
hold on
b2=plot(mean_A_other_pattern(pattern,:)','--','linewidth',1.5); 
b2(1).Color=b1(1).Color;
b2(2).Color=b1(2).Color;
legend({[pattern_names{pattern(1)},', SOZ'],[pattern_names{pattern(2)},', SOZ'],...
    [pattern_names{pattern(1)},': Other'],[pattern_names{pattern(2)},': Other']},'Location','northwest');
subplot(3,2,3)
title('Excitation/inhibition changes in "l" and "r" type seizures'); 
b1=errorbar(mean_B_soz_pattern(pattern,:)',std_B_soz_pattern(pattern,:)'./sqrt(n_pattern),'linewidth',2); 
ylabel('B (mean + SE)','FontWeight','bold'); xlim([0.5 4.5]); xticks(1:4); xticklabels({'intictal','preonset','onset','ictal'}); y=ylim;
hold on
b2=plot(mean_B_other_pattern(pattern,:)','--','linewidth',1.5); 
b2(1).Color=b1(1).Color;
b2(2).Color=b1(2).Color;
subplot(3,2,5)
title('Excitation/inhibition changes in "l" and "r" type seizures'); 
b1=errorbar(mean_G_soz_pattern(pattern,:)',std_G_soz_pattern(pattern,:)'./sqrt(n_pattern),'linewidth',2); 
ylabel('G (mean + SE)','FontWeight','bold'); xlim([0.5 4.5]); xticks(1:4); xticklabels({'intictal','preonset','onset','ictal'}); y=ylim;
hold on
b2=plot(mean_G_other_pattern(pattern,:)','--','linewidth',1.5); 
b2(1).Color=b1(1).Color;
b2(2).Color=b1(2).Color;

subplot(3,2,2)
b=bar(mean(mean_A_soz_pattern(pattern,:)-mean_A_other_pattern(pattern,:),2));
b.FaceColor=[0.6510 0.6510 0.6510]; hold on
er=errorbar(mean(mean_A_soz_pattern(pattern,:)-mean_A_other_pattern(pattern,:),2),...
    [nanstd(nanmean(mean_A_soz(files1,:)-mean_A_other(files1,:),2))./sqrt(n_files1);...
     nanstd(nanmean(mean_A_soz(files2,:)-mean_A_other(files2,:),2))./sqrt(n_files2)]);    
er.Color = [0 0 0]; er.LineStyle = 'none';  
xticklabels(pattern_names(pattern)); ylabel('A (mean + SE)','FontWeight','bold');
title('Excitation (SOZ - Other)'); 
subplot(3,2,4)
b=bar(mean(mean_B_soz_pattern(pattern,:)-mean_B_other_pattern(pattern,:),2));
b.FaceColor=[0.6510 0.6510 0.6510]; hold on
er=errorbar(mean(mean_B_soz_pattern(pattern,:)-mean_B_other_pattern(pattern,:),2),...
  [nanstd(nanmean(mean_B_soz(files1,:)-mean_B_other(files1,:),2))./sqrt(n_files1);...
     nanstd(nanmean(mean_B_soz(files2,:)-mean_B_other(files2,:),2))./sqrt(n_files2)]);  
er.Color = [0 0 0]; er.LineStyle = 'none';  
xticklabels(pattern_names(pattern)); ylabel('B (mean + SE)','FontWeight','bold');
title('Slow inhibition (SOZ - Other)'); 
subplot(3,2,6)
b=bar(mean(mean_G_soz_pattern(pattern,:)-mean_G_other_pattern(pattern,:),2));
b.FaceColor=[0.6510 0.6510 0.6510]; hold on
er=errorbar(mean(mean_G_soz_pattern(pattern,:)-mean_G_other_pattern(pattern,:),2),...
  [nanstd(nanmean(mean_G_soz(files1,:)-mean_G_other(files1,:),2))./sqrt(n_files1);...
     nanstd(nanmean(mean_G_soz(files2,:)-mean_G_other(files2,:),2))./sqrt(n_files2)]);  
er.Color = [0 0 0]; er.LineStyle = 'none';  
xticklabels(pattern_names(pattern)); ylabel('G (mean + SE)','FontWeight','bold');
title('Fast inhibition in (SOZ - Other)'); 


% test onset difference in ABG
disp('mean onset A for type "l" and "r":')
[nanmean(nanmean(mean_A_soz(files1,3))) nanmean(nanmean(mean_A_soz(files2,3)))]
[p_A_onset,~,stats_A_onset] = ranksum(nanmean(mean_A_soz(files1,3),2),...
    nanmean(mean_A_soz(files2,3),2),'tail','left')

disp('mean onset B for type "l" and "r":')
[nanmean(nanmean(mean_B_soz(files1,3))) nanmean(nanmean(mean_B_soz(files2,3)))]
[p_B_onset,~,stats_B_onset] = ranksum(nanmean(mean_B_soz(files1,3),2),...
    nanmean(mean_B_soz(files2,3),2),'tail', 'left')

disp('mean onset G for type "l" and "r":')
[nanmean(nanmean(mean_G_soz(files1,3))) nanmean(nanmean(mean_G_soz(files2,3)))]
[p_G_onset,~,stats_G_onset] = ranksum(nanmean(mean_G_soz(files1,3),2),...
    nanmean(mean_G_soz(files2,3),2),'tail', 'right')

[~, ~, ~, adj_p]=fdr_bh([p_A_onset;p_B_onset;p_G_onset],0.05,'dep','yes');
disp('adjusted p vlaues A, B, G for onset, then SOZloc:')
adj_p

% test mean difference in  SOZ-Other ABG (all time points) against zero
disp('mean SOZ-Other A for type "l" and "r":')
[mean(nanmean(mean_A_soz(files1,:)-mean_A_other(files1,:))) mean(nanmean(mean_A_soz(files2,:)-mean_A_other(files2,:)))]
[p_A_sozloc_r,~,stats_A_sozloc_r] = signrank(nanmean(mean_A_soz(files2,:)-mean_A_other(files2,:),2),0,'tail','right')
[p_A_sozloc_l,~,stats_A_sozloc_l] = signrank(nanmean(mean_A_soz(files1,:)-mean_A_other(files1,:),2),0,'tail','right')

disp('mean SOZ-Other B for type "l" and "r":')
[mean(nanmean(mean_B_soz(files1,:)-mean_B_other(files1,:))) mean(nanmean(mean_B_soz(files2,:)-mean_B_other(files2,:)))]
[p_B_sozloc_r,~,stats_B_sozloc_r] = signrank(nanmean(mean_B_soz(files2,:)-mean_B_other(files2,:),2),0,'tail','right')
[p_B_sozloc_l,~,stats_B_sozloc_l] = signrank(nanmean(mean_B_soz(files1,:)-mean_B_other(files1,:),2),0,'tail','right')

disp('mean SOZ-Other G for type "l" and "r":')
[mean(nanmean(mean_G_soz(files1,:)-mean_G_other(files1,:))) mean(nanmean(mean_G_soz(files2,:)-mean_G_other(files2,:)))]
[p_G_sozloc_r,~,stats_G_sozloc_r] = signrank(nanmean(mean_G_soz(files2,:)-mean_G_other(files2,:),2),0,'tail','right')
[p_G_sozloc_l,~,stats_G_sozloc_l] = signrank(nanmean(mean_G_soz(files1,:)-mean_G_other(files1,:),2),0,'tail','right')

[~, ~, ~, adj_p]=fdr_bh([p_A_sozloc_l;p_B_sozloc_l;p_G_sozloc_l;p_A_sozloc_r;p_B_sozloc_r;p_G_sozloc_r],0.05,'dep','yes');



%% Analysis 4: PATIENT VARIABLES 
run('make_patient_table'); %script that collects different patient varaibles from the database

change_A_all_pat=(mean_A_all_pat(:,end)-mean_A_all_pat(:,1))./((mean_A_all_pat(:,end) + mean_A_all_pat(:,1)));%./mean_A_all_pat(:,1);
change_B_all_pat=(mean_B_all_pat(:,end)-mean_B_all_pat(:,1))./((mean_B_all_pat(:,end) + mean_B_all_pat(:,1)));%./mean_B_all_pat(:,1);
change_G_all_pat=(mean_G_all_pat(:,end)-mean_G_all_pat(:,1))./((mean_G_all_pat(:,end) + mean_G_all_pat(:,1)));%./mean_G_all_pat(:,1);

change_A_soz_pat=(mean_A_soz_pat(:,end)-mean_A_soz_pat(:,1))./((mean_A_soz_pat(:,end) + mean_A_soz_pat(:,1)));%./mean_A_soz_pat(:,1);
change_B_soz_pat=(mean_B_soz_pat(:,end)-mean_B_soz_pat(:,1))./((mean_B_soz_pat(:,end) + mean_B_soz_pat(:,1)));%./mean_B_soz_pat(:,1);
change_G_soz_pat=(mean_G_soz_pat(:,end)-mean_G_soz_pat(:,1))./((mean_G_soz_pat(:,end) + mean_G_soz_pat(:,1)));%./mean_G_soz_pat(:,1);

change_A_other_pat=(mean_A_other_pat(:,end)-mean_A_other_pat(:,1))./((mean_A_other_pat(:,end) + mean_A_other_pat(:,1)));%./mean_A_other_pat(:,1);
change_B_other_pat=(mean_B_other_pat(:,end)-mean_B_other_pat(:,1))./((mean_B_other_pat(:,end) + mean_B_other_pat(:,1)));%./mean_B_other_pat(:,1);
change_G_other_pat=(mean_G_other_pat(:,end)-mean_G_other_pat(:,1))./((mean_G_other_pat(:,end) + mean_G_other_pat(:,1)));%./mean_G_other_pat(:,1);


% localisation
left=contains(Localisation,'L');
right=contains(Localisation,'R');


figure
disp('Localisation left vs right');
subplot(1,3,1)
p=ranksum(change_A_all_pat(left),change_A_all_pat(right))
boxchart(change_A_all_pat,'GroupByColor',Localisation); title('Change in A (All)'); xticklabels([]);
subplot(1,3,2)
p=ranksum(change_B_all_pat(left),change_B_all_pat(right))
boxchart(change_B_all_pat,'GroupByColor',Localisation); title('Change in B (All)'); xticklabels([]);
subplot(1,3,3)
p=ranksum(change_G_all_pat(left),change_G_all_pat(right))
boxchart(change_G_all_pat,'GroupByColor',Localisation); title('Change in G (All)'); xticklabels([]);
l=legend; l.Title.String='localisation';
set(gcf, 'Position', [347 538 763 235]);


% HS
disp('HS vs no HS');
figure
subplot(1,3,1)
p=ranksum(change_A_all_pat(HS),change_A_all_pat(~HS))
boxchart(change_A_all_pat,'GroupByColor',HS); title('Change in A (All)'); xticklabels([]);
subplot(1,3,2)
p=ranksum(change_B_all_pat(HS),change_B_all_pat(~HS))
boxchart(change_B_all_pat,'GroupByColor',HS); title('Change in B (All)'); xticklabels([]);
subplot(1,3,3)
p=ranksum(change_G_all_pat(HS),change_G_all_pat(~HS))
boxchart(change_G_all_pat,'GroupByColor',HS); title('Change in G (All)'); xticklabels([]);
l=legend; l.Title.String='Hippocampal Sclerosis';
set(gcf, 'Position', [347 538 763 235]);


% surgery outcome
disp('outcome II change (soz - other) vs outcome Ia change (soz - other)');
figure
subplot(1,3,1)
p=ranksum(change_A_soz_pat(contains(outcome,'II'))-change_A_other_pat(contains(outcome,'II')),...
    change_A_soz_pat(contains(outcome,'Ia'))-change_A_other_pat(contains(outcome,'Ia')), 'tail','left')
boxchart(change_A_soz_pat-change_A_other_pat,'GroupByColor',outcome); title('Change in A (SOZ - Other)'); xticklabels([]);
subplot(1,3,2)
p=ranksum(change_B_soz_pat(contains(outcome,'II'))-change_B_other_pat(contains(outcome,'II')),...
    change_B_soz_pat(contains(outcome,'Ia'))-change_B_other_pat(contains(outcome,'Ia')), 'tail','left')
boxchart(change_B_soz_pat-change_B_other_pat,'GroupByColor',outcome); title('Change in B (SOZ - Other)'); xticklabels([]);
subplot(1,3,3)
p=ranksum(change_G_soz_pat(contains(outcome,'II'))-change_G_other_pat(contains(outcome,'II')),...
    change_G_soz_pat(contains(outcome,'Ia'))-change_G_other_pat(contains(outcome,'Ia')), 'tail','left')
boxchart(change_G_soz_pat-change_G_other_pat,'GroupByColor',outcome); title('Change in G (SOZ - Other)'); xticklabels([]);
l=legend; l.Title.String='outcome';
set(gcf, 'Position', [347 538 763 235]);


