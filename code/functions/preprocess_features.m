function [X,N_features,featurenames,featureind,coeff]=...
 preprocess_features(X,pars)

N_features=length(X.datapars.featurenames);
featureind=1:X.datapars.N_features;
featurenames=X.datapars.featurenames;

if pars.oldonly
    featurenames=X.datapars.featurenames(X.datapars.oldfeatures);
    N_features=length(featurenames);
    featureind=X.datapars.oldfeatures;
    X.features=X.features(:,X.datapars.oldfeatures);
end

if pars.pcafeat==true
    [coeff,score,~,~,explained,~] = pca(X.features, 'NumComponents',pars.n_comp);

    figure;
    set(gcf,'Position',[777.0000   98.6000  727.2000  624.0000]);
    
    imagesc(coeff);
    set(gca,'YTick',0.5:1:length(featurenames)-0.5);
    yticklabels(featurenames);
    xlabel('PCA component');
    set(gca,'XTick',1:pars.n_comp,'XTick',1:N_features);
    set(gca,'XTick',1:pars.n_comp);
    
    c=colorbar;
    c.Title.String='PCA coeff';
    title(['explained var = ' num2str(sum(explained(1:pars.n_comp)))]);
    
    disp('explained variance');
    disp(explained');
    
    X.features=score;
    pcnames={'PC1','PC2','PC3','PC4','PC5','PC6'};
    featurenames=pcnames(1:pars.n_comp);
    X.datapars.featurenames=pcnames(1:pars.n_comp);
    N_features=pars.n_comp;
    clear pcnames
else
    coeff=1;
end
