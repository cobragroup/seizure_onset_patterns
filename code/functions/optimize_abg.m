function [minerrpars, minerr, minerrfeatures]= optimize_abg(featuresF0,featuresPARSWEEP,pars)

if size(featuresF0,2)~=size(featuresPARSWEEP,2)
    error('different amount of features in featuresF0 than in featuresPARSWEEP!');
end

errABG=zeros(size(featuresF0,1),size(featuresPARSWEEP,1),size(featuresF0,2));
for segm=1:size(featuresF0,1) 
    errABG(segm,:,:)=abs(featuresPARSWEEP-featuresF0(segm,:));
end
errsumABG=sum(errABG,3); 


%find parset with closest errsum
[minerr,minerrIND]=min(errsumABG,[],2,'omitnan'); 
minerrpars=pars(minerrIND,:);
minerrfeatures=featuresPARSWEEP(minerrIND,:);

end