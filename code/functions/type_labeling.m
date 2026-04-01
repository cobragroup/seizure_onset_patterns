function [centroids_labeled,unique_type_match,type_merged]=type_labeling(centroids,real_type,cluster_type)

%for each centroid, find most commonly assigned real type
for t=1:max(cluster_type)
    type_label(t)=mode(real_type(cluster_type==t));
end

%check whether several centrois received the same label
[found_types, ~, ia] = unique(type_label, 'stable');  % Stable keeps it in the same order
bincounts = accumarray(ia, 1);
unique_type_match=sort(found_types(bincounts==1));
unique_type_match=unique_type_match(~isnan(unique_type_match));
type_merged=found_types(bincounts>1);
check=find(bincounts>1);
if ~isempty(check)
    for l=1:length(check)
        same_label=find(type_label==found_types(check(l)));
        
        %if same label, omit the centroid with lower type overlap
        for i=1:length(same_label)
            [~,F(i)]=mode(real_type(cluster_type==same_label(i)));
        end
        [~,best]=max(F);
        remove_centroid=same_label(same_label~=same_label(best));
        type_label(remove_centroid)=NaN; %#ok<FNDSB>
    end
end
    
%if a centroid was never chosen as closest or it lost same label comparison, omit it!
centroids=centroids(~isnan(type_label),:);
type_label=type_label(~isnan(type_label));

% reorder centroids so that position 1 will be centroid type 1 and so on
centroids_labeled=NaN(max(real_type),size(centroids,2));
centroids_labeled(type_label,:)=centroids; %