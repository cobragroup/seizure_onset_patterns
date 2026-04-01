function TYPES = get_types(RESULT,first_ictal_ind,last_ictal_ind)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define TYPE onset indices
interictal=1;
preonset=round(first_ictal_ind/2);
onset=first_ictal_ind;
ictal=first_ictal_ind+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for channel_ind=1:size(RESULT,2)
    if isempty(RESULT(channel_ind).minerrparsABG)
        TYPES(channel_ind).filename= {'empty channel'};
        continue
    end
    
    %extract type-specific simulation results
    TYPES(channel_ind).filename={RESULT(channel_ind).fileID};
    TYPES(channel_ind).channel_ind=channel_ind;
    
    TYPES(channel_ind).types_intictal=RESULT(channel_ind).minerrtype(interictal:preonset-1)';
    TYPES(channel_ind).types_preonset=RESULT(channel_ind).minerrtype(preonset:onset-1)';
    TYPES(channel_ind).types_onset=RESULT(channel_ind).minerrtype(onset:ictal-1)';
    TYPES(channel_ind).types_ictal=RESULT(channel_ind).minerrtype(ictal:last_ictal_ind)';
    
    TYPES(channel_ind).a_intictal=RESULT(channel_ind).minerrparsABG(interictal:preonset-1,1)';
    TYPES(channel_ind).a_preonset=RESULT(channel_ind).minerrparsABG(preonset:onset-1,1)';
    TYPES(channel_ind).a_onset=RESULT(channel_ind).minerrparsABG(onset:ictal-1,1)';
    TYPES(channel_ind).a_ictal=RESULT(channel_ind).minerrparsABG(ictal:last_ictal_ind,1)';
    
    TYPES(channel_ind).b_intictal=RESULT(channel_ind).minerrparsABG(interictal:preonset-1,2)';
    TYPES(channel_ind).b_preonset=RESULT(channel_ind).minerrparsABG(preonset:onset-1,2)';
    TYPES(channel_ind).b_onset=RESULT(channel_ind).minerrparsABG(onset:ictal-1,2)';
    TYPES(channel_ind).b_ictal=RESULT(channel_ind).minerrparsABG(ictal:last_ictal_ind,2)';
    
    TYPES(channel_ind).g_intictal=RESULT(channel_ind).minerrparsABG(interictal:preonset-1,3)';
    TYPES(channel_ind).g_preonset=RESULT(channel_ind).minerrparsABG(preonset:onset-1,3)';
    TYPES(channel_ind).g_onset=RESULT(channel_ind).minerrparsABG(onset:ictal-1,3)';
    TYPES(channel_ind).g_ictal=RESULT(channel_ind).minerrparsABG(ictal:last_ictal_ind,3)';
    
end

end


