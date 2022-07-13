% function call_MIF_master


%%
subjects={...
    'al',... 
    'fp_h',... % subject fp, hand feedback    
    'fp_t',... % subject fp, tongue feedback
    'hh',... 
    'jc'... 
    };

%%
for k=1:length(subjects)
MIF_master(subjects{k})
end




