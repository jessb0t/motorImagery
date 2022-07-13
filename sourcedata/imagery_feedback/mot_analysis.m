function [r_HFB, r_LFB,p_HFB,p_LFB]=mot_analysis(fname,fb_cond)
% function mot_analysis
% this function calculates basic move-rest differences
% as in 
% note - this is rewritten code, and differs slightly from the original
% analysis in the following ways:
% 1) shift by 500ms - this may be commented out
% 2) compare vs paired rest following (bc of rebound)
% 3) variable length of blocks depending on if cues were 2 or 3 secs.
% kjm 2015

%% set defaults
    srate=1000;
    warning('off','signal:psd:PSDisObsolete'); %annoying
    BF_correct='n'; % to Bonferroni correct p-values (should be 'y' if so)

%% load data
    load(['data/' fname(1:2) '/' fname],'data','StimulusCode'); stim=double(StimulusCode);

%% re-reference data
    data=car(data);

%% reframe stimulus variable 
    stimtypes=unique(stim); stimtypes(stimtypes==0)=[];
    % find length of time for each, error if there is a cut.
    a=[0; (stim(2:end)-stim(1:(end-1)))]; blocksize=[]; 
    for k=1:length(stimtypes) 
        a1=find(a==-stimtypes(k)); a2=find(a==stimtypes(k)); 
            if length(a1)==length(a2), blocksize=mean(abs(a1-a2));end    
    end    
    if or(isempty(blocksize), floor(blocksize/1000)~=(blocksize/1000))
    error('can''t determine blocksize b/c task cut between between task block','can''t determine blocksize b/c task cut between between task block')
    end    

    % find rest conditions specific to behavior that preceded them
    stimtypes=unique(stim); stimtypes(stimtypes==0)=[];
    for k=1:length(stimtypes)    
        stim(find(and(...
            [zeros(blocksize,1); stim(1:(length(stim)-blocksize)) ]==stimtypes(k), ...
            stim==0 ...
            )))=10*stimtypes(k); %rest block condition following stimtype
    end
    

%% catalogue trials
    trialnr=0*stim; tr_sc=0; %initialize
    trtemp=1; trialnr(1)=trtemp; %start variable
    for n=2:length(stim)
        if stim(n)~=stim(n-1),
            trtemp=trtemp+1; tr_sc=[tr_sc stim(n)];
        end
        trialnr(n)=trtemp;
    end
    clear n trtemp
    trials=unique(trialnr);

%% calculate spectra for all trials
    num_chans=length(data(1,:));
    mean_PSD=zeros(200,num_chans); % to calculate the mean spectra over the entire period
    all_PSD=[];
    for cur_trial=1:max(trials), 
        % isolate relevant data 
        curr_data=data(find((trialnr == cur_trial)),:);
        % keep only .5 to block end of data, for each block
        if length(curr_data(:,1))>=blocksize
            curr_data=curr_data((floor(srate/2)+1):blocksize,:);
        end
        block_PSD=[];
        for p=1:num_chans
    %         temp_PSD=psd(curr_data(:,p),srate,srate,floor(srate*.25),floor(srate/10),'mean');
            temp_PSD=psd(curr_data(:,p),srate,srate,floor(srate*.25),floor(srate/10));
            block_PSD=[block_PSD temp_PSD];
        end
        block_PSD=block_PSD(1:200,:); % downsample - we only want to keep spectra up to 200 hz, not all of the others
        mean_PSD=mean_PSD+block_PSD;
        all_PSD=cat(3,all_PSD,block_PSD);
    end % session
    mean_PSD=mean_PSD./(max(trials));
    freq_bins=1:200;

    clear cur* *temp* p
%% Isolate HFB and LFB,
    all_PSD_temp=all_PSD./repmat(mean_PSD,[1 1 size(all_PSD,3)]);
    LFB_trials= squeeze(sum(all_PSD_temp( 8:32, :,:),1));
    HFB_trials=squeeze(sum(all_PSD_temp(76:100,:,:),1));
    clear all_PSD_temp

%% comparisons - signed r-squared and unpaired t-test 
% (note: this is done for the feedback condition only. change / expand code here if you would like it for all types of stimuli)
for m=1:num_chans      
        r_HFB(m)=rsa(HFB_trials(m,find(tr_sc==fb_cond)),HFB_trials(m,find(tr_sc==fb_cond*10)));
        r_LFB(m)= rsa(LFB_trials(m, find(tr_sc==fb_cond)), LFB_trials(m,find(tr_sc==fb_cond*10)));        
        [temp, p_HFB(m)]=ttest2(HFB_trials(m,find(tr_sc==fb_cond)),HFB_trials(m,find(tr_sc==fb_cond*10)));
        [temp, p_LFB(m)]= ttest2(LFB_trials(m,find(tr_sc==fb_cond)), LFB_trials(m, find(tr_sc==fb_cond*10))); 
        %
        % Bonferroni correct
        if BF_correct=='y'
        p_HFB(m)=p_HFB(m)*num_chans;
        p_LFB(m)=p_LFB(m)*num_chans;
        end
end

%% average spectra for different conditions
    mean_PSDmove=squeeze(mean(all_PSD(:,:,find(tr_sc==fb_cond)),3));
    mean_PSDrest=squeeze(mean(all_PSD(:,:,find(tr_sc==fb_cond*10)),3));
    

%% save analysis results
    save(['data/' fname(1:2) '/an_' fname],'*PSD*','*FB*','tr_sc','blocksize')
    disp('-------------------------------------------')

%%
                