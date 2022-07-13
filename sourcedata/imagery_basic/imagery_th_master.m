function imagery_th_master
% function mot_th_master
% this function calculates basic move-rest differences
% as in 
% note - this is rewritten code, and differs slightly from the original
% analysis in the following ways:
% 1) shift by 500ms, to account for behavioral delay - this may be commented out
% 2) compare vs paired rest following (bc of rebound)
% 3) variable length of blocks depending on if cues were 2 or 3 secs.
% kjm 2015

%% set defaults
    srate=1000;
    warning('off','signal:psd:PSDisObsolete'); %annoying
    BF_correct='n'; % to Bonferroni correct p-values (should be 'y' if so)

%%
subjects={'bp',...
%     'fp',...
%     'hh',...
%     'jc',...
%     'jm',...
%     'rh',...
    'rr'...
    };

%%
for k=1:length(subjects)
for reim={'mot','im'}    
%% load data and electrodes
    disp(['Subject ' subjects{k}])
    load(['data/' subjects{k} '_' reim{1} '_t_h'],'data','stim');
    load(['locs/' subjects{k} '_electrodes'],'electrodes'); 

%% re-reference data
    data=car(data);

%% reshape stimulus variable 
% note stim==11 is tongue movement; stim==12 is hand movement
% assign hand and tongue rest block conditions

    % find length of time for each, error if there is a cut.
    a=[0; (stim(2:end)-stim(1:(end-1)))]; a1=find(a==-12); a2=find(a==-11); a3=find(a==11); a4=find(a==12);
    if length(a1)==length(a4), blocksize=mean(abs(a1-a4));
    elseif length(a2)==length(a3), blocksize=mean(abs(a2-a3));
    else error('can''t determine blocksize b/c task cut between between task block','can''t determine blocksize b/c task cut between between task block')
    end
    if floor(blocksize/1000)~=(blocksize/1000), error('can''t determine blocksize b/c task cut between between task block','can''t determine blocksize b/c task cut between between task block'), end
    clear a*

    % find rest conditions specific to behavior that preceded them
    stim(find(and(...
        [zeros(blocksize,1); stim(1:(length(stim)-blocksize)) ]==12, ...
        stim==0 ...
        )))=120; %hand rest block condition

    stim(find(and(...
        [zeros(blocksize,1); stim(1:(length(stim)-blocksize)) ]==11, ...
        stim==0 ...
        )))=110; %tongue rest block condition
    
%     % shift stim by half-second to account for behavioral lag
%     stim=[zeros(500,1); stim(1:(end-500),1)];

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
for m=1:num_chans
    % hand (tr_sc==12)        
        r_hand_HFB(m)=rsa(HFB_trials(m,find(tr_sc==12)),HFB_trials(m,find(tr_sc==120)));
        r_hand_LFB(m)= rsa(LFB_trials(m, find(tr_sc==12)), LFB_trials(m,find(tr_sc==120)));        
        [temp, p_hand_HFB(m)]=ttest2(HFB_trials(m,find(tr_sc==12)),HFB_trials(m,find(tr_sc==120)));
        [temp, p_hand_LFB(m)]= ttest2(LFB_trials(m,find(tr_sc==12)), LFB_trials(m, find(tr_sc==120))); 
        %
        % Bonferroni correct
        if BF_correct=='y'
        p_hand_HFB(m)=p_hand_HFB(m)*num_chans;
        p_hand_LFB(m)=p_hand_LFB(m)*num_chans;
        end

    % tongue (tr_sc==11)        
        r_tongue_HFB(m)=rsa(HFB_trials(m,find(tr_sc==11)),HFB_trials(m,find(tr_sc==110)));
        r_tongue_LFB(m)= rsa(LFB_trials(m, find(tr_sc==11)), LFB_trials(m,find(tr_sc==110)));        
        [temp, p_tongue_HFB(m)]=ttest2(HFB_trials(m,find(tr_sc==11)),HFB_trials(m,find(tr_sc==110)));
        [temp, p_tongue_LFB(m)]= ttest2(LFB_trials(m,find(tr_sc==11)), LFB_trials(m, find(tr_sc==110)));  
        %
        % Bonferroni correct
        if BF_correct=='y'
        p_tongue_HFB(m)=p_tongue_HFB(m)*num_chans;
        p_tongue_LFB(m)=p_tongue_LFB(m)*num_chans;
        end
end

%% average spectra for different conditions
    mean_PSD_handmove=squeeze(mean(all_PSD(:,:,find(tr_sc==12)),3));
    mean_PSD_tonguemove=squeeze(mean(all_PSD(:,:,find(tr_sc==11)),3));
    mean_PSD_handrest=squeeze(mean(all_PSD(:,:,find(tr_sc==120)),3));
    mean_PSD_tonguerest=squeeze(mean(all_PSD(:,:,find(tr_sc==110)),3));
    
%% make figures
    for q={'hand','tongue'}
        for bnd={'HFB','LFB'}
            figure
            wts=eval(['r_' q{1} '_' bnd{1} '.*' '(p_' q{1} '_' bnd{1} '<.05)']); % weights to plot
            tail_gauss_plot_redux(electrodes,wts);
            hold on, plot3(electrodes(:,1)*1.01,electrodes(:,2),electrodes(:,3),'.','MarkerSize',4,'Color',[.99 .99 .99])
            title([subjects{k} ', ' q{1} '-' bnd{1} ', max r^2=' num2str(max(abs(wts)))])
            exportfig(gcf, ['figs/' subjects{k} '_' reim{1} '_GMap_' q{1} '_' bnd{1}], 'format', 'png', 'renderer', 'opengl', 'Color', 'cmyk', 'Resolution', 600, 'Width', 1.5, 'Height', 1);
            close    
        end
    end

%%    
    save(['data/' subjects{k} '_' reim{1} '_th_analyzed'],'*PSD*','*FB*','tr_sc','blocksize')
    disp('-------------------------------------------')
end    
end
%%
                