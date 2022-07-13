function MIF_master(subject)



%%
switch subject
    case 'al', 
        loc_file='al_locs_40.mat'; % electrode locations file
        mot_file='al_mot_ih_shrug.mat'; % motor movement file
        im_file='al_im_ih_shrug.mat'; % imagery file
        fb_file='al_fb_shrug.mat'; % target feedback file
        fb_cond=19; %shrug 
        fb_elec=36; % electrode used for feedback
    case 'fp_h', % subject fp, hand feedback
        loc_file='fp_locs_40.mat'; % electrode locations file
        mot_file='fp_mot_t_h_40.mat'; % motor movement file
        im_file='fp_im_t_h_40.mat'; % imagery file
        fb_file='fp_fbLR_hand.mat'; % left-right targets by hand feedback file
        fb_cond=12; %hand
        fb_elec=24; % electrode used for feedback
        sc_fbt=sc_fb; sc_fb(find(sc_fbt==1))=2; sc_fb(find(sc_fbt==2))=1; clear fbt %flip targets for appropriate sign on r2
    case 'fp_t', % subject fp, tongue feedback
        loc_file='fp_locs_40.mat'; % electrode locations file
        mot_file='fp_mot_t_h_40.mat'; % motor movement file
        im_file='fp_im_t_h_40.mat'; % imagery file
        fb_file='fp_fbUD_tongue.mat'; % up-down targets by tongue feedback file
        fb_cond=11; %tongue
        fb_elec=22; % electrode used for feedback
    case 'hh', 
        loc_file='hh_locs_48.mat'; % electrode locations file
        mot_file='hh_mot_t.mat'; % motor movement file
        im_file='hh_im_t.mat'; % imagery file
        fb_file='hh_fb_tongue.mat'; % target feedback file
        fb_cond=11; % tongue
        fb_elec=44; % electrode used for feedback
    case 'jc', 
        loc_file='jc_locs_48.mat'; % electrode locations file
        mot_file='jc_mot_l_mov.mat'; % motor movement file
        im_file='jc_mot_l_mov.mat'; % NOTE: NO IMAGERY FILE FOR THIS PATIENT -- motor movement file put in place to allow easy cycling through without error
        fb_file='jc_fb_mov.mat'; % target feedback file
        fb_cond=15; % saying the word "move"
        fb_elec=13; % electrode used for feedback
end

%% do basic analyses for each task type
    disp('analyzing motor data')
    [r_mot_HFB, r_mot_LFB,p_mot_HFB,p_mot_LFB]=mot_analysis(mot_file,fb_cond); % movement data analysis
    disp('analyzing imagery data')
    [r_im_HFB, r_im_LFB,p_im_HFB,p_im_LFB]=mot_analysis(im_file,fb_cond); % imagery data analysis
    disp('analyzing feedback data')
    [r_fb_HFB, r_fb_LFB,p_fb_HFB,p_fb_LFB]=fb_analysis(fb_file); % feedback data analysis
    %
    save(['data/' subject '/' subject '_MIF_rp_stats'],'r_*','p_*')


    
%% calculate overlap by resampling
% cycle through to calculate 
bnd={'HFB','LFB'};
q={'mot','im','fb'};
    for k=1:2
        for q1=1:2
        for q2=(q1+1):3
            d1=eval(['r_' q{q1} '_' bnd{k} '.*' '(p_' q{q1} '_' bnd{k} '<.05)']); % 
            d2=eval(['r_' q{q2} '_' bnd{k} '.*' '(p_' q{q2} '_' bnd{k} '<.05)']); % 
            [OLM(q1,q2,k), p_val(q1,q2,k), rs_kurt(q1,q2,k)]=spat_reshuffle(d1,d2, 'y'); figure
        end
        end
    end
    


%% make brain activity figures

%load locs
load(['locs/' loc_file])

% cycle through to make brain figures
    for q={'mot','im','fb'}
        for bnd={'HFB','LFB'}
            figure
            wts=eval(['r_' q{1} '_' bnd{1} '.*' '(p_' q{1} '_' bnd{1} '<.05)']); % weights to plot
            tail_gauss_plot_redux(electrodes,wts);
            hold on, plot3(electrodes(:,1)*1.01,electrodes(:,2),electrodes(:,3),'.','MarkerSize',4,'Color',[.99 .99 .99])
            title([subject ', ' q{1} '-' bnd{1} ', max r^2=' num2str(max(abs(wts)))])
            if subject=='al', loc_view(270,30), end
            exportfig(gcf, ['figs/' subject '_GMap_' q{1} '_' bnd{1}], 'format', 'png', 'renderer', 'opengl', 'Color', 'cmyk', 'Resolution', 600, 'Width', 1.5, 'Height', 1);
            close    
        end
    end



















