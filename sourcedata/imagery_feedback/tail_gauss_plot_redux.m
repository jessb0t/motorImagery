function [electrodes]=tail_gauss_plot_redux(electrodes,weights)
% function [electrodes]=tail_gauss_plot_redux(electrodes,weights)
% projects electrode locations onto their cortical spots using a gaussian kernel
% originally from "Location on Cortex" package, (Miller, KJ, et al NeuroImage, 2007) 
% altered for ECoG library by kjm 1/2016

%% load brain  
    load('plotting_fxns/halfbrains.mat')
    load('plotting_fxns/dg_colormap')
    if mean(electrodes(:,1))<0 % left 
        cortex=leftbrain;
    else % right
        cortex=rightbrain;
    end
    brain=cortex.vert;

%% fill out gaussian 
% gaussian "cortical" spreading parameter - in mm, so if set at 10, its 1 cm
%- distance between adjacent electrodes
gsp=50;

c=zeros(length(cortex(:,1)),1);

for k=1:length(electrodes(:,1)) %cycle through, adding gaussian shading
    b_z=abs(brain(:,3)-electrodes(k,3));
    b_y=abs(brain(:,2)-electrodes(k,2));
    b_x=abs(brain(:,1)-electrodes(k,1));
    d=weights(k)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian 
    c=c+d';
end
% c=(c/max(c));
a=tripatch(cortex, 'nofigure', c');
shading interp;
a=get(gca);
d=a.CLim;
set(gca,'CLim',[-max(abs(d)) max(abs(d))])
l=light;
colormap(cm)
lighting gouraud;
material dull;
axis off

% change view and lighting 
if mean(electrodes(:,1))<0 % left
view(270, 0); set(l,'Position',[-1 0 1]) 
else % right
view(90, 0); set(l,'Position',[1 0 1]) 
end



