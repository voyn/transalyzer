function [smoothed_fPA0,smoothed_SPA0] = smoothing(fPA0, PxxPA0)

% span = 5; % Size of the averaging window
% window = ones(span,1)/span; 
% smoothed_PxxPA0 = convn(PxxPA0,window,'same');
% smoothed_PxxPC0 = convn(PxxPC0,window,'same');
% smoothed_PxxPS0 = convn(PxxPS0,window,'same');
% 
% smoothed_PxxPA0 = convn(smoothed_PxxPA0,window,'same');
% smoothed_PxxPC0 = convn(smoothed_PxxPC0,window,'same');
% smoothed_PxxPS0 = convn(smoothed_PxxPS0,window,'same');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency averaging
% http://www.mathworks.nl/matlabcentral/newsreader/view_thread/133024
% http://www.mathworks.nl/matlabcentral/newsreader/view_thread/131742


% fPA0 = 1:1:100;
% PxxPA0 = abs(50*rand(size(fPA0))-(fPA0*3)+100*ones(size(fPA0)));
% figure()
% plot(fPA0,PxxPA0,'-b')
% hold on


nave=4; % You will need to play with this.
ns=length(PxxPA0); % number of smaples
nbins=fix(ns/nave); % round to zero, size of window
smoothed_SPA0=reshape(PxxPA0(1:nbins*nave),nave,nbins);
smoothed_fPA0=reshape(fPA0(1:nbins*nave),nave,nbins);
smoothed_SPA0=mean(smoothed_SPA0);
smoothed_fPA0=mean(smoothed_fPA0);

% smoothed_SPA0_m=mean(smoothed_SPA0);
% smoothed_fPA0_m=mean(smoothed_fPA0);

% plot(smoothed_fPA0,smoothed_SPA0,'-r')

