function  [em_thr,em_thr_behavior,P,Mean,Std,pdf_x,xx,pdf_xx,cdf_xx,P_0] = em_1dim(x,nr_groups,fig_nr,P_0,max_nr_it);

%  Syntax:  
%
%  [mixture_parameters,mixture_density] = em_1dim(x,nr_groups,fig_nr,P_0,max_nr_it);
%
%  EM algo for 1-dimensional Gaussian Mixture Model
%
%  INPUT:
% -------
%
% x =  data  with  [nr_pts,nr_dim] = size(x)
%
% nr_groups = number of Gaussians
%
% P_0  =  (OPTIONAL) initial soft assignment of data points to different clusters
%          (matrix of size = (nr_pts,k)
%          If omitted, we use either 
%              * k-means to divide group in approx equal parts
%              * or home-made function (see end-of-file) to divide RANGE in
%                equal parts (this should be beneficial for distributions
%                that have long tails).
%
%  OUTPUT:
%  -------
%
%    mixture_parameters = structure that contains the Gaussian Mixture parameters:
%       mixture_parameters.props  =  priors for each gaussian (size = [nr_gauss,1]);
%       mixture_parameters.means  =  means for each gaussian (size = [nr_gauss,1]);
%       mixture_parameters.stds   =  std for each gaussian (size = [nr_gauss,1]);
%
%    mixture_density = structure that contains various density components for plotting:
%
%       mixture_density.xplot :  x-values for plotting
%       mixture_density.pdf :    pdf values for plotting;
%       mixture_density.cdf :    cdf values for plotting;
%
%
%
% History:
%  THIS IS SPECIALIZATION OF THE GENERAL  "em_ndim.m" algorithm FOR 1-DIM
%  DATA   (CREATED:   26 Nov 2003, EP)
% 
%  8 feb 04:  eliminated potential "divide by zero" problem when
%             renormalizing P  (EP)
%  11 may 04:  streamlined output (more extensive original version:  em_1dim_extra.m)
%

em_thr = 0;
em_thr_behavior = 0;

plot_ind = 0;

% Define default values wherever needed.

if size(x,1) < size(x,2), x = x'; end   % turn into column if need be

nr_pts = length(x);

if nargin < 3,  fig_nr = 0;  end      %  default: no plotting 

if nargin < 4  % construct initial assignment -----(-5)
   
   if 0   %nr_groups > 1  % use K-means to divide dataset in equal GROUPS
     P_0 = zeros(nr_pts,nr_groups);  
     disp('Creating initial assignment based on k-means clustering')
     [cluster_thr,cluster_label] = k_means_cluster_1dim(x,nr_groups);  
     for j = 1:nr_groups
         zz = find(cluster_label==j);
         P_0(zz,j) = 1;
     end
   elseif nr_groups >1  % divide RANGE in equal parts
       
     P_0 = initialize_em_range_based(x,nr_groups);    % See end of this file
       
   else % only one group: trivial assignment
       P_0 = ones(nr_pts,1);   
   end
    
end %----------------------------------------(-5)

if nargin < 5
  max_nr_it = 1000;
end


P = P_0;
Mean = zeros(nr_groups,1);
Std = ones(nr_groups,1);

green_light = 1;

nr_it = 0;

while (green_light == 1) & (nr_it < max_nr_it)  %  begin---------------(100)

   nr_it  = nr_it + 1;

   P_new = zeros(size(P));

   for k = 1 : nr_groups  %--------------------------(50)

       PP = P(:,k);

       D = x.* PP;    %  Data weighted with P-matrix
           
       if sum(P(:,k)) ~=0   % there are datapoints assigned to this group
           
           mean_grp = sum(D)/sum(PP);
           var_grp = sum(((x - mean_grp).^2).*PP)/sum(PP);    %  should this be sum(PP)-1 ??
	       std_grp = sqrt(var_grp);
        else
           mean_grp = 0;
           std_grp = 1;
        end
        
       F =  normpdf(x,mean_grp,std_grp);

       Mean(k,:) = mean_grp;
       Std(k,:) = std_grp(:)';

       P_new(:,k) = F;

    end  %--------------------------(50)

    P_old = P;
    P = P_new;

    %  Renormalize
    
    P_sum = sum(P,2);  PP_sum = P_sum *ones(1,nr_groups);
    
    % Precautions to avoid "divide by zero"
    
    u_zero = find(P_sum < 10^(-200)); %
    
    if ~isempty(u_zero)
        % create uniform distribution
        Q = zeros(nr_pts,nr_groups);  Q(u_zero,:) = 1/nr_groups;
        N  = ones(nr_pts,nr_groups);  N(u_zero,:)=0;
        PP_sum(u_zero,:) = 1;
        P = (P./PP_sum).*N  + Q;
    else
          P = P./(sum(P,2)*ones(1,nr_groups));   % 
    end 

  if  sum(sum(abs(P-P_old))) < 0.001*nr_pts

   %    disp(['Convergence occurred after '  int2str(nr_it)  ' iterations'])
       green_light = 0;
   end

end   % end -------------------------------------------(100)

% Compute stuff to SHOW DENSITY if requested
%===========================================

r = range(x);   
dx = r/250;

xx = (min(x)-0.05*r:dx:max(x)+0.05*r);

gx = zeros(nr_groups,length(x));
g = zeros(nr_groups,length(xx));
G = zeros(nr_groups,length(xx));
s = sum(P,1);
s = s/sum(s);

for i = 1:nr_groups

   g(i,:) = s(i)*(normpdf(xx,Mean(i),Std(i)));
   G(i,:) = s(i)*(normcdf(xx,Mean(i),Std(i)));
   gx(i,:) = s(i)*(normpdf(x',Mean(i),Std(i)));
end


g_tot = sum(g,1);
gx_tot = sum(gx,1);
G_tot = sum(G,1)/nr_pts;

pdf_x = gx_tot;
pdf_xx = g_tot;
cdf_xx = G_tot;

if fig_nr > 0 
	figure(fig_nr), clf
	hold on
	for k = 1:nr_groups
        plot(xx,g(k,:),'b')
	end
	plot(xx,pdf_xx,'r'),title('EM-based Gaussian Mixture density')
	hold off
	drawnow
end 

% 
% PACKAGE RESULTS FOR EXPORT
%===========================

% P is matrix of size [nr_pts,nr_gauss] such that each row contains the
% soft assignment of the corresponding point to the different Gaussian components. 
% By summing over the rows, we get the proportional contribution of each
% Gaussian. 

props = sum(P,1);
props = props'/nr_pts;


mixture_parameters.props = props;
mixture_parameters.means = Mean;
mixture_parameters.stds = Std;

% 
% thresholds.vals = em_thr;
% thresholds.behavior = em_thr_behavior;
% 
% local_mins.xvals = xmin;
% local_mins.fvals = fmin;
% local_mins.quality = qmin;

mixture_density.xplot = xx;
mixture_density.pdf = pdf_xx;
mixture_density.cdf = cdf_xx;


%  PLOT RESULTING DENSITY ON REQUEST

if plot_ind == 1
    [nbin,xbin]=hist(x,30);   dxbin = xbin(2)-xbin(1);
    total_mass = sum(nbin).*dxbin;
    figure(fignr), 
    hist(x,30);
    hold on
    plot(xx,total_mass*pdf_xx,'r')
    title('Proposed GMM density')
    hold off
    drawnow
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function P_0 = initialize_em_range_based(x,nr_groups)

% divide the range up in equal parts (rather than the dataset itself, as is
% done by k-means).
% This should favour the creation of groups in long and slender tails

nr_pts = length(x);
P_raw = zeros(nr_pts,nr_groups);

xrange = range(x);
dx = xrange/(nr_groups-1);
ss = 0.5*xrange/(2*(nr_groups-1));

xc = min(x) + dx*(0:nr_groups-1);

for j = 1:nr_groups

    P_raw(:,j) = normpdf(x,xc(j),ss);  
end

%  make sure each row sums to 1:

P_0 = P_raw./(sum(P_raw,2)*ones(1,nr_groups));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TEST SUITE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 300;      % sample size
x1 = randn(n,1);
x2 = 0.5*randn(n,1)+3;

x = [x1; x2];   % mixture of two gaussians

figure(1), hist(x,30);
title('Input data')
drawnow

nr_gauss = 4

[mixture_parameters,mixture_density] = em_1dim(x,nr_gauss,77);

disp(['Extracted number of Gaussians: '  int2str(nr_gauss)])


