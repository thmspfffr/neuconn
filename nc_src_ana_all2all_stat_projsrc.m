%% COMPUTES WHOLE BRAIN STATISTICS ON INTERACTIONS
% nc_src_ana_all2all_stat_projsrc

% Projects significant differences in connecitity onto the brain.
% todo: implement fieldtrip sourceplot

clear all

m = 1;

% --------------------------------------------------------
% VERSION 2: thresh 0.05 - coarse grid
% --------------------------------------------------------
v       = 7;          % version of head model
perm_v  = 2;
NPERM   = 500;
NSUBJ   = 22;
perm    = 1;
thresh  = 0.05;
test    = 'ranksum';
tri     = 0;
% --------------------------------------------------------
% VERSION 3: thresh 0.01 - coarse grid
% --------------------------------------------------------
% v       = 7;          % version of head model
% perm_v  = 3;
% NPERM   = 1000;
% NSUBJ   = 24;
% perm    = 1;
% thresh  = 0.01;
% test    = 'ranksum';
% tri     = 0;
% --------------------------------------------------------
% VERSION 4: thresh 0.05 - cortex grid
% --------------------------------------------------------
% v       = 8;          % version of head model
% perm_v  = 4;
% NPERM   = 500;
% NSUBJ   = 22;
% perm    = 1;
% thresh  = 0.05;
% test    = 'ranksum';
% tri     = 0;
% --------------------------------------------------------
% VERSION 5: thresh 0.01 - cortex grid
% --------------------------------------------------------
% v       = 8;          % version of head model
% perm_v  = 5;
% NPERM   = 500;
% NSUBJ   = 22;
% perm    = 1;
% thresh  = 0.01;
% test    = 'ranksum';
% tri     = 0;
% --------------------------------------------------------
% VERSION 6: thresh 0.05 - coarse grid, 10000 perm
% --------------------------------------------------------
% v       = 7;          % version of head model
% perm_v  = 6;
% NPERM   = 10000;
% NSUBJ   = 22;
% perm    = 1;
% thresh  = 0.05;
% test    = 'ranksum';
% tri     = 0;
% --------------------------------------------------------
% VERSION 7: thresh 0.05 - cortex grid
% --------------------------------------------------------
% v       = 8;          % version of head model
% perm_v  = 7;
% NPERM   = 500;
% NSUBJ   = 22;
% perm    = 1;
% thresh  = 0.05;
% test    = 'ranksum';
% tri     = 0;

restoredefaultpath
foi_range = unique(round(2.^[1:.25:7]));

% addpaths
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest/

ft_defaults

% Define relevant input/output paths
indir = '/home/gnolte/neuconn/meg_out/rest/conn/';
mridir = '/home/gnolte/neuconn/mri_data/';
outdir = '/home/gnolte/neuconn/meg_out/rest/stats/';

% call subj info startup function
allsubj = nc_allsubj_start(m);
v_out = v;      % version of output

% rng('shuffle');

%%

for ifoi = 1 : 23

  rng('default')
  %   pause(ceil(rand*100));
%   
%   if ~exist([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d_processing.txt',ifoi,perm_v)])
%     system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d_processing.txt',ifoi,perm_v)]);
%   else
%     continue
%   end
  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  if ~tri
    allpow        = single(allpow(:,:,:,1:NSUBJ));
  else
    allpow        = single(allpow(:,:,1:NSUBJ));
  end
  

      a     = jh_ranksum(cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)))); 
      p     = 2*normcdf(-abs(a));
      t     = p<0.05;
      s(:,ifoi)     = sum(t);
%       cnt 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
      % >> cnt(:,3) => MS > HC, cnt(:,4) => HC > MS

end


%% PLOT

for ifoi = 1 : 23

  ifoi 
  
  load sa_meg_template;
  grid = sa_meg_template.grid_coarse;
  mri  = sa_meg_template.mri;

  para                  = [];
  para.mydotmarkersize  = 40;
  para.orientation      = 'axial';
  % para.colorlimits      = [0 1];
  para.colormaps        = {'jet'};

  h = figure; hold on
  set(h,'color','w');

  showmri_transp_v3(mri,para,[grid s(:,ifoi)],grid);

  saveas(gcf,sprintf('/home/gnolte/neuconn/meg_out/rest/plots/nc_src_num_conn_f%d.fig',ifoi))
 

end

