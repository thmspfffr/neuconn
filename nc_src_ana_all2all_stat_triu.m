% nc_src_ana_all2all_stat_triu

% Reduce nc_src_allpowcorr to just the upper triangle to save memory

clear all

m = 1;

% --------------------------------------------------------
% VERSION 2 all2all power correlations (SEGLEN too short for freq res)
% --------------------------------------------------------
v       = 5;          % version of head model
perm_v  = 5;
NSUBJ   = 18;
% --------------------------------------------------------

foi_range = unique(round(2.^[1:.25:7]));

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest/

ft_defaults

% Define relevant input/output paths
indir = '/home/gnolte/neuconn/meg_out/rest/conn/';
outdir = '/home/gnolte/neuconn/meg_out/rest/stats/';

% call subj info startup function
allsubj = nc_allsubj_start(m);
v_out = v;      % version of output

%%

for ifoi = 1 : length(foi_range)
  
  if ~exist([outdir sprintf('nc_src_ana_all2all_stat_triu_f%d_v%d_processing.txt',ifoi,perm_v)])
    system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_triu_f%d_v%d_processing.txt',ifoi,perm_v)]);
  else
    continue
  end

  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));

  allpow = single(allpow);
  
  allpow = triu(allpow);
  
  save([outdir sprintf('nc_src_ana_all2all_stat_triu_f%d_v%d.mat',ifoi,perm_v)],'cnt','cnt_perm');

end


