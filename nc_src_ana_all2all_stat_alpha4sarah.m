%% COMPUTES WHOLE BRAIN STATISTICS ON INTERACTIONS
% nc_src_ana_all2all_stat

% function nc_src_ana_all2all_stat(rand_num)

% rng(rand_num,'twister')

% Counts the number of significant connections and compares those
% across the two conditions using a permutation test. Does this in
% a frequency dependent manner.

clear all

m = 1;

% --------------------------------------------------------
% VERSION 2: thresh 0.05 - coarse grid
% --------------------------------------------------------
% v       = 7;          % version of head model
% perm_v  = 2;
% NPERM   = 500;
% NSUBJ   = 36;
% perm    = 1;
% thresh  = 0.05;
% test    = 'ranksum';
% tri     = 0;
% foi     = 7:9;
% --------------------------------------------------------
% VERSION 3: thresh 0.05 - coarse grid
% --------------------------------------------------------
% v       = 7;          % version of head model
% perm_v  = 3;
% NPERM   = 500;
% NSUBJ   = 36;
% perm    = 1;
% thresh  = 0.05;
% test    = 'ranksum';
% tri     = 0;
% foi     = 11:14;
% --------------------------------------------------------
% VERSION 3: thresh 0.05 - coarse grid
% --------------------------------------------------------
% v       = 7;          % version of head model
% perm_v  = 4;
% NPERM   = 500;
% NSUBJ   = 36;
% perm    = 1;
% thresh  = 0.05;
% test    = 'ranksum';
% tri     = 0;
% foi     = 15:18;
% --------------------------------------------------------
% VERSION 3: thresh 0.05 - coarse grid
% --------------------------------------------------------
v       = 7;          % version of head model
perm_v  = 5;
NPERM   = 500;
NSUBJ   = 36;
perm    = 1;
thresh  = 0.05;
test    = 'ranksum';
tri     = 0;
foi     = 4:6;
% --------------------------------------------------------


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
indir   = '/home/gnolte/neuconn/meg_out/rest/conn/';
mridir  = '/home/gnolte/neuconn/mri_data/';
outdir  = '/home/gnolte/neuconn/meg_out/rest/stats/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';
% call subj info startup function
allsubj = nc_allsubj_start(m);
v_out = v;      % version of output

% rng('shuffle');

%%
cnt = 0;
for ifoi = foi
  
  rng('default')
  %   pause(ceil(rand*100));

  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  cnt = cnt + 1;
  allallpow(:,:,:,:,cnt)        = single(allpow(:,:,:,1:NSUBJ));

  
end
  
clear allpow   
allpow = nanmean(allallpow,5);

load  /home/tpfeffer/pconn/matlab/aalmask_grid_coarse

aalpow = nan(91,91,2,NSUBJ);



for reg1 = 1 : 91
  reg1
  for reg2 = 1 : 91
    
    if reg1 == reg2
      continue
    end
    
    idx1 = find(aalgrid.mask==reg1);
    idx2 = find(aalgrid.mask==reg2);
    
    for ivox = 1 : length(idx1)
      tmp(ivox,:,:) = squeeze(nanmean(allpow(idx1(ivox),idx2,:,:),2));
    end
    
    %       tmp = nanmean(tmp);
    
    allcorr(reg1,reg2,:,:) = squeeze(nanmean(tmp,1));
    
    clear tmp
    
  end
end

save([outdir sprintf('nc_src_ana_all2all_stat_alpha4sarah_v%d.mat',perm_v)],'allcorr','allsubj','-v7.3');
    


