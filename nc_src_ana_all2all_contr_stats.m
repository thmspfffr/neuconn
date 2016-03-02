%% COMPUTES WHOLE BRAIN CONTRAST IN DIFFERENT FREQ BANDS
% nc_src_ana_all2all_contr_stats

% function nc_src_ana_all2all_stat(rand_num)

% rng(rand_num,'twister')

% Counts the number of significant connections and compares those 
% across the two conditions using a permutation test. Does this in 
% a frequency dependent manner.

clear all

m = 1;

% --------------------------------------------------------
% VERSION 5
% --------------------------------------------------------
v         = 5;         % version of interaction measure (powcorr / lph)
perm_v    = 5;
NSUBJ     = 24;
foi_range = unique(round(2.^[1:.25:7]));
thresh    = 0.05;
NPERM     = 1000;
perm      = 1;
% --------------------------------------------------------


restoredefaultpath

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
v_out   = v;      % version of output

%%

for ifoi = 1 : length(foi_range)
  
  %   pause(ceil(rand*100));
  
  if ~exist([outdir sprintf('nc_src_ana_all2all_contr_stats_f%d_v%d_processing.txt',ifoi,perm_v)])
    system(['touch ' outdir sprintf('nc_src_ana_all2all_contr_stats_f%d_v%d_processing.txt',ifoi,perm_v)]);
  else
    continue
  end
  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  allpow = single(allpow(:,:,:,1:NSUBJ));
  
  a     = jh_ranksum(cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:))));
  a     = 2*normcdf(-abs(a));
  
  if perm
    
    cnt     = nansum(a(:)<thresh)/2; clear a
    allpow	= cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)));
    
    for iperm = 1 : NPERM
      
      disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
      
      permidx         = randperm(NSUBJ*2);
      a_perm         	= jh_ranksum(cat(3,squeeze(allpow(:,:,permidx(1:NSUBJ))),squeeze(allpow(:,:,permidx(NSUBJ+1:end)))));
      a_perm          = 2*normcdf(-abs(a_perm));
      cnt_perm(iperm) = nansum(a_perm(:)<thresh)/2;
      
    end
    
    save([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d.mat',ifoi,perm_v)],'cnt','cnt_perm');
  else
    cnt(ifoi) = nansum(a(:)); clear a
  end
end

%% PLOT STATISTICS

% for ifoi = 1 : 23
% %   try
%     load([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d.mat',ifoi,perm_v)]);
% %   catch me
% %   end
%   
%   p(ifoi) = sum(cnt_perm>cnt)/length(cnt_perm);
%   c(ifoi) = cnt;
%   
%   clear cnt cnt_perm
% end
% 
% figure;
% 
% plot(p,'r','LineWidth',5)
% 
% line([1 23],[.05 0.05],'LineStyle','--','color','k','LineWidth',2)
% 
% ylabel('p-value');
% xlabel('Carrier frequency (in Hz)');
% title('Differences in connectivity');
% 
% set(gca,'TickDir','out','XTick',[1 4 7 11 15 19 23],'XTickLabel',[2 4 8 16 32 64 128])
% 
%   


