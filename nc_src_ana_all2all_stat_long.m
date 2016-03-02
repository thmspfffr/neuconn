%% COUNTS THE NUMBER OF ALTERED S/L CONNECTIONS

% Computes strength of connections from one hemisphere to the
% other, i.e., taking only long range connections into account

% Takes input from the following function(s):
% (1) nc_src_ana.m (v_pow!)
% (2) nc_src_powcorr_diff.m

% last update: 23-03-2015, tpfeffer

% to be implemented:
% - MC correction

% execute as
% --------------------------------------------------------
% nc_src_ana_all2all_stat_long
% --------------------------------------------------------
% function nc_src_ana_all2all_stat_long(rand_num)

clear all
%
% VERSION 1:
% --------------------------------------------------------
v_out = 1;
v_pow = 7; % 7:  + powercorrelations
tri = 0;
NSUBJ = 36;
NPERM = 2000;
thresh = 0.05;
% --------------------------------------------------------

if exist('rand_num','var')
  rng(rand_num,'twister')
end

m = 1;

if NPERM > 1
  perm  = 1;
end

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

% rng('shuffle');

%%

for ifoi = 1 : 23
  
  % Waiting time for torque queue
  if exist('rand_num','var')
    w = ceil(rand*150);
    disp(sprintf('Waiting time: %d sec ...',w));
    pause(ceil(rand*150));
    rng(1)
  end
  
  %
  if ~exist([outdir sprintf('nc_src_ana_all2all_stat_long_f%d_v%d_processing.txt',ifoi,v_out)])
    system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_long_f%d_v%d_processing.txt',ifoi,v_out)]);
  else
    continue
  end
  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v_pow)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  if ~tri
    allpow        = single(allpow(:,:,:,1:NSUBJ));
  else
    allpow        = single(allpow(:,:,1:NSUBJ));
  end
  
  % WITHIN:
%   load  /home/tpfeffer/pconn/matlab/aalmask_grid_coarse

  % extract right (1st dim) and left (2nd dim) hemisphere
  allpow_lr = allpow(1:ceil(2113/2)-150,floor(2113/2)+151:2113,:,1:NSUBJ); clear allpow

  a     = jh_ranksum(cat(3,squeeze(allpow_lr(:,:,1,:)),squeeze(allpow_lr(:,:,2,:))));
  b     = a(logical(triu(ones(size(a)),1))); clear a
  pn    = [b>0 b<0];
  p     = 2*normcdf(-abs(b));
  cnt 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a

  clear fc x_ms m_hc s_hc tmp
  
  if perm
    
    cntperm = nan(NPERM,4);

    % if statistical test is unpaired ttest
    allpow_lr	= cat(3,squeeze(allpow_lr(:,:,1,:)),squeeze(allpow_lr(:,:,2,:)));    
    
    for iperm = 1 : NPERM
      
      clear a_perm
      disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
      
      permidx           = randperm(NSUBJ*2);
      a_perm            = jh_ranksum(cat(3,squeeze(allpow_lr(:,:,permidx(1:NSUBJ))),squeeze(allpow_lr(:,:,permidx(NSUBJ+1:end)))));
      b_perm            = a_perm(logical(triu(ones(size(a_perm)),1))); clear a
      pn_perm           = [b_perm>0 b_perm<0];
      p_perm            = 2*normcdf(-abs(b_perm));
      cnt_perm(iperm,:)	= [sum(p_perm<thresh) sum(p_perm<thresh)/length(p_perm) sum(p_perm(pn_perm(:,1))<thresh) sum(p_perm(pn_perm(:,2))<thresh)];
      
    end
  end
    save([outdir sprintf('nc_src_ana_all2all_stat_long_f%d_v%d.mat',ifoi,v_out)],'cnt','cnt_perm','-v7.3');
    
    clear cnt allpow
end





%%
v_out = 1;
NFREQ = 23;
clear all_cnt acnt idx_D idx_R idx_D_max

it = 0;
  for ifoi = 1 : 2 : NFREQ*2
    it = it + 1;
    load([outdir sprintf('nc_src_ana_all2all_stat_long_f%d_v%d.mat',it,v_out)]);
    all_cnt(:,ifoi:ifoi+1)  = cnt_perm(:,3:4);
    acnt(:,ifoi:ifoi+1)     = cnt(:,3:4);
  end
  

 for ifreq = 1 : NFREQ*2
    [idx_D(:,ifreq)] = floor(tiedrank(all_cnt(:,ifreq)));
 end
%   
  % get maximum rank across frequencies and directions
  idx_R     = max(idx_D,[],2);
  idx_D_max = sort(all_cnt,'ascend');
  idx_D_max = idx_D_max(idx_R,:);
  
%
  figure; set(gcf,'color','white'); hold on; box on
  % ----------------------------------
  % PLOT CORRECTED P-VALUES
  % ----------------------------------
  % even freqs: HC > MS
  cnt = 0;
  for ifoi = 2 : 2 : 46
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
    
  end
  
  plot(log10(foi_range),-log10(corp),'b','LineWidth',3); hold on; 
  
  
 
  clear corp
  % odd freqs: MS > HC
  cnt = 0;
  for ifoi = 1 : 2 : 45
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
    
  end
  
  plot(log10(foi_range),-log10(corp),'r','LineWidth',3); hold on
  line(log10([foi_range(1) foi_range(end)]),-log10([.05 .05]),'LineStyle','--','color','k','LineWidth',2)
  
  set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
  set(gca,'YTick',[0 -log10(0.5) -log10(0.25) 1 2 3],'YTickLabel',[1 0.5 0.25 0.1 0.01 0.001])
  set(gca,'TickDir','out')

  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')
set(gca,'YLim',[0 2.5]);
saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_long_corr_v%d.fig'],v_out),'fig')


% ----------------------------------
  % UNCORRECTED P-VALUES
  % ----------------------------------

  figure; set(gcf,'color','white');

  cnt = 0;
  
  for ifoi = 2 : 2 : 46
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<all_cnt(:,ifoi))/NPERM;
    
  end
  
%   plot(log10(foi_range),corp,'b','LineWidth',3); hold on; 
  plot(log10(foi_range),-log10(corp),'b','LineWidth',3); hold on; 

  % odd freqs: MS > HC
   cnt = 0;
  for ifoi = 1 : 2 : 45
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<all_cnt(:,ifoi))/NPERM;
    
  end
  
%   plot(log10(foi_range),corp,'r','LineWidth',3); hold on
  plot(log10(foi_range),-log10(corp),'r','LineWidth',3); hold on; 

  line(log10([foi_range(1) foi_range(end)]),[-log10(.05) -log10(.05)],'LineStyle','--','color','k','LineWidth',2)
  
  set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
  set(gca,'YTick',[1 -log10(0.05) 2 3],'YTickLabel',[0.1 0.05 0.01 0.001])
  set(gca,'TickDir','out')
set(gca,'YLim',[0 2.5]);

  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')

  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_long_uncorr_v%d.fig'],v_out),'fig')





