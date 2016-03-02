%% COUNTS THE NUMBER OF ALTERED LOCAL CONNECTIONS

% Computes connectivity strength WITHIN modules of the whole brain.
% Each region is inferred from the AAL atlas and local connectivity
% is defined as connections within a AAL region but not between regions

% Takes input from the following function(s):
% (1) nc_src_ana.m (v_pow!)
% (2) nc_src_powcorr_diff.m

% last update: 23-03-2015, tpfeffer

% to be implemented:
% - MC correction

% execute as
% --------------------------------------------------------
% nc_src_ana_all2all_stat_short
% --------------------------------------------------------
% function nc_src_ana_all2all_stat_short(rand_num)

clear all
%
% --------------------------------------------------------
% VERSION 1:
% --------------------------------------------------------
v_out = 1;
v_pow = 7; % 7: eloreta + powercorrelations
tri = 0;
NSUBJ = 22;
NPERM = 500;
thresh = 0.05;
% --------------------------------------------------------

if exist('rand_num','var')
  rng(rand_num,'twister')
end

m = 1;

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
  if ~exist([outdir sprintf('nc_src_ana_all2all_stat_short_f%d_v%d_processing.txt',ifoi,v_out)])
    system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_short_f%d_v%d_processing.txt',ifoi,v_out)]);
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
  load  /home/tpfeffer/pconn/matlab/aalmask_grid_coarse
  
  for ireg = 1 : 91
    
    idx = find(aalgrid.mask==ireg);
    allpow_redux{ireg} = allpow(idx,idx,:,:); clear idx
    
    a     = jh_ranksum(cat(3,squeeze(allpow_redux{ireg}(:,:,1,:)),squeeze(allpow_redux{ireg}(:,:,2,:))));
    b     = a(logical(triu(ones(size(a)),1))); clear a
    pn    = [b>0 b<0];
    p     = 2*normcdf(-abs(b));
    cnt(ireg,:) 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
    
  end
  
  clear fc x_ms m_hc s_hc tmp
  
  % ---------------------------------------------------
  % PERMUTATION TEST
  % ---------------------------------------------------
  
  allpow	= squeeze(cat(4,allpow(:,:,1,:),allpow(:,:,2,:)));
  tic
  
  for iperm = 1 : NPERM
    
    disp(sprintf('Processing freq #%d / perm #%d ...',ifoi,iperm))
    
    permidx	= randperm(NSUBJ*2);
    
    allpow_perm(:,:,1,:) = allpow(:,:,permidx(1:NSUBJ));
    allpow_perm(:,:,2,:) = allpow(:,:,permidx(NSUBJ+1:end));   
    
    for ireg = 1 : 91
      
      idx = find(aalgrid.mask==ireg);
      allpow_redux = allpow_perm(idx,idx,:,:); clear idx
      
      a     = jh_ranksum(cat(3,squeeze(allpow_redux(:,:,1,:)),squeeze(allpow_redux(:,:,2,:))));
      b     = a(logical(triu(ones(size(a)),1))); clear a
      pn    = [b>0 b<0];
      p     = 2*normcdf(-abs(b));
      cnt_perm(ireg,:,iperm) 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
      
      clear allpow_redux
      
    end
    
    
  end
  
  save([outdir sprintf('nc_src_ana_all2all_stat_short_f%d_v%d.mat',ifoi,v_out)],'cnt','cnt_perm','-v7.3');
  
  clear cnt allpow
end

%%
clear all_cnt acnt idx_R idx_D idx_D_max

v_out = 1;
NFREQ = 23;



it = 0;
  for ifoi = 1 : 2 : NFREQ*2
    it = it + 1;
    load([outdir sprintf('nc_src_ana_all2all_stat_short_f%d_v%d.mat',it,v_out)]);
    all_cnt(:,ifoi:ifoi+1)  = squeeze(sum(cnt_short_perm(:,3:4,:),1))';
    acnt(:,ifoi:ifoi+1)     = sum(cnt_short(:,3:4),1);
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
saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_short_corr_v%d.fig'],v_out),'fig')


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

  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_short_uncorr_v%d.fig'],v_out),'fig')


