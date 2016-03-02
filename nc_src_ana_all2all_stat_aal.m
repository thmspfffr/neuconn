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
v       = 1;          % version of head model
perm_v  = 1;
NPERM   = 500;
NSUBJ   = 22;
perm    = 1;
thresh  = 0.05;
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
load([outdir sprintf('nc_src_ana_all2all_stat_long_aal_fc_v%d.mat',v_out)]);

for ifoi = 1 : 23
  
  rng('default')
  %   pause(ceil(rand*100));
  
  if ~exist([outdir sprintf('nc_src_ana_all2all_stat_aal_f%d_v%d_processing.txt',ifoi,perm_v)])
    system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_aal_f%d_v%d_processing.txt',ifoi,perm_v)]);
  else
    continue
  end

      a     = jh_ranksum(cat(3,squeeze(allcorr(:,:,1,1:NSUBJ,ifoi)),squeeze(allcorr(:,:,2,1:NSUBJ,ifoi))));
      b     = a(logical(triu(ones(size(a)),1))); clear a
      pn    = [b>0 b<0];
      p     = 2*normcdf(-abs(b));
      cnt 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
      % >> cnt(:,3) => MS > HC, cnt(:,4) => HC > MS

  if perm
    
    cntperm = nan(NPERM,4);
    
        allpow	= cat(3,squeeze(allcorr(:,:,1,:,ifoi)),squeeze(allcorr(:,:,2,:,ifoi)));
        
        for iperm = 1 : NPERM
          
          clear a_perm
          disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
          
          permidx           = randperm(NSUBJ*2);
          a_perm            = jh_ranksum(cat(3,squeeze(allpow(:,:,permidx(1:NSUBJ))),squeeze(allpow(:,:,permidx(NSUBJ+1:end)))));
          b_perm            = a_perm(logical(triu(ones(size(a_perm)),1))); clear a
          pn_perm           = [b_perm>0 b_perm<0];
          p_perm            = 2*normcdf(-abs(b_perm));
          cnt_perm(iperm,:)	= [sum(p_perm<thresh) sum(p_perm<thresh)/length(p_perm) sum(p_perm(pn_perm(:,1))<thresh) sum(p_perm(pn_perm(:,2))<thresh)];
          

          
    end
    
    save([outdir sprintf('nc_src_ana_all2all_stat_aal_f%d_v%d.mat',ifoi,perm_v)],'cnt','cnt_perm','-v7.3');
    
  end
end


%% MULTIPLE COMPARISONS CORRECTION AND PLOTTING

v  = 1;

plt = 1;
NFREQ = 23;

if plt == 1
  
  clear corp thresh acnt all_cnt srt
  
  it = 0;
  
  for ifoi = 1 : 2 : NFREQ*2
    it = it + 1;
    load([outdir sprintf('nc_src_ana_all2all_stat_aal_f%d_v%d.mat',it,v)]);
    all_cnt(:,ifoi:ifoi+1)  = cnt_perm(:,3:4);
    acnt(:,ifoi:ifoi+1)     = cnt(:,3:4);
  end
  
%   ranks, where 1 = min, NPERM = max
%   for ifreq = 1 : NFREQ*2
%     [~,~,idx_D(:,ifreq)] = unique(all_cnt(:,ifreq));
%   end
%   
 for ifreq = 1 : NFREQ*2
    [idx_D(:,ifreq)] = floor(tiedrank(all_cnt(:,ifreq)));
 end
%   
  % get maximum rank across frequencies and directions
  idx_R     = max(idx_D,[],2);
  idx_D_max = sort(all_cnt,'ascend');
  idx_D_max = idx_D_max(idx_R,:);
  
  % not sure about this part
  % ----------------------------------
  % thresh	= prctile(idx_R,95);
  % all_srt = sort(all_cnt,'ascend');
  % threshs = all_srt(thresh,:);
  % ----------------------------------
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
saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_corr_v%d.fig'],v),'fig')


%%
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

  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')

  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_uncorr_v%d.fig'],v),'fig')


end
