%% COMPUTES WHOLE BRAIN STATISTICS ON INTERACTIONS
% nc_src_ana_all2all_stat_redux

% function nc_src_ana_all2all_stat(rand_num)

% rng(rand_num,'twister')

% LESS FREQUENCIES

clear all

m = 1;

% --------------------------------------------------------
% VERSION 1:
% --------------------------------------------------------
% v_out = 2;
% v_pow = 7; % 7:  + powercorrelations
% tri = 0;
% NSUBJ = 36;
% NPERM = 2000;
% thresh = 0.05;
% perm_v = 2;
% --------------------------------------------------------
% VERSION 3: best redux verion
% --------------------------------------------------------
v_out = 3;
v_pow = 9; % powercorrelations
tri = 0;
NSUBJ = 36;
NPERM = 1000;
thresh = 0.05;
perm_v = 3;
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

% rng('shuffle');

%%

for iff = 1 : 6
  

  if ~exist([outdir sprintf('nc_src_allpowcorr_redux_f%d_v%d_processing.txt',iff,perm_v)])
    system(['touch ' outdir sprintf('nc_src_allpowcorr_redux_f%d_v%d_processing.txt',iff,perm_v)]);
  else
    continue
  end
  
  if v_pow ~= 9
    if iff == 1
      %     f = [2 3 4]; % FREQ
      f = [1 2 3];   % INDEX
    elseif iff == 2
      %     f = [5 6 7];
      f = [4 5 6];
    elseif iff == 3
      f = [7 8 9 10];
    elseif iff == 4
      f = [11 12 13 14];
    elseif iff == 5
      f = [15 16 17 18];
    elseif iff == 6
      f = [19 20 21 22 23];
    end
  else
    if iff == 1
      %     f = [2 3 4]; % FREQ
      f = [1 2 3];   % INDEX
    elseif iff == 2
      %     f = [5 6 7];
      f = [4 5 6];
    elseif iff == 3
      f = [7:13];
    elseif iff == 4
      f = [14:22];
    elseif iff == 5
      f = [23:33];
    elseif iff == 6
      f = [35:46];
    end
  end
  
 	for isubj = 1 : NSUBJ

    for ifoi = f
      
      rng('default')
      %   pause(ceil(rand*100));
      %
      disp(sprintf('Processing s%d f%d ... Loading data ...',isubj,ifoi));

      load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v_pow)])

      disp(sprintf('Processing s%d f%d ... Data loaded ...',isubj,ifoi));

      if ~tri
        allpow        = single(allpow(:,:,:,isubj));
      else
        allpow        = single(allpow(:,:,isubj));
      end
         
      allpowfreq(:,:,:,ifoi) = allpow; clear allpow

    end
    
    allpows(:,:,:,isubj) = squeeze(nanmean(allpowfreq,4)); clear allpowfreq
    
  end
    
  save([outdir sprintf('nc_src_allpowcorr_redux_f%d_v%d.mat',iff,v_out)],'allpows','-v7.3');
  
end

% error('DONE!')

%%
%
% perm_v = 2;

for ifoi = 1 : 6
  
  if ~exist([outdir sprintf('nc_src_ana_all2all_stat_redux_f%d_v%d_processing.txt',ifoi,perm_v)])
    system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_redux_f%d_v%d_processing.txt',ifoi,perm_v)]);
  else
    continue
  end
  
  load([outdir sprintf('nc_src_allpowcorr_redux_f%d_v%d.mat',ifoi,perm_v)])
  
  a     = jh_ranksum(cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:))));
  b     = a(logical(triu(ones(size(a)),1))); clear a
  pn    = [b>0 b<0];
  p     = 2*normcdf(-abs(b));
  cnt 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
  % >> cnt(:,3) => MS > HC, cnt(:,4) => HC > MS
  
  
  cntperm = nan(NPERM,4);
  allpow	= cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)));
  
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
  
  save([outdir sprintf('nc_src_ana_all2all_stat_redux_f%d_v%d.mat',ifoi,perm_v)],'cnt','cnt_perm','-v7.3');
  
  
end

error('!')

%% MULTIPLE COMPARISONS CORRECTION AND PLOTTING

v = 2;
clear all_cnt acnt

plt = 1;
NFREQ = 6;

if plt == 1
  
  clear corp thresh acnt all_cnt srt idx_D idx_R idx_D_max
  
  it = 0;
  for ifoi = 1 : 2 : NFREQ*2
    
    it = it + 1;
    %     try
    load([outdir sprintf('nc_src_ana_all2all_redux_f%d_v%d.mat',it,v)]);
    all_cnt(:,ifoi:ifoi+1)  = cnt_perm(:,3:4);
    acnt(:,ifoi:ifoi+1)     = cnt(:,3:4);
    %     catch me
    %     end
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
  for ifoi = 2 : 2 : NFREQ*2
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
    
  end
  
  plot(log10(1:6),-log10(corp),'b','LineWidth',3); hold on;
  
  % odd freqs: MS > HC
  cnt = 0;
  for ifoi = 1 : 2 : NFREQ*2
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
    
  end
  
  plot(log10(1:6),-log10(corp),'r','LineWidth',3); hold on
  line(log10([1 6]),-log10([.05 .05]),'LineStyle','--','color','k','LineWidth',2)
  
  set(gca,'XTick',[log10([1:6])],'XTickLabel',[1:6])
  set(gca,'YTick',[0 -log10(0.5) -log10(0.25) 1 2 3],'YTickLabel',[1 0.5 0.25 0.1 0.01 0.001])
  set(gca,'TickDir','out')
  
  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')
  set(gca,'YLim',[0 2.5]);
  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_redux_corr_v%d.fig'],v),'fig')
  
  
  %%
  % ----------------------------------
  % UNCORRECTED P-VALUES
  % ----------------------------------
  
  figure; set(gcf,'color','white');
  
  cnt = 0;
  
  for ifoi = 2 : 2 : NFREQ*2
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<all_cnt(:,ifoi))/NPERM;
    
  end
  
  %   plot(log10(foi_range),corp,'b','LineWidth',3); hold on;
  plot(log10(1:6),-log10(corp),'b','LineWidth',3); hold on;
  
  % odd freqs: MS > HC
  cnt = 0;
  for ifoi = 1 : 2 : NFREQ*2-1
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<all_cnt(:,ifoi))/NPERM;
    
  end
  
  %   plot(log10(foi_range),corp,'r','LineWidth',3); hold on
  plot(log10(1:6),-log10(corp),'r','LineWidth',3); hold on;
  
  line(log10([1:6]),[-log10(.05) -log10(.05)],'LineStyle','--','color','k','LineWidth',2)
  
  set(gca,'XTick',[log10([1:6])],'XTickLabel',[1:6])
  set(gca,'YTick',[1 -log10(0.05) 2 3],'YTickLabel',[0.1 0.05 0.01 0.001])
  set(gca,'TickDir','out')
  
  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')
  
  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_uncorr_redux_v%d.fig'],v),'fig')
  
  
end
