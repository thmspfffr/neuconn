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
% nc_src_ana_all2all_stat_shortvslong
% --------------------------------------------------------
% function nc_src_ana_all2all_stat_shortvslong(rand_num)

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
    cnt_short(ireg,:) 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
    
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
      cnt_short_perm(ireg,:,iperm) 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
      
      clear allpow_redux
      
    end
    
    
  end
  
  save([outdir sprintf('nc_src_ana_all2all_stat_short_f%d_v%d.mat',ifoi,v_out)],'cnt_short','cnt_short_perm','-v7.3');
  
  clear cnt allpow
end
%% LONG RANGE

clear allpow_redux allpow 
for ifoi = 1 : 23
  
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
  
  load  /home/tpfeffer/pconn/matlab/aalmask_grid_coarse
  
  for ireg = 1 : 91
    
    nvox = 1 : size(allpow,1);
    
    idx = find(aalgrid.mask==ireg);
    
    nvox(idx) = [];
    
    allpow_redux = allpow(nvox,nvox,:,:);
    
    a     = jh_ranksum(cat(3,squeeze(allpow_redux(:,:,1,:)),squeeze(allpow_redux(:,:,2,:))));
    b     = a(logical(triu(ones(size(a)),1))); clear a
    pn    = [b>0 b<0];
    p     = 2*normcdf(-abs(b));
    cnt_long(ireg,:) 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
    
  end
  
  allpow	= squeeze(cat(4,allpow(:,:,1,:),allpow(:,:,2,:)));
  tic
  
  for iperm = 1 : NPERM
    
    disp(sprintf('Processing freq #%d / perm #%d ...',ifoi,iperm))
    
    permidx	= randperm(NSUBJ*2);
    
    allpow_perm(:,:,1,:) = allpow(:,:,permidx(1:NSUBJ));
    allpow_perm(:,:,2,:) = allpow(:,:,permidx(NSUBJ+1:end));
    
    for ireg = 1 : 91
      
      nvox = 1 : size(allpow,1);
      
      idx = find(aalgrid.mask==ireg);
      
      nvox(idx) = [];
      
      allpow_redux = allpow(nvox,nvox,:,:);
      
      a     = jh_ranksum(cat(3,squeeze(allpow_redux(:,:,1,:)),squeeze(allpow_redux(:,:,2,:))));
      b     = a(logical(triu(ones(size(a)),1))); clear a
      pn    = [b>0 b<0];
      p     = 2*normcdf(-abs(b));
      cnt_long_perm(ireg,:) 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
      
      clear allpow_redux
    end
  end
  save([outdir sprintf('nc_src_ana_all2all_stat_long_f%d_v%d.mat',ifoi,v_out)],'cnt_long','cnt_long_perm','-v7.3');
  
  clear cnt allpow
end

































% %
% v_out = 2;
% %
% for ifoi = 1 : 23
%   %   ifoi
%   %   try
%   load([outdir sprintf('nc_src_ana_all2all_stat_svsl_f%d_v%d.mat',ifoi,v_out)])
%   
%   p_abs(ifoi) = sum((sum(cnt.emp_abs)>sum(cnt.perm_abs,2)))/cnt.nperm;
%   pos(ifoi)   = 1-sum(sum(cnt.emp_pos)>sum(cnt.perm_pos,2))/cnt.nperm;
%   neg(ifoi)   = 1-sum(sum(cnt.emp_neg)>sum(cnt.perm_neg,2))/cnt.nperm;
%   
% end
% 
% figure; hold on; set(gcf,'color','white');
% 
% plot(log10(foi_range),pos,'LineWidth',5,'color','r')
% plot(log10(foi_range),neg,'LineWidth',5,'color','b')
% 
% line([log10(2) log10(128)],[0.05 0.05],'color','k','LineWidth',2,'LineStyle','--')
% 
% set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
% set(gca,'TickDir','out')
% 
% title('Connectivity differences')
% xlabel('CARRIER FREQUENCY (HZ)')
% ylabel('P-VALUE')
% 
% saveas(gcf,sprintf([plotdir 'nc_src_all2all_stat_indiv_v%d.fig'],v_out),'fig')

% close

%   for isubj = 1 : NSUBJ
%
%    tmp = triu(a(:,:,isubj));
%    cnt(ifoi).abs(isubj) = sum(abs(tmp(:)));
%    cnt(ifoi).rel(isubj) = sum(abs(tmp(:)))/length(tmp(:));
%    cnt(ifoi).pos(isubj) = sum(tmp(:)>0)/length(tmp(:));
%    cnt(ifoi).neg(isubj) = sum(tmp(:)<0)/length(tmp(:));
%
%   end
% %     cnt(ifoi).p = nanmean(cnt(ifoi).rel);
%   p(ifoi) = nanmean(cnt(ifoi).pos);
%     n(ifoi) = nanmean(cnt(ifoi).neg);
%
%   catch me
%   end
%
% end
%
% h=figure; set(h,'color','w'); hold on;
% plot(n,'b','LineWidth',5);
% plot(p,'r','LineWidth',5);
%
% set(gca,'XTick',[1:2:23],'XTickLabel',[foi_range(1:2:end)])
% line([0 23],[0.05 0.05],'LineStyle','--','color','k')
%
%
%   para = [];
% para.mydotmarkersize=20;
% para.orientation='coronal';
% % para.colorlimits = [0 0.2];
%
% h=figure; hold on
% set(h,'color','w');
%
% showmri_transp_v3(mri,para,[grid nanmean(a(:,:,1),1)'],grid(1,:));




%









%
%   switch test
%
%     case 'ttest'
%
%       [a,~,~,t]     = ttest2(squeeze(allpow(:,1,:)),squeeze(allpow(:,2,:)),'dim',2);
%       t             = t.tstat;
%       cnt           = [nansum(a) sum(t(logical(a))>0) sum(t(logical(a))<0)]; clear a
%
%     case 'ranksum'
%
%       a     = jh_ranksum(cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:))));
%       b     = a(logical(triu(ones(size(a)),1))); clear a
%       pn    = [b>0 b<0];
%       p     = 2*normcdf(-abs(b));
%       cnt 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
%       % >> cnt(:,3) => MS > HC, cnt(:,4) => HC > MS
%   end
%
%   if perm
%
%     cntperm = nan(NPERM,4);
%
%     switch test
%
%       % if statistical test is unpaired ttest
%       case 'ttest'
%         allpow	= squeeze(cat(3,allpow(:,1,:),allpow(:,2,:)));
%
%         for iperm = 1 : NPERM
%
%           clear a_perm
%           disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
%
%           permidx           = randperm(NSUBJ*2);
%           [a_perm,~,~,t]    = ttest2(allpow(:,permidx(1:NSUBJ)),allpow(:,permidx(NSUBJ+1:end)),'dim',2);
%           t                 = t.tstat;
%           cnt_perm(iperm,:)	= [nansum(a_perm(:)) sum(t(logical(a_perm))>0) sum(t(logical(a_perm))<0)];
%
%         end
%
%         % if statistical test is ranksum test
%       case 'ranksum'
%
%         allpow	= cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)));
%
%         for iperm = 1 : NPERM
%
%           clear a_perm
%           disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
%
%           permidx           = randperm(NSUBJ*2);
%           a_perm            = jh_ranksum(cat(3,squeeze(allpow(:,:,permidx(1:NSUBJ))),squeeze(allpow(:,:,permidx(NSUBJ+1:end)))));
%           b_perm            = a_perm(logical(triu(ones(size(a_perm)),1))); clear a
%           pn_perm           = [b_perm>0 b_perm<0];
%           p_perm            = 2*normcdf(-abs(b_perm));
%           cnt_perm(iperm,:)	= [sum(p_perm<thresh) sum(p_perm<thresh)/length(p_perm) sum(p_perm(pn_perm(:,1))<thresh) sum(p_perm(pn_perm(:,2))<thresh)];
%
%         end
%     end
%
%     save([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d.mat',ifoi,perm_v)],'cnt','cnt_perm','-v7.3');
%
%   end
% end
%
%
% %% MULTIPLE COMPARISONS CORRECTION AND PLOTTING
%
% v  = 7;
%
% plt = 1;
% NFREQ = 23;
%
% if plt == 1
%
%   clear corp thresh acnt all_cnt srt
%
%   it = 0;
%   for ifoi = 1 : 2 : NFREQ*2
%     it = it + 1;
%     load([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d.mat',it,perm_v)]);
%     all_cnt(:,ifoi:ifoi+1)  = cnt_perm(:,3:4);
%     acnt(:,ifoi:ifoi+1)     = cnt(:,3:4);
%   end
%
% %   ranks, where 1 = min, NPERM = max
% %   for ifreq = 1 : NFREQ*2
% %     [~,~,idx_D(:,ifreq)] = unique(all_cnt(:,ifreq));
% %   end
% %
%  for ifreq = 1 : NFREQ*2
%     [idx_D(:,ifreq)] = floor(tiedrank(all_cnt(:,ifreq)));
%  end
% %
%   % get maximum rank across frequencies and directions
%   idx_R     = max(idx_D,[],2);
%   idx_D_max = sort(all_cnt,'ascend');
%   idx_D_max = idx_D_max(idx_R,:);
%
%   % not sure about this part
%   % ----------------------------------
%   % thresh	= prctile(idx_R,95);
%   % all_srt = sort(all_cnt,'ascend');
%   % threshs = all_srt(thresh,:);
%   % ----------------------------------
%
%   % ----------------------------------
%   % PLOT CORRECTED P-VALUES
%   % ----------------------------------
%   % even freqs: HC > MS
%   for ifoi = 2 : 2 : 46
%
%     corp(ifoi)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
%
%   end
%
%   plot(corp(2 : 2 : 46),'b','LineWidth',3); hold on;
%
%   % odd freqs: MS > HC
%   for ifoi = 1 : 2 : 45
%
%     corp(ifoi)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
%
%   end
%
%   plot(corp(1 : 2 : 45),'r','LineWidth',3); hold on
%   line([1 23],[.05 .05],'LineStyle','--','color','k','LineWidth',2)
%
%   ylabel('Corrected p-value');
%   xlabel('Carrier frequency (in Hz)');
%   title('Differences in connectivity');
%
%   set(gca,'TickDir','out','XTick',[1 4 7 11 15 19 23],'XTickLabel',[2 4 8 16 32 64 128])
%   saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_corr_v%d.fig'],v),'fig')
%
%   % ----------------------------------
%   % UNCORRECTED P-VALUES
%   % ----------------------------------
%
%   figure
%
%   for ifoi = 2 : 2 : 46
%
%     corp(ifoi)  = sum(acnt(ifoi)<all_cnt(:,ifoi))/NPERM;
%
%   end
%
%   plot(corp(2 : 2 : 46),'b','LineWidth',3); hold on;
%
%   % odd freqs: MS > HC
%   for ifoi = 1 : 2 : 45
%
%     corp(ifoi)  = sum(acnt(ifoi)<all_cnt(:,ifoi))/NPERM;
%
%   end
%
%   plot(corp(1 : 2 : 45),'r','LineWidth',3); hold on
%   line([1 23],[.05 .05],'LineStyle','--','color','k','LineWidth',2)
%
%   ylabel('Uncorrected p-value');
%   xlabel('Carrier frequency (in Hz)');
%   title('Differences in connectivity');
%
%   set(gca,'TickDir','out','XTick',[1 4 7 11 15 19 23],'XTickLabel',[2 4 8 16 32 64 128])
%
%   saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_uncorr_v%d.fig'],v),'fig')
%
%
% end
