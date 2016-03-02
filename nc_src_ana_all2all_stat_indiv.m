%% COUNTS THE NUMBER OF ALTERED CONNECTIONS INDIVIDUALLY FOR EACH SUBJECT

% Computes strength of connection from each voxel to the rest of the brain.
% Takes this as statistical threshold to compare each voxel's connectivity
% in the diseased brain. Then counts the number of altered connections for
% each subject individually and performs statistical comparison as
% implemented in nc_src_ana_all2all_stat.m

% Takes input from the following function(s):
% (1) nc_src_ana.m (v_pow!)
% (2) nc_src_powcorr_diff.m

% last update: 22-02-2015, tpfeffer

% to be implemented:
% - statistical test
% - paired test instead of taking average across all HCs

% execute as
% --------------------------------------------------------
% nc_src_ana_all2all_stat_indiv
% --------------------------------------------------------
% function nc_src_ana_all2all_stat_indiv(rand_num)

clear all
%
% --------------------------------------------------------
% VERSION 1:
% --------------------------------------------------------
% v_out = 1;
% v_pow = 10; % 10: beamformer + powercorrelations
% tri = 0;
% NSUBJ = 22;
% NPERM = 1000;
% --------------------------------------------------------
% VERSION 2:
% --------------------------------------------------------
% v_out = 2;
% v_pow = 7; % 7: eloreta + powercorrelations
% tri = 0;
% NSUBJ = 22;
% NPERM = 1000;
% --------------------------------------------------------
% VERSION 3:
% --------------------------------------------------------
% v_out = 3;
% v_pow = 7; % 7: eloreta + powercorrelations
% tri = 0;
% NSUBJ = 22;
% NPERM = 10000;
% --------------------------------------------------------
% VERSION 4: maybe-bug fixed (or introduced)
% --------------------------------------------------------
v_out = 4;
v_pow = 7; % 7: eloreta + powercorrelations
tri = 0;
NSUBJ = 22;
NPERM = 500;
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
%   if ~exist([outdir sprintf('nc_src_ana_all2all_stat_indiv_f%d_v%d_processing.txt',ifoi,v_out)])
%     system(['touch ' outdir sprintf('nc_src_ana_all2all_stat_indiv_f%d_v%d_processing.txt',ifoi,v_out)]);
%   else
%     continue
%   end
  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load([indir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v_pow)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  if ~tri
    allpow        = single(allpow(:,:,:,1:NSUBJ));
  else
    allpow        = single(allpow(:,:,1:NSUBJ));
  end
  
  % ---------------------------------------------------
  % START WITH EMPIRICAL DATA
  % ---------------------------------------------------
    
  fc = zeros(size(allpow,1),NSUBJ);
  
  % Compute average correlation of ivox + fisher z-transform
  m_hc = squeeze(nanmean(nanmean(atanh(allpow(:,:,2,:)),2),4));
  s_hc = squeeze(nanmean(nanstd(atanh(allpow(:,:,2,:)),[],2),4));

%   x_ms = squeeze(atanh(allpow(:,:,1,:)));
  x_ms = squeeze(nanmean(atanh(allpow(:,:,1,:)),2));

  for i = 1 :  size(allpow,1)
    disp(sprintf('Processing freq #%d perm #0 / voxel #%d ...',ifoi,i))
    
    % has been z-scored previously, now use t-score (divide by sqrt(N))
%     fc(i,:,:) = sign(x_ms(i,:,:) - m_hc(i)).*(abs(((x_ms(i,:,:) - m_hc(i)) ./ (s_hc(i)/sqrt(NSUBJ)))) > 1.96);
    fc(i,:) = sign(x_ms(i,:) - m_hc(i)).*(abs(((x_ms(i,:) - m_hc(i)) ./ (s_hc(i)/sqrt(NSUBJ)))) > 1.96);
  
  end
  %       end
  save([outdir sprintf('nc_src_ana_all2all_stat_indiv_fc_f%d_v%d.mat',ifoi,v_out)],'fc','-v7.3');
  
  for isubj = 1 : NSUBJ
    
%     tmp = triu(fc(:,:,isubj),1);
        tmp = fc(:,isubj);

    cnt.emp_abs(isubj) = sum(abs(tmp(:)));
    cnt.emp_rel(isubj) = sum(abs(tmp(:)))/length(tmp(:));
    cnt.emp_pos(isubj) = sum(tmp(:)>0);
    cnt.emp_neg(isubj) = sum(tmp(:)<0);
    
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
     
    fc = zeros(size(allpow,1),NSUBJ);
      
    % Compute average correlation of ivox + fisher z-transform
    m_hc = squeeze(nanmean(nanmean(atanh(allpow_perm(:,:,2,:)),2),4));
    s_hc = squeeze(nanmean(nanstd(atanh(allpow_perm(:,:,2,:)),[],2),4));
    
    % it's not clear to me, whether threshold should be fixed or 
    % whether it's rather the empirical value (or maybe nothing fixed?)
    
%     x_ms = squeeze(atanh(allpow_perm(:,:,1,:)));
    x_ms = squeeze(nanmean(atanh(allpow_perm(:,:,1,:)),2));

  for i = 1 :  size(allpow,1)
    disp(sprintf('Processing freq #%d perm #0 / voxel #%d ...',ifoi,i))
    
    % has been z-scored previously, now use t-score (divide by sqrt(N))
%     fc(i,:,:) = sign(x_ms(i,:,:) - m_hc(i)).*(abs(((x_ms(i,:,:) - m_hc(i)) ./ (s_hc(i)/sqrt(NSUBJ)))) > 1.96);
    fc(i,:) = sign(x_ms(i,:) - m_hc(i)).*(abs(((x_ms(i,:) - m_hc(i)) ./ (s_hc(i)/sqrt(NSUBJ)))) > 1.96);
  
  end
    for isubj = 1 : NSUBJ
      
%       tmp = triu(fc(:,:,isubj),1);
            tmp = fc(:,isubj);

      cnt.perm_abs(iperm,isubj) = sum(abs(tmp(:)));
      cnt.perm_rel(iperm,isubj) = sum(abs(tmp(:)))/length(tmp(:));
      cnt.perm_pos(iperm,isubj) = sum(tmp(:)>0);
      cnt.perm_neg(iperm,isubj) = sum(tmp(:)<0);
      
    end
    
    cnt.v_pow = v_pow;
    cnt.nperm = NPERM;
    
    clear tmp  x_ms m_hc s_hc permidx
    
  end
  
  save([outdir sprintf('nc_src_ana_all2all_stat_indiv_f%d_v%d.mat',ifoi,v_out)],'cnt','-v7.3');
  
  clear cnt allpow
end

%% PLOT
if plt
    str = {'Beamformer';'eLoreta'};

  v_out = 1;
  %
  NFREQ = 23;
  
  for ifoi = 1 : NFREQ
  %   ifoi
%     try
    load([outdir sprintf('nc_src_ana_all2all_stat_indiv_f%d_v%d.mat',ifoi,v_out)])

    p_abs(ifoi) = sum((sum(cnt.emp_abs)>sum(cnt.perm_abs,2)))/cnt.nperm;
    pos(ifoi)   = 1-sum(sum(cnt.emp_pos)>sum(cnt.perm_pos,2))/cnt.nperm;
    neg(ifoi)   = 1-sum(sum(cnt.emp_neg)>sum(cnt.perm_neg,2))/cnt.nperm;
%     catch me
%     end
  end

  figure; hold on; set(gcf,'color','white');

  plot(log10(foi_range(1:NFREQ)),-log10(pos),'LineWidth',5,'color','r')
  plot(log10(foi_range(1:NFREQ)),-log10(neg),'LineWidth',5,'color','b')

  line([log10(2) log10(128)],-log10([0.05 0.05]),'color','k','LineWidth',2,'LineStyle','--')

  set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
  set(gca,'YTick',[0 1 2 3],'YTickLabel',[1 0.1 0.01 0.001])

  set(gca,'TickDir','out')
  set(gca,'YLim',[0 2.5]);

  title(sprintf('Connectivity differences: %s',str{v_out}))
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')

  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stat_indiv_conn_v%d.fig'],v_out),'fig')

%   close
  
end

%% PLOT INDIV
if plt
    str = {'Beamformer';'eLoreta'};

  v_out = 2;
  %
  NFREQ = 23;
  
  for ifoi = 1 : NFREQ
  %   ifoi
%     try
    load([outdir sprintf('nc_src_ana_all2all_stat_indiv_f%d_v%d.mat',ifoi,v_out)])

    pos(:,ifoi) = cnt.emp_pos;
    neg(:,ifoi) = cnt.emp_neg;
%     catch me
%     end
  end

  figure; hold on; set(gcf,'color','white');

  plot(log10(foi_range(1:NFREQ)),-log10(pos),'LineWidth',5,'color','r')
  plot(log10(foi_range(1:NFREQ)),-log10(neg),'LineWidth',5,'color','b')

  line([log10(2) log10(128)],-log10([0.05 0.05]),'color','k','LineWidth',2,'LineStyle','--')

  set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
  set(gca,'YTick',[0 1 2 3],'YTickLabel',[1 0.1 0.01 0.001])

  set(gca,'TickDir','out')
  set(gca,'YLim',[0 2.5]);

  title(sprintf('Connectivity differences: %s',str{v_out}))
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')

  saveas(gcf,sprintf([plotdir 'nc_src_all2all_stat_indiv_conn_v%d.fig'],v_out),'fig')

%   close
  
end