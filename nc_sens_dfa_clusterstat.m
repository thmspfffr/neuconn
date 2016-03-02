%% COMPUTE DFA SENSOR LEVEL TOPO STATISTICS
% pconn_sens_dfa_clusterstat

% Implement statistical test as described in Nichols & Holmes, 2001

% (1) Single threshold permutation test
% (2) Cluster-based permutation test

% tpfeffer | thms.pfffr@gmail.com | 05-05-15

clear

% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
v         = 7;
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 17 19 20 21 22 23 24];
% --------------------------------------------------------

addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/

outdir   = '/home/gnolte/neuconn/meg_out/rest/dfa/';

subsel = 1:length(SUBJLIST);
clear dfa_all idx


%% CLUSTER-BASED PERMUTATION TEST

load('/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/nc_postproc_c1_s49_b1_v6.mat');
%%
clear data_hi
% -------------------------------------------------------
% GET RID OF MISSING SENSORS
% -------------------------------------------------------
cfg             = [];
cfg.method      = 'template';
cfg.layout      = 'CTF275';
n       = ft_prepare_neighbours(cfg);
load ~/pconn/matlab/pconn_idxlab.mat
clear dfa_all

for ifoi = 1 : 6
  
  
  for icond = 1 : 2
    
    subj_cnt{icond} = [];
    
    for isubj = 1:49
      
      d = dir(sprintf([outdir 'nc_sens_dfa_s%d_b*_c%d_f%d_v%d.mat'],isubj,icond,ifoi,v));
      
      if length(d)==0
        continue
      elseif length(d)==1
        if isubj < 10
          blocks = str2num(d.name(17));
        else
          blocks = str2num(d.name(18));
        end
        subj_cnt{icond} = [subj_cnt{icond} isubj];
      elseif length(d)==2
        if isubj < 10
          blocks = [str2num(d(1).name(17)) str2num(d(2).name(17))];
        else
          blocks = [str2num(d(1).name(18)) str2num(d(2).name(18))];
        end
        subj_cnt{icond} = [subj_cnt{icond} isubj];
      end
      
      for iblock = blocks
        
        load(sprintf([outdir 'nc_sens_dfa_s%d_b%d_c%d_f%d_v%d.mat'],isubj,iblock,icond,ifoi,v));
        
        dfa_all(:,isubj,icond,ifoi,iblock) = par.dfa(idx_lab); clear par
        
      end
    end
  end
end

%%

dfa_all = squeeze(dfa_all(:,intersect(subj_cnt{1},subj_cnt{2}),:,:,:));

dfa_f = squeeze(nanmean(nanmean(dfa_all(:,:,:,3,:),4),5));

dat     = [dfa_f(:,:,2) dfa_f(:,:,1)];

data_low.dimord                  = 'subj_chan_freq_time';
data_low.time                   = [1 2];
data_low.freq                   = [1 2];
data_low.powspctrm(:,:,1:2,1:2) = repmat(permute(dat,[2 1]),[1 1 2 2]);

cfg                  = [];
cfg.channel          = 'all';
cfg.latency          = [1 1];
cfg.frequency        = [1 1];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.computeprob      = 'yes';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025;
cfg.clusterstatistic = 'maxsum';
cfg.clustertail      = 0; %1 = right
cfg.tail             = 0; %1 = right
cfg.alpha            = 0.025;
cfg.minnbchan        = 2;
cfg.numrandomization = 10000;
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';

% specifies with which sensors other sensors can form clusters
cfg_neighb.method           = 'template';
cfg_neighb.template         = 'CTF275_neighb.mat';
cfg_neighb.feedback         = 'no';
cfg.neighbours = n;

n_subj = size(dfa_all,2);
design = zeros(2,2*n_subj);
design(1,:) = repmat(1:n_subj,1,2);
design(2,:) = mod(floor([0:(2*n_subj-1)]/(n_subj/1)),2)+1;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stats] = ft_freqstatistics(cfg, data_low);

%%

dfa_f = squeeze(nanmean(nanmean(dfa_all,2),5));
par = dfa_f;

par1 = squeeze(nanmean(dfa_all(:,:,1,:,:),5));
par2 = squeeze(nanmean(dfa_all(:,:,2,:,:),5));

cnt_pos = zeros(268,38,16);
cnt_neg = zeros(268,38,16);

for ifoi = 1 : 16
  for isens = 1 : 268
    
    x_ms = par1(isens,:,ifoi);
    hc_thresh = [prctile(par2(isens,:,ifoi),5) prctile(par2(isens,:,ifoi),95)];
    
    cnt_pos(isens,x_ms>hc_thresh(2),ifoi) = cnt_pos(isens,x_ms>hc_thresh(2),ifoi) + 1;
    cnt_neg(isens,x_ms<hc_thresh(1),ifoi) = cnt_neg(isens,x_ms<hc_thresh(1),ifoi) + 1;
    
  end
  
end

pos = squeeze(sum(sum(cnt_pos,2),1))
neg = squeeze(sum(sum(cnt_neg,2),1))


cnt_pos_perm = zeros(268,38,16,NPERM);
cnt_neg_perm = zeros(268,38,16,NPERM);

par_all	= squeeze(cat(2,par1,par2));


for iperm = 1 : NPERM
  
  
  permidx	= randperm(NSUBJ*2);
  
  par1_perm = par_all(:,permidx(1:NSUBJ),:);
  par2_perm = par_all(:,permidx(NSUBJ+1:end),:);
  
  for ifoi = 1 : 16
    disp(sprintf('Processing freq #%d / perm #%d ...',ifoi,iperm))
    for isens = 1 : 268
      
      x_ms = par1_perm(isens,:,ifoi);
      hc_thresh = [prctile(par2_perm(isens,:,ifoi),5) prctile(par2_perm(isens,:,ifoi),95)];
      
      cnt_pos_perm(isens,x_ms>hc_thresh(2),ifoi,iperm) = cnt_pos_perm(isens,x_ms>hc_thresh(2),ifoi,iperm) + 1;
      cnt_neg_perm(isens,x_ms<hc_thresh(1),ifoi,iperm) = cnt_neg_perm(isens,x_ms<hc_thresh(1),ifoi,iperm) + 1;
    end
  end
end

pos_perm = squeeze(sum(sum(cnt_pos_perm,2),1))
neg_perm = squeeze(sum(sum(cnt_neg_perm,2),1))



%       par = squeeze(nanmean(dfa_all(:,isubj,:),2));
%

%     par = abs(log10(tpdf(stats.stat,37)));

%



%% PLOT DFA AVERAGED OVER SPACE

NFOI = 6;


s = squeeze(nanstd(squeeze(nanmean(nanmean(dfa_all,5),1)),[],1));
m = squeeze(nanmean(squeeze(nanmean(nanmean(dfa_all,5),1)),1));

figure; set(gcf,'color','white'); hold on
bar(2:2:12,m(1,:),'r','barwidth',0.4,'edgecolor','w')
bar(2.5:2:12.5,m(2,:),'b','barwidth',0.4,'edgecolor','w')

ll = m-(s./sqrt(38));
ul = m+(s./sqrt(38));

cnt = 0;
for i = 2 : 2 : 2*NFOI
  cnt = cnt + 1;
  line([i i],[ll(1,cnt) ul(1,cnt)],'color','r','linewidth',3)
  line([i+0.5 i+0.5],[ll(2,cnt) ul(2,cnt)],'color','b','linewidth',3)
end

xlabel('Frequency [Hz]');ylabel('DFA Exponent'); axis([0 14 0.5 0.8]); title('Detrended Fluctuation Analysis')
set(gca,'TickDir','out','XTick',[1:2:12],'XTickLabel',{'2-4';'4-7';'8-13';'13-30';'30-60';'60-120'});

par = squeeze(squeeze(nanmean(nanmean(dfa_all,5),1)));

for i = 1 : NFOI
  
  [~,p(i)]=ttest(par(:,1,i),par(:,2,i));
  
end

print(gcf,'-depsc2',sprintf('/home/gnolte/neuconn/meg_out/rest/plots/nc_sens_dfa_wholebrain_v%d.eps',v))


%% PLOT TOPOGRAPHIES
scale = [0.5 0.7; 0.5 0.75; 0.6 0.8; 0.6 0.75; 0.5 0.65; 0.5 0.6];
ifoi = 1;

par = squeeze(squeeze(nanmean(nanmean(dfa_all,5),2)));

pars.cbar=0;
pars.linewidth = 9;
cfg.markersize = 0;
pars.scale = scale(ifoi,:);

pars.resolution = 300;
pars.markersize = 0;

%     subplot(1,2,1)
showfield_colormap(par(:,1,ifoi),sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg100',sprintf('/home/gnolte/neuconn/meg_out/rest/plots/nc_sens_dfa_topo_c%d_f%d_v%d.jpg',1,ifoi,v))

figure; set(gcf,'color','white');

showfield_colormap(par(:,2,ifoi),sa.locs_2D,pars);
colormap(hot)

print(gcf,'-djpeg100',sprintf('/home/gnolte/neuconn/meg_out/rest/plots/nc_sens_dfa_topo_c%d_f%d_v%d.jpg',2,ifoi,v))

figure; set(gcf,'color','white');
pars.scale = [-0.05 0.05];
showfield_colormap(par(:,2,ifoi)-par(:,1,ifoi),sa.locs_2D,pars);
%   colormap(hot)

print(gcf,'-djpeg100',sprintf('/home/gnolte/neuconn/meg_out/rest/plots/nc_sens_dfa_topo_diff_f%d_v%d.jpg',ifoi,v))

close all
%        subplot(1,2,2)
%     showfield_colormap(par(:,1),sa.locs_2D,pars);
%       colormap(hot)