%% PLOT SOURCE SPACE SEED RESULTS
% nc_src_plot_powcorr_seeds

clear all

m = 1;

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v_in = 7; % 3: powcorr / 4: lpc
foi_range = unique(round(2.^[1:.25:7]));
gridsize = 'coarse';
NSUBJ = 22;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest/
addpath /home/gnolte/neuconn/meg_out/rest/dti/

ft_defaults

indir   = '/home/gnolte/neuconn/meg_out/rest/src/';
outdir = '/home/gnolte/neuconn/meg_out/rest/conn/';
mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';

% call subj info startup function
allsubj = nc_allsubj_start(m);

%%
load sa_meg_template;
if strcmp(gridsize,'medium')
  grid = sa_meg_template.grid_medium;
elseif strcmp(gridsize,'coarse')
  grid = sa_meg_template.grid_coarse;
elseif strcmp(gridsize,'cortex')
  grid = sa_meg_template.grid_cortex3000;
end

mri  = sa_meg_template.mri; 
vc = sa_meg_template.vc;

clear sa_meg_template;
% allpow = zeros(5003,5003,2,16);

ifoi = 11;

load([outdir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v_in)])

% somatosensory
r(1,:) = [-42 -26 54]/10;
% visual cortex
r(2,:) = [-20 -86 18]/10;
% auditory
r(3,:) = [-54 -22 10]/10;
% MPFC
r(4,:) = [-3 39 -2]/10;
% SMA
r(5,:) = [-2 1 51]/10;

for p = 1 : length(r)
  i(p)=findgridpoint(r(p,:),grid);
end

allpow = allpow(i,:,:,:);

m = nanmean(allpow,4);
% s = nanstd(allpow,[],4);

%% PLOT

% ROI: 1 = somatosens, 2: visual, 3: aud, 4: MPFC, 5: SMA
iroi  = 2;
icond = 1;
diff = 0;

if ~diff
  xmean = squeeze(nanmean(allpow(iroi,:,icond,1:NSUBJ),2));
  s     = squeeze(nanstd(allpow(iroi,:,icond,1:NSUBJ),[],2));

  for iloc = 1 : size(grid,1)

    iloc
    x = squeeze(allpow(iroi,iloc,icond,1:NSUBJ));

    z = (x - xmean)./s;

    if ~isnan(z)
      [~,p(iloc)] = ttest(z,0,'tail','right');
    else
      p(iloc)  = 1;
    end
  end

%   th = p<0.05;
  th = fdr(p,0.05);

  dat = squeeze(m(iroi,:,icond));
  dat(~th) = 0;
  
end

if diff == 1
  dat = squeeze(m(iroi,:,1))-squeeze(m(iroi,:,2));
elseif diff == 2
  [dat, p] = ttest(squeeze(allpow(iroi,:,1,:)),squeeze(allpow(iroi,:,2,:)),'dim',2)';
end

para                  = [];
para.mydotmarkersize  = 40;
para.orientation      = 'axial';
para.colorlimits      = [0 0.08];
para.colormaps        = {'jet'};

h = figure; hold on
set(h,'color','w');

showmri_transp_v3(mri,para,[grid dat'],grid(i(iroi),:));


%% PLOT
v = 5;

tc = 0;
if tc
clear i pow
% ROI: 1 = somatosens, 2: visual, 3: aud, 4: MPFC, 5: SMA
iroi  = 5;
NSUBJ = 24;

if ~exist([outdir sprintf('nc_src_interhem_r%d_v%d.mat',iroi,v)])

  for ifoi = 1 : 23

    ifoi
    load([outdir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
    allpow = single(allpow);

    if iroi == 1
      % somatosensory
      r(1,:) = [-42 -26 54]/10;
      r(2,:) = [42 -26 54]/10;
    elseif iroi == 2
      % visual cortex
      r(1,:) = [-20 -86 18]/10;
      r(2,:) = [20 -86 18]/10;
    else
      % auditory
      r(1,:) = [-54 -22 10]/10;
      r(2,:) = [54 -22 10]/10;
    end

    for p = 1 : size(r,1)
      i(p)=findgridpoint(r(p,:),grid);
    end

    pow(:,:,ifoi) = squeeze(allpow(i(1),i(2),:,1:NSUBJ)); clear allpow

  end
  
  save([outdir sprintf('nc_src_interhem_r%d_v%d.mat',iroi,v)],'pow')
                       
end
%% PLOT
% PLOT INTERHEMISPHERIC CORRELATION BETWEEN TWO HOMOLOGUOUS AREAS
% FOR BOTH HCs AND PATs. 


titles = {'Somatosensory';'Visual';'Auditory'};

iroi = 1;

% figure 

hold on

if exist([outdir sprintf('nc_src_interhem_r%d_v%d.mat',iroi,v)])
  load([outdir sprintf('nc_src_interhem_r%d_v%d.mat',iroi,v)])

  col = {'r';'b'};

  for icond = 1 : 2
    plot(squeeze(nanmean(pow(icond,:,:),2)),col{icond},'LineWidth',4,'LineStyle','-')
    hold on
  end

  ylabel('Correlation');
  xlabel('Carrier frequency (in Hz)');
  title(sprintf('Interhemispheric connectivity: %s',titles{iroi}));

  % set(gca,'TickDir','out','XTick',foi_range,'XTickLabel',foi_range)
  set(gca,'TickDir','out','XTick',[1 3 7 11 15 19 23],'XTickLabel',[2 4 8 16 32 64 128])
else
  error('File not found.');
end





%%
xmean = squeeze(nanmean(allpow(iroi,:,icond,:),2));
s     = squeeze(nanstd(allpow(iroi,:,icond,:),[],2));

for iloc = 1 : 5003
  
  iloc
  x = squeeze(allpow(iroi,iloc,icond,:));
  
  z = (x - xmean)./s;
  
  if ~isnan(z)
    [~,p(iloc)] = ttest(z);
  else
    p(iloc)  = 1;
  end
end

th = p<0.05;
% th = fdr(p,0.05);

dat = squeeze(m(iroi,:,icond));
dat(~th) = 0;

para                  = [];
para.mydotmarkersize  = 80;
para.orientation      = 'coronal';
para.colorlimits      = [min(dat) max(dat)];

h=figure; hold on
set(h,'color','w');

showmri_transp_v3(mri,para,[grid dat'],grid(i(1:2),:));


end



