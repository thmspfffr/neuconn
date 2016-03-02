%% PLOT SOURCE SPACE SEED RESULTS
% nc_src_plot_powcorr_whole

clear all

m = 1;

% --------------------------------------------------------
% VERSION 7
% --------------------------------------------------------
v_in      = 10; 
v         = 10;
foi_range = unique(round(2.^[1:.25:7]));
gridsize  ='coarse';
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

ifoi = 16;

load([outdir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])

a = ttest(squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)),'dim',3); 

b = nansum(a);

b = b/max(b);
%% PLOT

% ROI: 1 = somatosens, 2: visual, 3: aud, 4: MPFC, 5: SMA
iroi  = 1;
icond = 2;
diff = 0;

if ~diff
  
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
para.colorlimits      = [0 1];
para.colormaps        = {'jet'};

h = figure; hold on
set(h,'color','w');

showmri_transp_v3(mri,para,[grid dat'],grid);


%% PLOT
v = 5;

tc = 0;
if tc
clear i pow
% ROI: 1 = somatosens, 2: visual, 3: aud, 4: MPFC, 5: SMA
iroi  = 1;
NSUBJ = 22;

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


titles = {'Somatosensory';'Visual';'Auditory';'MPFC';'SMA'};

iroi = 5;

% figure 

hold on

if exist([outdir sprintf('nc_src_interhem_r%d_v%d.mat',iroi,v)])
  
  load([outdir sprintf('nc_src_interhem_r%d_v%d.mat',iroi,v)])
  s = squeeze(nanstd(pow,[],2));
  e = s./sqrt(22);
  
  col = {'r';'b'};
  
m = squeeze(nanmean(pow(:,:,:),2));


  for icond = 1 : 2
    
    
    plot(m(icond,:),col{icond},'LineWidth',4,'LineStyle','-')
    hold on
    
    X = [(1:23),fliplr(1:23)];
    Y = [m(1,:)-e(1,:),fliplr(m(1,:)+e(1,:))];
    b = fill(X',Y,[1 .55 0.2],'EdgeColor','none'); alpha(b,0.2)

    X = [(1:23),fliplr(1:23)];
    Y = [m(2,:)-e(2,:),fliplr(m(2,:)+e(2,:))];
    c = fill(X',Y,[.11 .56 1],'EdgeColor','none'); alpha(c,0.2)

%     set(gca,'xlim',[0 25.101])
    set(gca,'ylim',[-0.0101 0.1002])
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



