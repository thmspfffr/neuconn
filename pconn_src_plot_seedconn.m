%% PLOT SOURCE SPACE SEED RESULTS
% nc_src_plot_powcorr_seeds

clear all

m = 1;

% --------------------------------------------------------
% VERSION 1 - ELORETA
% --------------------------------------------------------
% v_in = 1; % 3: powcorr / 4: lpc
% v = 1;
% foi_range = unique(round(2.^[1:.25:7]));
% gridsize = 'coarse';
% --------------------------------------------------------
% VERSION 2 - BEAMFORMER
% --------------------------------------------------------
v_in = 3; % 3: powcorr / 4: lpc
v = 3;
foi_range = unique(round(2.^[1:.25:7]));
gridsize = 'coarse';
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath ~/pconn/matlab/

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir = '/home/tpfeffer/pconn/proc/conn/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

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

ifoi = 11;

load([outdir sprintf('pconn_src_allpowcorr_f%d_v%d.mat',ifoi,v)])

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

% s = nanstd(allpow,[],4);

%% PLOT

clear m th dat p xmean s x z para

thresh        = 1;
collapse_all  = 0;
diff = 1;
im = 3;


if collapse_all && ~diff
  m = nanmean(allpow,4);
  clear p
elseif ~collapse_all && ~diff
  m = nanmean(nanmean(allpow,4),3);
  clear p
else
  m = nanmean(allpow,4);
  clear p
end

if ~diff
for iroi  = 1 : 5


  if collapse_all
    
    xmean = squeeze(nanmean(allpow(iroi,:,im,:),2));
    s     = squeeze(nanstd(allpow(iroi,:,im,:),[],2));
  else
    xmean = squeeze(nanmean(nanmean(allpow(iroi,:,:,:),3),2));
    s     = squeeze(nanstd(nanmean(allpow(iroi,:,:,:),3),[],2));
  end
  for iloc = 1 : size(grid,1)

    iloc
    if collapse_all
      x = squeeze(allpow(iroi,iloc,im,:));
    else
      x = squeeze(nanmean(allpow(iroi,iloc,:,:),3));
    end
    z = (x - xmean)./s;

    if ~isnan(z)
      [~,p(iloc)] = ttest(z);
    else
      p(iloc)  = 1;
    end
  end
  
  if thresh == 1 
%     th = p<0.05;
    th = fdr(p',0.05);
  end
  
  if collapse_all
    dat = squeeze(m(iroi,:,im));
  else
    dat = squeeze(m(iroi,:));
  end
  
  if thresh == 1
    dat(~th) = 0;
  end
  
end

elseif diff == 1
  iroi = 1;
  dat = squeeze(nanmean(m(iroi,:,2),2))-squeeze(nanmean(m(iroi,:,1),2));
elseif diff == 2
  [dat, p] = ttest(squeeze(allpow(iroi,:,1,:)),squeeze(allpow(iroi,:,2,:)),'dim',2)';
end

para                  = [];
para.mydotmarkersize  = 40;
para.orientation      = 'axial';
para.colorlimits      = [0 max(dat)];
para.colormaps        = {'jet'};

h = figure; hold on
set(h,'color','w');

showmri_transp_v3(mri,para,[grid abs(dat')],grid(i(iroi),:));

saveas(gcf,sprintf([plotdir 'pconn_src_seedconn_f%d_r%d_t%d_v%d.fig'],ifoi,iroi,thresh,v),'fig')

% close(gcf)
close all
%% PLOT
v = 1;


clear i pow
% ROI: 1 = somatosens, 2: visual, 3: aud, 4: MPFC, 5: SMA
iroi  = 3;
NSUBJ = 12;

if ~exist([outdir sprintf('pconn_src_interhem_r%d_v%d.mat',iroi,v)])

  for ifoi = 1 : 23

    ifoi
    load([outdir sprintf('pconn_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
    allpow = single(allpow);

    if iroi == 1
      % somatosensory
      r(1,:) = [-42 -26 54]/10;
      r(2,:) = [42 -26 54]/10;
    elseif iroi == 2
      % visual cortex
      r(1,:) = [-20 -86 18]/10;
      r(2,:) = [20 -86 18]/10;
    elseif iroi == 3
      % auditory
      r(1,:) = [-54 -22 10]/10;
      r(2,:) = [54 -22 10]/10;   
    end

    for p = 1 : size(r,1)
      i(p)=findgridpoint(r(p,:),grid);
    end

    pow(:,:,ifoi) = squeeze(allpow(i(1),i(2),:,1:NSUBJ)); clear allpow

  end
  
  save([outdir sprintf('pconn_src_interhem_r%d_v%d.mat',iroi,v)],'pow')
end                      

%% PLOT
% PLOT INTERHEMISPHERIC CORRELATION BETWEEN TWO HOMOLOGUOUS AREAS
% FOR BOTH HCs AND PATs. 

p_ord = 1;
iroi = 3;

ord   = pconn_randomization;
ord   = ord([4 5 6 7 8 9 11 16 17 19 20 23],:);

titles  = {'Somatosensory';'Visual';'Auditory'};
h       = figure; hold on; set(h,'color','w')
col     = {[1 0 0];[0.1 0.6 0.8];[1 0.8 0.2]};

if ~p_ord
  for iroi = 1 : 3

    if exist([outdir sprintf('pconn_src_interhem_r%d_v%d.mat',iroi,v)])

      load([outdir sprintf('pconn_src_interhem_r%d_v%d.mat',iroi,v)])

      m = squeeze(nanmean(nanmean(pow,2)));
      s = squeeze(nanstd(nanmean(pow,1),[],2));
      s = s./sqrt(NSUBJ);

      u = m + s;
      l = m - s;

      X = [(1:23),fliplr(1:23)];
      Y = [u',fliplr(l')];

      hold on
      a{iroi} = plot(m,'color',col{iroi},'LineWidth',4,'LineStyle','-');
      b = fill(X,Y,col{iroi},'EdgeColor','none'); alpha(b,0.2)

      ylabel('Correlation');
      xlabel('Carrier frequency (in Hz)');
      title(sprintf('Interhemispheric connectivity'));

      % set(gca,'TickDir','out','XTick',foi_range,'XTickLabel',foi_range)
      set(gca,'TickDir','out','XTick',[1 3 7 11 15 19 23],'XTickLabel',[2 4 8 16 32 64 128])
    else
      error('File not found.');
    end

  end

  legend([a{1},a{2},a{3}],'Somatosensory','Visual','Auditory');   
  ylim([-0.01 0.25])
  saveas(gcf,sprintf([plotdir 'pconn_src_interhem_v%d.fig'],v),'fig')

else
  
  
  load([outdir sprintf('pconn_src_interhem_r%d_v%d.mat',iroi,v)])

  for im = 1 : 3       
    for isubj = 1 : length(ord)
      
      ord_pow(im,isubj,:) = pow(find(ord(isubj,:)==im),isubj,:);
      
    end      
 
    m = squeeze(nanmean(ord_pow(im,:,:),2));
    s = squeeze(nanstd(ord_pow(im,:,:),[],2));
    s = s./sqrt(NSUBJ);

    u = m + s;
    l = m - s;

    X = [(1:23),fliplr(1:23)];
    Y = [u',fliplr(l')];

    hold on
    a{im} = plot(m,'color',col{im},'LineWidth',4,'LineStyle','-');
    b = fill(X,Y,col{im},'EdgeColor','none'); alpha(b,0.2)
    
  end
	ylabel('Correlation');
	xlabel('Carrier frequency (in Hz)');
	title(sprintf('Interhemispheric connectivity -- %s',titles{iroi}));

	% set(gca,'TickDir','out','XTick',foi_range,'XTickLabel',foi_range)
	set(gca,'TickDir','out','XTick',[1 3 7 11 15 19 23],'XTickLabel',[2 4 8 16 32 64 128])
  legend([a{1},a{2},a{3}],'Placebo','Atomoxetine','Donepezil');   

end

