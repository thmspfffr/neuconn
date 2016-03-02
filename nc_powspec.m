%% COMPUTE POWER SPECTRUM FOR ALL SENSORS
% fits a straight line to the powerspectrum in order
% to identify peaks in the spectrum for subsequent
% DFA analysis.

% last update: 18-02-2015, tpfeffer

% to be implemented:
% ***nothing left to implement***

% execute as
% --------------------------------------------------------
% pconn_powspec
% --------------------------------------------------------

clear all

% --------------------------------------------------------
% VERSION 1 - all freq, freqoi 2-200
% --------------------------------------------------------
v         = 1;
d         = 'data_low.trial{1}+data_hi.trial{1};';
v_rawdata = 3;
fsample   = 300;
FOI       = 2:1:200;
T         = 0.5;
% --------------------------------------------------------

SUBJLIST{1} = [1 2 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 24 25 28 30 32 33 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49];
SUBJLIST{2} = [2 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 28 29 30 32 33 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49];

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/fieldtrip-20150215/
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/

ft_defaults

indir    = '/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/';
outdir   = '/home/gnolte/neuconn/meg_out/rest/sens/';
plotdir  = '/home/gnolte/neuconn/meg_out/rest/plots/';

%%

for icond = 1 : 2
  
  if ~exist(sprintf([outdir 'nc_powspec_c%d_v%d_processing.txt'],icond,v))
    system(['touch ' outdir sprintf('nc_powspec_c%d_v%d_processing.txt',icond,v)]);
  else
    continue
  end
  
  if ~exist([outdir sprintf('nc_powspec_m%d_v%d.mat',icond,v)])
    cnt = 0;
    
    for isubj = SUBJLIST{icond}
      
      disp(sprintf('Processing s%d c%d ...', isubj,icond))
      
      cnt = cnt + 1;
      
      for iblock = 1 : 2
        
        disp(sprintf('Loading MEG data ...'));
                
        load(sprintf('/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6));
        data_low.trial{1} = eval(d);
        
        cfg = [];
        cfg.length = 1/T;
        cfg.overlap = 0;
        
        data      = ft_redefinetrial(cfg,data_low); clear data_low


        cfg = [];
        
        cfg.trials  = 1:size(data.trial,2)-1;         
        data        = ft_redefinetrial(cfg,data);
          
        FOI = 1 : 1 : 40;
%         
%         cfg             = [];
%         cfg.method      = 'template';
%         cfg.layout      = 'CTF275';
%         neighbours       = ft_prepare_neighbours(cfg);
% 
%       	cfg              = [];
%         cfg.feedback     = 'no';
%         cfg.method       = 'template';
%       	cfg.planarmethod = 'sincos';
%        	cfg.channel      = {'MEG'};
%       	cfg.neighbours   = neighbours;
%       	data             = ft_megplanar(cfg, data);
%         
      	cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.output      = 'pow';
        cfg.taper       = 'hanning';
        cfg.channel     = {'MEG'};
        cfg.foi         = 0.1:0.1:100;
        cfg.keeptrials  = 'yes';
        dat             = ft_freqanalysis(cfg, data); % run freqanalysis here
%             
%         cfg = [];
%         dat = ft_combineplanar(cfg,dat);
%         
        tp(:,:,iblock) = squeeze(nanmean(dat.powspctrm,1));
        
      end
      if ndims(tp) == 3
        power(:,:,cnt) = nanmean(tp,3); clear tp
      elseif ndims(tp) == 2
        power(:,:,cnt) = tp; clear tp
      end
      
    end
    save([outdir sprintf('nc_powspec_c%d_v%d.mat',icond,v)],'power');
  end
  % save preliminary result
end

error('Stop here!')

%% Timescale: Pupil dilation analysis
% Loads preprocessed .mat files and plots processed data

clear
% -------------------------------------------------------------------------
% VERSION 01 - STIMULUS LOCKED
% -------------------------------------------------------------------------
v = 1;
fsample     = 400; % original sampling rate
lpfilt      = 5;
comp = 'dfa';
% -------------------------------------------------------------------------

restoredefaultpath

addpath('/home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/')
addpath /home/tpfeffer/pconn/matlab/

outdir      = '/home/tpfeffer/pconn/proc/pup/';
sampledir   = '/home/tpfeffer/pconn/proc/pup/';
eventdir    = '/home/tpfeffer/pconn/proc/pup/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

%% POWER SPECTRUM

% if strcmp(comp,'pow')
powspctrm = nan(20001,2,24,3);

for m = 1 : 3
  for isubj = 1 : 24
    
    d=dir([outdir sprintf('pconn_pup_preproc_s%d_m%d_b*_v%d.mat',isubj,m,v)]);
    
    if length(d)<2
      
      fprintf('\nNo s%d m%d ... \n',isubj,m)
      powspctrm(:,:,isubj,m) = NaN;
      powspec_subjlist(isubj,m) = 0;
      continue
      
    else
      
      powspec_subjlist(isubj,m) = 1;
      for iblock = 1 : 2
        
        fprintf('\nProcessing s%d b%d m%d ... \n',isubj,iblock,m)
        
        load([outdir d(iblock).name])
        
        if sum(isnan(pup.dil))>0
          warning('Still NaNs in the Data!')
          idx = find(isnan(pup.dil));
          warning(sprintf('Deleting NaNs at pos %d of %d',idx(1),size(pup.dil,1)))
          pup.dil(idx) = [];
        end
        
        nwin = 40000;
        w = hanning(nwin);
        
%         if fsample/size(w,1)>0.04; error('Error!'); end
       
        [p(:,iblock),f(:,iblock)]=pwelch(pup.dil,w,[],nwin,fsample,'power');
        
      end
      
      powspctrm(:,:,isubj,m) = p(1:end,:); clear p
      
      
    end
  end
end

save(sprintf([outdir 'pconn_pup_powspctrm_v%d.mat'],v),'powspctrm','powspec_subjlist','f','-v7.3')

error('Done.')

%%
v = 4;

col = {'r';'b';'k'}
ord       = pconn_randomization;
load(sprintf([outdir 'pconn_pup_powspctrm_v%d.mat'],v))

NCOND = 1;
% SUBJLIST = find(powspec_subjlist(:,1:NCOND));
SUBJLIST  = [4 5 6 7 9 10 11 12 13 15 16 19 21 24];

ses_col = {'k';'r';'b'};
blo_col = {'r','b'};

if ~exist(sprintf([outdir 'pconn_pup_powspctrm_v%d.mat'],v))
  error('No input data found!')
end

h=figure; set(gcf,'color','white'); hold on
 cnt = 0; 
for subj = 1 : length(SUBJLIST)
  
  isubj = SUBJLIST(subj);
  cnt = cnt + 1;
  subplot(4,4,cnt); hold on
  for m = 1 : 3
      
  	im = find(ord(isubj,:)==m);
    idx   = find(f(:,1)>1,1,'first');
    idx10 = find(f(:,1)==1,1,'first');
  
    powspctrm1 = squeeze(nanmean(powspctrm,2));
    
    
    
      X = [ones(idx10-2,1) log10(f(3:idx10,1))];
      Y = log10(squeeze(powspctrm1(3:idx10,isubj,im)));

      t = X\Y;
      reg(m,cnt) = t(2);

      k=plot(log10(f(3:idx10,1)),log10(squeeze(powspctrm1(3:idx10,isubj,im))),'color',col{m},'LineWidth',1.5);
      plot(log10(f(3:idx10,1)),t(2)*log10(f(3:idx10,1))+t(1),'LineWidth',3,'color',col{m})

      box off; 
        


   
    if any(cnt == 1 : 4 : 20); 
      ylabel('power'); 
    else
      set(gca,'YTick',[]);
    end 
    cnt
    if any(cnt == 17 : 1 : 20); 
      xlabel('Frequency (Hz)'); 
    else
%       set(gca,'XTick',[]);
    end  
    box on
  end
         title(sprintf('s%d: %0.2f / %0.2f / %0.2f',isubj,reg(1,cnt),reg(2,cnt),reg(3,cnt)));

end

saveas(gcf,sprintf('~/pconn/proc/plots/pconn_pup_powspctrm_allsubj_v%d.fig',v),'fig')
 
ms1   = nanmean(reg,2);
ss1   = nanstd(reg,[],2)/sqrt(size(reg,2));

[~,p] = ttest(reg(1,:),reg(2,:));

subplot(4,4,cnt+1); set(gca,'visible','off');
title('Average Slopes:')
plot(ms1,'o','MarkerSize',7,'MarkerFace','b');
line([1 1],[ms1(1)-ss1(1) ms1(1)+ss1(1)],'LineWidth',3)
line([2 2],[ms1(2)-ss1(2) ms1(2)+ss1(2)],'LineWidth',3)
line([3 3],[ms1(3)-ss1(3) ms1(3)+ss1(3)],'LineWidth',3)
axis([0.8 3.2 min(ms1)-0.1 max(ms1)+0.1]); grid on;
set(gca,'XTick',[1 2 3],'XTickLabel',{'Pla';'Ato';'Don'})




 
clear reg
subplot(4,4,cnt+2);  hold on

for i = 1 : 3
  
  X        = [ones(idx10-2,1) log10(f(3:idx10,1))];
  Y        = log10(squeeze(nanmean(powspctrm1(3:idx10,:,i),2)));
  t(:,i)   = X\Y; 
  reg(i)   = t(2,i);
  
  plot(log10(f(3:idx10,1)),log10(squeeze(nanmean(powspctrm1(3:idx10,:,i),2))),'LineWidth',3,'color',col{i});
  plot(log10(f(3:idx10,1)),t(2,i)*log10(f(3:idx10,1))+t(1,i),'LineWidth',3,'color',col{i});

end

set(gca,'XTick',log10([0.01 0.1 1]),'XTickLabel',[0.01 0.1 1]);

% set(gca,'visible','off');
% 
% tx = sprintf('Mean values:\nPlacebo: %.2f \nAtomoxetin: %.2f \nDonepezil: %.2f',ms1(1),ms1(2),ms1(3))
% tx = [tx sprintf('\n\nPlac vs. Atomox: p = %.3f',p)];
% text(0,0.5,tx)
% set(h,'Position',[40 78 1060 720]);
% 
% saveas(gcf,sprintf('~/pconn/proc/plots/pconn_pup_powspec_allsubj_v%d.fig',v),'fig')


% % elseif strcmp(comp,'dfa')





