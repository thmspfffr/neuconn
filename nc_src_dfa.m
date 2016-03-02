%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_src_dfa

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 2;
% --------------------------------------------------------

% restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/gnolte/neuconn/matlab/rest/
addpath /home/gnolte/neuconn/meg_out/rest/dti/
addpath(genpath('/home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/'))


indir   = '/home/gnolte/neuconn/meg_out/rest/src/';
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';

% call subj info startup function
% allsubj = ed;


v_filt = 2;

siginfo = nbt_Info;
siginfo.converted_sample_frequency = 300;

foi = [1 4; 4 8; 8 13; 13 30];

%%

for ifoi = 1:length(foi)
  for icond = 1 : 2
    for isubj = 1 : 49
      
      if ~exist(sprintf([outdir 'nc_src_dfa_f%d_c%d_s%d_v%d_processing.txt'],ifoi,icond,isubj,v))
        system(['touch ' outdir sprintf('nc_src_dfa_f%d_c%d_s%d_v%d_processing.txt',ifoi,icond,isubj,v)]);
      else
        continue
      end
      
      if ifoi <= 30
        freq = 1;
      else
        freq = 2;
      end
      
      try
        % Load eloreta filter, leadfield and grid
        load([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v_filt)],'A','L','grid');
        
        clear L grid
        
        
        for iblock = 1 : 2
          
          
          disp(sprintf('Computing dipole for f%d c%d s%d b%d ...',ifoi,icond,isubj,iblock))
          
          disp(sprintf('Loading MEG data ...'));
          load(sprintf('/home/gnolte/neuconn/meg_out/rest/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6));
          
          if freq == 1
            clear data_hi
            [mydata,epleng] = megdata2mydata(data_low);
            clear data_low
          elseif freq == 2
            clear data_low
            [mydata,epleng] = megdata2mydata(data_hi);
            clear data_hi
          else
            error('Missing information on frequency!')
          end
          
          % load cross spectrum
          load([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_filt)]);
          
          F = getdipdir(cs(:,:,ifoi),A); clear cs
          
          for iloc = 1 : size(F,2)
          
            dat = mydata*F(:,iloc);
              
         	% computing amplitude envelopes
            ampenv = nbt_GetAmplitudeEnvelope(dat,siginfo,min(foi(ifoi,:)),max(foi(ifoi,:)),2/min(foi(ifoi,:)));
            
        	% compute DFA
            [tmp,expo(iblock,iloc)] = nbt_doDFA(ampenv, siginfo, [3 100],[1 150],0.5,0,0,[]);
            
          end
          
          
        end
        
        expo = nanmean(expo,1);
        
        save(sprintf([outdir 'nc_src_dfa_f%d_c%d_s%d_v%d.mat'],ifoi,icond,isubj,v),'expo','-v7.3');
        
        clear A
      catch me
        save(sprintf([outdir 'nc_src_dfa_f%d_c%d_s%d_v%d_error.mat'],ifoi,icond,isubj,v),'me','-v7.3');
      end
      
    end
  end
end

error('!')


%% PLOT DFA EXPONENTS
NSUBJ = 34;

plt = 1;
if plt
  
  clear expo_all
  
  load sa_meg_template;
  grid = sa_meg_template.grid_medium;
  mri  = sa_meg_template.mri;
  vc = sa_meg_template.vc;
  clear sa_meg_template
  
  load('/home/gnolte/neuconn/meg_out/rest/pasat/perf_ms.mat');
  
  for ifoi = 3 : 3
    for icond = 1 : 2
      cnt = 0;
      for isubj = 1 : 49
        
        
        try
          load(sprintf([outdir 'nc_src_dfa_f%d_c%d_s%d_v%d.mat'],ifoi,icond,isubj,v));
          cnt = cnt + 1;
          expo_all(:,cnt,icond) = expo;
          if icond == 1
            %           perf_all(:,isubj) = perf(:,isubj);
          end
        end
        
      end
    end
  end
  
  % perf_all(perf_all==0)=nan;
  
  % mask = zeros(5003,1);
  % p = zeros(5003,1);
  %
  % for iloc = 1 : 5003
  %
  %   iloc
  %   [mask(iloc) p(iloc)] = ttest(squeeze(expo_all(iloc,:,1)),squeeze(expo_all(iloc,:,2)));
  %
  % end
  %
  % plt = squeeze(nanmean(expo_all(:,:,1),2))-squeeze(nanmean(expo_all(:,:,1),2));
  % plt(~mask)=0;
  % plt(find(mask))=1;
  
  % vout = spatfiltergauss(plt,grid,.5,sa_meg_template.cortex10K.vc);
  %
  expo_all = expo_all(:,1:NSUBJ,:);
  %%
  para = [];
  para.mydotmarkersize =20;
  para.orientation     ='coronal';
  para.colormaps       = {'jet'};
  para.showcenterline  = 0;
  para.nsubplotrows     = 1;
  
  for icond = 1  : 2
    
    para.plt(:,icond) = squeeze(nanmean(expo_all(:,:,icond),2));
    
    %   para.colorlimits = [0.6 .9];
    % para.colormaps = 'cool';
    
    h=figure; hold on
    set(h,'color','w');
    
    subplot(2,1,icond)
    % showsurface(sa_meg_template.cortex10K,para,vout);
    
    % degree = squeeze(c(:,ifreq,2))/5003;
    % subplot(2,1,2)
    showmri_transp_v3(mri,para,[grid para.plt(:,icond)]);
  end
  
  h=figure; hold on
  set(h,'color','w');
  showmri_transp_v3(mri,para,[grid para.plt(:,1)-para.plt(:,2)]);
  
  % plot
  
  a = squeeze(nanmean(expo_all,1));
  e = nanstd(a)/sqrt(16);
  
  figure; hold on;
  
  plot([1 1.1],nanmean(a),'k.','MarkerSize',10)
  plot([1],nanmean(a(:,1))+e(1),'.')
  plot([1],nanmean(a(:,1))-e(1),'.')
  plot([1.1],nanmean(a(:,2))+e(2),'.')
  plot([1.1],nanmean(a(:,2))-e(2),'.')
  
  plot([1 1.1],a,'kx')
  axis([0.95 1.12 0.5 1])
  
end

%% STATISTICS
% permutation test

expo_mean = squeeze(nanmean(expo_all,1));
allexps   = expo_mean(:);

NPERM = 10000;

for iperm = 1 : NPERM
  
  disp(sprintf('Permutation %d ...',iperm));
  
  permdat         = allexps(randperm(16*2));
  permexp         = [permdat(1:16) permdat(17:end)];
  permdiff(iperm) = diff(nanmean(permexp));
  
end

figure; hold on;
hist(permdiff);

prc = [prctile(permdiff,5) prctile(permdiff,95)];
line([diff(nanmean(expo_mean)) diff(nanmean(expo_mean))],[0 3000],'LineWidth',5,'color','r','LineStyle','--')
line([prc(1) prc(1)],[0 3000],'LineWidth',2,'color','y','LineStyle','-')
line([prc(2) prc(2)],[0 3000],'LineWidth',2,'color','y','LineStyle','-')
