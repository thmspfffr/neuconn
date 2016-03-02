%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_src_powcorr_degree_surface

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 2;
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
allsubj = nc_allsubj_start;

%%
% load sa_meg_template;
% grid = sa_meg_template.grid_medium;

allpow = zeros(5003,5003,2,16);

for ifoi = 1:64
  for icond = 1 : 2
    
    if ~exist(sprintf([outdir 'nc_powcorr_threshold_c%d_f%d_v%d_processing.txt'],icond,ifoi,v))
      system(['touch ' outdir sprintf('nc_powcorr_threshold_c%d_f%d_v%d_processing.txt',icond,ifoi,v)]);
    else
      continue
    end
    
    disp(sprintf('Computing f%d c%d ...',ifoi,icond))
    
    cnt = 0;
    for isubj = 1 : 40
%       disp(sprintf('Processing subject %d ...',isubj))
      % exclude excluded subjects
      if cell2mat(allsubj{icond}(isubj,3)) == 0
        continue
      end
      % average the two blocks
      for iblock = 1 : 2
        if exist([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
          load([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
          br = 0;
          % correction factor for correlations and fisher z-transform
          res(:,:,iblock) = (resout+resout')./2; clear resout
        else
          br = 1;
          break
        end
      end
      % break & continue with next subject
      if br
        continue
      end
      disp(sprintf('Processig s%df%dc%d',isubj,ifoi,icond));
      % count only if subject data exists
      cnt = cnt + 1;
      allpow(:,:,icond,cnt)	= nanmean(res,3); clear res
    end
    
    %% Statistical thresholding
    j = 1 : size(allpow,2);
    m = squeeze(nanmean(allpow(:,:,icond,:),2));
    s = squeeze(nanstd(allpow(:,:,icond,:),[],2));
    th  = zeros(5003,5003);

    for il = 1 : size(allpow,1)
      
      disp(sprintf('Computing location %d ...',il));
      
      th1 = zeros(5003,1);
      th2 = zeros(5003,1);
      
      jj = j(j~=il);
      
      % compute average correlation per subject
        
      for jl = jj
        
        xmean = squeeze(allpow(il,jl,icond,:));
     
        z     = (xmean-m(jl,:)')./s(jl,:)';  
        [~, p]   = ttest(z);
        th1(jl)  = p < (0.01/2);
        
        z     = (xmean-m(il,:)')./s(il,:)';      
        [~, p]   = ttest(z);
        th2(jl)  = p < (0.01/2);
        
      end
    	% any connection significant?
     	th(il,:) = (th1+th2) > 0; clear th1 th2        

    end
   	save([outdir sprintf('nc_powcorr_threshold_c%d_f%d_v%d.mat',icond,ifoi,v)],'th','-v7.3')

  end
end
  

%% plot degree (with statistical threshold) 
plt = 1;

if plt

load sa_meg_template;
grid = sa_meg_template.grid_medium;
mri  = sa_meg_template.mri; 
clear sa_meg_template

for ifoi = 1 : 25
  ifoi
  for icond = 1 : 2
    
   	load([outdir sprintf('nc_powcorr_threshold_c%d_f%d_v%d.mat',icond,ifoi,v)])
    c(:,ifoi,icond) = sum(th,2);
    
  end
end
   
%%
ifreq = 17;

degree = (squeeze(c(:,ifreq,1))/5003)-(squeeze(c(:,ifreq,2))/5003);

para = [];
para.mydotmarkersize=20;
para.orientation='axial';
% para.colorlimits = [-0.1 0.1];



h=figure; hold on
set(h,'color','w');

subplot(2,1,1)
showmri_transp_v3(mri,para,[grid degree],grid(1,:));

% degree = squeeze(c(:,ifreq,2))/5003;
% subplot(2,1,2)
% showmri_transp_v3(mri,para,[grid degree],grid(1,:));

% saveas(h,sprintf([plotdir 'nc_powcorr_diff_f%d_loc%d_v%d.fig'],ifoi,idx,1));
% close

% lim(1) = min([min(allpow_plt(1,:)) min(allpow_plt(2,:))]);
% lim(2) = max([max(allpow_plt(1,:)) max(allpow_plt(2,:))]);

% para.colorlimits = [lim(1) lim(2)];

end