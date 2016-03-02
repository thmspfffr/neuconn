%% COMPUTE POWCORR FOR ALL SUBJECTS / CONDITIONS
% nc_src_powcorr_diff

% v3: amplitude correlations
% v4: lagged phase coherence
% v5: amplitude correlations, with fake log-scaling
% v6: lagged phase coherence, with fake log-scaling

% implement seed-based?

% this step computes matrices containing all voxel values
% for all subjects and both conditons. This is needed for
% statistics (e.g. nc_src_ana_all2all_stats.m).
% The actual difference between conditions is computed in
% a next step at the end of this script.
% (09-16-2014)

% nc_src_powcorr_diff

clear all


m = 1; 

% --------------------------------------------------------
% VERSION 7
% --------------------------------------------------------
v = 7;
foi_range = unique(round(2.^[1:.25:7]));
grid_size = 'cortex3000';
NSUBJ = 49;
usetriu = 0;
method = 'powcorr';
% --------------------------------------------------------
% VERSION 9
% --------------------------------------------------------
% v = 9;
% foi_range = unique(round(2.^[1:.25:7]));
% grid_size = 'cortex';
% NSUBJ = 36;
% usetriu = 0;
% method = 'powcorr';
% foi_range = [1:13 14:2:30 31:3:60 61:5:128];
% -----------------------------------------------------
% VERSION 10
% --------------------------------------------------------
% v = 10;
% foi_range = unique(round(2.^[1:.25:7]));
% grid_size = 'coarse';
% NSUBJ = 22;
% usetriu = 0;
% method = 'powcorr';
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

ft_defaults

indir   = '/home/gnolte/neuconn/meg_out/rest/src/';
outdir = '/home/gnolte/neuconn/meg_out/rest/conn/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';


%%
% First create one big matrix for all subjects
% individually for each frequency of interest
% At the end of this block, this matrix is saved

load sa_meg_template;

if usetriu
  if strcmp(grid_size,'coarse')
    allpow = zeros((2113*2113-2113)/2,2,NSUBJ);
  elseif strcmp(grid_size,'medium')
    allpow = zeros((5003*5003-5003)/2,2,NSUBJ);
  elseif strcmp(grid_size,'cortex')
    allpow = zeros((3000*3000-3000)/2,2,NSUBJ);
  end 
else
  if strcmp(grid_size,'coarse')
    allpow = zeros(2113,2113,2,NSUBJ);
  elseif strcmp(grid_size,'medium')
    allpow = zeros(5003,5003,2,NSUBJ);
  elseif strcmp(grid_size,'cortex')
    allpow = zeros(3000,3000,2,NSUBJ);
  end
end

for ifoi = 1 : length(foi_range)
  
  if ~exist(sprintf([outdir 'nc_src_allpowcorr_f%d_v%d_processing.txt'],ifoi,v))
    system(['touch ' outdir sprintf('nc_src_allpowcorr_f%d_v%d_processing.txt',ifoi,v)]);
  else
    continue
  end
  
  
  for icond = 1 : 2
    
    disp(sprintf('Computing f%d c%d ...',ifoi,icond))
    
     allsubj{icond} = [];
    
    cnt = 0;
    for isubj = 1 : 49

      % average the two blocks
      for iblock = 1 : 2
        if exist([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
          
          clear resout
          
          load([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
          
          resout = single(resout);
          
          br = 0;
          % correction factor for correlations and fisher z-transform
         switch method
           case 'powcorr'
           	res(:,:,iblock) = (resout+resout')./2; clear resout
%             res(:,:,iblock) = (atanh(resout)+atanh(resout)')./2; clear resout
           case 'lpc'
            res(:,:,iblock) = resout; clear resout
         end        
        else
        	system(['touch ' outdir sprintf('nc_src_allpowcorr_s%d_c%d_b%d_f%d_v%d_missing.txt',isubj,icond,iblock,ifoi,v)]);
          br = 1;
          break
        end
      end
              
      
      % break & continue with next subject
      if br
        continue
        br = 0;
      end
      
      disp(sprintf('Processig s%df%dc%d',isubj,ifoi,icond));
      cnt = cnt + 1;
      
      if usetriu
        dummy = triu(ones(size(res,1),size(res,2)),1); 
        res = single(res);
        res = res(find(dummy)); clear dummy
        allpow(:,icond,cnt)	= single(nanmean(res,3)); clear res
      else
        allpow(:,:,icond,cnt)	= single(nanmean(res,3)); clear res
        allsubj{icond}(cnt) = isubj;
      end
      
      
      if cnt == NSUBJ
        disp(sprintf('%d subjects computed ...',cnt));
        break
      end
      
    end
  end
  
  save([outdir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)],'allpow','allsubj','-v7.3')
  clear allpow allsubj
end

%% PLOT AND COMPUTE DIFFERENCES
% This actually computes and plots raw differences between
% the two conditions.

% v = 1;
% 
% subj = [4 5 6 7 8 9 11 16 17 19 20 23];
% 
% ord = pconn_randomization;
% 
% for im = 1 : 3
% 
%   [idx(:,1,im) idx(:,2,im)] = find(ord==im);
% 
% end
% 
% for ifoi = 10 : 10
%   
%   load([outdir sprintf('nc_src_allpowcorr_f%d_v%d.mat',ifoi,v)])
%   
%   m1 = squeeze(nanmean(allpow(:,:,1,:),4));
%   m2 = squeeze(nanmean(allpow(:,:,2,:),4));
%   
%   m = m1-m2;
%   
%   clear m1 m2
% 
% end
% % 
% m = nanmean(m,1);
% 
% load sa_meg_template;
% grid = sa_meg_template.grid_coarse;
% mri  = sa_meg_template.mri;
% clear sa_meg_template
% 
% para = [];
% para.mydotmarkersize=20;
% para.orientation='coronal';
% % para.colorlimits = [0 0.2];
% 
% h=figure; hold on
% set(h,'color','w');
% 
% showmri_transp_v3(mri,para,[grid m'],grid(1,:));

% 
% %% Statistical thresholding
% 
% % X_D: average difference between subjects
% % m = squeeze(nanmean(allpow(:,:,1,:),4))-squeeze(nanmean(allpow(:,:,2,:),4));
% % % S_D: standard deviaton of difference values
% % s = squeeze(nanstd(allpow(:,:,1,:)-allpow(:,:,2,:),[],4));
% 
% j = 1 : size(allpow,1);
% 
% for il = 1 : size(allpow,1)
%   
%   disp(sprintf('Computing location %d ...',il));
%   
%   th1 = zeros(5003,1);
%   th2 = zeros(5003,1);
%   th  = zeros(5003,5003);
%   
%   jj = j(j~=il);
%   
%   % compute average correlation per subject
%       
%   th1   = ttest(allpow(il,jj,1,:),allpow(il,jj,2,:),'dim',4);
%   
% 
% end
% save([outdir sprintf('nc_powcorr_diff_threshold_c%d_f%d_v%d.mat',icond,ifoi,v)],'th','-v7.3')
% 



