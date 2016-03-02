%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_src_ana

clear all

% --------------------------------------------------------
% VERSION 2 (all2all)
% --------------------------------------------------------
v_out = 9;
m = 1;
% --------------------------------------------------------

restoredefaultpath

% addpaths
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest
ft_defaults

% Define relevant input/output paths
indir  = '/home/gnolte/neuconn/meg_data/m1/';
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/';

% call subj info startup function
allsubj = nc_allsubj_start(m);

% Define frequency range of interest
foi_range = 1:50;


%% clean nc_powcorr
v_out = 7;
cnt = 0;
for icond = 1 : 2
  for isubj = 1 : 49
    
    % Only process selected datasets
    if cell2mat(allsubj{icond}(isubj,3)) == 0
      continue
    end
%     profile on
    for ifoi = foi_range(1) : foi_range(end)
      ifoi
      for iblock = 1 : 2
         
        if     exist(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat'],isubj,icond,iblock,ifoi,v_out)) ...
            && ~exist(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v_out))
          
          system(['touch ' outdir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,ifoi,v_out)]);
          
        elseif exist(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v_out)) ...
            && ~exist(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat'],isubj,icond,iblock,ifoi,v_out))
          warning('deleting something')  
          delete(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v_out));
          cnt = cnt + 1;
        end
      end      
    end
%     profile viewer
  end
end

%%
error('Be careful!')
%% clean nc_src_powcorr_all2all (nc_threshold)

outdir = '/home/gnolte/neuconn/meg_out/rest/conn/';
v = 2;

for icond = 1 : 2
  
  for ifoi = 1 : 60
    
    if     exist(sprintf([outdir 'nc_powcorr_threshold_c%d_f%d_v%d.mat'],icond,ifoi,v)) ...
        && ~exist(sprintf([outdir 'nc_powcorr_threshold_c%d_f%d_v%d_processing.txt'],icond,ifoi,v))
      
      system(['touch ' outdir sprintf('nc_powcorr_threshold_c%d_f%d_v%d_processing.txt',icond,ifoi,v)]);
      
    elseif exist(sprintf([outdir 'nc_powcorr_threshold_c%d_f%d_v%d_processing.txt'],icond,ifoi,v)) ...
        && ~exist(sprintf([outdir 'nc_powcorr_threshold_c%d_f%d_v%d.mat'],icond,ifoi,v)) ...
        
      delete(sprintf([outdir 'nc_powcorr_threshold_c%d_f%d_v%d_processing.txt'],icond,ifoi,v));
      
      
    end
  end
end
%%
error('Be careful!')
%% clean nc_prep_src_ana (nc_threshold)

outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
v = 1;

for icond = 1 : 2
  
  for isubj = 1 : 50
    isubj
    for iblock = 1 : 2 
      
      if     exist(sprintf([outdir 'nc_sa_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v)) ...
          && ~exist(sprintf([outdir 'nc_sa_c%d_s%d_v%d_processing.txt'],icond,isubj,v))

        system(['touch ' outdir sprintf('nc_sa_c%d_s%d_v%d_processing.txt',icond,isubj,v)]);

      elseif exist(sprintf([outdir 'nc_sa_c%d_s%d_v%d_processing.txt'],icond,isubj,v)) ...
          && ~exist(sprintf([outdir 'nc_sa_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v)) ...

        delete(sprintf([outdir 'nc_sa_c%d_s%d_v%d_processing.txt'],icond,isubj,v));

      end
    end
  end
end

error('Be careful!')
%% clean nc_prep_src_ana (nc_threshold)

outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
v = 2;

for icond = 1 : 2
  cnt = 0;
  for isubj = 1 : 50
    isubj
    for iblock = 1 : 2
      for ifoi = 1 : 23
        
        if     exist(sprintf([outdir 'nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat'],isubj,icond,iblock,ifoi,v)) ...
            && ~exist(sprintf([outdir 'nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v))
          
          system(['touch ' outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,ifoi,v)]);
          
        elseif exist(sprintf([outdir 'nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v)) ...
            && ~exist(sprintf([outdir 'nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat'],isubj,icond,iblock,ifoi,v)) ...
            cnt = cnt + 1;
          delete(sprintf([outdir 'nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v));
          
        end
      end
    end
  end
end


%%
error('Be careful!')
%% clean nc_prep_src_ana (nc_threshold)

outdir = sprintf('/home/gnolte/neuconn/meg_out/rest/raw/m1/');
v = 4;

for icond = 1 : 2
  
  for isubj = 1 : 50
    isubj
    for iblock = 1 : 2 
      
      if     exist(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v)) ...
          && ~exist(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt'],icond,isubj,iblock,v))

        system(['touch ' outdir sprintf('nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,iblock,v)]);

      elseif exist(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt'],icond,isubj,iblock,v)) ...
          && ~exist(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v)) ...

        delete(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt'],icond,isubj,iblock,v));

      end
    end
  end
end


