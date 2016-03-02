%% COMPUTE AMPLITUDE ENVELOPE CORRELATIONS
% nc_src_ana

% function nc_src_ana(rand_num)

% rng(rand_num,'twister')

clear all

% --------------------------------------------------------
% VERSION 2 all2all power correlations (SEGLEN too short for freq res)
% --------------------------------------------------------
% v = 1;          % version of head model
% v_cs = 1; 
% version of cross spectrum
% v_out = 2;      % version of output
% SEGLENG = 200;
% gridinfo = 'grid = sa.grid_medium;';
% --------------------------------------------------------
% VERSION 5 all2all power correlations 
% --------------------------------------------------------
% v         = 1;          % version of head model
% v_cs      = 2;          % version of cross spectrum
% v_out     = 5;          % version of output
% lphc      = 0;          % compute lagged phase coherence?
% powcorr   = 1;          % compute power correlations?
% SEGLENG   = 300;        % segment length
% foi_range = unique(round(2.^[1:.25:7]));
% gridinfo = 'grid = sa.grid_medium;';
% % difference to v_out = 3: foi range log instead of lin
% --------------------------------------------------------
% VERSION 6 all2all lagged phase coherence
% --------------------------------------------------------
% v = 1;            % version of head model
% v_cs = 2;         % version of cross spectrum
% v_out = 6;        % version of output
% lphc = 1;         % compute lagged phase coherence?
% powcorr = 0;      % compute power correlations?
% SEGLENG = 300;
% foi_range = unique(round(2.^[1:.25:7]));
% gridinfo = 'grid = sa.grid_medium;';
% difference to v_out = 4: foi range log instead of lin
% --------------------------------------------------------
% VERSION 7 all2all power correlations SMALL GRID
% --------------------------------------------------------
v         = 1;          % version of head model
v_cs      = 2;          % version of cross spectrum
v_out     = 7;          % version of output
lphc      = 0;          % compute lagged phase coherence?
powcorr   = 1;          % compute power correlations?
SEGLENG   = 300;        % segment length
foi_range = unique(round(2.^[1:.25:7]));
gridinfo = 'grid = sa.grid_coarse_indi;';
% --------------------------------------------------------


restoredefaultpath

% addpaths
addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest/
ft_defaults

% Define relevant input/output paths
indir  = '/home/gnolte/neuconn/meg_data/m1/';
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/';

% call subj info startup function
allsubj = nc_allsubj_start;

% Define frequency range of interest


%%
for icond = 2 : 2
  for isubj = 1 : 40
    
    % Only process selected datasets
    if cell2mat(allsubj{icond}(isubj,3)) == 0
      continue
    end
    
   	% -----------------------------------------------------------------
    % Loop through frequencies of interest
    % -----------------------------------------------------------------
    for foi = foi_range
      
      ifoi = find(foi==foi_range);
      
      if foi > 35
        freq = 2;
      else
        freq = 1;
      end
      
      for iblock = 1 : 2
        
        clear A L grid cs nave coh
        
        try
          
          % in case, job is stuck in queue, pause:
          pause(ceil(rand*500));
          % -----------------------------------------------------------------
          % Check whether file has already been processed
          % -----------------------------------------------------------------
          
          if     ~exist(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat'],isubj,icond,iblock,ifoi,v_out)) ...
              && ~exist(sprintf([outdir 'nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v_out))
            
            system(['touch ' outdir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,ifoi,v_out)]);
            
            disp(sprintf('Processing s%d c%d f%d ...',isubj,icond,ifoi));
            disp(sprintf('Processing block %d ...',iblock));
            
            % -----------------------------------------------------------------
            % Load data
            % -----------------------------------------------------------------
            disp(sprintf('Loading MEG data ...'));
            load(sprintf('/home/gnolte/neuconn/meg_out/rest/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6));
            % -----------------------------------------------------------------
            
            if freq == 1
              clear data_hi
              [mydata,epleng] = megdata2mydata(data_low);
            elseif freq == 2
              clear data_low
              [mydata,epleng] = megdata2mydata(data_hi);
            else
              error('Missing information on frequency!')
            end
            
            % -----------------------------------------------------------------
            % Check whether leadfield, grid and A are loaded already
            % -----------------------------------------------------------------
            
            if ~exist('A','var') && ~exist('grid','var')
              if ~exist([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)]) ...
               && exist([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)])
                
                disp(sprintf('Computing spatial filter ...'));
                
                % LOAD OUTPUT FROM NC_PREP_SRC_ANA.M
                load([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)]);
                
                eval(sprintf('%s',gridinfo)); 
                
                A = mkfilt_eloreta_v2(sa.L_medium);
%                 A = mkfilt_eloreta_v2(L);

                save([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','grid');
                
              elseif exist([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)])
                disp(sprintf('Loading forward model ...'));
                load([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','grid');
              end
            end
            
            % -----------------------------------------------------------------
            % Compute cross-spectrum if not yet computed
            % -----------------------------------------------------------------
            
            if ~exist('cs','var')
              
              % if no output yet and file not currently being processed 
              if     ~exist([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]) ...
                  && ~exist([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,freq,v_cs)])

                % tell system that file is now being processed
                system(['touch ' outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,freq,v_cs)]);

                disp(sprintf('Computing cross spectrum ...'));

                segleng         = SEGLENG;
                segshift        = segleng/2;
                
                if freq == 1
                  maxfreqbin = max(foi_range);
                else
                  maxfreqbin = max(foi_range);
                end

                [cs, coh, nave] = data2cs_event(mydata,segleng,segshift,epleng,maxfreqbin);
                
                disp(sprintf('Saving cross spectrum ...'));
                save([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)],'cs');
                save([outdir sprintf('nc_cs_fm_coh_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)],'coh','nave');
              
              % no output yet, but file is already being processed
              elseif ~exist([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]) ...
                  && exist([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,freq,v_cs)])
    
                disp(sprintf('Waiting for input ...'));
                
                while ~exist([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]) ...
                  pause(60)
                end
                
                % if output, load
                disp(sprintf('Loading cross spectrum ...'));
                load([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]);         
              
              % output exists already
              else
                disp(sprintf('Loading cross spectrum ...'));
                load([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]);
              end
            end
            
            % -----------------------------------------------------------------
            % Compute power correlations per carrier-foi
            % -----------------------------------------------------------------
            
            disp(sprintf('Processing frequency %d ...',ifoi));
            
            [F1,p] = getdipdir(cs(:,:,foi),A);
            
            segleng   = SEGLENG;
            segshift  = segleng/2;
            fsample   = 300;
            
            disp(sprintf('Computing amplitude correlations ...'));
            
            if powcorr 
              resout = nc_orthopowercorr(mydata,segleng,segshift,epleng,foi,fsample,F1,F1);
            elseif lphc
              resout = nc_lphcoh(mydata,segleng,segshift,epleng,foi,fsample,F1,F1);
            end
              
            save([outdir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v_out)],'resout');
            clear resout F1 p
            % -----------------------------------------------------------------
          else
            continue
          end
      	% -----------------------------------------------------------------
        % Catch error message and save
        % -----------------------------------------------------------------    
        catch me
          save([outdir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing_err.mat',isubj,icond,iblock,ifoi,v_out)],'me');
          system(['touch ' outdir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d_processing_err.txt',isubj,icond,iblock,ifoi,v_out)]);
          warning('Something went wrong!')
        end
        clear data_hi data_low
      end
     
    end
  end
end

