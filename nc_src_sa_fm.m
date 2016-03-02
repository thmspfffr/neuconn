%% COMPUTE AMPLITUDE ENVELOPE CORRELATIONS
% nc_src_sa_fm

% function nc_src_ana(rand_num)

% rng(rand_num,'twister')

clear all
m = 1;

% --------------------------------------------------------
% VERSION 1 - ALL WAVELETS, POWCORR, VARIABLE SEGS
% --------------------------------------------------------
% v = 3;            % head model: cortex
% v_cs = 20;        % new cs for wavelets
% v_out = 1;        % version of output
% lphc = 0;         % compute lagged phase coherence?
% powcorr = 1;      % compute power correlations?
% foi_range = 2.^(0:.5:7);
% SEGLENG   = (5.8./foi_range)*300;
% filt = 'mkfilt_eloreta_v2';
% gridinfo  = 'grid = sa.grid_cortex3000;';
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
% VERSION 7 all2all power correlations COARSE GRID
% --------------------------------------------------------
v         = 2;          % version of head model
v_cs      = 2;          % version of cross spectrum
v_out     = 7;          % version of output
lphc      = 0;          % compute lagged phase coherence?
powcorr   = 1;          % compute power correlations?
SEGLENG   = 300;        % segment length
foi_range = unique(round(2.^[1:.25:7]));
gridinfo  = 'grid = sa.grid_coarse;';
filt = 'mkfilt_eloreta_v2';
% --------------------------------------------------------
% VERSION 8 all2all power correlations CORTEX GRID
% --------------------------------------------------------
% v         = 3;          % version of head model
% v_cs      = 2;          % version of cross spectrum
% v_out     = 8;          % version of output
% lphc      = 0;          % compute lagged phase coherence?
% powcorr   = 1;          % compute power correlations?
% SEGLENG   = 300;        % segment length
% foi_range = unique(round(2.^[1:.25:7]));
% gridinfo  = 'grid = sa.grid_cortex3000;';
% filt = 'mkfilt_eloreta_v2';
% --------------------------------------------------------
% VERSION 9 all2all powcorr/cortex
% --------------------------------------------------------
% v         = 3;          % version of head model
% v_cs      = 2;          % version of cross spectrum
% v_out     = 9;          % version of output
% lphc      = 0;          % compute lagged phase coherence?
% powcorr   = 1;          % compute power correlations?
% foi_range = [1:13 14:2:30 31:3:60 61:5:128];
% SEGLENG   = ones(1,length(foi_range))*300;        % segment length

% seg_len   = (5.8./foi_range)*300;
% round(ans*10)/10
% gridinfo  = 'grid = sa.grid_cortex3000;';
% filt = 'mkfilt_eloreta_v2';
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
allsubj = nc_allsubj_start(m);

% Define frequency range of interest


%%
for icond = 1 : 2
  for isubj = 8 : 8
    
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
%            sprintf('/home/gnolte/neuconn/meg_out/rest/m%d/proc/ica/',m)
          % in case, job is stuck in queue, pause:
%           pause(ceil(rand*500));
          % -----------------------------------------------------------------
          % Check whether file has already been processed
          % -----------------------------------------------------------------
          
          if  ~exist(sprintf([outdir 'nc_sa_fm_s%d_c%d_b%d_v%d_processing.txt'],isubj,icond,iblock,ifoi,v_out)) ...
              && exist(sprintf('/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6)) ...
              && exist(sprintf([outdir 'nc_sa_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v))
            
            system(['touch ' outdir sprintf('nc_sa_fm_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,ifoi,v_out)]);
            
            disp(sprintf('Processing s%d c%d f%d ...',isubj,icond,ifoi));
            disp(sprintf('Processing block %d ...',iblock));
            
            % -----------------------------------------------------------------
            % Load data
            % -----------------------------------------------------------------
            disp(sprintf('Loading MEG data ...'));
            load(sprintf('/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6));
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
            % Compute cross-spectrum if not yet computed
            % -----------------------------------------------------------------
            
            if ~exist('cs','var')
              
              % if no output yet and file not currently being processed 
              if     ~exist([outdir sprintf('nc_cs_fm_done_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]) ...
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
                
                if v ~= 9
                  [cs, coh, ss] = data2cs_wavelet(mydata,segleng,segshift,epleng,ifoi,300);
                elseif v == 9
                 [cs, coh, ss] = data2cs_event(mydata,segleng,segshift,epleng,ifoi,300);
                end
                
                disp(sprintf('Saving cross spectrum ...'));
                save([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)],'cs','-v7.3');
                save([outdir sprintf('nc_cs_fm_coh_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)],'coh','-v7.3');
               	
                %  Dummy file: When file is saved, other processes know
                %  that cross spectrum is complete and saved!
                save([outdir sprintf('nc_cs_fm_done_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)],'segleng','-v7.3');

              % no output yet, but file is already being processed
              elseif ~exist([outdir sprintf('nc_cs_fm_done_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]) ...
                  && exist([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,freq,v_cs)])
%                 
                disp(sprintf('Waiting for input ...'));
                
                while ~exist([outdir sprintf('nc_cs_fm_done_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]) ...
                  pause(60)
                end

                pause(randi(85)+25)
                
                % if output, load
                disp(sprintf('Loading cross spectrum ...'));
                load([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]);         
%               
              % output exists already
              else
                disp(sprintf('Loading cross spectrum ...'));
                load([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]);
              end
            end
             
            % -----------------------------------------------------------------
            % COMPUTE SPATIAL FILTER
            % -----------------------------------------------------------------
            
            if ~exist('A','var') && ~exist('grid','var')
              if ~exist([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)]) ...
               && exist([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)])
                
                disp(sprintf('Computing spatial filter ...'));
                
                % LOAD OUTPUT FROM NC_PREP_SRC_ANA.M
                load([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)]);             
                
                if strcmp(filt,'mkfilt_eloreta_v2')
                  if      strcmp(sa.leadfield,'medium')
                    A = mkfilt_eloreta_v2(sa.L_medium);
                  elseif  strcmp(sa.leadfield,'coarse')
                    A = mkfilt_eloreta_v2(sa.L_coarse);
                  elseif strcmp(sa.leadfield,'cortex')
                    % accidentely called 'L_coarse' in nc_prep_src_ana.m
                    A = mkfilt_eloreta_v2(sa.L_coarse);
                  end
                elseif strcmp(filt,'dicspow')
                  if      strcmp(sa.leadfield,'medium')
                    [~,A] = dicspow(sa.L_medium,real(cs(:,:,ifoi)),0.05);
                  elseif  strcmp(sa.leadfield,'coarse')
                    [~,A] = dicspow(sa.L_coarse,real(cs(:,:,ifoi)),0.05);
                  elseif strcmp(sa.leadfield,'cortex')
                    [~,A] = dicspow(sa.L_cortex,real(cs(:,:,ifoi)),0.05);
                  end
                elseif strcmp(filt,'pconn_beamformer')
                  if      strcmp(sa.leadfield,'medium')
                    [A] = pconn_beamformer(cs(:,:,ifoi),sa.L_medium);
                  elseif  strcmp(sa.leadfield,'coarse')
                    [A] = pconn_beamformer(cs(:,:,ifoi),sa.L_coarse);
                  elseif strcmp(sa.leadfield,'cortex')
                    [A] = pconn_beamformer(cs(:,:,ifoi),sa.L_cortex);
                  end
                end
                
                eval(sprintf('%s',gridinfo))

                save([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','grid','-v7.3');
                
              elseif exist([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)])
                disp(sprintf('Loading filter ...'));
                load([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','grid');
              end
            end
            % -----------------------------------------------------------------
          else
            continue
          end
      	% -----------------------------------------------------------------
        % Catch error message and save
        % -----------------------------------------------------------------    
        catch me
          save([outdir sprintf('nc_sa_fm_s%d_c%d_b%d_f%d_v%d_processing_err.mat',isubj,icond,iblock,ifoi,v_out)],'me');
          warning('Something went wrong!')
        end
        clear data_hi data_low
      end 
    end
  end
end

