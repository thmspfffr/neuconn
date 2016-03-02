%% nc_preproc_proc.m
% nc_preproc_proc

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% M1 or M2?
% -------------------------------------------------------------------------
m = 1;
% -------------------------------------------------------------------------

indir  = sprintf('/home/gnolte/neuconn/meg_data/m%d/',m);
outdir = sprintf('/home/gnolte/neuconn/meg_out/rest/raw/m%d/',m);

addpath /home/gnolte/neuconn/matlab/rest/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults

if m == 1 
  allsubj = nc_allsubj_start(m);
elseif m == 2
  allsubj = nc_allsubj_start(m);
end


% -------------------------------------------------------------------------
% VERSION 1
% ------------------------------------------------------------------------- 
% v = 1;
% pad = 1200/2; % 500 ms
% -------------------------------------------------------------------------
% VERSION 2
% ------------------------------------------------------------------------- 
% v_in = 1;
% v_out = 2;
% pad = 1200/2; % 500 ms
% CFG.dftfreq     = [];
% CFG.dftfilter   = 'no';
% ------------------------------------------------------------------------- 
% VERSION 3
% ------------------------------------------------------------------------- 
% % change compared to version 2: cut out first 20 seconds of data
% v_in = 1;
% v_out = 3;
% pad = 1200/2; % 500 ms
% CFG.dftfreq     = [];
% CFG.dftfilter   = 'no';
% ------------------------------------------------------------------------- 
% VERSION 4
% ------------------------------------------------------------------------- 
% diff to version 3: delete train segments and last segment of data
v_in = 1;
v_out = 4;
pad = 1200/2; % 500 ms
CFG.dftfreq     = [];
CFG.dftfilter   = 'no';
% ------------------------------------------------------------------------- 


addpath /home/gnolte/neuconn/matlab/rest/
%%
for icond = 1 : 1

for isubj = 7 : 7
  for ibl = 1 : 2
%     try
      
%       if icond == 1
%         a(1) = exist([outdir sprintf('nc_preproc_artvec_s%d_b%d_v%d.mat',isubj,ibl,v_in)]);
%       elseif icond == 2
%         a(1) = exist([outdir sprintf('nc_preproc_artvec_hc_s%d_b%d_v%d.mat',isubj,ibl,v_in)]);
%       end
%       
      a(1) = exist([outdir sprintf('nc_preproc_artvec_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v_in)]);
      a(2) = ~exist(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat'],icond, isubj,ibl,v_out));
      a(3) = ~exist(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt'],icond, isubj,ibl,v_out));  
      a = [1 1 1] 
      if a(1) && a(2) && a(3)
        
      system(['touch ' outdir sprintf('nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt',icond, isubj,ibl,v_out)]);
      
      load([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v_in)])

      % NEW VERSION
      load([outdir sprintf('nc_preproc_artvec_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v_in)])
%       load([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v_in)])

      disp(sprintf('Processing subject %d block %d cond %d',isubj,ibl,icond));
      
      cfg = [];

      cfg = cfgs;
%       cfg = rmfield(cfg,'trl');

      cfg.channel     = {'MEG'};
      cfg.padding     = 10;
      cfg.continuous  = 'yes';
      system(['touch ' outdir sprintf('nc_preproc_artifacts_c%d_s%d_b%d_v%d_processing.txt',icond, isubj,ibl,v_out)]);

      cfg.bsfilter    = 'yes';
      cfg.hpfilter    = 'yes';
      cfg.lpfilter    = 'yes';
      cfg.dftfilter   = CFG.dftfilter;
      cfg.dftfreq     = CFG.dftfreq;

      cfg.bsfreq      = [49 51; 99 101; 149 151];
      cfg.hpfreq      = 30;
      cfg.lpfreq      = 150;

      data.Fs         = 1200;
      data_hi         = ft_preprocessing(cfg);

      cfg.hpfiltord    = 4;
      cfg.lpfiltord    = 4;

      cfg.lpfreq      = 40;
      cfg.hpfreq      = 0.5;
      data_low        = ft_preprocessing(cfg);
      
      %%
      % -------------------------------------------------------------------------
      % REJECT ARTIFACTS
      % -------------------------------------------------------------------------
      
      art(art>data.sampleinfo(end))=data.sampleinfo(end);   
      art(art<1) = 1;
      
      for iart = 1 : size(art,1) 
        data_low.trial{1}(:,art(iart,1):art(iart,2)) = NaN;
      end

      
      
      
      %% PROBLEM SOLUTION
      
      % CUT OUT LAST BIT OF IDX. it does not matter anyway
      % then data and idx are of same length
      % add idx as second sensor for ft_resample
      % test
      
      
      idx = isnan(data_low.trial{1}(1,:));
      
      data_low.trial{1}(:,idx)=[];
      data_low.time{1}(:,idx) =[];
      if v_out > 2
        data_low.trial{1}(:,1:20*1200) =[];
        data_low.time{1}(:,1:20*1200)  =[];
      end
      
      data_low.sampleinfo = [1 length(data_low.time{1})];

      for iart = 1 : size(art,1) 
        data_hi.trial{1}(:,art(iart,1):art(iart,2)) = NaN;
      end

      data_hi.trial{1}(:,idx)=[];
      data_hi.time{1}(:,idx) =[];
      
      if v_out > 2
        data_hi.trial{1}(:,1:20*1200) =[];
        data_hi.time{1}(:,1:20*1200)  =[];
      end
      
      data_hi.sampleinfo = [1 length(data_hi.time{1})];

      %% RESAMPLE DATA

      cfg2            = [];
      cfg2.resamplefs = 300;
      cfg2.detrend    = 'yes'; % resampling with NaNs does not work if 'yes'
      cfg2.demean     = 'no'; % resampling with NaNs does not work if 'yes'
      data_hi         = ft_resampledata(cfg2, data_hi);
      data_low        = ft_resampledata(cfg2, data_low);
      data_low.idx    = resample(double(idx),300,1200);
      data_hi.idx     = resample(double(idx),300,1200);

      %% SAVE PREPROCESSED FILES

      save([outdir sprintf('nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v_out)], 'data_hi','data_low','idx', '-v7.3')
      disp(sprintf('Saved.'));
      
      else 
        system(['touch ' outdir sprintf('nc_preproc_artifacts_c%d_s%d_b%d_v%d_FILENOTFOUND.txt',icond, isubj,ibl,v_out)]);
        continue
      end
%     catch me
%       system(['touch ' outdir sprintf('nc_preproc_artifacts_c%d_s%d_b%d_v%d_ERROR.txt',icond, isubj,ibl,v_out)]);
%       continue
%     end
  end
end
end