%% READ IN DATA AND SAVE AS MAT
% nc_preproc_read_data

clear all
restoredefaultpath

% -------------------------------------------------------------------------
% M1 or M2?
% -------------------------------------------------------------------------
m = 	1;
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
v = 1;

pad = 1200/2; % 500 ms

icond = 1;

subj = allsubj{icond};
  
numdir = [1 2; 3 4];

if icond == 1 
  s = '';
else
  s = 'hc_';
end
%%

for isubj = 7 : 7
  for ibl = 1 : 2
    try

    % not exist: both
%     if ~exist([outdir sprintf('nc_preproc_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,ibl,v)]) ...
%     && ~exist([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)])
%       
%       iblock = ibl;   
%       system(['touch ' outdir sprintf('nc_preproc_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,ibl,v)]);
%       disp(sprintf('Processing subject c%d s%d b%d',icond,isubj,ibl));
% 
%     % no proc file, but output: create proc file
%     elseif ~exist([outdir sprintf('nc_preproc_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,ibl,v)]) ...
%     && exist([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)])
%   
%      	system(['touch ' outdir sprintf('nc_preproc_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,ibl,v)]);
%       continue
%       
%     elseif exist([outdir sprintf('nc_preproc_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,ibl,v)]) ...
%     && ~exist([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)])
%  
%       iblock = ibl;   
%       disp(sprintf('Processing subject c%d s%d b%d',icond,isubj,ibl));
%   
%     else
%       continue
%     end
      
    iblock = ibl;   

    data_dir = [indir subj{isubj}];
    cont_dir = dir(data_dir);

    % READ DATA
    % ---------------------------------------------------------------
    idir = numdir(iblock,subj{isubj,2});

    for i = 1 : length(cont_dir)
      tmp(i) = ~isempty(regexp(cont_dir(i).name,'(.ds)','once'));
    end

    ind = find(tmp>0);
    idir = ind(idir);

    tmp_dataset      = [cont_dir(idir).name];
    cfgs.datafile    = [data_dir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'meg4'];
    cfgs.headerfile  = [data_dir '/' tmp_dataset '/' cont_dir(idir).name(1:end-2) 'res4'];    

    cfg = [];
    % cfg.continuous = 'yes';
    cfg.dataset    = [data_dir '/' cont_dir(idir).name];
    
    a=ft_read_event(cfg.dataset); 

    if size(cell2mat({a.value}),2)>10
      error('Too many triggers!');
    end

    cfg.channel   = {'MEG'};
    cfg.demean    = 'yes';
    cfg.precision = 'single';
    data          = ft_preprocessing(cfg);

    cfg = [];
    cfg.trialfun              = 'ft_trialfun_general';
    cfg.trialdef.triallength  = 1;
    cfg.dataset               = [data_dir '/' cont_dir(idir).name];
    cfg1 = ft_definetrial(cfg);

    save([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)], 'data','cfg1','cfgs')

    catch me
      sprintf('ERROR! s%d b%d', isubj,iblock);
      save([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d_err.mat',icond,isubj,iblock,v)], 'me')
    end 
  end
end