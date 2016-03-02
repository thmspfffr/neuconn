%% nc_preproc_artvec.m
% nc_preproc_artvec

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
v = 1;
pad = 1200/2; % 500 ms
% -------------------------------------------------------------------------
%%
for icond = 1 : 1

for isubj = 13 : 13
  try
  subj = allsubj{icond};
  if subj{isubj,3} == 0
    continue
  end
  for ibl = 1 : 2
   
    clear cfgs
    
%     if ~exist([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)])      
%       continue
%     else
%       
%       if exist([outdir sprintf('nc_preproc_artvec_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)])
%         continue
%       else
        
        disp(sprintf('Processing c%d s%d b%d ...',icond,isubj,ibl))
        
       load([outdir sprintf('nc_preproc_data_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)])
        
        cfg1.trl = cfg1.trl;
        cfgs.trl = cfg1.trl;

        % CUT OUT 10 MINUTES
        strt = find(data.trial{1}(strcmp(data.label,'UPPT001'),:)==200,50,'first');

        % VISUAL ARTIFACT REJECTION
 
        cfg = [];
        cfg.method    = 'summary';
        cfg.channel   = {'MEG'};
        tmp_data      = ft_rejectvisual(cfg,data);
        artifact_vis  = tmp_data.cfg.artifact; tmp_data

        cfg = [];
        cfg.artfctdef.rejvis.artifact = artifact_vis;
        cfg.viewmode = 'vertical';
        cfg = ft_databrowser(cfg,data);
        artifact_vis = cfg.artfctdef.rejvis.artifact;


        % -------------------------------------------------------------------------
        % DETECT SQUID JUMPS 
        % -------------------------------------------------------------------------
%         cfg = cfgs;
%         cfg.continuous  = 'yes';
%         cfg.artfctdef.zvalue.channel       = 'MEG';
%         cfg.artfctdef.zvalue.cutoff        = 150;
%         cfg.artfctdef.zvalue.trlpadding    = 0;
%         cfg.artfctdef.zvalue.artpadding    = 0.3;
%         cfg.artfctdef.zvalue.fltpadding    = 0;
%         cfg.artfctdef.zvalue.cumulative    = 'yes';
%         cfg.artfctdef.zvalue.medianfilter  = 'yes';
%         cfg.artfctdef.zvalue.medianfiltord = 9;
%         cfg.artfctdef.zvalue.absdiff       = 'yes';
%         cfg.artfctdef.zvalue.interactive   = 'yes';
% 
%         [cfg, artifact_jump] = ft_artifact_zvalue(cfg);
% 
%         % -------------------------------------------------------------------------
%         % DETECT CAR ARTIFACTS
%         % -------------------------------------------------------------------------
%         cfg = cfgs;
%         cfg.continuous  = 'yes';
%         cfg.artfctdef.zvalue.channel       = 'MEG';
%         cfg.artfctdef.zvalue.cutoff        = 40;
%         cfg.artfctdef.zvalue.trlpadding    = 0;
%         cfg.artfctdef.zvalue.artpadding    = 0.3;
%         cfg.artfctdef.zvalue.fltpadding    = 0;
%         cfg.artfctdef.zvalue.lpfilter    = 'yes';
%         cfg.artfctdef.zvalue.lpfreq      = 0.3;
%         cfg.artfctdef.zvalue.lpfiltord   = 2;
%         cfg.artfctdef.zvalue.lpfilttype  = 'but';
%         cfg.artfctdef.zvalue.cumulative    = 'yes';
%         cfg.artfctdef.zvalue.medianfilter  = 'yes';
%         cfg.artfctdef.zvalue.medianfiltord = 9;
%         cfg.artfctdef.zvalue.absdiff       = 'yes';
%         cfg.artfctdef.zvalue.interactive   = 'yes';
% 
%         [cfg, artifact_car] = ft_artifact_zvalue(cfg);

        % -------------------------------------------------------------------------
        % DETECT MUSCLE ARTIFACTS
        % -------------------------------------------------------------------------
        cfg = cfgs;
        cfg.continuous  = 'yes';
        cfg.artfctdef.zvalue.channel     = {'MRT','MLT'};
        cfg.artfctdef.zvalue.cutoff      = 8;
        cfg.artfctdef.zvalue.trlpadding  = 0;
        cfg.artfctdef.zvalue.artpadding  = 0.3;
        cfg.artfctdef.zvalue.fltpadding  = 0;
        cfg.artfctdef.zvalue.bpfilter    = 'yes';
        cfg.artfctdef.zvalue.bpfreq      = [110 140];
        cfg.artfctdef.zvalue.bpfiltord   = 7;
        cfg.artfctdef.zvalue.bpfilttype  = 'but';
        cfg.artfctdef.zvalue.hilbert     = 'yes';
        cfg.artfctdef.zvalue.boxcar      = 0.2;
        cfg.artfctdef.zvalue.interactive = 'yes';

        [cfg, artifact_muscle] = ft_artifact_zvalue(cfg);

        % SAVE ALL ARTIFACTS

        art = [artifact_muscle; artifact_vis];
        art(:,1) = art(:,1)-pad; art(:,2)=art(:,2)+pad;
        art(art>data.sampleinfo(end))=data.sampleinfo(end);
        art(art<0) = 1;

        save([outdir sprintf('nc_preproc_artvec_c%d_s%d_b%d_v%d.mat',icond,isubj,ibl,v)], 'art')
%       end
%     end
  end
end
end
end