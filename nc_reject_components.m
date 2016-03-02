%% NEUCONN REJECT ICA COMPONENTS
% after nc_viewcomps_v2.m

clear all
restoredefaultpath
m = 1;

% =========================================================================
% DEFINITIONS
% -------------------------------------------------------------------------

datadir  = sprintf('/home/gnolte/neuconn/meg_data/m%d/',m);


% outdir = sprintf('/home/gnolte/neuconn/meg_out/rest/raw/m%d/',m);
outdir = sprintf('/home/gnolte/neuconn/meg_out/rest/m%d/proc/ica/',m);

indir  = sprintf('/home/gnolte/neuconn/meg_out/rest/m%d/proc/ica/',m);

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults

% -------------------------------------------------------------------------
% VERSION 4
% -------------------------------------------------------------------------
v     = 6; % ICA version
v_art = 4; % preproc version
split = 0;
% -------------------------------------------------------------------------

if split == 0
  nseg = 1;
else
  nseg = 2;
end
%%
for icond = 1 : 2
  for isubj = 1 : 49
    for iblock = 1 : 2
      
      if exist([outdir sprintf('nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)])
        continue
      else
      end
      
      for iseg = 1: nseg
        
        ex(1) = exist([outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)]);
        
        if ex(1)
          
          load([outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)])
          try
            for ifr = 1 : 2
              
              ex(2) = exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,ifr,isubj,iblock,iseg,v));
              
              if ex(2)
                disp(sprintf('Processing c%d s%d b%d ...',icond,isubj,iblock))
                load(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,ifr,isubj,iblock,iseg,v));
                
                cfg = [];
                cfg.component = find(rej_comp);
                
                if ifr == 1
                  data_low	= ft_rejectcomponent(cfg,comp_low);
                elseif ifr == 2
                  data_hi   = ft_rejectcomponent(cfg,comp_hi);
                end
              else
                error('Both output files needed!')
              end
            end
            
            save([outdir sprintf('nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)],'data_low','data_hi','-v7.3')
          catch me
            continue
          end
        else
          continue
        end
      end
    end
  end
end







