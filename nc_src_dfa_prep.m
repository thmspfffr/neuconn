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
addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/

ft_defaults

indir   = '/home/gnolte/neuconn/meg_out/rest/src/';
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';

% call subj info startup function
allsubj = nc_allsubj_start;

%%

for ifoi = 10 : 10
  for icond = 1 : 2
    for isubj = 1 : 47
      
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
        v_filt = 1;
        load([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v_filt)],'A','L','grid');
        
        clear L grid
                
        for iblock = 1 : 2
          
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
          F = getdipdir(cs(:,:,ifoi),A); clear A cs
          
          for iloc = 1 : 5003
            
            disp(sprintf('Computing DFA for location %d ...',iloc));
            dat = mydata*F(:,iloc);
            
            siginfo = nbt_Info;
            siginfo.converted_sample_frequency = 300;
            
            % computing amplitude envelopes
            ampenv = nbt_GetAmplitudeEnvelope(dat,siginfo,8,13,2/8);
            
            % compute DFA
            [tmp,expo(iblock,iloc)] = nbt_doDFA(ampenv, siginfo, [1 15],[.8 25],0.5,1,0,[]);
            
          end
          
          
        end
        
        expo = nanmean(expo,1);
        
        save(sprintf([outdir 'nc_src_dfa_f%d_c%d_s%d_v%d.mat'],ifoi,icond,isubj,v),'expo','-v7.3');
        
        
      catch me
        save(sprintf([outdir 'nc_src_dfa_f%d_c%d_s%d_v%d_error.mat'],ifoi,icond,isubj,v),'me','-v7.3');
      end
      
    end
  end
end


