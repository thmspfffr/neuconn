%% COMPUTE EXCITATION/INHIBITION BALANCE ACROSS SENSOR SPACE
% pconn_sens_dfa

% Method by R Hardstone, discussed during meeting on 4th of March
% Takes band-pass filtered amplitude and time-resolved DFA as input
% and estimates E/I balance from that.


% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
% v         = 1;
% v_rawdata = 6;
% is_src    = 0;
% fsample   = 400;
% SUBJLIST  = [4 5 6 7 8 9 11 16 17 19 20 23];
% foi       = [2 4; 4 8; 8 12; 12 24];
% i_fit     = [3 28];
% i_calc    = [.8 40];
% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
% v         = 2;
% v_rawdata = 6;
% is_src    = 0;
% fsample   = 300;
% SUBJLIST  = [4 5 6 7 8 9 11 16 17 19 20 23];
% foi       = [2 4; 4 8; 8 12; 12 24; 24 48; 48 96; 96 192];
% i_fit     = [1 100];
% i_calc    = [0.5 120];
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v         = 3;
% v_rawdata = 6;
% is_src    = 0;
% fsample   = 300;
% foi       = [unique(round(2.^[1:.25:7]))-1; unique(round(2.^[1:.25:7]))+1]';
% i_fit     = [1 30];
% i_calc    = [0.5 50];
% --------------------------------------------------------
% VERSION 4
% --------------------------------------------------------
% v         = 4;
% v_rawdata = 6;
% is_src    = 0;
% fsample   = 300;
% foi       = [7 13];
% i_fit     = [3 50];
% i_calc    = [1 70];
% --------------------------------------------------------
% VERSION 5
% --------------------------------------------------------
% v         = 5;
% v_rawdata = 6;
% is_src    = 0;
% fsample   = 300;
% foi       = [2 4; 4 7; 7 13; 13 30; 30 60; 60 120; 100 180];
% i_fit     = [1 15];
% i_calc    = [.8 25];
% --------------------------------------------------------
% VERSION 6
% --------------------------------------------------------
% v         = 6;
% v_rawdata = 6;
% is_src    = 0;
% fsample   = 300;
% foi       = [2 4; 4 7; 7 13; 13 30; 30 60; 60 120; 100 180];
% i_fit     = [3 30];
% i_calc    = [2 50];
% --------------------------------------------------------
% VERSION 7
% --------------------------------------------------------
v         = 7;
v_rawdata = 6;
is_src    = 0;
fsample   = 300;
foi       = [2 4; 4 7; 7 13; 13 30; 30 60; 60 120; 100 180];
i_fit     = [3 70];
i_calc    = [2 100];
% --------------------------------------------------------

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
outdir = '/home/gnolte/neuconn/meg_out/rest/dfa/';
mridir = '/home/gnolte/neuconn/mri_data/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%
for ifoi = 1 : length(foi)
  
  for icond = 1 : 2
    for isubj = 1 : 49
      
      if ~exist(sprintf([outdir 'nc_sens_dfa_s%d_c%d_f%d_v%d_processing.txt'],isubj,icond,ifoi,v))
        system(['touch ' outdir sprintf('nc_sens_dfa_s%d_c%d_f%d_v%d_processing.txt',isubj,icond,ifoi,v)]);
      else
        continue
      end
      
      for iblock = 1 : 2
        
        disp(sprintf('Processing subj %d block %d cond %d freq %d ...',isubj,iblock,icond,ifoi));
        
        % -----------------------------------------------------------------
        % Load data
        % -----------------------------------------------------------------
        if exist(sprintf('/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v_rawdata))
          disp(sprintf('Loading MEG data ...'));
          load(sprintf('/home/gnolte/neuconn/meg_out/rest/m1/proc/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v_rawdata));
        else
          continue
        end
        % -----------------------------------------------------------------
        
        if foi(ifoi,1)<30
          clear data_hi
          [mydata,epleng] = megdata2mydata(data_low);
        else
          clear data_low
          [mydata,epleng] = megdata2mydata(data_hi);
        end
        
        siginfo = nbt_Info;
        siginfo.converted_sample_frequency = fsample;
        
        % compute bp-filtered signal
        tmp    = nbt_filter_fir(mydata,foi(ifoi,1),foi(ifoi,2),siginfo.converted_sample_frequency,4/foi(ifoi,1));
        ampenv = abs(hilbert(tmp)); clear tmp mydata
        
        dfa               = nbt_doDFA(ampenv, siginfo, i_fit,i_calc,0.5,0,0,[]);
        
        par.dfa	= dfa.MarkerValues;
        
        %         cfg             = [];
        %         cfg.method      = 'mtmfft';
        %         cfg.output      = 'pow';
        %         cfg.taper       = 'hanning';
        %         cfg.channel     = {'MEG'};
        %         cfg.foilim      = [8 12];
        %         cfg.keeptrials  = 'yes';
        %         dat             = ft_freqanalysis(cfg, data_low);
        %
%         par.pow(:,iblock) = squeeze(nanmean(dat.powspctrm,3));
        
        save(sprintf([outdir 'nc_sens_dfa_s%d_b%d_c%d_f%d_v%d.mat'],isubj,iblock,icond,ifoi,v),'par','-v7.3');
        
        clear par dfa data_low dat ampenv
        
      end
    end
  end
end

error('STOP')

