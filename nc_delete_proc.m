%% delete processing files without output

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 1;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

ft_defaults

outdir = '/home/gnolte/neuconn/meg_out/rest/src/';

%%

for icond = 1 : 2
  for isubj = 1 : 40
    for iblock = 1 : 2
      
      ifwin = 7;
      
      name_proc = sprintf('nc_src_ana_s%d_c%d_v%d_processing.txt',isubj,icond,v);
      name_out  = sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifwin,v);

      if exist([outdir name_proc]) && ~exist([outdir name_out])

        delete([outdir name_proc])

      end
    
    end   
  end 
end

cd(outdir)
delete *error*

