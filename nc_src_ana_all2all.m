%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_src_ana

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v = 1;
v_out = 2;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

ft_defaults

indir  = '/home/gnolte/neuconn/meg_data/m1/';
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/';

% call subj info startup function
allsubj = nc_allsubj_start;

strseg  = {'ms';'hc'};
strseg2 = {'pat';'con'};

%%
for icond = 1 : 2
  for isubj = 1 : 40
    
    if cell2mat(allsubj{icond}(isubj,3)) == 0
      continue
    end
    
    try
      
      if ~exist(sprintf([outdir 'nc_src_ana_s%d_c%d_v%d_processing.txt'],isubj,icond,v_out))
        system(['touch ' outdir sprintf('nc_src_ana_s%d_c%d_v%d_processing.txt',isubj,icond,v_out)]);
      else
        continue
      end
      
      disp(sprintf('Processing s%d c%d ...',isubj,icond));

      for iblock = 1 : 2
        
        disp(sprintf('Processing block %d ...',iblock));
        
        if ~exist([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)]) && exist([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)])
          
          disp(sprintf('Computing forward model ...'));
          
          % LOAD OUTPUT FROM NC_PREP_SRC_ANA.M
          load([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)]);
          grid=sa.grid_medium;
          
          A=mkfilt_eloreta_v2(L);
          save([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','L','grid');
          
        elseif exist([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)])
          disp(sprintf('Loading forward model ...'));
          load([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','L','grid');
        else
          continue
        end
        
        disp(sprintf('Loading MEG data ...'));
        load(sprintf('/home/gnolte/neuconn/meg_out/rest/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6));
        clear data_hi
        
        [mydata,epleng]=megdata2mydata(data_low);
        
        if ~exist([outdir sprintf('cs_fm_s%d_c%d_b%d_v%d.mat',isubj,icond,iblock)])
          segleng=200;
          segshift=segleng/2;
          maxfreqbin=30;
          [cs, coh, nave]=data2cs_event(mydata,segleng,segshift,epleng,maxfreqbin);
          %%
          disp(sprintf('Saving cross spectrum ...'));
          save([outdir sprintf('cs_fm_s%d_c%d_b%d_v%d.mat',isubj,icond,iblock,v)],'cs','coh','nave');
        else
          disp(sprintf('Loading cross spectrum ...'));
          load([outdir sprintf('cs_fm_s%d_c%d_b%d_v%d.mat',isubj,icond,iblock,v)]);
        end
        
        for ifoi = 1 : maxfreqbin
          
          disp(sprintf('Processing frequency %d ...',ifoi));
          
          [F1,p]=getdipdir(cs(:,:,ifoi),A);
          
%           r=[-40 -40 60]/10;
          % somatosensory
%           r(1,:) = [-42 -26 54]/10;
%           % visual cortex
%           r(2,:) = [-20 -86 18]/10;
%           % auditory
%           r(3,:) = [-54 -22 10]/10;
%           % MPFC
%           r(4,:) = [-3 39 -2]/10;
%           % SMA
%           r(5,:) = [-2 1 51]/10;
%           
%           for p = 1 : length(r)
%             i(p)=findgridpoint(r(p,:),grid);
%           end
          
          segleng=200;
          segshift=segleng/2;
          fsample=300;
          f0=ifoi;
          
          disp(sprintf('Computing amplitude correlations ...'));
          
          [resout] = nc_orthopowercorr(mydata,segleng,segshift,epleng,f0,fsample,F1,F1);
          
          save([outdir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v_out)],'resout');
        end
      end
    catch me
      save([outdir sprintf('nc_src_ana_s%d_c%d_v%d_error.mat',isubj,icond,v_out)],'me');
      system(['touch ' outdir sprintf('nc_src_ana_s%d_c%d_v%d_error.txt',isubj,icond,v)]);
      
    end
  end
end

%   load sa_meg_template;
%
%   foi = 15;
%
%   for icond = 1 : 2
%     cond = allcond{icond};
%     for iblock = 1 : 2
%       load([indir sprintf('nc_powcorr_%s_s%d_b%d_f%d.mat',cond,isubj,iblock,foi)]);
%       tmp_allpow(icond,iblock,:)=resout;
%     end
%     allpow(icond,:) = squeeze(nanmean(tmp_allpow(icond,:,:),2));
%   end
%
%   diffpow = allpow(1,:)-allpow(2,:);
%
%   para = [];
%   para.mydotmarkersize=50;
%   para.orientation='coronal';
%
%   figure;showmri_transp_v3(sa_meg_template.mri,para,[grid resout'],grid(i,:));
%
