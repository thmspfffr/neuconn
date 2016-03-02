%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS

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
addpath /home/gnolte/neuconn/matlab/rest/
addpath /home/gnolte/neuconn/meg_out/rest/dti/

ft_defaults

indir   = '/home/gnolte/neuconn/meg_out/rest/src/';
outdir = '/home/gnolte/neuconn/meg_out/rest/conn/';
mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';

% call subj info startup function
allsubj = nc_allsubj_start;

% somatosensory
r_all(1,:) = [-42 -26 54 42 -26 54]/10;
% visual cortex
r_all(2,:) = [-20 -86 18 20 -86 18]/10;
% auditory
r_all(3,:) = [-54 -22 10 54 -22 10]/10;
% MPFC
r_all(4,:) = [-3 39 -2 3 39 -2]/10;
% SMA
r_all(5,:) = [-2 1 51 2 1 51]/10;

%%
load sa_meg_template;
grid = sa_meg_template.grid_medium;

for idx = 1 : 1%5
  
  for ifoi = 17
    for icond = 1 : 2
      cnt = 0;
      for isubj = 1 : 40
        % exclude excluded subjects
        if cell2mat(allsubj{icond}(isubj,3)) == 0
          continue
        end 
        % average the two blocks
        for iblock = 1 : 2
          if exist([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
            load([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])          
            br = 0;
            % correction factor for correlations and fisher z-transform
            res(:,iblock) = atanh(resout(idx,:)./0.577); clear resout
          else
            br = 1;
            break
          end        
        end
        % break & continue with next subject
        if br
          continue
        end
      	disp(sprintf('Processig s%df%dc%didx%d',isubj,ifoi,icond,idx));
        % count only if subject data exists
        
       	cnt = cnt + 1;
        subj(icond,cnt) = isubj;
        allpow(icond,:,cnt)	= nanmean(res,2); clear res
      end
      allpow_std(icond,:) = nanstd(allpow(icond,:,:),[],3);  
      clear powcorr
    end
     
    pow = nanmean(allpow(1,:,:),3) - nanmean(allpow(2,:,:),3);
    
    i(1)=findgridpoint(r_all(idx,1:3),grid);
    i(2)=findgridpoint(r_all(idx,4:6),grid);

    para = [];
    para.mydotmarkersize=20;
    para.orientation='sagittal';
    %   para.colorlimits = [-0.075 0.075];
  
    M = [nanmean(nanmean(allpow(1,:,:),3)) nanmean(nanmean(allpow(2,:,:),3))];
    Mall = nanmean(M);
    
% =========================================================================
% PROCESS DTI FA VALUES
% =========================================================================
    
    load z_FA_values_06052014
    load subjects_06052014
    
    strseg = {'NCMS';'NCHC';};
    
    num = { '01','02','03','04','05','06','07','08','09','10', ...
            '11','12','13','14','15','16','17','18','19','20', ...
            '21','22','23','24','25','26','27','28','29','30', ...
            '31','32','33','34','35','36','37','38','39','40',};
      
    clear subj_num
    for icond = 1 : 2
      subjcnt = 0;
      for itmp = 1 : size(subj,2)
        
        if subj(icond,itmp) < 10
          arg = '_0';
        else
          arg = '_';
        end
        
        tmp_str = [strseg{icond} arg num2str(subj(icond,itmp))];
        
        a = regexp(subjects,tmp_str,'start');
        a = find(~cellfun(@isempty,a));
        
        if isempty(a)
          subjcnt = subjcnt + 1;
          subj_num(icond,subj(icond,itmp)) = 0;
        else
          subjcnt = subjcnt + 1;
          subj_num(icond,subj(icond,itmp)) = a;
        end
      end
    end
% =========================================================================
       
    corr(:,:) = allpow(:,i(2),:);
   
    for j = 1 : 2
    	for is = 1 : cnt
        z(is,:,j) = (allpow(j,:,is) - nanmean(allpow(j,:,is),2))./nanstd(allpow(j,:,is),[],2);
      end
      
      t  = (nanmean(z(:,:,j))-0)./(nanstd(z(:,:,j))/sqrt(is));
      p  = tpdf(t,15);
      th = fdr(p,0.05);
      
      allpow_plt(j,:) = nanmean(allpow(j,:,:),3);
      allpow_plt(j,~th | allpow_plt(j,:)<M(j)) = 0;
            
    end
    
    h = figure;showmri_transp_v3(sa_meg_template.mri,para,[grid pow'],grid(i,:));
%     saveas(h,sprintf([plotdir 'nc_powcorr_diff_f%d_loc%d_v%d.fig'],ifoi,idx,1));
%     close
    
    lim(1) = min([min(allpow_plt(1,:)) min(allpow_plt(2,:))]);
    lim(2) = max([max(allpow_plt(1,:)) max(allpow_plt(2,:))]);
    
    para.colorlimits = [lim(1) lim(2)];
    
    h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid allpow_plt(1,:)'],grid(i,:));
%     saveas(h,sprintf([plotdir 'nc_powcorr_c%d_f%d_loc%d_v%d.fig'],1,ifoi,idx,1));
%     close
    
    h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid allpow_plt(2,:)'],grid(i,:));
%     saveas(h,sprintf([plotdir 'nc_powcorr_c%d_f%d_loc%d_v%d.fig'],2,ifoi,idx,1),'fig');
%     close
    
  end
end














