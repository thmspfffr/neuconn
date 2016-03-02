%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_plot_src_powcorr_all2all_mean

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

ft_defaults

indir   = '/home/gnolte/neuconn/meg_out/rest/src/';
outdir = '/home/gnolte/neuconn/meg_out/rest/conn/';
mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';

% call subj info startup function
allsubj = nc_allsubj_start;

%%


for ifoi = 1:60
  for icond = 1 : 2
    cnt = 0;
    
    for isubj = 1 : 40
      if ~exist(sprintf([outdir 'nc_src_powcorr_all2all_mean_s%d_c%d_f%d_v%d_processing.txt'],isubj,icond,ifoi,v));
        system(['touch ' outdir sprintf('nc_src_powcorr_all2all_mean_s%d_c%d_f%d_v%d_processing.txt',isubj,icond,ifoi,v)]);
      else
        continue
      end
      
      % exclude excluded subjects
      if cell2mat(allsubj{icond}(isubj,3)) == 0
        continue
      end
      
      % average the two blocks
      for iblock = 1 : 2
        if exist([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
          
          load([indir sprintf('nc_powcorr_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifoi,v)])
          br = 0;
          res(:,:,iblock) = (resout+resout')./2; clear resout
          
        else
          br = 1;
          break
        end
      end
      
      % break & continue with next subject
      if br
        continue
      end
      disp(sprintf('Processig s%df%dc%d',isubj,ifoi,icond));
      % count only if subject data exists
      cnt = cnt + 1;
      %       allpow(:,:,icond,cnt)	= nanmean(res,3); clear res
      c = nanmean(res(:)); clear res
      save(sprintf([outdir 'nc_src_powcorr_all2all_mean_s%d_c%d_f%d_v%d.mat'],isubj,icond,ifoi,v),'c','-v7.3');
    end
  end
end

%% load correlations and plot

for ifoi = 1:60
  disp(sprintf('Processing f%d',ifoi))
  for icond = 1 : 2
    cnt = 0;
    for isubj = 1 : 40
      if exist(sprintf([outdir 'nc_src_powcorr_all2all_mean_s%d_c%d_f%d_v%d.mat'],isubj,icond,ifoi,v))
        cnt = cnt + 1;
        load(sprintf([outdir 'nc_src_powcorr_all2all_mean_s%d_c%d_f%d_v%d.mat'],isubj,icond,ifoi,v));
        all_corr(ifoi,icond,isubj) = c;
      else
        continue
      end      
    end
  end
end

%%
figure; hold on;
set(gcf,'color','w')

sem(1,:) = nanstd(all_corr(:,1,:),[],3)/sqrt(16);
sem(2,:) = nanstd(all_corr(:,2,:),[],3)/sqrt(16);

X1 = [log2(1:60) fliplr(log2(1:60))];
Y1 = [nanmean(all_corr(:,1,:),3)'-sem(1,:) fliplr(nanmean(all_corr(:,1,:),3)'+sem(1,:))];

X2 = [log2(1:60) fliplr(log2(1:60))];
Y2 = [nanmean(all_corr(:,2,:),3)'-sem(2,:) fliplr(nanmean(all_corr(:,2,:),3)'+sem(2,:))];

a = plot(log2(1:60),nanmean(all_corr(:,1,:),3),'color',[1 .55 0.2],'LineWidth',4);
b = patch(X1,Y1,[1 .55 0.2],'EdgeColor','none'); alpha(b,.4);

c = plot(log2(1:60),nanmean(all_corr(:,2,:),3),'color',[.11 .56 1],'LineWidth',4);
d = patch(X2,Y2,[.11 .56 1],'EdgeColor','none'); alpha(d,.4);

set(gca,'TickDir','out','XTick',log2([2,4,8,16,32,64]),'XTickLabel',[2,4,8,16,32,64])

xlabel('Carrier Frequency (in Hz)');
ylabel('Correlation');
title('Mean Correlation');

legend([a,c],'MS','HC')

set(gcf,'renderer','painters')


