%% nc_src_powcorr_degree

% compute degree across frequencies of interest.

clear all

m = 1;

% --------------------------------------------------------
% VERSION 2 - initial version
% --------------------------------------------------------
% v = 2;
% --------------------------------------------------------
% VERSION 3 - same thresholding for both HC and MS
% --------------------------------------------------------
% v_in          = 8; % amplitude corr with log-scaling
% v             = 4;
% foi_range     = unique(round(2.^[1:.25:7]));
% fixed_thresh  = 0;
% thresh        = 0.15;
% NSUBJ         = 22;
% gridsize      = 'cortex';
% --------------------------------------------------------
% VERSION 3 - same thresholding for both HC and MS
% --------------------------------------------------------
v_in          = 7; % amplitude corr with log-scaling
v             = 5;
foi_range     = unique(round(2.^[1:.25:7]));
fixed_thresh  = 1;
% thresh        = 0.15;
NSUBJ         = 22;
gridsize      = 'coarse';
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

outdir = '/home/gnolte/neuconn/meg_out/rest/conn/';
mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';
% call subj info startup function
allsubj = nc_allsubj_start(m);

%% DEGREE

% load sa_meg_template;
% grid = sa_meg_template.grid_medium;

for ifoi = 1 : length(foi_range)
  
  
  if ~exist(sprintf([outdir 'nc_powcorr_degree_f%d_v%d_processing.txt'],ifoi,v))
    system(['touch ' outdir sprintf('nc_powcorr_degree_f%d_v%d_processing.txt',ifoi,v)]);
  else
    continue
  end
  
  disp(sprintf('Computing f%d ...',ifoi))
  
  load(sprintf([outdir 'nc_src_allpowcorr_f%d_v%d.mat'],ifoi,v_in));
  
  
  if strcmp(gridsize,'coarse')
    th = zeros(2113,2113,2);
  elseif strcmp(gridsize,'medium')
    th = zeros(5003,5003,2);
  elseif strcmp(gridsize,'medium')
    th = zeros(3000,3000,2);
  end
  
  j = 1 : size(th,1);

  m1 = squeeze(nanmean(allpow(j,:,:,1:NSUBJ),2));
  s1 = squeeze(nanstd(allpow(j,:,:,1:NSUBJ),[],2));
  
  m2 = squeeze(nanmean(allpow(:,j,:,1:NSUBJ),1));
  s2 = squeeze(nanstd(allpow(:,j,:,1:NSUBJ),[],1));

  
 
  for il = 1 : 100%size(allpow,1)
        
    jj = j(j~=il);
    
  
    for jl = 1:2113
      
    	disp(sprintf('Computing location %d %d ...',il,jl));

      x = squeeze(allpow(il,jl,:,1:NSUBJ));
      
      z1 = (x - squeeze(m1(jl,:,:)))./squeeze(s1(jl,:,:));
      z2 = (x - squeeze(m2(jl,:,:)))./squeeze(s2(jl,:,:));
      
      [~,p1] = ttest(z1');
      [~,p2] = ttest(z2');
      
      th1(jl,:) = p1 < 0.01/2;
      th2(jl,:) = p2 < 0.01/2;     
      
    end
    
    th1 = th1(jj,:);
    th2 = th2(jj,:);

    % any connection significant?
    th(il,jj,:) = (th1 + th2) > 0; clear th1 th2
    
  end
  
  save([outdir sprintf('nc_powcorr_degree_f%d_v%d.mat',ifoi,v)],'th','-v7.3')

end

%% AVERAGE CORRELATION

% load sa_meg_template;
% grid = sa_meg_template.grid_medium;

for ifoi = 1 : length(foi_range)
  
  
  if ~exist(sprintf([outdir 'nc_powcorr_avgcorr_f%d_v%d_processing.txt'],ifoi,v))
    system(['touch ' outdir sprintf('nc_powcorr_avgcorr_f%d_v%d_processing.txt',ifoi,v)]);
  else
    continue
  end
  
  disp(sprintf('Computing f%d ...',ifoi))
  
  load(sprintf([outdir 'nc_src_allpowcorr_f%d_v%d.mat'],ifoi,v_in));
  
  
  if strcmp(gridsize,'coarse')
    th = zeros(2113,2113,2);
  elseif strcmp(gridsize,'medium')
    th = zeros(5003,5003,2);
  elseif strcmp(gridsize,'medium')
    th = zeros(3000,3000,2);
  end
  
  j = 1 : size(th,1);

  m1 = squeeze(nanmean(allpow(j,:,:,1:NSUBJ),2));
  s1 = squeeze(nanstd(allpow(j,:,:,1:NSUBJ),[],2));
  
  m2 = squeeze(nanmean(allpow(:,j,:,1:NSUBJ),1));
  s2 = squeeze(nanstd(allpow(:,j,:,1:NSUBJ),[],1));

  
 
  for il = 1 : 100%size(allpow,1)
        
    jj = j(j~=il);
    
  
    for jl = 1:2113
      
    	disp(sprintf('Computing location %d %d ...',il,jl));

      x = squeeze(allpow(il,jl,:,1:NSUBJ));
      
      z1 = (x - squeeze(m1(jl,:,:)))./squeeze(s1(jl,:,:));
      z2 = (x - squeeze(m2(jl,:,:)))./squeeze(s2(jl,:,:));
      
      [~,p1] = ttest(z1');
      [~,p2] = ttest(z2');
      
      th1(jl,:) = p1 < 0.01/2;
      th2(jl,:) = p2 < 0.01/2;     
      
    end
    
    th1 = th1(jj,:);
    th2 = th2(jj,:);

    % any connection significant?
    th(il,jj,:) = (th1 + th2) > 0; clear th1 th2
    
  end
  
  save([outdir sprintf('nc_powcorr_avgcorr_f%d_v%d.mat',ifoi,v)],'th','-v7.3')

end
%% plot degree (with statistical threshold)
%
plt = 0;

if plt == 1
  NFREQ = 23;
end

if plt
  if v == 3
    
    for ifoi = 1 : NFREQ
      ifoi
      for icond = 1 : 2
        
        load([outdir sprintf('nc_powcorr_threshold_c%d_f%d_thresh%d_v%d.mat',icond,ifoi,fixed_thresh,v)])
        c(ifoi,icond) = sum(th(:))/(size(th,1)*size(th,1));
        
      end
    end
    
    figure; hold on
    set(gcf,'color','w')
    
    plot(log2(1:NFREQ),c(1:end,1),'color',[1 .55 0.2],'LineWidth',4)
    plot(log2(1:NFREQ),c(1:end,2),'color',[.11 .56 1],'LineWidth',4)
    
    ylabel('Number of connections (in %)');
    xlabel('Carrier frequency (in Hz)');
    title('Degree Distribution');
    
    set(gca,'TickDir','out','XTick',log2([2,4,8,16,32,64]),'XTickLabel',[2,4,8,16,32,64])
    
    b = area(log2(1:NFREQ),c(1:end,2),'FaceColor',[.11 .56 1]);
    a = area(log2(1:NFREQ),c(1:end,1),'FaceColor',[1 .55 0.2]);
    
  elseif v == 4
    
    load([outdir sprintf('nc_powcorr_allpow_thresh%d_v%d.mat',fixed_thresh,v)])
    
    m = nanmean(allpow,3);
    s = nanstd(allpow,[],3)./sqrt(size(allpow,3));
    
    figure; hold on
    set(gcf,'color','w')
    
    plot(1:NFREQ,m(1,:)','color',[1 .55 0.2],'LineWidth',4)
    plot(1:NFREQ,m(2,:)','color',[.11 .56 1],'LineWidth',4)
    
    ylabel('Correlation');
    xlabel('Carrier frequency (Hz)');
    title('Average Correlation');
    
    set(gca,'TickDir','out','XTick',[1,3,7,11,15,19,23],'XTickLabel',[2,4,8,16,32,64 128])
    
    X = [(1:NFREQ),fliplr(1:NFREQ)];
    Y = [m(1,:)-s(1,:),fliplr(m(1,:)+s(1,:))];
    b = fill(X',Y,[1 .55 0.2],'EdgeColor','none'); alpha(b,0.2)
    
    X = [(1:NFREQ),fliplr(1:NFREQ)];
    Y = [m(2,:)-s(2,:),fliplr(m(2,:)+s(2,:))];
    c = fill(X',Y,[.11 .56 1],'EdgeColor','none'); alpha(c,0.2)
    
    set(gca,'xlim',[0 25.101])
    set(gca,'ylim',[-0.0101 0.0502])
    
  end
  
  for ifoi = 1 : NFREQ
    [~, p(ifoi)] = ttest(allpow(1,ifoi,:),allpow(2,ifoi,:),'tail');
  end
  
  
  
  
  %
  %     sem(1,:) = nanstd(c(5:30,1,:),[],3)./sqrt(isubj);
  %     sem(2,:) = nanstd(c(5:30,2,:),[],3)./sqrt(isubj);
  %
  % X1 = [5:30,fliplr(5:30)];
  % Y1 = [nanmean(c(5:30,1,:),3)'-sem(1,:),fliplr(nanmean(c(5:30,1,:),3)'+sem(1,:))];
  %
  % h = figure; hold on;
  % set(h,'color','w')
  %
  % a = plot(5:30,nanmean(c(5:30,1,:),3),'LineWidth',3); hold on
  % b = patch(X1,Y1,'b','EdgeColor','none'); alpha(b,0.2);
  %
  % X1 = [5:30,fliplr(5:30)];
  % Y1 = [nanmean(c(5:30,2,:),3)'-sem(2,:),fliplr(nanmean(c(5:30,2,:),3)'+sem(2,:))];
  %
  % d = plot(5:30,nanmean(c(5:30,2,:),3),'color','r','LineWidth',3); hold on
  % e = patch(X1,Y1,'r','EdgeColor','none'); alpha(e,0.2);
  %
  % legend([a,d],'MS','HC')
  %
  % title('Mean correlation across frequency bands')
  % xlabel('Frequency (in Hz)')
  % ylabel('Mean correlation')
  
  
  
  %   pow = nanmean(allpow(1,:,:),3) - nanmean(allpow(2,:,:),3);
  
  
  %   para = [];
  %   para.mydotmarkersize=20;
  %   para.orientation='sagittal';
  %     para.colorlimits = [-0.075 0.075];
  
  %   M = [nanmean(nanmean(allpow(1,:,:),3)) nanmean(nanmean(allpow(2,:,:),3))];
  %   Mall = nanmean(M);
  
  
  %   for j = 1 : 2
  %     for is = 1 : cnt
  %       z(is,:,j) = (allpow(j,:,is) - nanmean(allpow(j,:,is),2))./nanstd(allpow(j,:,is),[],2);
  %     end
  %
  %     t = (nanmean(z(:,:,j))-0)./(nanstd(z(:,:,j))/sqrt(is));
  %     p = tpdf(t,15);
  %     th = fdr(p,0.05);
  %
  %     allpow_plt(j,:) = nanmean(allpow(j,:,:),3);
  %     allpow_plt(j,~th | allpow_plt(j,:)<M(j)) = 0;
  %
  %   end
  
  
  %   h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid pow'],grid(1,:));
  %     saveas(h,sprintf([plotdir 'nc_powcorr_diff_f%d_loc%d_v%d.fig'],ifoi,idx,1));
  %     close
  
  %   lim(1) = min([min(allpow_plt(1,:)) min(allpow_plt(2,:))]);
  %   lim(2) = max([max(allpow_plt(1,:)) max(allpow_plt(2,:))]);
  
  %   para.colorlimits = [lim(1) lim(2)];
  
  %   h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid allpow_plt(1,:)'],grid(i,:));
  %     saveas(h,sprintf([plotdir 'nc_powcorr_c%d_f%d_loc%d_v%d.fig'],1,ifoi,idx,1));
  %     close
  
  %   h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid allpow_plt(2,:)'],grid(i,:));
  %     saveas(h,sprintf([plotdir 'nc_powcorr_c%d_f%d_loc%d_v%d.fig'],2,ifoi,idx,1),'fig');
  %     close
  
  % end
  
  
end












