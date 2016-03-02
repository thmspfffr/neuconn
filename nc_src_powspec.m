%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_src_powspec

clear all

% --------------------------------------------------------
% VERSION 2 (all2all)
% --------------------------------------------------------
v = 3;
m = 1;
% -------------------------------------------------------

restoredefaultpath

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
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/';

% call subj info startup function
allsubj = nc_allsubj_start(m);

% Define frequency range of interest
foi_range = 1:23;


%%
for icond = 1 : 2
  for isubj = 1 : 47
    
    % Only process selected datasets
    if cell2mat(allsubj{icond}(isubj,3)) == 0
      continue
    end
    
    % -----------------------------------------------------------------
    % Loop through frequencies of interest
    % -----------------------------------------------------------------
    for ifreq = 1 : 2
      
      for iblock = 1 : 2
        
        try
          
          % -----------------------------------------------------------------
          % Check whether file has already been processed
          % -----------------------------------------------------------------
          
          if     ~exist(sprintf([outdir 'nc_src_powspec_s%d_c%d_b%d_f%d_v%d_processing.txt'],isubj,icond,iblock,ifreq,v))
            
            system(['touch ' outdir sprintf('nc_src_powspec_s%d_c%d_b%d_f%d_v%d_processing.txt',isubj,icond,iblock,ifreq,v)]);
            
            disp(sprintf('Processing s%d c%d f%d ...',isubj,icond,ifreq));
            disp(sprintf('Processing block %d ...',iblock));
            
            % -----------------------------------------------------------------
            % Load data
            % -----------------------------------------------------------------
            disp(sprintf('Loading MEG data ...'));
            load(sprintf('/home/gnolte/neuconn/meg_out/rest/ica/nc_postproc_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,6));
            % -----------------------------------------------------------------
            
            if ifreq == 1
              clear data_hi;
              data = data_low;
            else
              clear data_low
              data = data_hi;
            end
            
            for ichan = 1 : size(data.trial{1},1)
              
              disp(sprintf('Processing channel %d ...',ichan))
              px(ichan,:) = pwelch(data.trial{1}(ichan,:),[],[],[1:150],300);
              
            end
            
            px = nanmean(px);
          	save([outdir sprintf('nc_src_powspec_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifreq,v)],'px');

            
          else 
            continue
          end
          % -----------------------------------------------------------------
          % Catch error message and save
          % -----------------------------------------------------------------
        catch me
          save([outdir sprintf('nc_src_powspec_s%d_c%d_b%d_f%d_v%d_processing_err.mat',isubj,icond,iblock,ifreq,v)],'me');
          system(['touch ' outdir sprintf('nc_src_powspec_s%d_c%d_b%d_f%d_v%d_processing_err.txt',isubj,icond,iblock,ifreq,v)]);
          warning('Something went wrong!')
          clear data_hi data_low
        end
        
      end
    end
  end
end
%%
clear px_all

for icond = 1 : 2
  for isubj = 1 : 47
    
    isubj
    % Only process selected datasets
    if cell2mat(allsubj{icond}(isubj,3)) == 0
      continue
    end
    
    % -----------------------------------------------------------------
    % Loop through frequencies of interest
    % -----------------------------------------------------------------
    for ifreq = 1 : 2
      
      for iblock = 1 : 2
        if exist(([outdir sprintf('nc_src_powspec_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifreq,v)]))
          
          load([outdir sprintf('nc_src_powspec_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,ifreq,v)]);
          tmp(iblock,:) = px;
          clear px
        else
          continue
        end
      end
      
      if exist('tmp','var')
        px_all(:,isubj,icond,ifreq) = nanmean(tmp)'; clear tmp
      else
      	px_all(:,isubj,icond,ifreq) = nan; 
      end
    end
  end
end

px_all(px_all==0) = nan;
px_all = px_all*10e27;

% px_all(:,6,:,:) = [];

% plot results
%%
h=figure; hold on;
set(h,'color','w')

% low frequencies
clear sem X1 Y1 X2 Y2 

sem(1,:) = nanstd(log10(px_all(1:32,:,1,1)),[],2)/sqrt(16);
sem(2,:) = nanstd(log10(px_all(1:32,:,2,1)),[],2)/sqrt(16);

X1 = [log2(1:32),fliplr(log2(1:32))];
Y1 = [log10(nanmean(px_all(1:32,:,1,1),2))'-sem(1,:) fliplr(log10(nanmean(px_all(1:32,:,1,1),2))'+sem(1,:))];

X2 = [log2(1:32),fliplr(log2(1:32))];
Y2 = [log10(nanmean(px_all(1:32,:,2,1),2))'-sem(2,:) fliplr(log10(nanmean(px_all(1:32,:,2,1),2))'+sem(2,:))];

a = plot(log2(1:32),log10(nanmean(px_all(1:32,:,1,1),2)),'color',[1 .55 0.2],'LineWidth',4);
b = patch(X1,Y1,[1 .55 0.2],'EdgeColor','none'); 

c = plot(log2(1:32),log10(nanmean(px_all(1:32,:,2,1),2)),'color',[.11 .56 1],'LineWidth',4);
d = patch(X2,Y2,[.11 .56 1],'EdgeColor','none'); 

% high frequencies
clear sem X1 Y1 X2 Y2

sem(1,:) = nanstd(log10(px_all(36:64,:,1,2)),[],2)/sqrt(16);
sem(2,:) = nanstd(log10(px_all(36:64,:,2,2)),[],2)/sqrt(16);

X1 = [log2(36:64),fliplr(log2(36:64))];
Y1 = [log10(nanmean(px_all(36:64,:,1,2),2))'-sem(1,:) fliplr(log10(nanmean(px_all(36:64,:,1,2),2))'+sem(1,:))];

X2 = [log2(36:64),fliplr(log2(36:64))];
Y2 = [log10(nanmean(px_all(36:64,:,2,2),2))'-sem(2,:) fliplr(log10(nanmean(px_all(36:64,:,2,2),2))'+sem(2,:))];

a = plot(log2(36:64),log10(nanmean(px_all(36:64,:,1,2),2)),'color',[1 .55 0.2],'LineWidth',4);
b = patch(X1,Y1,[1 .55 0.2],'EdgeColor','none'); 
hold on
c = plot(log2(36:64),log10(nanmean(px_all(36:64,:,2,2),2)),'color',[.11 .56 1],'LineWidth',4);
d = patch(X2,Y2,[.11 .56 1],'EdgeColor','none'); 

set(gca,'TickDir','out','XTick',log2([2,4,8,16,32,64]),'XTickLabel',[2,4,8,16,32,64])

camva(gca,camva+(camva/100000000));

set(gca,'TickDir','out')

legend([a,c],'MS','HC');
title('Power Spectrum')
xlabel('Frequency (in Hz)'); ylabel('log-power');grid off

axis([log2(1) log2(64) log10(.05) log10(120)])


