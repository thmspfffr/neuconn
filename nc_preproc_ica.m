% =========================================================================
% INDEPENDENT COMPONENT ANALYSIS
% =========================================================================

% nc_preproc_ica


clear all
restoredefaultpath

% -------------------------------------------------------------------------
% VERSION 6
% -------------------------------------------------------------------------
v = 4;
v_out = 6;
CFG.method  = 'runica';
split       = 0;
% -------------------------------------------------------------------------


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



% SUBJECT CODES + TASK-REST SEQUENCE + INCLUDE
% 1 = RTRT; 2 = TRTR
%%

for icond = 1 : 1
  
  subj = allsubj{icond};
  
  for isubj = 7 : 7
    %       isubj
    for iblock = 1 : 2
      
      if subj{isubj,3}~=1; continue; end % check if subj is excluded
      %
%       if exist(sprintf([outdir 'nc_preproc_ica_c%d_s%d_b%d_v%d_processing.txt'],icond,isubj,iblock,v_out)) ...
%           && exist(sprintf([outdir 'nc_preproc_ica_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v_out))
%         continue
%         
%       elseif ~exist(sprintf([outdir 'nc_preproc_ica_c%d_s%d_b%d_v%d_processing.txt'],icond,isubj,iblock,v_out)) ...
%           && exist(sprintf([outdir 'nc_preproc_ica_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v_out))
%         system(['touch ' outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,iblock,v_out)]);
%         continue
%         
%       end
      
      try
      disp(sprintf('Processing subject %d, block %d, cond %d',isubj,iblock,icond))
      %         try
      % READ DATA
      load(sprintf([outdir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v));
      system(['touch ' outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d_processing.txt',icond,isubj,iblock,v_out)]);
      catch me
        continue
      end
      %         catch me
      %           continue
      %         end
      % SPLIT DATA IN TWO HALVES
      if split == 1
        split = floor(size(data_low.trial{1},2)/2);
        data_low_1.trial = data_low.trial{1}(:,1:split);
        data_low_1.time  = data_low.time{1}(:,1:split);
        data_low_2.trial = data_low.trial{1}(:,split+1:end);
        data_low_2.time  = data_low.time{1}(:,split+1:end);
        
        data_hi_1.trial = data_hi.trial{1}(:,1:split);
        data_hi_1.time  = data_hi.time{1}(:,1:split);
        data_hi_2.trial = data_hi.trial{1}(:,split+1:end);
        data_hi_2.time  = data_hi.time{1}(:,split+1:end);
        
        % -----------------------------------------------------
        % ICA SETTINGS
        % -----------------------------------------------------
        cfg = [];
        cfg.method = 'runica';
        cfg.numcomponent = 80;
        cfg.runica.pca = 80;
        % -----------------------------------------------------
        % COMPUTE ICA FOR LOW-FREQ
        % -----------------------------------------------------
        % first half
        data_low.trial{1} = data_low_1.trial;
        data_low.time{1}  = data_low_1.time;
        comp_low_1 = ft_componentanalysis(cfg,data_low);
        % second half
        data_low.trial{1} = data_low_2.trial;
        data_low.time{1}  = data_low_2.time;
        comp_low_2 = ft_componentanalysis(cfg,data_low);
        % -----------------------------------------------------
        % COMPUTE ICA FOR HIGH-FREQ
        % -----------------------------------------------------
        % first half
        data_hi.trial{1}  = data_hi_1.trial;
        data_hi.time{1}   = data_hi_1.time;
        comp_hi_1  = ft_componentanalysis(cfg,data_hi);
        % second half
        data_hi.trial{1}  = data_hi_2.trial;
        data_hi.time{1}   = data_hi_2.time;
        comp_hi_2  = ft_componentanalysis(cfg,data_hi);
        
        % -----------------------------------------------------
        % SAVE COMPONENTS
        % -----------------------------------------------------
        save([outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)],'comp_low_1','comp_low_2','comp_hi_1','comp_hi_2','-v7.3')
        
      elseif split == 0
        
        % -----------------------------------------------------
        % ICA SETTINGS
        % -----------------------------------------------------
        cfg = [];
        cfg = CFG;
        % -----------------------------------------------------
        % COMPUTE ICA FOR LOW-FREQ
        % -----------------------------------------------------
        switch cfg.method
          case 'fastica'
            
            [A1,W1] = fastica(data_low.trial{1}.*1e12,'g','pow3','approach','defl','verbose','on','stabilization','on','displayMode','off');
            [dummy,index] = sort(-sum(A1.^2));
            A1   = A1(:,index);
            W1   = W1(index,:);
            
            [A2,W2] = fastica(data_hi.trial{1}.*1e12,'g','pow3','approach','defl','verbose','on','stabilization','on','displayMode','off');
            [dummy,index] = sort(-sum(A2.^2));
            A2 = A2(:,index);
            W2 = W2(index,:);
            
          case 'runica'
            
            % -----------------------------------------------------
            % ICA SETTINGS
            % -----------------------------------------------------
            cfg = [];
            cfg.method = 'runica';
            cfg.numcomponent = 80;
            cfg.runica.pca = 80;
            % -----------------------------------------------------
            % COMPUTE ICA FOR LOW-FREQ
            % -----------------------------------------------------
            % first half
            comp_low = ft_componentanalysis(cfg,data_low);
            % -----------------------------------------------------
            % COMPUTE ICA FOR HIGH-FREQ
            % -----------------------------------------------------
            % first half
            comp_hi = ft_componentanalysis(cfg,data_hi);
            
            %         [A1,W1] = runica(data_low.trial{1}.*1e12,'pca',100);
            %         [dummy,index] = sort(-sum(A1.^2));
            %         A1   = A1(:,index);
            %         W1   = W1(index,:);
            %
            %         [A2,W2] = runica(data_hi.trial{1}.*1e12,'pca',100);
            %         [dummy,index] = sort(-sum(A2.^2));
            %         A2   = A2(:,index);
            %         W2   = W2(index,:);
            
        end
        
        
        
        % -----------------------------------------------------
        % SAVE COMPONENTS
        % -----------------------------------------------------
        save([outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v_out)],'comp_low','comp_hi','-v7.3')
        
      end
    end
  end
end








