%% PREPARE HEAD MODEL FOR SOURCE ANALYSIS
% nc_prep_src_ana

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
% v         = 1;
% grid_size = 'medium';
% m         = 1;
% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
v         = 2;
grid_size = 'coarse';
m         = 1;
% --------------------------------------------------------
% VERSION 3
% --------------------------------------------------------
% v         = 3;
% grid_size = 'cortex';
% m         = 1;
% --------------------------------------------------------

restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest/

ft_defaults

indir  = '/home/gnolte/neuconn/meg_data/m1/';
outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
mridir = '/home/gnolte/neuconn/mri_data/'; 

% call subj info startup function
allsubj = nc_allsubj_start(m);

strseg  = {'ms';'hc'};
strseg2 = {'pat';'con'};

%%
for icond = 2 : 2
  for isubj = 1:50
    
    if cell2mat(allsubj{icond}(isubj,3)) == 0
      continue
    else
      
%     try
    if ~exist(sprintf([outdir 'nc_sa_c%d_s%d_v%d_processing.txt'],icond,isubj,v))    
    	system(['touch ' outdir sprintf('nc_sa_c%d_s%d_v%d_processing.txt',icond,isubj,v)]);
    else
      continue
    end
    
    disp(sprintf('Processing c%d s%d ...',icond,isubj));
    
    % ----- ---------------------------------------------------
    % SELECT MRI DATA
    % --------------------------------------------------------

    load sa_meg_template;

    if ~exist(sprintf('/home/gnolte/neuconn/mri_data/%s/nc_mri_s%d_v2.mat',lower(strseg{icond}),isubj))
      
      mri_content = dir(sprintf([mridir '%s/'],strseg{icond})); 
      mri_content = mri_content(3:end);
      dat = 0;

      disp(sprintf('Looking for MRI ...'));

      for imri = 1 : length(mri_content)

%         suffix1 = mri_content(imri).name(regexp(mri_content(3).name,'[.]')-2:regexp(mri_content(3).name,'[.]')-1);
        suffix1 = mri_content(imri).name(11:12);
%         suffix2 = mri_content(imri).name(regexp(mri_content(3).name,'[.]')-5:regexp(mri_content(3).name,'[.]')-4);
        suffix2 = mri_content(imri).name(8:9);
        if strcmp(suffix1,'V2') && isequal(isubj,str2num(suffix2))

          dat = 1;
          mri_data = sprintf([mridir '%s/%s'],strseg{icond},mri_content(imri).name);      

          disp(sprintf('Looking for MRI ... Found!'));

          break

        end  
      end

      if ~dat
        disp(sprintf('Looking for MRI ... Not Found!'));
        continue
      end
    else
      load(sprintf('/home/gnolte/neuconn/mri_data/%s/nc_mri_s%d_v2.mat',lower(strseg{icond}),isubj));
    end
    
    if exist('mri_data','var')
      sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri_data);
    elseif exist('mri','var')
      sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri);
    elseif exist('mri_realigned','var')
      sa_meg1 = nc_mk_sa_meg_mri(sa_meg_template,mri_realigned);
    end
    % --------------------------------------------------------
    % SELECT SUBJECT & BLOCK
    % --------------------------------------------------------
    
    disp(sprintf('Searching MEG-Data ...'));
    clear subj
    
    cont = dir(indir);
    
    for i = 1 : length(cont)-3
      if isubj < 10
        if ~isempty(regexp(cont(i+3).name,sprintf('%s%d%d',strseg2{icond},0,isubj)))
          subj = cont(i+3).name;
          break
        end
      else
        if ~isempty(regexp(cont(i+3).name,sprintf('%s%d',strseg2{icond},isubj)))
          subj = cont(i+3).name;
          break
        end
      end
    end

    meg_data = [indir subj '/'];
    cont     = dir(meg_data);

    cnt  = 0;
    for i = 3 : length(cont)
      if isempty(findstr(cont(i).name,'.ds'))
        continue
      end
      a = ft_read_event([meg_data cont(i).name]);
      a = cell2mat({a.value});
      if length(a)>10
        warning(sprintf('Too many events in block %s',cont(i).name(end-4:end-3)));
      else
        cnt = cnt + 1;
        blocks{cnt} = cont(i).name;
      end
    end

    % --------------------------------------------------------

    for iblock = 1 : 2
      
      clear meg_data
      
      disp(sprintf('Processing block %d MEG-Data ...',iblock));
      
      meg_data    = [indir subj '/' blocks{iblock}];
      sa          = mk_sa_meg_forward(sa_meg1, meg_data);
      if      strcmp(grid_size,'medium')
        
        L            = grid2L(sa.grid_medium_indi,sa.fp_indi);
        sa.L_medium  = L;
        sa.leadfield = 'medium';
        
      elseif  strcmp(grid_size,'coarse')
        
        L            = grid2L(sa.grid_coarse_indi,sa.fp_indi);
        sa.L_coarse  = L;
        sa.leadfield = 'coarse';
        
      elseif strcmp(grid_size,'cortex')
        
        L             = grid2L(sa.grid_cortex3000_indi,sa.fp_indi);
        sa.L_coarse   = L;
        sa.leadfield  = 'cortex';
        
      end

      save([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)],'sa','-v7.3');
      close all
    end
%     catch me
%       save([outdir sprintf('nc_sa_c%d_s%d_v%d_ERROR.mat',icond,isubj,v)],'me');
%       continue
%     end
    end 
    
  end
end
