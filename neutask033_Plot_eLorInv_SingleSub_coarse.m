% prepare sa structure
% (standard mni slices, lead fields, etc).
startup_NeuconnTask_ae

% task subjects complete
subs=1:49;
nsubs=size(subs,2); % number of subjects

% script settings
pfigs=0;    % print figures?  (0:no, 1:yes);
efigs=1;    % export figures? (0:no, 1:yes);
if efigs==1
    pfigs=1;
end
sv=1;       % save results?   (0:no, 1:yes);
% vb=1;       % verbose script? (0:no, 1:yes);
% dr=1;       % save diary?     (0:no, 1:yes);
eform=0;    % export format   (0:jpg, 1: pdf);
resol='-r600';  % resolution

%% plot options

close all

opl={'FontName'   , 'AvantGarde', ...
    'FontSize', 12, ...
    'Box', 'off',  ...
    };
opgca={'FontName', 'AvantGarde', ...
    'FontSize', 12, ...
    'FontWeight', 'Bold' ...
    'Box', 'off',  ...
    'XColor' , [.3 .3 .3], ...
    'YColor' , [.3 .3 .3]};
oplabel={'FontName', 'AvantGarde', ...
    'FontSize', 12, ...
    'FontWeight', 'Bold', ...
    'Color', [.3 .3 .3]};
optitle={'FontName', 'AvantGarde', ...
    'FontSize', 16, ...
    'FontWeight', 'Bold', ...
    'Color', [.3 .3 .3]};

opcbar={'FontName', 'AvantGarde', ...
    'FontSize', 10, ...
    'FontWeight', 'Bold', ...
    'Color', [.3 .3 .3]};

load('mri_template.mat');

%%


% bands in Hz
bd0=[1 3]; % delta
bd1=[4 7]; % theta
bd2=[8 13]; % alpha
bd3=[14 30]; % beta

% Hz to bin (-> +1)
bands=[bd0;bd1;bd2;bd3]+1;
[nbands, ~]=size(bands);

flag=0;

%for band=1:nbands
for band=3
    
    bd=bands(band,:); % frequency of interest in bin
    
    
    %% Directories
    
    
    % where are the cross-spectra?
    csmatdir=[res_mat_dir '20141210_m2_cs' sh];
    
    % where are the sa files?
    sadir=[datapref 'meg_out' sh 'res_mat_ae' sh 'task' sh '20150408_sa_m2_coarse' sh];
    
    resdir_postfix=['20150409_eLor_coarse_m2_' num2str(bd(1)-1) '-' num2str(bd(2)-1) 'Hz'];
    
    rpicsdir=[res_pics_dir resdir_postfix sh];
    
    % create directories if appropiate
    if efigs==1
        if ~exist(rpicsdir, 'dir')
            mkdir(rpicsdir);
        end
    end
    
    
    % where to save results?
    rmatdir=[res_mat_dir resdir_postfix sh];
    
    
    if sv==1
        if ~exist(rmatdir, 'dir')
            mkdir(rmatdir);
        end
    end
    
    
    
    %% load files
    
    s=2; % subject number
    c=1; % MS patient=1; healthy control=2
    b=1; % block
    m=2; % which measurement time point
    f=1; % low frequency=1; high freuquency = 2;
    
    sid=0;
    for s=subs % subjects
        sid=sid+1;
        for c=1:2   % condition (1:MS, 2:HC)
            
            if c==1     % condition
                pc='pat';
            else
                pc='con';
            end
            
            if s<10
                sub=['0' num2str(s)];
            else
                sub=num2str(s);
            end
            
            %         % sbjInfo
            %         sbmfile=[subInfoDir pc num2str(sub) '_m' num2str(m) '.m'];
            %         if exist(sbmfile, 'file')
            %             run(sbmfile);
            %         else
            %             warning([ sbmfile ' does not exist.' ...
            %                 'Proceeding with next subject!']);
            %             continue
            %         end
            %
            
            for b=1:2   % block
                
                fprintf(['** Processing Subject #: %d, C: %d, B: %d, M: %d , f: %d, bd %d-%d \n'], ...
                    s, c, b, m, f, bd(1), bd(2));
                
                
                safile=['nc_sa' ...
                    '_c' num2str(c) ...
                    '_s' num2str(s) ...
                    '_b' num2str(b) ...
                    '_v4.mat'];
                
                csfile=['nc_cs' ...
                    '_s' num2str(s) ...
                    '_c' num2str(c) ...
                    '_b' num2str(b) ...
                    '_f' num2str(f) '.mat'];
                
                if exist([csmatdir csfile], 'file')
                    load([csmatdir csfile]);
                else
                    warning(['CS for "' csfile '" does not exist']);
                    flag=1;
                    continue
                    
                end
                
                if exist([sadir safile], 'file')
                    load([sadir safile]);
                else
                    warning(['SA for "' safile '" does not exist']);
                    flag=1;
                    continue
                    
                end
                
                
                
                
                [nch,~,nf]=size(cs);
                df=1;
                xf=(0:nf-1)*df;
                
                
                mygrid=sa.grid_coarse_indi;
                L=sa.L_coarse;
                [~, ngrid, ndim]=size(L);
                
                A=sa.A;
                
                csmean=mean((cs(:,:,bd(1):bd(2))), 3);
                
                po_sens=diag(csmean);
                
                po_sour2=zeros(ngrid,ndim);
                % Apply eLoreta to sensor power
                for dim=1:ndim
                    po_sour2(:, dim)=squeeze(A(:,:,dim))'*po_sens;
                end
                
                % fix dipole direction
                po_source=zeros(ngrid,1);
                for v=1:ngrid
                    po_source(v)=norm(squeeze(po_sour2(v,:)));
                end
                
                if b==1
                    pow1=po_source;
                else
                    pow2=po_source;
                end
                
            end % block
            %%
            
            if flag==0
                
                powm=(pow1'+pow2')/2;
                powm=powm./norm(powm);
             
                para=[];
                %para.colormaps={'cool'};
                
                fig=figure;
                showmri_transp_v31(mri, para, [sa.grid_coarse powm']);
                
%                 vout=spatfiltergauss(powm',sa.grid_cortex3000,3,sa.cortex10K.vc);
%                 
%                 fig=figure;
%                 subplot(2,2,1);
%                 para=[];
%                 para.myviewdir=[1 0 .5];
%                 showsurface(sa.cortex10K, para, vout);
%                 
%                 colorbar off
%                 %set((get(gca,'title')), 'String', 'text')
%                 title('Right', optitle{:});
%                 subplot(2,2,2);
%                 para=[];
%                 para.myviewdir=[-1 0 .5];
%                 showsurface(sa.cortex10K, para, vout);
%                 colorbar off
%                 title('Left', optitle{:});
%                 subplot(2,2,3);
%                 para=[];
%                 para.myviewdir=[0,1,0.5];
%                 showsurface(sa.cortex10K, para, vout);
%                 colorbar off
%                 title('Front', optitle{:});
%                 subplot(2,2,4);
%                 para=[];
%                 para.myviewdir=[0,-1,0.5];
%                 showsurface(sa.cortex10K, para, vout);
%                 %colorbar off
%                 title('Back', optitle{:});
                
                
                if efigs==1
                    set(gcf,'PaperUnits','centimeters')
                    xSize = 10; ySize = 10;
                    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
                    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
                    set(gcf,'Position',[100 100 xSize*50 ySize*50])
                    
                    figname=[rpicsdir ...
                        'nc' ...
                        '_s' num2str(s) ...
                        '_c' num2str(c) ...
                        'gcoarse_eLorPow_' num2str(bd(1)-1) '-' num2str(bd(2)-1) 'Hz_M2_Task'];
                    
                    
                    if eform==0
                        ex='.jpg';
                        print(fig, '-djpeg', resol, [figname ex]);
                    else
                        set(gcf,'Units','centimeters');
                        screenposition = get(gcf,'Position');
                        set(gcf,...
                            'PaperPosition',[0 0 screenposition(3:4)],...
                            'PaperSize',[screenposition(3:4)]);
                        ex='.pdf';
                        print(fig, '-dpdf', resol, [figname ex]);
                        %saveas(fig, [figname ex]);
                    end
                    fprintf('\nFigure saved: as\n %s \n\n', [figname ex]);
                end
                
                
                %%
            end
            
            flag=0;
            close all
            
            
            %%
            
        end  % condition
    end   % subjects
    
end % bands