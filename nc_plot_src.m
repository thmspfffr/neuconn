%% PLOT POWER IN SOURCE SPACE

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
addpath /home/gnolte/neuconn/matlab/rest/
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';

outdir = '/home/gnolte/neuconn/meg_out/rest/src/';
v_cs      = 2;          % version of cross spectrum
v         = 2;
freq = 1;
% for isubj = 1 

pow = nan(128,1,2);

for isubj = 1 : 49

  
  for icond = 1 : 2
    
    try
    
    for iblock = 1 : 2

%   freq = 1;
  
        fprintf('Loading Data s%d b%d c%d...\n',isubj,iblock,icond)

        load([outdir sprintf('nc_cs_fm_s%d_c%d_b%d_f%d_v%d.mat',isubj,icond,iblock,freq,v_cs)]);

        load([outdir sprintf('nc_sa_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)]);             

        load([outdir sprintf('nc_sa_fm_c%d_s%d_v%d.mat',icond,isubj,v)],'A','grid');

        fprintf('Done.\n',isubj,icond)

        [nch,~,nf]=size(cs);
        df=1;
        xf=(0:nf-1)*df;


        mygrid=sa.grid_coarse_indi;
        L=sa.L_coarse;
        [~, ngrid, ndim]=size(L);

        % A=sa.A;

        pfreq = 10;

        csmean=mean((cs(:,:,pfreq)), 3);


        for i = 1 : size(cs,3)
          pow(i,:,iblock) = mean(diag(cs(:,:,i)));
        end
      end

      powall = squeeze(nanmean(pow,3));

    	fprintf('Plotting...\n',isubj,icond)

      figure; 
      set(gcf,'color','white');

      plot(powall,'LineWidth',5); hold on
      line([pfreq pfreq],[0 max(powall)+0.1*max(powall)],'LineStyle','--','LineWidth',3,'color','r')

      set(gca,'TickDir','out'); box off
      xlabel('FREQUENCY (in Hz)'); ylabel('Power (a.U.)'); title(sprintf('pwrspctrm s%d c%d',isubj,icond))

      saveas(gcf,sprintf([plotdir 'nc_pwrspctrum_s%d_c%d.fig'],isubj,icond),'fig')

      fprintf('Done.\n',isubj,icond)

      clear pow
    
    catch me
    end
  end
end

% 
% po_sens=diag(csmean);
% 
% po_sour2=zeros(ngrid,ndim);
% % Apply eLoreta to sensor power
% for dim=1:ndim
%   po_sour2(:, dim)=squeeze(A(:,:,dim))'*po_sens;
% end

% fix dipole direction
% po_source=zeros(ngrid,1);
% for v=1:ngrid
%   po_source(v)=norm(squeeze(po_sour2(v,:)));
% end
% 
% [F1,p,dipori]=getdipdir(csmean,A);
% 
% load sa_meg_template;
% 
% mri  = sa_meg_template.mri; 
% vc = sa_meg_template.vc;
% 
% % for i = 1 : size(mri.data,3)
% %   mri.data(:,:,i) = rot90(rot90(mri.data(:,:,i)));
% % end
% 
% % powm=po_source;
% % powm=powm./norm(powm);
% 
% para=[];
% %para.colormaps={'cool'};
% 
% fig=figure;
% showmri_transp_v3(mri, para, [grid p]);
%                 