%% NC_PLOT_FC_AAL

foi_range = unique(round(2.^[1:.25:7]));

outdir  = '/home/gnolte/neuconn/meg_out/rest/stats/';
plotdir = '/home/gnolte/neuconn/meg_out/rest/plots/';


v = 5;

lim = [nan nan; 0 0.15; 0 0.08; 0 0.01; 0 0.075];

freq = {'Nan - Nan'; '8 - 11 Hz' ; '16 - 27 Hz'; '32 - 54 Hz'; '5 - 7 Hz'};

%%

load(sprintf([outdir 'nc_src_ana_all2all_stat_alpha4sarah_v%d.mat'],v))

pow = nanmean(allcorr,4);

s = [22 36:42, 71:79];

pow(s,:,:) = [];
pow(:,s,:) = [];


figure; set(gcf,'color','white');

subplot(2,2,1); title('FC MS') 
imagesc(tril(pow(:,:,1)),lim(v,:));
axis equal; axis tight
xlabel('AAL REGION'); 
colorbar
title('MS');

subplot(2,2,2);
imagesc(tril(pow(:,:,2)),lim(v,:))
axis equal; axis tight
xlabel('AAL REGION'); 
colorbar
title('HC');
 
subplot(2,2,4); title('MS - HC');
imagesc(tril(pow(:,:,1))-tril(pow(:,:,2)))
axis equal; axis tight
xlabel('AAL REGION'); 
colorbar

title(sprintf('FREQ: %s',freq{v}))

saveas(gcf,sprintf([plotdir 'nc_src_fc_aal_v%d.fig'],v),'fig');
disp('Image saved.')
close


