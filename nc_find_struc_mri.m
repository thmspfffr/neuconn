%% fafe

% nc_find_struc_mri
icond = 1;

if icond == 1
  cond = 'MS';
else
  cond = 'HC';
end
  
indir   = sprintf('/home/sbuetof/NeuConn/rawdata/latest/Neuconn/%s/',cond);
outdir  = sprintf('/home/gnolte/neuconn/mri_data/Raw/%s/',lower(cond));

% addpath('~/Documents/MATLAB/toolboxes/fieldtrip-20130101/')
addpath /home/gnolte/neuconn/matlab/rest/
ft_defaults

%%

for isubj = 40 : 55

  

  f = dir(indir);
  cd(indir)

  for i = 7 : size(f,1)

    if isubj < 10
      strfind(sprintf('0%d',isubj),f(i).name(1:2)) & length(f(i).name)>3;
    else
      strfind(sprintf('%d',isubj),f(i).name(1:2)) & length(f(i).name)>3;
    end

    if ans == 1 

       subjfold = f(i).name;
       mkdir([outdir subjfold])

       d    = dir([f(i).name '/T0/MRT/DICOM/']);
        if ~isempty(d) 
         
       date = d(end).name;
              d    = dir([f(i).name '/T0/MRT/DICOM/' d(end).name '/']);

        end

       if isempty(d) 
         
          disp(sprintf('Trying to copy files to common neuconn folder ...'))
          system(['scp ' f(i).name '/T0/MRT/DICOM/* ' outdir subjfold '/']);

       else
         for id = 1 : 3

           disp(sprintf('Trying to copy files to common neuconn folder ...'))
           system(['scp ' f(i).name '/T0/MRT/DICOM/' date '/' d(end-(id-1)).name '/* ' outdir subjfold '/']);

         end
       end


       break    
    end
  end

  warning(sprintf('No MRIs for subject %d',isubj));
  
end

addpath('~/Documents/MATLAB/toolboxes/fieldtrip-20130101/')
ft_defaults

%%
clear all

for isubj = 2 : 50
  try
% isubj = 1;
icond = 1;
m = 1;
subj = nc_allsubj_start(m);

if icond == 1
  cond = 'MS';
else
  cond = 'HC';
end

outdir  = sprintf('/home/gnolte/neuconn/mri_data/Raw/%s/',lower(cond));

ls(outdir)
% subjnum = subj{icond}{isubj,1}(4:5);
subjnum = isubj;

if isubj < 10
  d = dir([outdir sprintf('0%d',subjnum) '*']);
else
  d = dir([outdir sprintf('%d',subjnum) '*']);
end
n = d.name;
d = dir([outdir n]);

clear bin
for id = 3 : length(d)
  try
  disp(sprintf('Processing image %d ...',id))
  file  = d(id).name;
  hdr   = dicominfo([outdir n '/' file]);
  hdr   = hdr.ProtocolName;
  bin(id-2) = strncmp(hdr,'t1_mprage',9);
  catch me
    continue
  end
end

length(find(bin))  
mkdir([outdir n],'struc')

disp(sprintf('Copying structural images ...'))
for i = find(bin)
  copyfile([outdir n '/' d(i+2).name],[outdir n '/struc/'])
end

disp(sprintf('Compressing DICOMs ...'))
zip(sprintf([outdir 's%d_c%d_struc.zip'],isubj,icond),[outdir n '/struc/'])
  catch me
    isubj
  end
end
% ERROR MESSAGES

% s10c1: the specified file is not in DICOM format
% s17c1: no files
% s18c1: no files
% s19c1: not enough files
% s26c1: no images
% s27c1: no images
% s28c1: no images
% s30c1: empty
% s31c1: empty
% s32c1: 512 images as output
% s33c1: 512 images as output
% s34c1: empty
% s35c1: 512 images as output
% s36c1: 512 images as output
% s37c1: 512 images as output
% s38c1: empty
% s39c1: 512 images as output
% s40c1: 512 images as output
% s42c1: 512 images as output
% s44c1: 512 images as output
  %%
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/
ft_defaults

icond   = 1;
cond    = {'MS';'HC'};
outdir  = sprintf('/home/gnolte/neuconn/mri_data/V2/%s/',cond{icond});
isubj   = 47;

datadir = sprintf([outdir 's%d_c%d_struc/'],isubj,icond);
d       = dir(datadir);
d       = d(3).name;


mri = ft_read_mri([datadir d]);

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
% cfg.parameter = 'rotanatomy';
mri_data = ft_volumerealign(cfg,mri);

save([outdir sprintf('nc_mri_s%d_v2.mat',isubj)],'mri_data');



 