%% NEUCONN RESTING-STATE VIEW ICA COMPONENTS
% v2: difference to older version: ica was computed w/o the use of fieldtrip. 
% therefore, data needs to be party restructured.

% views ica components in gui
% components can be marked for later rejection
% thomas pfeffer, february 2014

clear all
restoredefaultpath

% =========================================================================
% DEFINITIONS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% M1 or M2?


% icond = 1, s48, b2 - whatup???
% icond = 1, s43 ???
% -------------------------------------------------------------------------
% VERSION 4
% -------------------------------------------------------------------------
v     = 6; % ICA version
v_art = 4; % preproc version
split = 0;
hc    = 0;
cfg.method = 'runica';
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
m = 1;
% -------------------------------------------------------------------------

indir   = sprintf('/home/gnolte/neuconn/meg_out/rest/m%d/proc/ica/',m);
datadir = sprintf('/home/gnolte/neuconn/meg_data/m%d/',m);
outdir  = sprintf('/home/gnolte/neuconn/meg_out/rest/m%d/proc/ica/',m);

addpath /home/gnolte/neuconn/matlab/rest/
addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults

if m == 1 
  allsubj = nc_allsubj_start(m);
elseif m == 2
  allsubj = nc_allsubj_start(m);
end

datadir = '/home/tpfeffer/NEUCONN/out_ana/M1/';

addpath('/home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/')

ft_defaults


NSUBJ = 60;

% =========================================================================
% =========================================================================
% DO NOT MODIFY FOLLOWING SEGMENTS
% =========================================================================
% =========================================================================

if hc == 0
  % PATIENTS
  icond = 1;
  strseg = '';
else
  % CONTROLS
  icond = 2;
  strseg = '_hc';
end

total_err = 0;

for iblock = 1 : 2
  for is = 1 : NSUBJ
    ex(iblock,is) = exist([outdir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d.mat',icond,is,iblock,v)])>0;
  end
end

subj_col = {'r';'g'};

while total_err == 0 
  back = 1;
  while back == 1
  back = 0;
  err = 0;
  % INTERFACE
  figure('units','normalized','outerposition',[0.3 0.3 0.5 0.5]);
  
  for is = 1 : NSUBJ
    subj_cb{is} = '';
  end
  
  if split == 0
    nseg = 1;
  else
    nseg = 2;
  end
  
  for is = 1 : 20
    cnt = 0;
    for ibl = 1 : 2; for ifr = 1 : 2; for ise = 1 : nseg
      cnt = cnt + 1;
      a(cnt) = exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,ifr,is,ibl,ise,v));    
    end; end; end; 
    if ~any(a<1)
      subj_col = {'r';'g'};
    else
      subj_col = {'r';'w'};
    end
    subj{is} = eval(sprintf('uicontrol(''Units'',''normalized'',''Position'',[-0.05+is/20 0.90 0.05 0.05],''Style'',''pushbutton'',''Backgroundcolor'',subj_col{ex(iblock,is)+1},''String'',''S%d'',''Callback'',''isubj = %d'')',is,is));
  end
  for is = 21 : 40
    cnt = 0;
    for ibl = 1 : 2; for ifr = 1 : 2; for ise = 1 : nseg
      cnt = cnt + 1;
      a(cnt) = exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,ifr,is,ibl,ise,v));    
      end; end; end;
    if ~any(a<1)
      subj_col = {'r';'g'};
    else
      subj_col = {'r';'w'};
    end
    subj{is} = eval(sprintf('uicontrol(''Units'',''normalized'',''Position'',[-0.05+(is-20)/20 0.85 0.05 0.05],''Style'',''pushbutton'',''Backgroundcolor'',subj_col{ex(iblock,is)+1},''String'',''S%d'',''Callback'',''isubj = %d; clear comp_var'')',is,is));
  end 
  if NSUBJ > 40 
    for is = 41 : NSUBJ
      cnt = 0;
      for ibl = 1 : 2; for ifr = 1 : 2; for ise = 1 : nseg
        cnt = cnt + 1;
        a(cnt) = exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,ifr,is,ibl,ise,v));    
        end; end; end;
      if ~any(a<1)
        subj_col = {'r';'g'};
      else
        subj_col = {'r';'w'};
      end
      subj{is} = eval(sprintf('uicontrol(''Units'',''normalized'',''Position'',[-0.05+(is-40)/20 0.80 0.05 0.05],''Style'',''pushbutton'',''Backgroundcolor'',subj_col{ex(iblock,is)+1},''String'',''S%d'',''Callback'',''isubj = %d; clear comp_var'')',is,is));
    end 
  end
  bl1 = uicontrol('Units','normalized','Position',[0.38 0.67 0.1 0.1],'Style','pushbutton','String','Block 1','Callback','iblock = 1;');
  bl2 = uicontrol('Units','normalized','Position',[0.52 0.67 0.1 0.1],'Style','pushbutton','String','Block 2','Callback','iblock = 2;');
  
  seg1 = uicontrol('Units','normalized','Position',[0.1 0.3 0.3 0.2],'Style','pushbutton','String','Low Frequencies','Callback','which_freq = 1; close');
  seg2 = uicontrol('Units','normalized','Position',[0.6 0.3 0.3 0.2],'Style','pushbutton','String','High Frequencies','Callback','which_freq = 2; close');
  
  abort = uicontrol('Units','normalized','Position',[0.4 0.01 0.2 0.1],'Style','pushbutton','String','QUIT','Callback','close; err = 1;');

  uiwait

  if err == 1
    error('Aborted!');
  end

  % LOAD DATA
  load([indir sprintf('nc_preproc_ica_c%d_s%d_b%d_v%d.mat',icond,isubj,iblock,v)])
  if hc == 0 && split == 0 && v< 4
    load(sprintf([indir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v_art));
  elseif hc == 1 && split == 0 && v< 4
    load(sprintf([indir 'nc_preproc_artifacts_c%d_s%d_b%d_v%d.mat'],icond,isubj,iblock,v_art));
  end
  
  % IF DATA IN ONE SEGMENT, CREATE COMP VECTOR FIRST
  if split == 0 && v<4
    comp_low.fsample    = data_low.fsample;  comp_hi.fsample   = data_hi.fsample;
    if v ~= 4
    comp_low.unmixing   = W1;                comp_hi.unmixing  = W2;
    comp_low.topo       = A1;                comp_hi.topo      = A2;
    else
    comp_low.unmixing   = A1*W1;             comp_hi.unmixing  = A2*W2;
    comp_low.topo       = pinv(A1*W1);       comp_hi.topo      = pinv(A1*W1);
    end
    comp_low.topolabel  = data_low.label;    comp_hi.topolabel = data_hi.label;
    comp_low.time       = data_low.time;     comp_hi.time      = data_hi.time;
    comp_low.grad       = data_low.grad;     comp_hi.grad      = data_hi.grad;
    
    comp_low.trial{1} = W1*data_low.trial{1};
    comp_hi.trial{1}  = W2*data_hi.trial{1};
  end

  if split == 1
    for iseg = 1 : 2
      if exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,which_freq,isubj,iblock,iseg,v))   
        col(iseg) = 'r';
      else
        col(iseg) = 'g';
      end
    end
  elseif split == 0
    if exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_v%d.mat'],icond,which_freq,isubj,iblock,v))   
      col = 'r';
    else
      col = 'g';
    end
  end
  figure('units','normalized','outerposition',[0.3 0.3 0.5 0.5]);

  if split == 1
    seg1 = uicontrol('Units','normalized','Position',[0.1 0.3 0.3 0.2],'Style','pushbutton','Backgroundcolor',col(1),'String','Segment 1','Callback','iseg = 1; close;');
    seg2 = uicontrol('Units','normalized','Position',[0.6 0.3 0.3 0.2],'Style','pushbutton','Backgroundcolor',col(2),'String','Segment 2','Callback','iseg = 2; close;');
  else
    seg1 = uicontrol('Units','normalized','Position',[0.1 0.3 0.3 0.2],'Style','pushbutton','Backgroundcolor',col(1),'String','Segment 1','Callback','iseg = 1; close;');
  end    
  abort= uicontrol('Units','normalized','Position',[0.4 0.01 0.2 0.1],'Style','pushbutton','String','Quit','Callback','close; err = 1;');
  back = uicontrol('Units','normalized','Position',[0.4 0.85 0.2 0.1],'Style','pushbutton','String','Back','Callback','close; back = 1;'); 
  
  uiwait
  if err == 1
    error('Aborted!');
  end
  end

  if which_freq == 1
    hilow = 'low';
  else
    hilow = 'hi';
  end

  subpl = 4;
  rej_comp = zeros(80,1);
  
  if exist(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,which_freq,isubj,iblock,iseg,v))   
    load(sprintf([outdir 'nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat'],icond,which_freq,isubj,iblock,iseg,v));
  end
  
  l = 0;
  if split == 0 
    mean_spec = eval(sprintf('nanmean(comp_%s.trial{1},1);',hilow));
  elseif split == 1
    mean_spec = eval(sprintf('nanmean(comp_%s_%d.trial{1},1);',hilow,iseg));
  end
  
  cnt = 1;
  il = 0;

while il < subpl
    
  il = il + 1;
  i = (cnt-1)*subpl+il;
  
  % COMPUTE VARIANCE FOR 2s-WINDOWS
  if split == 1
    smax =  eval(sprintf('floor(size(comp_%s_%d.trial{1},2)/100);',hilow,iseg));
    for s = 1 : smax
      eval(sprintf('comp_var_%d(i,s)=var(comp_%s_%d.trial{1}(i,(s-1)*100+1:s*100));',iseg,hilow,iseg));
    end
  elseif split == 0
    smax =  eval(sprintf('floor(size(comp_%s.trial{1},2)/100);',hilow));
    for s = 1 : smax
      eval(sprintf('comp_var(i,s)=var(comp_%s.trial{1}(i,(s-1)*100+1:s*100));',hilow));
    end
  end
  % -----------------
  
  if mod(i-1,subpl)==0
    if exist('manpos')
      figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    else
      figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    end
    l = l + 1;
  end
  
  smo = 100;
  steps = 20;
  Fs = 300;
  if split == 1
    N = eval(sprintf('length(comp_%s_%d.trial{1}(i,:));',hilow,iseg));
    xdft = eval(sprintf('fft(comp_%s_%d.trial{1}(i,:));',hilow,iseg));
  elseif split == 0
    N = eval(sprintf('length(comp_%s.trial{1}(i,:));',hilow));
    xdft = eval(sprintf('fft(comp_%s.trial{1}(i,:));',hilow));
  end
  xdft = xdft(1:N/2+1);
  psdx = (1/(Fs*N)).*abs(xdft).^2;
  psdx(2:end-1) = 2*psdx(2:end-1);
  
  j = 1;
  k = 1;
  while j < length(psdx)-smo
    smoothed(k)=mean(psdx(j:j+smo));
    j = j + steps;
    k = k + 1;
  end
  
  freq = linspace(0,Fs/2,size(smoothed,2));
  strt = find(freq > 4,1,'first');
  
  % PLOT POWER SPECTRUM
  subcomp{1}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-2);
  plot(log10(freq(strt:end)),log10(smoothed(strt:end))); grid on;
  set(gca,'TickDir','out','XTickLabel',round([10.^(0.8:0.2:2.2)]))
  xlabel('Frequency (Hz)'); ylabel('(dB/Hz)');
  axis tight
  
  % PLOT VARIANCE OVER TIME
  subcomp{2}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3-1);
  if split == 1
    eval(sprintf('scatter(1:smax,comp_var_%d(i,:),''k.'');',iseg))
  elseif split == 0
    eval(sprintf('scatter(1:smax,comp_var(i,:),''k.'');'))
  end    
  xlabel('Time'); ylabel('Variance');
  axis tight
  
  % PLOT COMPONENT TOPOGRAPHY
  subcomp{3}{il} = subplot(subpl,3,(i-(l-1)*subpl)*3);
  cfg = [];
  cfg.component = [i];       % specify the component(s) that should be plotted
  cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
  cfg.comment   = 'no';
  cfg.highlight = 'off';
  cfg.marker    = 'off';
  if split == 1
    eval(sprintf('ft_topoplotIC(cfg, comp_%s_%d);',hilow,iseg));
  elseif split == 0
    eval(sprintf('ft_topoplotIC(cfg, comp_%s);',hilow));
  end
  
  if mod(i,subpl)==0 || i == 80
    
    pos = [0.76 0.73 0.075 0.035; ...
      0.76 0.51 0.075 0.035; ...
      0.76 0.29 0.075 0.035; ...
      0.76 0.07 0.075 0.035];
    
    rej_callback1 = ['if (rej_comp(i-3) == 0), set(rej1,''Backgroundcolor'',''r''),rej_comp(i-3)=1;' ...
      'else set(rej1,''Backgroundcolor'',''g''), rej_comp(i-3)=0; end'];
    rej_callback2 = ['if (rej_comp(i-2) == 0), set(rej2,''Backgroundcolor'',''r''),rej_comp(i-2)=1;' ...
      'else set(rej2,''Backgroundcolor'',''g''), rej_comp(i-2)=0; end'];
    rej_callback3 = ['if (rej_comp(i-1) == 0), set(rej3,''Backgroundcolor'',''r''),rej_comp(i-1)=1;' ...
      'else set(rej3,''Backgroundcolor'',''g''), rej_comp(i-1)=0; end'];
    rej_callback4 = ['if (rej_comp(i-0) == 0), set(rej4,''Backgroundcolor'',''r''),rej_comp(i-0)=1;' ...
      'else set(rej4,''Backgroundcolor'',''g''), rej_comp(i-0)=0; end'];
    
    for ibgc = 1 : 4
      if rej_comp(i+(ibgc-4)) == 1
        bgc{ibgc} = 'r';
      else
        bgc{ibgc} = 'g';
      end
    end
    
    rej1 = uicontrol('Units','normalized','Position',pos(1,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{1},'Callback',rej_callback1);
    rej2 = uicontrol('Units','normalized','Position',pos(2,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{2},'Callback',rej_callback2);
    rej3 = uicontrol('Units','normalized','Position',pos(3,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{3},'Callback',rej_callback3);
    rej4 = uicontrol('Units','normalized','Position',pos(4,:),'Style','pushbutton','String','Keep','Backgroundcolor',bgc{4},'Callback',rej_callback4);
    
    if split == 1 
      tc1_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-3]; ft_databrowser(cfg, comp_%s_%d);',hilow,iseg);
      tc2_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-2]; ft_databrowser(cfg, comp_%s_%d);',hilow,iseg);
      tc3_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-1]; ft_databrowser(cfg, comp_%s_%d);',hilow,iseg);
      tc4_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-0]; ft_databrowser(cfg, comp_%s_%d);',hilow,iseg);
    elseif split == 0
      tc1_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-3]; ft_databrowser(cfg, comp_%s);',hilow);
      tc2_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-2]; ft_databrowser(cfg, comp_%s);',hilow);
      tc3_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-1]; ft_databrowser(cfg, comp_%s);',hilow);
      tc4_cb = sprintf('cfg = [];  cfg.layout = ''CTF275.lay''; cfg.viewmode = ''vertical''; cfg.channel = [i-0]; ft_databrowser(cfg, comp_%s);',hilow);
    end
    
    tc1 = uicontrol('Units','normalized','Position',[0.86 0.73 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc1_cb);
    tc2 = uicontrol('Units','normalized','Position',[0.86 0.51 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc2_cb);
    tc3 = uicontrol('Units','normalized','Position',[0.86 0.29 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc3_cb);
    tc4 = uicontrol('Units','normalized','Position',[0.86 0.07 0.075 0.035],'Style','pushbutton','String','Timecourse','Callback',tc4_cb);
    
    % SAVE COMPONENTS
    % ------------------------------------------------
    sc_cb1 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{1},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{1},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{1},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots%s_f%d_s%d_b%d_seg%d_comp%d_v%d.pdf'']); close(h)',strseg,which_freq,isubj,iblock,iseg,1,v);
    sc_cb2 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{2},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{2},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{2},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots%s_f%d_s%d_b%d_seg%d_comp%d_v%d.pdf'']); close(h)',strseg,which_freq,isubj,iblock,iseg,2,v);
    sc_cb3 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{3},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{3},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{3},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots%s_f%d_s%d_b%d_seg%d_comp%d_v%d.pdf'']); close(h)',strseg,which_freq,isubj,iblock,iseg,3,v);
    sc_cb4 = sprintf('h = figure; set(h,''Position'',[200 200 1000 300]); set(h,''Units'',''inches''); screenposition = get(h,''Position''); set(h, ''PaperPosition'',[0 0 screenposition(3:4)],''PaperSize'',[screenposition(3:4)]); new = copyobj(subcomp{1}{4},h); set(new,''Position'',[.05 .1 0.25 0.85]); new = copyobj(subcomp{2}{4},h);set(new,''Position'',[.35 .1 0.25 0.85]); set(new,''LineWidth'',2); new = copyobj(subcomp{3}{4},h); set(new,''Position'',[.55 .05 0.5 0.95]); print(h,''-dpdf'',[outdir ''/ica_plots/nc_icaplots%s_f%d_s%d_b%d_seg%d_comp%d_v%d.pdf'']); close(h)',strseg,which_freq,isubj,iblock,iseg,4,v);
    
    savecomp1 = uicontrol('Units','normalized','Position',[0.86 0.78 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb1);
    savecomp2 = uicontrol('Units','normalized','Position',[0.86 0.56 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb2);
    savecomp3 = uicontrol('Units','normalized','Position',[0.86 0.34 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb3);
    savecomp4 = uicontrol('Units','normalized','Position',[0.86 0.12 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',sc_cb4);
    % ------------------------------------------------
  
    if i > 4
      prev = uicontrol('Units','normalized','Position',[0.1 0.01 0.075 0.05],'Style','pushbutton','String','Prev','Callback','cnt = cnt - 1; il = 0; l = l - 2; manpos = get(gcf,''Position''); close');
    else
      prev = uicontrol('Units','normalized','Position',[0.1 0.01 0.075 0.05],'Style','pushbutton','String','');
    end
    if i < 77
      next = uicontrol('Units','normalized','Position',[0.2 0.01 0.075 0.05],'Style','pushbutton','String','Next','Callback','cnt = cnt + 1; il = 0; close');
    else
      next = uicontrol('Units','normalized','Position',[0.2 0.01 0.075 0.05],'Style','pushbutton','String','');
    end
    
    save_callback = ['idx = find(rej_comp==1); save(sprintf([outdir ''nc_rejected_comps_c%d_f%d_s%d_b%d_seg%d_v%d.mat''],icond,which_freq,isubj,iblock,iseg,v),''idx'',''rej_comp''); close;'];    
    
    s = uicontrol('Units','normalized','Position',[0.90 0.01 0.075 0.05],'Style','pushbutton','String','Save','Callback',save_callback);
    err = 0;
    quit = uicontrol('Units','normalized','Position',[0.80 0.01 0.075 0.05],'Style','pushbutton','String','Quit','Callback','close; err = 1;');
    orient landscape
    uiwait
    
    if err == 1
      error('Abort!')
    end  
  end
end
end


