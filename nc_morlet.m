function y = nc_morlet(cfg)


dt = 1/cfg.fsample;

foi = cfg.foi;
gwidth = cfg.gwidth;
endnsample = cfg.endnsample;
width = cfg.width;

  sf = foi / width;
  st = 1/(2*pi*sf);
  toi2 = -gwidth*st:dt:gwidth*st;
  
  A = 1/sqrt(st*sqrt(pi));
  
  tap = (A*exp(-toi2.^2/(2*st^2)))';
  
  acttapnumsmp = size(tap,1);
  taplen(ifoi) = acttapnumsmp;
  ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* foi(ifoi));
  
  % create wavelet and fft it
  yw = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));

  