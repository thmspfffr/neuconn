function [cs,coh,nave]=nc_cscoh_wavelet(data,f,para)

% usage: [cs,coh,nave]=data2cs_event(data,segleng,segshift,epleng,maxfreqbin,para)
%
% calculates cross-spectra and coherency from data for event-related measurement
% input:
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% f: frequency of interest
% para: optional structure:
%       para.segave=0  -> no averaging across segments
%       para.segave neq 0 -> averaging across segments (default is 1)% \
%       para.subave =1 subtracts the average across epochs,
%       para.subave ~= 1 -> no subtraction (default is 1)
%       IMPORTANT: if you just one epoch (e.g. for continuous data)
%         set para.subave=0
%
%       -> averaging across segments (default is 0)
%       para.proj must be a set of vector in channel space,
%       if it exists then the output raw contains the single trial
%       Fourier-transform in that channel
%       para.zeropad=n  adds  n zeros at the end of each segment and at the end
%                       of the window. default n=0
%       para.mydetrend=1 (detrends linear trends in all segments; default=0 (no detrending))
%
% output:
% cs: nchan by chan by maxfreqbin by nseg tensor cs(:,:,f,i) contains
%     the cross-spectrum at frequency f and segment i
% coh: complex coherency calculated from cs
% nave: number of averages

subave=0;

if nargin<6
  para=[];
end

segave=1;
mydetrend=0;
proj=[];
zeropad=0;
if isfield(para,'segave')
  segave=para.segave;
end
if isfield(para,'detrend')
  mydetrend=para.detrend;
end
if isfield(para,'proj')
  proj=para.proj;
end
if isfield(para,'subave')
  subave=para.subave;
end
if isfield(para,'zeropad')
  zeropad=para.zeropad;
end
[ndum,npat]=size(proj);

[eplen,nchan]=size(data);
if npat>0
  data=data*proj;
  nchan=npat;
end


if segave==0
  cs=zeros(nchan,nchan,maxfreqbin);
  av=zeros(nchan,maxfreqbin);
else
  cs=zeros(nchan,nchan,maxfreqbin);
  av=zeros(nchan,maxfreqbin);
end

if npat>0
  if segave==0
    cs=zeros(nchan,nchan,maxfreqbin,nep,nseg);
    av=zeros(nchan,maxfreqbin,nep,nseg);
  else
    cs=zeros(nchan,nchan,maxfreqbin,nep);
    av=zeros(nchan,maxfreqbin,nep);
  end
end

if ishanning
  mywindow=repmat(hanning(segleng),1,nchan);
elseif iswavelet
  
  % wavelet analysis
  
  
  
  fd          = 2^0.25;
  f_l         = f/fd;
  f_h         = f*fd;
  s_f         = f/(f/(f_h-f_l)*2);
  s_t         = 1/(2*pi*s_f);
  seglen      = floor(300*2*3*s_t/2)*2;
  segshift    = seglen/2;
  
  nseg = (eplen-seglen)/segshift+1;

  
  A           = 1/sqrt(s_t*sqrt(pi));
  
  wavelet = A * exp(2*1i*pi*f.*time) .* exp(-time.^2 ./(2*s_t^2));
  
  datfft = zeros(size(data,2),n_conv);
  
  for ichan = 1 : size(data,2)
    ichan
    datfft(ichan,:) = conv(wavelet,data(:,ichan));
  end
  
  datfft = datfft(:,half_of_wl+1:end-half_of_wl);
  

    for i = 1 : nseg; %average over all segments;
      
      dataloc = data((i-1)*segshift+1:(i-1)*segshift+seglen,:);
      
      if mydetrend==1;
        dataloc=detrend(dataloc);
      end
      
      time        = -seglen/2:seglen/2-1;
      n_wavelet   = length(time);
      n_conv      = n_wavelet + n_wavelet - 1;
      n_conv_pow2 = pow2(nextpow2(n_conv));
      half_of_wl  = (n_wavelet-1)/2;
      
      for ichan = 1 : nchan
        ichan
        dataconv(ichan,:) = mean(abs(conv(dataloc(:,ichan),wavelet)).^2);
      end
      
%       
%       
%       if zeropad>0
%         dataloc=[dataloc;zeros(zeropad,nchan)];
%       end
%       
%       datalocfft=fft(dataloc.*mywindow);
%       
      
      
      for f=1:maxfreqbin % for all frequencies
        if npat==0
          if segave==0
            cs(:,:,f)=cs(:,:,f)+conj(datfft'*datfft);
            av(:,f)=av(:,f)+conj(datalocfft(f,:)');
          else
            %disp([i,f,size(datalocfft)])
            cs(:,:,f)=cs(:,:,f)+conj(datalocfft(f,:)'*datalocfft(f,:));
            av(:,f)=av(:,f)+conj(datalocfft(f,:)');
          end
        else
          if segave==0
            cs(:,:,f,j,i)=conj(datalocfft(f,:)'*datalocfft(f,:));
            av(:,f,j,i)=conj(datalocfft(f,:)');
          else
            %disp([i,f,size(datalocfft)])
            cs(:,:,f,j)=cs(:,:,f,j)+conj(datalocfft(f,:)'*datalocfft(f,:));
            av(:,f,j)=av(:,f,j)+conj(datalocfft(f,:)');
          end
        end
        
      end
    end
    nave=nave+1;
  end
  
  if segave==0
    cs=cs/nave;
    av=av/nave;
  else
    nave=nave*nseg;
    cs=cs/nave;
    av=av/nave;
  end
  
  for f=1:maxfreqbin
    if subave==1
      if npat==0
        if segave==0
          for i=1:nseg;cs(:,:,f,i)=cs(:,:,f,i)-av(:,f,i)*av(:,f,i)';end;
        else
          cs(:,:,f)=cs(:,:,f)-av(:,f)*av(:,f)';
        end
      else
        if segave==0
          for i=1:nseg;for j=1:nep;
              cs(:,:,f,j,i)=cs(:,:,f,j,i)-av(:,f,j,i)*av(:,f,j,i)';
            end;end;
        else
          for j=1:nep;cs(:,:,f,j)=cs(:,:,f,j)-av(:,f,j)*av(:,f,j)';end
        end
      end
    end
  end
  
  ndim=length(size(cs));
  if ndim==3;
    [n1,n2,n3]=size(cs);
    coh=cs;
    for i=1:n3;
      c=squeeze(cs(:,:,i));
      coh(:,:,i)=c./sqrt(diag(c)*diag(c)');
    end
  elseif ndim==4
    [n1,n2,n3,n4]=size(cs);
    coh=cs;
    for i=1:n3;for j=1:n4;
        c=squeeze(cs(:,:,i,j));
        coh(:,:,i,j)=c./sqrt(diag(c)*diag(c)');
      end;end;
  end
  
  
  
  
  
  
  
  
  return;
  
  
