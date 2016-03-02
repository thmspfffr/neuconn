function [A A1 po]=mkfilt_dics(L,C,alpha0);
% Calculates spacial filters, and  power over all grid points using lcmv
% Input:
% L: Lead field matrix, NxMx3 matrix for N channels, M voxels and 3 dipole
%    directions
% C:  NxN matrix for N channels covariance matrix or real part of cross-spectrum
%    (The program uses only the real of C)
% alpha0: (relative) regularization parameter. In the algorthm C+alpha*eye(N) is
%        inverted with alpha=alpha0*trace(C)/N
%
% Output
% A  :3-dimensional filter
% A1 : 1-dimensional filter along direction with strongest power
% % po: Mx1 vector for M voxels, po(i) is the power at the i.th voxel along
%          strongest direction

C=real(C);

if nargin<3
    alpha0=.05;
end
 alpha=alpha0*trace(C)/length(C);
 
[nchan nchan]=size(C);
[nchan ns ndum]=size(L);
Cr=C+alpha*eye(nchan);

Crinv=inv(Cr);

A=zeros(nchan,ns,3);
for i=1:ns
    Lloc=squeeze(L(:,i,:));
    A(:,i,:)=reshape((inv((Lloc'*Crinv*Lloc))*Lloc'*Crinv)',nchan,3);
end


po=zeros(ns,1);
A1=zeros(nchan,ns);
for i=1:ns
    Aloc=transpose(squeeze(A(:,i,:)));
    Ploc=Aloc*C*Aloc';
    [u s v]=svd(Ploc);
    po(i)=s(1,1);
    A1(:,i)=Aloc'*u(:,1);
end


    
return;
