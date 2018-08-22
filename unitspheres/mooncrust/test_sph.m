% tests
clear all; clc; tic;
addpath('../../packages/distmesh/');
addpath('../../packages/Spherical-Harmonic-Transform/');
fsvtk = './vtk/Moon08CMI60';

R = 1737.0;
%fd=@(p) dsphere(p,0,0,0,R);
%[p,t]=distmeshsurface(fd,@huniform,0.2*R,R*1.1*[-1,-1,-1;1,1,1]);
%[p,t]=distmeshsurface(@sphdist,@huniform,0.2*R,R*1.1*[-1,-1,-1;1,1,1]);
%[p,t]=distmeshsurface(@sphdistMd,@huniform,0.03*R,R*1.2*[-1,-1,-1;1,1,1]);
[p,t]=distmeshsurface(@sphdistCMI,@huniform,0.08*R,R*1.2*[-1,-1,-1;1,1,1]);


if 0
save('Moon03CMI60_40k.mat','p','t');
F = inversesphH(p)*sqrt(4*pi)+34;
F = max(0.6,F);
Fs = inversesphH_surf(p)*sqrt(4*pi);
vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs-F, 'depth (km)', fsvtk);
end

toc; 