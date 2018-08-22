clear all; clc;
addpath('../../packages/distmesh/');
addpath('../../packages/Spherical-Harmonic-Transform/');
fsvtk = './vtk/MoonICB';

R0 = 1737.; Rc = 1497;
R = 1737.0 - Rc;
%fd=@(p) dsphere(p,0,0,0,R);
%[p,t]=distmeshsurface(fd,@huniform,0.2*R,R*1.1*[-1,-1,-1;1,1,1]);
[p,t]=distmeshsurface(@sphdist,@huniform,0.08*R,R*1.1*[-1,-1,-1;1,1,1]);
%[p,t]=distmeshsurface(@sphdistMd,@huniform,0.08*R,R*1.2*[-1,-1,-1;1,1,1]);

if 0
F = inversesphH_surf(p);
vtktrisurf(t, p(:,1), p(:,2), p(:,3), 'ICB (km)', fsvtk);
end