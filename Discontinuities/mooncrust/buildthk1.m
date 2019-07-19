% build surface mesh
clear all; clc; tic;
addpath('../../packages/distmesh/');
global thk1map srfmap
tic 
load workdata/moon1thick.mat

load workdata/moonsurf.mat

R = 1737.15;
[p,t]=distmeshsurface(@cmi1dist,@huniform,0.3*R,R*1.2*[-1,-1,-1;1,1,1]);

if 0
fsvtk = './vtk/M1cmi364k';
save('Mcmi1_364k.mat','p','t');
dpth = cmi1dist0(p); 
vtktrisurf(t, p(:,1), p(:,2), p(:,3), dpth, 'deep (km)', fsvtk);
end