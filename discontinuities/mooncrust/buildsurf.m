% build surface mesh
clear all; clc; tic;
addpath('../../packages/distmesh/');
global srfmap 
tic 
load workdata/moonsurf.mat

R = 1737.15;
[p,t]=distmeshsurface(@surfdist,@huniform,0.3*R,R*1.2*[-1,-1,-1;1,1,1]);


if 0
fsvtk = './vtk/Msurf168k_renormal';
%save('Msurf_168k.mat','p','t');
Fs = surfdist0(p);
%vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs, 'depth (km)', fsvtk);
p = p/R; 
vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs, 'depth (km)', fsvtk);
end