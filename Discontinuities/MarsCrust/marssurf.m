% build surface mesh
clear all; clc; tic;
addpath('../../packages/distmesh/');
global topo Re SR
tic 
load data0/mars_crust.mat
SR = 4; % 720/180

%Re = 3396.0;
Re = 3389.5;
srf = 0.009; 

R = Re; 
[p,t]=distmeshsurface(@surfdist,@huniform,srf*R,R*1.2*[-1,-1,-1;1,1,1]);
save(['workdata/Msurf_',num2str(size(t,1)),'.mat'],'p','t');

if 0
fsvtk = 'vtk/Marstopo42k_unit';
Fs = surfdist0(p);
p = p/R;
vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs, 'Topo (km)', fsvtk);
end