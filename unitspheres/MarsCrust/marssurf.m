% build surface mesh
clear all; clc; tic;
addpath('../../packages/distmesh/');
global topo Re SR
tic 
load data0/mars_crust.mat
SR = 4; % 720/180

Re = 3396.0;
R = Re; 
[p,t]=distmeshsurface(@surfdist,@huniform,0.03*R,R*1.2*[-1,-1,-1;1,1,1]);


if 0
fsvtk = 'vtk/Marstopo42k_unit';
%save('workdata/Msurf_377k.mat','p','t');
Fs = surfdist0(p);
p = p/R;
vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs, 'Topo (km)', fsvtk);
end