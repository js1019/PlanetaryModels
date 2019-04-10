% build surface mesh
clear all; clc; tic;
addpath('../../packages/distmesh/');
global cmi0 Re SR
tic 
load data0/mars_crust.mat
SR = 4; % 720/180

Re = 3396.0;
cmi0 = topo-depth;
R = Re; 
[p,t]=distmeshsurface(@cmidist,@huniform,0.03*R,R*1.2*[-1,-1,-1;1,1,1]);

if 0
fsvtk = 'vtk/Marscmi40k_unit';
%save('workdata/Marscmi_377k.mat','p','t');
Fs = cmidist0(p); 
p = p/R;
vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs, 'Depth (km)', fsvtk);
end