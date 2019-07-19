% build surface mesh
clear all; clc; tic;
addpath('../../packages/distmesh/');
global cmi0 Re SR
tic 
load data0/mars_crust.mat
SR = 4; % 720/180

Re = 3389.5;
cmi0 = topo-depth;

srf = 0.015; 

R = Re; 
[p,t]=distmeshsurface(@cmidist,@huniform,srf*R,R*1.2*[-1,-1,-1;1,1,1]);
save(['workdata/Mcmi_',num2str(size(t,1)),'.mat'],'p','t');



if 0
fsvtk = 'vtk/Marscmi40k_unit';
%save('workdata/Marscmi_377k.mat','p','t');
Fs = cmidist0(p); 
p = p/R;
vtktrisurf(t, p(:,1), p(:,2), p(:,3), Fs, 'Depth (km)', fsvtk);
end