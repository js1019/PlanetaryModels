% build up an earth model
clear all; clc;
 
%fmesh = '/local/js116/NM_models/CONST/models/CONST256M/CONST_1L_256M';
fmesh  = '/jia/PNM/CONST/CONST3k/CONST_1L_3k_p1';
%tetgen = '/home/js116/Documents/tetgen1.5.0/tetgen';
tetgen = '../tetgen1.5.0/tetgen';
tic
% load radial information
load ../deal_prem/prem3L_noocean.mat

% radius 
R1 = RD(1,1); 

% load unit spheres
load ../unitspheres/data/Sph392.mat
p1 = R1*p;
np1 = size(p1,1); t1 = t; nt1 = size(t1,1);

istart = 1; iend = np1; 
pin(istart:iend,:) = p1;

istart = 1; iend = nt1; 
tin(istart:iend,:) = t1;

% generate internal surfaces 
trisurf2poly(fmesh,pin,tin);
toc 
% generate the mesh
a = 5e8; %392 3k
%a = 3e8; %392 5k 
%a = 1.5e8; % 956 10k
%a = 2e6; % 42k 1M
%a = 1e6; % 94k 2M
%a = 5e5; % 167k 4M 
%a = 2.5e5; % 167k 8M 
%a = 2.6e5; % 377k 8M
%a = 1.3e5; % 589k 16M
%a = 6.5e4; % 1047k 32M
%a = 3.2e4; % 1507k 64M
%a = 1.6e4; % 2353k 128M
%a = 8.0e3; % 3077k 256M
unix([tetgen,' -pq1.5nYVFAa',num2str(a,'%f'),' ',fmesh,'.poly']);
toc
[pout,tout,~,~,neigh] = read_mesh3d([fmesh,'.1']);

%vtk_write_general([fmesh,'_face.vtk'],'test',pin,tin);
toc
