% build up an earth model
clear all; clc;
addpath('../../modelbuilder/');  

fmesh  = 'output/PREM3k/prem_3L_3k';
%fmesh = '/pylon2/ac4s8pp/js116/NMmodels/PREM512M/prem_3L_512M';
tetgen = '../../packages/tetgen1.5.0/tetgen'; 

% finite element order (choose 1 or 2)
pOrder  = 1;

% set the value to control the degrees of freedom
a = 8e10; % 396 3k
%a = 4e10; % 3k6 20k
%a = 6e7; % 15k 100k
%a = 2e6; % 42k 1M
%a = 1.15e6; % 94k 2M
%a = 5.8e5; % 167k 4M
%a = 3.1e5; % 377k 8M
%a = 1.5e5; % 589k 16M
%a = 7.5e4; % 1047k 32M
%a = 3.6e4; % 1507k 64M
%a = 1.75e4; % 2353k 128M
%a = 8.0e3; % 3077k 256M
%a = 5.0e3; % 3077k 400M
%a = 4.0e3; % 3077k 512M
%a = 3.0e3; % 3077k 690M


tic
% load radial information
load ../../radialmodels/PREM/prem3L_noocean.mat

% radius 
R1 = RD(1,1); R2 = RD(2,1); R3 = RD(3,1);

% load unit spheres
load ../../discontinuities/data/Sph392.mat
p1 = R1*p;
np1 = size(p1,1); t1 = t; nt1 = size(t1,1);

load ../../discontinuities/data/Sph260.mat
p2 = p*R2; 
np2 = size(p2,1); t2 = t + np1; nt2 = size(t2,1);

load ../../discontinuities/data/Sph260.mat
p3 = p*R3; 
np3 = size(p3,1); t3 = t + np1 + np2;  nt3 = size(t3,1);

istart = 1; iend = np1; 
pin(istart:iend,:) = p1;
istart = iend + 1; iend = iend + np2; 
pin(istart:iend,:) = p2;
istart = iend + 1; iend = iend + np3; 
pin(istart:iend,:) = p3;

istart = 1; iend = nt1; 
tin(istart:iend,:) = t1;
istart = iend + 1; iend = iend + nt2; 
tin(istart:iend,:) = t2;
istart = iend + 1; iend = iend + nt3; 
tin(istart:iend,:) = t3;


% generate internal surfaces 
trisurf2poly(fmesh,pin,tin);
toc 

fhed = [fmesh,'.1_mesh.header'];
fele = [fmesh,'.1_ele.dat'];
fngh = [fmesh,'.1_neigh.dat'];
fnde = [fmesh,'.1_node.dat'];

% generate the mesh
unix([tetgen,' -pq1.5nYVFAa',num2str(a,'%f'),' ',fmesh,'.poly']);
toc
[pout,tout,~,at,neigh] = read_mesh3d([fmesh,'.1']);

dh =[size(tout,1) size(pout,1)];
fid = fopen(fhed,'w');
fprintf(fid,'%d %d',dh);
fclose(fid);

fid=fopen(fele,'w');
fwrite(fid,tout','int');
fclose(fid);
fid=fopen(fngh,'w');
fwrite(fid,neigh','int');
fclose(fid);
fid=fopen(fnde,'w');
fwrite(fid,pout','float64');
fclose(fid);

%vtk_write_general([fmesh,'_face.vtk'],'test',pin,tin);
size(tout)
toc


run PREM_models
