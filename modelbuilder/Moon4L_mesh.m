% build up a Moon model
clear all; clc;
 
%fmesh = '/local/js116/NM_models/CONST/models/CONST256M/CONST_1L_256M';
%fmesh  = '/jia/PNM/CONST/trueG/CONST3k/CONST_1L_3k';
fmesh = '../../models/Moon4L_5k';
%tetgen = '/home/js116/Documents/tetgen1.5.0/tetgen';
tetgen = '../packages/tetgen1.5.0/tetgen';

tic
% load radial information
load ../deal_prem/Moon4L.mat;

fhed = [fmesh,'.1_mesh.header'];
fele = [fmesh,'.1_ele.dat'];
fngh = [fmesh,'.1_neigh.dat'];
fnde = [fmesh,'.1_node.dat'];

pOrder  = 1;

R1 = RD(2,1); R2 = RD(3,1); R3 = RD(4,1);
R4 = RD(5,1);

load ../unitspheres/data/Sph392.mat
p1 = R1*p;
np1 = size(p1,1); t1 = t; nt1 = size(t1,1);

load ../unitspheres/data/Sph392.mat
p2 = p*R2; 
np2 = size(p2,1); t2 = t + np1; nt2 = size(t2,1);

load ../unitspheres/data/Sph392.mat
p3 = p*R3; 
np3 = size(p3,1); t3 = t + np1 + np2;  nt3 = size(t3,1);

load ../unitspheres/data/Sph956.mat
p4 = p*R4; 
np4 = size(p4,1); t4 = t + np1 + np2 +np3; nt4 = size(t4,1);

istart = 1; iend = np1; 
pin(istart:iend,:) = p1;
istart = iend + 1; iend = iend + np2; 
pin(istart:iend,:) = p2;
istart = iend + 1; iend = iend + np3; 
pin(istart:iend,:) = p3;
istart = iend + 1; iend = iend + np4; 
pin(istart:iend,:) = p4;

istart = 1; iend = nt1; 
tin(istart:iend,:) = t1;
istart = iend + 1; iend = iend + nt2; 
tin(istart:iend,:) = t2;
istart = iend + 1; iend = iend + nt3; 
tin(istart:iend,:) = t3;
istart = iend + 1; iend = iend + nt4; 
tin(istart:iend,:) = t4;

% generate internal surfaces 
trisurf2poly(fmesh,pin,tin);
toc
% generate the mesh
a = 4e10; % 956 5k

unix([tetgen,' -pq1.5nYVFAa',num2str(a,'%f'),' ',fmesh,'.poly']);
toc
[pout,tout,~,at] = read_mesh3d([fmesh,'.1']);
toc
vtk_write_general([fmesh,'_face.vtk'],'test',pin,tin);
toc
