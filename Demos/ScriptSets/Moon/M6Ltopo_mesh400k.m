% build up a Moon model
clear all; clc;

%-------------------------------------------------------------
fmesh = '/local/js116/NMmodels/GitTests/Moon/Mtopo400k/Mtopo_6L_400k';
pOrder  = 1;
% play with a for different number of elements
a = 3.3e5; % 42k 400k
%-------------------------------------------------------------

tetgen = '../../../packages/tetgen1.5.0/tetgen';
tic
% load radial information
load ../../../radialmodels/moon/Moon6L.mat; MI(:,2:4) = MI(:,2:4)/1e3;

fhed = [fmesh,'.1_mesh.header'];
fele = [fmesh,'.1_ele.dat'];
fngh = [fmesh,'.1_neigh.dat'];
fnde = [fmesh,'.1_node.dat'];


R1 = RD(2,1); R2 = RD(3,1); R3 = RD(4,1);
R4 = RD(5,1); R5 = RD(6,1); R6 = RD(7,1);

load ../../../discontinuities/data/Sph6k.mat
p1 = R1*p;
np1 = size(p1,1); t1 = t; nt1 = size(t1,1);

load ../../../discontinuities/data/Sph6k.mat
p2 = p*R2; 
np2 = size(p2,1); t2 = t + np1; nt2 = size(t2,1);

load ../../../discontinuities/data/Sph15k.mat
p3 = p*R3; 
np3 = size(p3,1); t3 = t + np1 + np2;  nt3 = size(t3,1);

load ../../../discontinuities/data/Sph15k.mat
p4 = p*R4; 
np4 = size(p4,1); t4 = t + np1 + np2 +np3; nt4 = size(t4,1);

load ../../../discontinuities/mooncrust/workdata/Mcmi1_14k.mat
p5 = p; 
np5 = size(p5,1); t5 = t + np1 + np2 +np3 + np4; 
nt5 = size(t5,1);

load ../../../discontinuities/mooncrust/workdata/Msurf_42k.mat
p6 = p;
np6 = size(p6,1); t6 = t + np1 + np2 +np3 + np4 + np5; 
nt6 = size(t6,1);

istart = 1; iend = np1; 
pin(istart:iend,:) = p1;
istart = iend + 1; iend = iend + np2; 
pin(istart:iend,:) = p2;
istart = iend + 1; iend = iend + np3; 
pin(istart:iend,:) = p3;
istart = iend + 1; iend = iend + np4; 
pin(istart:iend,:) = p4;
istart = iend + 1; iend = iend + np5; 
pin(istart:iend,:) = p5;
istart = iend + 1; iend = iend + np6; 
pin(istart:iend,:) = p6;

istart = 1; iend = nt1; 
tin(istart:iend,:) = t1;
istart = iend + 1; iend = iend + nt2; 
tin(istart:iend,:) = t2;
istart = iend + 1; iend = iend + nt3; 
tin(istart:iend,:) = t3;
istart = iend + 1; iend = iend + nt4; 
tin(istart:iend,:) = t4;
istart = iend + 1; iend = iend + nt5; 
tin(istart:iend,:) = t5;
istart = iend + 1; iend = iend + nt6; 
tin(istart:iend,:) = t6;

% generate internal surfaces 
trisurf2poly(fmesh,pin,tin);
toc


unix([tetgen,' -pq1.5nYVFAa',num2str(a,'%f'),' ',fmesh,'.poly']);

toc
vtk_write_general([fmesh,'_face.vtk'],'test',pin,tin);
toc

[pout,tout,~,at,neigh] = read_mesh3d([fmesh,'.1']);


dh =[size(tout,1) size(pout,1)];
fid = fopen(fhed,'w');
fprintf(fid,'%d %d',dh);

fid=fopen(fele,'w');
fwrite(fid,tout','int');
fid=fopen(fngh,'w');
fwrite(fid,neigh','int');
fid=fopen(fnde,'w');
fwrite(fid,pout','float64');

size(tout)
toc

run MoonTopo_models 
run Gravity
