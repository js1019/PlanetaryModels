% visualize the eigenfunctions
clear all;clc;
addpath('../modelbuilder/'); 

fmesh  = '/jia/PNM/CONST/trueG/CONST10k/';
fout   = '/jia/PNM/CONST/output/trueG/CONST10k/datan16/';
fbase  = 'CONST_1L_10k.1';
fdtail = '0.0000000E+00_1.000000';

JOB = 1; pOrder = 2; nproc = 16; nth = 12; 
Radial = 6.371E3;

fmeshorg = [fmesh,fbase];
fdat =  [fout,fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_',int2str(nth),'.dat'];

fvtk = [fout,'vtk/',fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_',int2str(nth),'.vtk'];

fvlist = [fout,fbase,'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_vlist.dat'];

[pxyz,tet,~,~,~] = read_mesh3d(fmeshorg);

nvtx = size(pxyz,1);

fid = fopen(fvlist,'r'); 
vlist = fread(fid,'int'); 
fclose(fid);

fid = fopen(fdat,'r');
vsol = fread(fid,'float64');
fclose(fid);
eigm0 = reshape(vsol,3,length(vsol)/3);
eigm1(:,vlist) = eigm0;
eigm  = eigm1(:,1:nvtx);

filename = fvtk;
data_title = 'eigenmodes';
% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct(1).type = 'vector';
data_struct(1).name = 'modes';
data_struct(1).data = eigm(:);

flipped = false;

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pxyz/Radial,...
    tet,data_struct,flipped);
