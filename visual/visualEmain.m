% visualize the eigenfunctions
clear all;clc;
addpath('../modelbuilder/'); 

fmesh  = '../Demos/PREM/output/PREM3k/';
fout   = 'demos/PREM3k/';
fbase  = 'prem_3L_3k.1';
fdtail = '0.1000000_1.000000';

JOB = 2; pOrder = 2; nproc = 48; nth = 4; 
Radial = 6.371E3;

fmeshorg = [fmesh,fbase];
fdat =  [fout,fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_',int2str(nth),'.dat'];

fvtk = [fout,'vtk/',fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_',int2str(nth),'.vtk'];

fvlist = [fout,fbase,'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_vlist.dat'];

fvstat = [fout,fbase,'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_vstat.dat'];

[pxyz,tet,~,~,~] = read_mesh3d(fmeshorg);

nvtx = size(pxyz,1);

fid = fopen(fvlist,'r'); 
vlist = fread(fid,'int'); 
fclose(fid);


fid = fopen(fvstat,'r'); 
vstat = fread(fid,'int'); 
fclose(fid);


fid = fopen(fdat,'r');
vsol = fread(fid,'float64');
fclose(fid);

eigm0 = reshape(vsol,3,length(vsol)/3);
ll = sum(vstat==2 | vstat==5); lsiz = length(vstat);
chk = true(ll+lsiz,1); k = 0; 
for i = 1:lsiz
    if (vstat(i)==2 || vstat(i) == 5) 
        chk(k+2) = false;
        k = k + 2;
    else
        k = k + 1;
    end
end
eigm00 = eigm0(:,chk);


eigm1(:,vlist) = eigm00;
eigm  = eigm1(:,1:nvtx);


rad = sqrt(sum(pxyz'.*pxyz'));
Rcom = sum(eigm.*pxyz')./rad;


filename = fvtk;
data_title = 'eigenmodes';
% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct(1).type = 'vector';
data_struct(1).name = 'modes';
data_struct(1).data = eigm(:);

data_struct(2).type = 'scalar';
data_struct(2).name = 'Radial';
data_struct(2).data = Rcom(:);

flipped = false;

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pxyz/Radial,...
    tet,data_struct,flipped);
