% mesh file names
fmeshorg = [fmesh,fbase];
fdat =  [fout,fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_',int2str(nth),'.dat'];

fvtk = [fout,'vtk/',fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_PG_',int2str(nth),'.vtk'];

fvtkeig = [fout,'vtk/',fbase,'_JOB',int2str(JOB),'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_',fdtail,'_eigv_',int2str(nth),'.vtk'];


fvlist = [fout,fbase,'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_vlist.dat'];

fvstat = [fout,fbase,'_pod',int2str(pOrder),...
    '_np',int2str(nproc),'_vstat.dat'];

% obtain meshes info
[pxyz,EToV,~,~,neigh] = read_mesh3d(fmeshorg);

% number of vertices
nvtx = size(pxyz,1); 
VX = pxyz(:,1)'; VY = pxyz(:,2)'; VZ = pxyz(:,3)'; 
% number of elements
K = size(EToV,1);


% read info 
fid = fopen(fvlist,'r'); 
vlist = fread(fid,'int'); 
fclose(fid);

fid = fopen(fvstat,'r'); 
vstat = fread(fid,'int'); 
fclose(fid);

% obtain the values of the eigenvector on the vertices
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
% end of obtaining the eigenvector 

% rescale eigm to meters
eigm = eigm*1.e3;

% clear 
clear eigm00 eigm0 eigm1 vsol 
