% create the model files 
fname = [fmesh,'.1']; 
tic
pOrder  = 1;

fmid    = ['_pod_',int2str(pOrder),'_'];
% true model filenames
ftail   = 'true.dat';
fvp     = [fname,'_vp', fmid,ftail];
fvs     = [fname,'_vs', fmid,ftail];
frho    = [fname,'_rho',fmid,ftail];
%fvpvtk  = [fname,'_vp', fmid,'true.vtk'];
%fvsvtk  = [fname,'_vs', fmid,'true.vtk'];
%frhovtk = [fname,'_rho',fmid,'true.vtk'];
fvtk = [fname,fmid,'model.vtk'];

accry = 'float64';
%pNp = (pOrder+1)*(pOrder+2)*(pOrder+3)/6;
%[x,y,z,tet] = construct(fname,pOrder);

vp  = zeros(size(tout,2),size(tout,1)); 
vp(:) = 10.0;
fid=fopen(fvp,'w');
fwrite(fid,vp(:),accry);

clear vp

vs  = zeros(size(tout,2),size(tout,1)); 
vs(:) = 10.0/sqrt(3.0);
fid=fopen(fvs,'w');
fwrite(fid,vs(:),accry);

clear vs

rho = zeros(size(tout,2),size(tout,1)); 
rho(:) = 5.51;
fid=fopen(frho,'w');
fwrite(fid,rho(:),accry);
fclose(fid);

clear rho

toc
if 0
% setup filename
filename = fvtk;
data_title = 'CONST';
% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct(1).type = 'scalar';
data_struct(1).name = 'Vp';
data_struct(1).data = vp(:);

data_struct(2).type = 'scalar';
data_struct(2).name = 'Vs';
data_struct(2).data = vs(:);

data_struct(3).type = 'scalar';
data_struct(3).name = 'Density';
data_struct(3).data = rho(:);

flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*
tnew = reshape(1:size(tout,1)*size(tout,2),size(tout,2),size(tout,1));
pnew = pout(tout',:);  


% write the file

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pnew,...
    tnew',data_struct,flipped);
toc
end

