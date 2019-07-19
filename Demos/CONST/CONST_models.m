% create the model files 
fname = [fmesh,'.1']; 
tic


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

pNp = (pOrder+1)*(pOrder+2)*(pOrder+3)/6;
[x,y,z,tet] = construct(fname,pOrder);



vp0  = zeros(pNp,size(x,2)); 
vp0(:) = 10.000;
fid=fopen(fvp,'w');
fwrite(fid,vp0(:),accry);

vp = vp0(tet');
%clear vp

vs0  = zeros(pNp,size(x,2)); 
%vs0(:) = 10.000/sqrt(3.0);
vs0(:) = 5.7735;
fid=fopen(fvs,'w');
fwrite(fid,vs0(:),accry);

vs = vs0(tet');
%clear vs

rho0 = zeros(pNp,size(x,2)); 
rho0(:) = 5.5100;
fid=fopen(frho,'w');
fwrite(fid,rho0(:),accry);
fclose(fid);

rho = rho0(tet');
%clear rho

toc
if 1
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
tnew = reshape(1:size(tet,1)*size(tet,2),size(tet,2),size(tet,1));
psiz = max(tet(:));
pnew0 = zeros(psiz,3);
pnew0(:,1) = x(:); 
pnew0(:,2) = y(:); 
pnew0(:,3) = z(:);  

pnew = pnew0(tet',:);
 
% write the file
stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pnew,...
    tnew',data_struct,flipped);
toc
end

