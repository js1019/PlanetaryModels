% create the model files 
fname = [fmesh,'.1']; 
%tic
%pOrder  = 2;

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
vs0  = zeros(pNp,size(x,2)); 
rho0 = zeros(pNp,size(x,2)); 


for i = 1:nlayer
   tid = find(at==i); 
   rvd = sqrt(x(:,tid).^2+y(:,tid).^2+z(:,tid).^2);
   
   ravg = sum(rvd(:))/length(tid)/pNp; 
   for j = 1:nlayer
   if (ravg > RD(j,1) && ravg < RD(j+1,1))
      tmpr   = MI(RD(j,2)+1:RD(j+1,2),1);
      tmpvp  = MI(RD(j,2)+1:RD(j+1,2),3);
      vptmp  = interp1(tmpr,tmpvp,rvd,'pchip');
      vp0(:,tid) = reshape(vptmp,pNp,length(tid)); 
      tmpvs  = MI(RD(j,2)+1:RD(j+1,2),4);
      vstmp  = interp1(tmpr,tmpvs,rvd,'pchip');
      vs0(:,tid) = reshape(vstmp,pNp,length(tid)); 
      tmprho  = MI(RD(j,2)+1:RD(j+1,2),2);
      rhotmp  = interp1(tmpr,tmprho,rvd,'pchip');
      rho0(:,tid) = reshape(rhotmp,pNp,length(tid)); 
   end
   end
end
toc

%clear ttrs vstmp vptmp rhotmp rvd

% write the true model 
fid=fopen(fvp,'w');
fwrite(fid,vp0(:),accry);
fid=fopen(fvs,'w');
fwrite(fid,vs0(:),accry);
fid=fopen(frho,'w');
fwrite(fid,rho0(:),accry);
fclose(fid);


vp  = vp0(tet');
vs  = vs0(tet');
rho = rho0(tet');

if 1
% setup filename
filename = fvtk;
data_title = 'Mars_DWAK';
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
%tnew = reshape(1:size(tout,1)*size(tout,2),size(tout,2),size(tout,1));
%pnew = pout(tout',:);  

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
