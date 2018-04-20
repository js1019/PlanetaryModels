% construct the gravity fields
clear all; clc; tic
addpath('../packages/fmmlib3d-1.2/matlab/');
%fmesh = '/local/js116/NM_models/Earth/models/PREM2M/prem_3L_2M';
%fmesh = '/pylon2/ac4s8pp/js116/NMmodels/PREM32M/prem_3L_32M';
fmesh = '/jia/PNM/CONST/trueG/CONST3k/CONST_1L_3k';


scaling = 6.371*10^3;

[pout,tout,~,at] = read_mesh3d([fmesh,'.1']);

pOrder  = 1;
fmid    = ['_pod_',int2str(pOrder),'_'];
% true model filenames
fname = [fmesh,'.1']; 
ftail   = 'true.dat';
frho    = [fname,'_rho',fmid,ftail];
fgfld = [fname,fmid,'potential_accceleration_',ftail];
fvtk = [fname,fmid,'gravity.vtk'];

accry = 'float64';
G = 6.6723*10^-5; % gravitational constant

Ne = size(tout,1); pNp = (pOrder+1)*(pOrder+2)*(pOrder+3)/6;

% read density model
fid=fopen(frho);
rho = fread(fid,Ne*pNp,accry);
fclose(fid);

rho0 = reshape(rho,pNp,Ne);

% compute the detJ
N = pOrder; 
[x,y,z] = Nodes3D(N); 
[r,s,t] = xyztorst(x,y,z);
V = Vandermonde3D(N,r,s,t);
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);
Mass = inv(V)'*inv(V);
va = tout(:,1)'; vb = tout(:,2)'; vc = tout(:,3)'; vd = tout(:,4)';
x = 0.5*(-(1+r+s+t)*pout(va,1)'+(1+r)*pout(vb,1)'+(1+s)*pout(vc,1)'+(1+t)*pout(vd,1)');
y = 0.5*(-(1+r+s+t)*pout(va,2)'+(1+r)*pout(vb,2)'+(1+s)*pout(vc,2)'+(1+t)*pout(vd,2)');
z = 0.5*(-(1+r+s+t)*pout(va,3)'+(1+r)*pout(vb,3)'+(1+s)*pout(vc,3)'+(1+t)*pout(vd,3)');
[~,~,~,~,~,~,~,~,~,J] = GeometricFactors3D(x,y,z,Dr,Ds,Dt);

% prepare the sources
iprec= 5; %5; 
nsource = Ne; 

source(1,:) = ones(1,pNp)*x/pNp;
source(2,:) = ones(1,pNp)*y/pNp;
source(3,:) = ones(1,pNp)*z/pNp;

ifcharge = 1; 
charge = J(1,:).*sum(Mass*rho0);

ifdipole = 0; dipstr = zeros(1,Ne); dipvec = rand(3,Ne);
ifpot = 0; iffld = 0;
%ntarget = Ne*pNp; 
%target(1,:) = reshape(x,ntarget,1);
%target(2,:) = reshape(y,ntarget,1);
%target(3,:) = reshape(z,ntarget,1);

ntarget = size(pout,1); 
target = pout'; 


ifpottarg = 1; iffldtarg = 1;
toc 
[U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,...
    ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
toc
rnrm = sqrt(sum(target.*target));
gnrm = sqrt(sum(real(U.fldtarg).*real(U.fldtarg)))*G;

max(gnrm(:))
min(real(U.pottarg(:))*G)

U.ier

if (max(gnrm(:))>1.E-5) 
gfld = -real(U.fldtarg(:,tout'))*G;
gpot = -real(U.pottarg(:))*G;

if 0
% save the data
fid=fopen(fgfld,'w');
fwrite(fid,gfld',accry);
fclose(fid);
end 


% visual
filename = fvtk;
data_title = 'Gravity';
% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct(1).type = 'scalar';
data_struct(1).name = 'potential';
data_struct(1).data = -real(U.pottarg(:))*G;

data_struct(2).type = 'vector';
data_struct(2).name = 'field';
data_struct(2).data = -real(U.fldtarg(:))*G;
flipped = false;

% otherwise, if you want to transpose the data, then set this to *true*
%tnew = reshape(1:size(tout,1)*size(tout,2),size(tout,2),size(tout,1));
%pnew = pout(tout',:)/scaling;  
% write the file

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pout/scaling,...
    tout,data_struct,flipped);
toc
end

