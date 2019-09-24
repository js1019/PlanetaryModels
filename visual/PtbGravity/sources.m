% compute values on each element

% read density 
fname1 = [fmesh,fbase];
accry = 'float64';
fmid    = ['_pod_',int2str(N),'_'];
frho    = [fname1,'_rho',fmid,'true.dat'];
fid=fopen(frho,'r');
%Np0 = (pOrder+1)*(pOrder+2)*(pOrder+3)/6; 
rho0=fread(fid,K*Np,accry);
fclose(fid);
rho1 = reshape(rho0,Np,K); 
rho = sum(rho1)/Np; % elementwise 

% values of each elements: first 3: location, last 1 values
v0elm = zeros(K,4); 

% values of each surfaces
% 1:3 locations; 4:5 elem ids; 6 values; 
v0srf = zeros(K*4,6);


for k1 = 1:K
   % Build local operators  
   Derv(:,:,1) = rx(1,k1)*Dr + sx(1,k1)*Ds + tx(1,k1)*Dt;   
   Derv(:,:,2) = ry(1,k1)*Dr + sy(1,k1)*Ds + ty(1,k1)*Dt;
   Derv(:,:,3) = rz(1,k1)*Dr + sz(1,k1)*Ds + tz(1,k1)*Dt; 
   
   nk = EToV(k1,:);
   for j = 1:3
      v0elm(k1,4) = v0elm(k1,4) + ...
          sum(MassMatrix*Derv(:,:,j)*(rho(k1)*eigm(j,nk)'))*J(1,k1);   
   end
   v0elm(k1,1) = sum(x(:,k1))/4;
   v0elm(k1,2) = sum(y(:,k1))/4;
   v0elm(k1,3) = sum(z(:,k1))/4;
   
   for i = 1:4
      ntsf = EToV(k1,Fmask(:,i));
      vsrf = nx(i*3,k1)*eigm(1,ntsf) + ny(i*3,k1)*eigm(2,ntsf) +...
             nz(i*3,k1)*eigm(3,ntsf);
      % values
      v0srf((k1-1)*4+i,6) = v0srf((k1-1)*4+i,6) - ...
              sum(Mass2D*vsrf')*sum(rho1(Fmask(:,i),k1))/3.0*sJ(i*3,k1); 
      % neigh
      v0srf((k1-1)*4+i,4) = min(k1,neigh(k1,i));
      v0srf((k1-1)*4+i,5) = max(k1,neigh(k1,i)); 
      % location 
      %v0srf((k1-1)*4+i,1) = sum(x(Fmask(:,i),k1))/3;
      %v0srf((k1-1)*4+i,2) = sum(y(Fmask(:,i),k1))/3;
      %v0srf((k1-1)*4+i,3) = sum(z(Fmask(:,i),k1))/3;
      v0srf((k1-1)*4+i,1) = sum(pxyz(ntsf,1))/3;
      v0srf((k1-1)*4+i,2) = sum(pxyz(ntsf,2))/3;
      v0srf((k1-1)*4+i,3) = sum(pxyz(ntsf,3))/3;
   end
   
end

[~,ids] = sort(v0srf(:,4)); 
v1srf = v0srf; 
for i = 1:6
    v1srf(:,i) = v0srf(ids,i);
end

Eid0 = zeros(K+1,1); j = 1;  
for i = 2:length(v1srf(:,4))
    if (v1srf(i-1,4)~=v1srf(i,4))
       j = j + 1;
       Eid0(j) = i-1; 
    else
        
    end
end

for i = 1:K
   if (Eid0(i+1)-Eid0(i)==0) 
   else
       tmpl = Eid0(i)+1:Eid0(i+1);
       %tmpl = find(v1srf(:,4)==i);
       [~,id0] = sort(v1srf(tmpl,5));
       for j = 1:6
           v1srf(tmpl,j) = v1srf(tmpl(id0),j);
       end
   end
end

lsrf0 = length(find(v1srf(:,4)==-1)); 

% reorganize
vtsrf = v1srf(1:lsrf0,:); 
v2srf = v1srf(lsrf0+1:end,:);
v3srf = (v2srf(1:2:end-1,:)+v2srf(2:2:end,:))/2.0;
v4srf = [vtsrf;v3srf];

id5 = find(abs(v4srf(:,6))>1e-10); 
v5srf = v4srf(id5,:);
lsrf  = size(v5srf,1); 
%v5srf(:,6) = 0.0;

clear vtsrf v0srf v1srf v2srf v3srf v4srf
