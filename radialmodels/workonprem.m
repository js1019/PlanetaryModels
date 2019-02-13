% deal with PREM data
clear all; clc; 

fname = 'PREM_1s.csv';
% Columns from left to right represent 
% radius,depth,density,Vpv,Vph,Vsv,Vsh,eta,Q-mu,Q-kappa
M = csvread(fname);

% construct an isotropic file with land only
nwater = sum((M(:,3)<1.5)); 
siz1   = size(M,1); 

MI = zeros(siz1-nwater,4); 
MI(:,1) = M(nwater+1:end,1); MI(1,1) = M(1,1); % radius 
MI(:,2) = M(nwater+1:end,3); % density
MI(:,3) = (M(nwater+1:end,4)+M(nwater+1:end,5))/2.0; % Vp
MI(:,4) = (M(nwater+1:end,6)+M(nwater+1:end,7))/2.0; % Vs

% find out intersections
nuniq = length(unique(MI(:,1)));
nitsc = length(MI(:,1)) - nuniq; 

insc = zeros(nitsc+1,2); 
j = 0; insc(1,1) = M(1,1); 
for i = 2:length(MI(:,1))
    if (abs(MI(i-1,1)-MI(i,1))<1.E-3)
       j = j + 1;
       insc(j+1,1) = MI(i,1);
       insc(j+1,2) = i - 1; 
    end
end

% build a 3L model for 3D interpolation 
if 1 
   % keep 1 11 & 12
   nlayer = 3; 
   keep = [1 11 12];
   RD = zeros(nlayer+1,2); 
   RD(nlayer+1,1) = 0.0; % origin 
   for i = 1:nlayer
      RD(i,1) = insc(keep(i),1); 
   end
   
   l = 1;
   for i = 1:size(insc,1)
      if (keep(l)==i) 
          RD(l,2) = insc(i,2);
          l = l + 1;
      else
          MI(insc(i,2),1) = MI(insc(i,2),1) + .01;
      end
   end
   RD(nlayer+1,2) = size(MI,1); 
end

clear i j l keep insc fname M nitsc nuniq nwater siz1
