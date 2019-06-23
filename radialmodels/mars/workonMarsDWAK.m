% work on a Mars model
clear all; clc; 
fname = '../models1D/Mars_DWAK_noheader.bm';
M0 = load(fname); 

% construct an isotropic file with land only
nwater = sum((M0(:,3)<1.5)); 
siz1   = size(M0,1); 

MI = zeros(siz1-nwater,4); 
MI(:,1) = M0(nwater+1:end,1); MI(1,1) = M0(1,1); % radius 
MI(:,2) = M0(nwater+1:end,2); % density
MI(:,3) = (M0(nwater+1:end,3)+M0(nwater+1:end,7))/2.0; % Vp
MI(:,4) = (M0(nwater+1:end,4)+M0(nwater+1:end,8))/2.0; % Vs

MI(:,:) = MI(:,:)/1e3; % from meter to km 

% find out intersections
nuniq = length(unique(MI(:,1)));
nitsc = length(MI(:,1)) - nuniq; 

% edit a few jumps
insc = zeros(nitsc+1,2); 
j = 0; insc(1,1) = M0(1,1); 
for i = 2:length(MI(:,1))
    if (abs(MI(i-1,1)-MI(i,1))<1.E-3)
       j = j + 1;
       insc(j+1,1) = MI(i,1);
       insc(j+1,2) = i - 1; 
    end
end

% build a 2L model for 3D interpolation 
if 1 
   % keep 1 11 & 12
   nlayer = 3; 
   keep = [1 2 4];
   RD = zeros(nlayer,2); 
   RD(nlayer+1,1) = MI(end,1); % origin 
   for i = 1:nlayer
      RD(i,1) = insc(keep(i),1); 
   end

   l = 1;
   for i = 1:size(insc,1)
      if (keep(l)==i) 
          RD(l,2) = insc(i,2);
          if (l<nlayer) 
              l = l + 1;
          end
      else
          MI(insc(i,2),1) = MI(insc(i,2),1) - .1;
      end
   end
   RD(nlayer+1,2) = size(MI,1); 
end

run semianalyticMars

clear mas siz1 nitsc nuniq 
