% build ellpsoid with Radau approximation 
clear all; clc;
load('../mars/marsDWAK_gravity.mat');
global RI epsl;

% average radius
Re = 3389.5; % km
Tp = 24.6229; % hours
g0 = 3.71; % m/s^2 surface gravity
% check
% plot(RI,G*gref,'+'); 

dsy = MI(end:-1:1,2);
lz = length(RI); 
RI0 = MI(end:-1:1,1); 
eta = zeros(lz,1); 

% build eta via Dahlen & Tromp (14.19)
d0 = 0; d1 = 0;
for i = 1:lz-1
   j = lz-i+1; 
   d0 = d0 + dsy(j-1)/3*(RI0(j-1)^3-RI0(j)^3); 
   d1 = d1 + dsy(j-1)/5*(RI0(j-1)^5-RI0(j)^5);
   eta(j-1) = 25/4*(1-d1/d0/RI0(j-1)^2)^2-1; 
end

%plot(RI0,eta);
%epa = 5/2*Re*1e3/g0/(eta(1)+2)*(2*pi/(Tp*3600))^2;
% change !! epa to the observed value
epa = 589e-5;

epsl = zeros(lz,1); epsl(:) = epa; 

d0 = 1; 
for i = 2:lz
    d0 = d0/(RI0(i-1)/RI0(i))^(eta(i));
    epsl(i) = epa*d0; 
end
%plot(RI0,epsl)

% compute the ellpicity of three layers: surface, CMB, ICB
% RE = zeros(nlayer,1); 
% RE(1) = epsl(1); 
% need to check!!
RE = epsl(lz-RD(2:3,2)+1); 


%% build map (r, cos \theta) with a 
rmax = Re*(1+epa/3); nr = 200; ny = 50;
xr = linspace(0,rmax,nr); co2 = linspace(0,1,ny);

axr = zeros(nr,ny); a0 = 0; a1 = Re;
for i = 2:nr
   for j = 1:ny
       b0 = Re*(1-epa*(co2(j)-1/3));
       if (xr(i) <= b0)
          axr(i,j) = bisection(@fundistant,a0,a1,xr(i),co2(j)); 
       else
          axr(i,j) = Re;
       end
   end
end

%imagesc(axr)
clear d0 d1 i j lz RI0 a0 a1 