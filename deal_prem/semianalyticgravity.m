% calculate gravity
clear all; clc;
load prem3L_noocean.mat
G = 6.6723*10^-5; % gravitational constant

dsy = MI(:,2); 
RI = MI(:,1); 

phi = zeros(size(RI)); gref = zeros(size(RI));

% for gravitational acceleration
l = size(RI,1); mas = 0.0; 
for i = 1:l-1
    j = l - i + 1;
    mas = mas + 4*pi/3*(RI(j-1)^3-RI(j)^3)*dsy(j-1); 
    %phi(j-1) = - mas/RI(j-1);
    gref(j-1) = mas/(RI(j-1)^2);
end

% for potential
phi(1) = - mas/RI(1);
for i = 2:l
   phi(i) = phi(i-1) - (gref(i-1)+gref(i))/2*(RI(i-1)-RI(i)); 
end

%phi = phi - mas/RI(1);
%plot(RI,G*phi,'*')
plot(RI,G*gref,'+');
clear i j l 
