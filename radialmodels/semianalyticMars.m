G = 6.67408*10^-5; % gravitational constant

dsy = MI(end:-1:1,2); 
RI  = MI(end:-1:1,1); 

% extended samplings?
ll = length(MI(:,1)); 
for i = 2:ll
    if (abs(MI(i-1,1)-MI(i,1))<1.E-6)
       %i
       RI(ll-i+2) = RI(ll-i+2) - 1.E-6; 
    end
end

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

phi = phi - mas/RI(1);
%plot(RI,G*phi,'-*')
%plot(RI,G*gref,'-+');
clear i j l ll