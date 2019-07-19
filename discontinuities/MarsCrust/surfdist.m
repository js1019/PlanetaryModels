function d = surfdist(p)
global topo Re SR
% transfer to spherical 
[azimuth,elevation,r] = cart2sph(p(:,1),p(:,2),p(:,3));

azimuth = azimuth+pi;
elevation = elevation+pi/2; 

fact = SR; 
az = azimuth/pi*180*fact + 1; 
elv = elevation/pi*180*fact + 1; 

nx = 180*fact*2+1; ny = 180*fact+1;
[X,Y] = meshgrid(1:nx,1:ny);


dint = interp2(X,Y,topo,az,elv);

%max(dint)
%min(dint)

R = Re;
d = dint + R -r;
toc
return