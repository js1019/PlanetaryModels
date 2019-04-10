function d = surfdist0(p)
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


d = interp2(X,Y,topo,az,elv);
%d = d/1e3;

toc
return
