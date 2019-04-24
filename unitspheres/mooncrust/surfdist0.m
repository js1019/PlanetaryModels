function dint = surfdist0(p)
global srfmap 
% transfer to spherical 
[azimuth,elevation,~] = cart2sph(p(:,1),p(:,2),p(:,3));

azimuth = azimuth+pi;
elevation = elevation+pi/2; 

az = azimuth/pi*180*4 + 1; 
elv = elevation/pi*180*4 + 1; 

[X,Y] = meshgrid(1:1441,1:721);


dint = interp2(X,Y,srfmap,az,elv);


return