function d = cmi1dist(p)
global srfmap thk1map
% transfer to spherical 
[azimuth,elevation,r] = cart2sph(p(:,1),p(:,2),p(:,3));

azimuth = azimuth+pi;
elevation = elevation+pi/2; 

az = azimuth/pi*180*4 + 1; 
elv = elevation/pi*180*4 + 1; 

[X,Y] = meshgrid(1:1441,1:721);


dint0 = interp2(X,Y,srfmap,az,elv,'spline');
dint1 = interp2(X,Y,thk1map,az,elv,'spline');

%max(dint)

R = 1737.15;
d = dint0 - dint1 + R -r;
toc
return