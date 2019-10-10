function d = cmidist0(p)

global cmi0 Re SR
% transfer to spherical 
[azimuth,elevation,r] = cart2sph(p(:,1),p(:,2),p(:,3));

azimuth = azimuth+pi;%+0.5/180*pi;
elevation = elevation+pi/2;%+0.5/180*pi; 

fact = SR; 
az = azimuth/pi*180*fact+1; 
elv = elevation/pi*180*fact+1; 

%az(find(az<0.5)) = 360 + az(find(az<0.5)); 
%elv(find(elv<0.5)) = 180 + elv(find(elv<0.5)); 


nx = 180*fact*2+1; ny = 180*fact+1;
[X,Y] = meshgrid(1:nx,1:ny);

d = interp2(X,Y,cmi0,az,elv,'spline');

%max(dint)

%R = Re;
%d = dint + R -r;


return
