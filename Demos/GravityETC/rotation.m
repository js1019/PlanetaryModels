% construct the rotation and centrifugal force

% set a few parameters
period = 24.6229; %24*3600+37*60+22; hours
aglar = 2*pi/(period*3600); 
Omega = [0 0 1]*aglar;

factor = 1e3; % km 

Ome0 = ones(size(pnew0,1),1)*Omega;
rot0 = cross(Ome0,pnew0); rot1 = rot0'*factor; 
fce0 = - cross(Ome0,rot0); fce1 = fce0'*factor; 

% visualize the force
if 1
filename = [fname,fmid,'Rotation.vtk'];
data_title = 'Rotation';
data_str(1).type = 'vector';
data_str(1).name = 'Rotation';
data_str(1).data = rot1(:);

data_str(2).type = 'vector';
data_str(2).name = 'Force';
data_str(2).data = fce1(:);
flipped = false;

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pnew0/scaling,...
    tet,data_str,flipped);

end
