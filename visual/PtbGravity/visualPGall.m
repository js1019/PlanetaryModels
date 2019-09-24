% visualize the results
filename = fvtk;
data_title = 'Perturbed Gravity';
% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct(1).type = 'scalar';
data_struct(1).name = 'PPot';
data_struct(1).data = real(U.pottarg(:))*G;

data_struct(2).type = 'vector';
data_struct(2).name = 'PField';
data_struct(2).data = real(U.fldtarg(:))*G;
flipped = false;

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pxyz/Radial,...
    EToV,data_struct,flipped);

clear data_struct

% visual eigenvectors
rad = sqrt(sum(pxyz'.*pxyz'));
Rcom = sum(eigm.*pxyz')./rad;


filename = fvtkeig;
data_title = 'solution';
% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct(1).type = 'vector';
data_struct(1).name = 'modes';
data_struct(1).data = eigm(:);

data_struct(2).type = 'scalar';
data_struct(2).name = 'Radial';
data_struct(2).data = Rcom(:);

flipped = false;

stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pxyz/Radial,...
    EToV,data_struct,flipped);
