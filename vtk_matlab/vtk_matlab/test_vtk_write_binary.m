%
% test_vtk_write_binary
clc;
clear all;

% generate data
% n = 20;
% [tet_connectivity, vtx_coord] = regular_tetrahedral_mesh(n,n,n);
MM = load('Tet_Mesh.mat','tet_connectivity','vtx_coord');
tet_connectivity = MM.tet_connectivity;
vtx_coord = MM.vtx_coord;
clear MM;

x = vtx_coord(:,1);
y = vtx_coord(:,2);
z = vtx_coord(:,3);

% make up a "pressure"
Pr = sin(2*pi*x) .* cos(2*pi*z);

% make up a "velocity" field
U = x.^2 + y.^3;
V = cos(2*pi*(y + z));
W = sin(2*pi*(x - y));

% setup filename
filename = 'Vel_Pressure_Data.vtk';
data_title = 'example velocity and pressure data';

% organize data
% .type field must be 'scalar' or 'vector'
% number of vector components must equal dimension of the mesh
data_struct.type = 'scalar';
data_struct.name = 'pressure';
data_struct.data = Pr;

data_struct(2).type = 'vector';
data_struct(2).name = 'velocity';
data_struct(2).data = [U,V,W];

% if *all data* is oriented column-wise, i.e. lots of rows, few columns,
% then set this to *false*
flipped = false;
% otherwise, if you want to transpose the data, then set this to *true*

% write the file
tic
stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,vtx_coord,tet_connectivity,data_struct,flipped);
toc

% END %