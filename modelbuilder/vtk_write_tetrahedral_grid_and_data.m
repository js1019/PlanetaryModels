function stat = vtk_write_tetrahedral_grid_and_data(filename,vtk_title,grid_X,grid_TET,data_struct,flipped)
%vtk_write_tetrahedral_grid_and_data
%
%   This writes a .vtk file (in legacy-binary format) that stores a given
%   tetrahedral mesh and associated data consisting of scalars and vectors.
%
%   Note: for now, we assume the tetrahedral grid is linear (i.e. each tet
%   is defined by 4 vertices only).
%
%   stat = vtk_write_tetrahedral_grid_and_data(filename,vtk_title,...
%                                grid_X,grid_TET,data_struct,flipped);
%
%   filename  = string containing the filename to write to.
%   vtk_title = string description of the data (written in the file).
%   grid_X    = MxD matrix of vertex coordinates, where D is the geometric
%               dimension of the grid (must be 3!); M is the number of
%               vertices.  Note: if flipped == true, then this data is
%               assumed to be transposed.
%   grid_TET  = RxK matrix of connectivity info, where K is the number
%               vertices (nodes) needed to define each tetrahedron; R is
%               the number of tetrahedra.  Note: that this is assumed to
%               use 1-based indexing.  This routine internally subtracts
%               "1" so that the vtk data uses 0-based indexing. Note: if
%               flipped == true, then this data is assumed to be transposed.
%   data_struct = 1-D array of structs with the following sub-fields:
%                 .type = string == 'scalar' or 'vector'.
%                 .name = name string .
%                 .data = MxS numerical data, where M is the number of
%                         vertices in grid_X, and S == 1 for 'scalar' or
%                         S == 3 for 'vector'.
%                 Note: if flipped == true, then this data is assumed to be
%                 transposed.
%   flipped   = true/false indicating if the data has already been
%               transposed; false means not transposed.

% Copyright (c) 07-01-2016,  Shawn W. Walker

% define write commands
write_ascii = @(fid,str) [fprintf(fid,str); fprintf(fid,'\n');];
write_data  = @(fid,dat,prec) [fwrite(fid,dat,prec); fprintf(fid,'\n');];

% check the filename
[P1, N1, E1] = fileparts(filename);
if isempty(E1)
    E1 = '.vtk'; % add the extension
elseif ~strcmp(E1,'.vtk')
    disp('warning: file extension is *not* ".vtk"!');
end
filename = fullfile(P1, [N1, E1]);

% open the file in binary mode
fopen_opts = {'wb','ieee-be'};
fid = fopen(filename,fopen_opts{:});
if fid == -1
    error('Unable to write file %s: permission denied.',filename);
end

% write the initial header
write_ascii(fid,'# vtk DataFile Version 2.0');
vtk_title_reduce = vtk_title(1:min(length(vtk_title),256));
write_ascii(fid,vtk_title_reduce);
write_ascii(fid,'BINARY\n'); % give extra line-feed

% write the vertex coordinates of the mesh
write_ascii(fid,['DATASET ', 'UNSTRUCTURED_GRID']);
if ~flipped
    grid_X = grid_X'; % need to transpose for fwrite
end
GD = size(grid_X,1);
if (GD~=3)
    stat = fclose(fid);
    error('Grid vertex coordinate data is not 3-D!');
end
Num_Vtx = size(grid_X,2);
write_ascii(fid, ['POINTS ', num2str(Num_Vtx), ' float']);
write_data(fid, grid_X, 'float32');

% write the tetrahedral grid connectivity
write_ascii(fid, ''); % write a blank line
if ~flipped
    Tet_Order = size(grid_TET,2);
    Num_Tet   = size(grid_TET,1);
else
    Tet_Order = size(grid_TET,1);
    Num_Tet   = size(grid_TET,2);
end
if (Tet_Order~=4)
    stat = fclose(fid);
    error('Grid tetrahedra connectivity does not have 4 nodes per tet!');
end
Cell_Size = Num_Tet * ( Tet_Order + 1 ); % vtk needs this
write_ascii(fid,['CELLS  ', num2str(Num_Tet), '  ', num2str(Cell_Size)]);

% 1-based to 0-based indexing
min_index = min(grid_TET(:));
if min_index==0
    disp('Tetrahedral connectivity is already 0-based!');
elseif min_index > 0
    grid_TET = grid_TET - 1;
else
    fclose(fid);
    error('Tetrahedral connectivity has negative indices!');
end

% modify and append extra data
if ~flipped
    DATA = uint32([(Tet_Order + 0*grid_TET(:,1))';
                    grid_TET']);
else
    DATA = uint32([(Tet_Order + 0*grid_TET(1,:));
                    grid_TET]);
end
write_data(fid,DATA, 'uint32');

% must write the cell types
%   VTK has a cell type 24 for quadratic tetrahedrons.  For standard linear
%   tetrahedrons, the cell type is 10.
write_ascii(fid,''); % blank line
write_ascii(fid,['CELL_TYPES ', num2str(Num_Tet)]);
if ( Tet_Order == 4 )
    TET_LABEL = 10; % cell type label == 10
elseif ( Tet_Order == 10 )
    TET_LABEL = 24; % cell type label == 24
else
    error('Invalid tetrahedral order!');
end
DATA = uint32(TET_LABEL*ones(1,Num_Tet));
write_data(fid,DATA, 'uint32');

% write the POINT_DATA
if ~isempty(data_struct)
    write_ascii(fid,''); % blank line
    write_ascii(fid,['POINT_DATA ', num2str(Num_Vtx)]);
    
    Num_DS = length(data_struct);
    for ii = 1:Num_DS
        if strcmpi(data_struct(ii).type,'scalar')
            write_ascii(fid,''); % blank line
            write_ascii(fid,['SCALARS ', data_struct(ii).name, ' float']);
            write_ascii(fid,'LOOKUP_TABLE default');
            write_data(fid,data_struct(ii).data(:)','float32'); % scalar data easy to write
        elseif strcmpi(data_struct(ii).type,'vector')
            write_ascii(fid,''); % blank line
            write_ascii(fid,['VECTORS ', data_struct(ii).name, ' float']);
            if ~flipped
                D1 = data_struct(ii).data';
            else
                D1 = data_struct(ii).data;
            end
            write_data(fid,D1,'float32');
        else
            stat = fclose(fid);
            error('Invalid data type!');
        end
    end
end

% Note: CELL_DATA is not implemented!

% end of file!
stat = fclose(fid);

end