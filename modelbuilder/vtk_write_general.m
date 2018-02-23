function vtk_write_general(fname,title,xyz,t,p,datatype)

fid=fopen(fname,'w');
Nnode=size(xyz,1);
[Nele,ele_order]=size(t);
fprintf(fid,'# vtk DataFile Version 2.0  \n');
fprintf(fid,'%s\n',title);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %d double\n',Nnode);
for node=1:Nnode
    fprintf(fid,'%f\t%f\t%f\n',xyz(node,:));
end
cell_size=Nele*(ele_order+1);
fprintf(fid,'CELLS %d %d\n',Nele,cell_size);
for ele=1:Nele
    fprintf(fid,'\t%d',ele_order);
    for order=1:ele_order
        fprintf(fid,'\t%d',t(ele,order)-1);
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'CELL_TYPES %d\n',Nele);
cell_type=[1,3,5,10];
cell_type=cell_type(ele_order);

for ele=1:Nele
    fprintf(fid,'%d\n',cell_type);
end
if(nargin>=5)
    if(datatype==0)
        fprintf(fid,'POINT_DATA %d\n',Nnode);
        fprintf(fid,'SCALARS value double\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        for node=1:Nnode
            fprintf(fid,'%f\n',p(node));
        end
    else
        fprintf(fid,'CELL_DATA %d\n',Nele);
        fprintf(fid,'SCALARS value double\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        for ele=1:Nele
            fprintf(fid,'%f\n',p(ele));
        end
    end
end
fclose(fid);
return;
end

