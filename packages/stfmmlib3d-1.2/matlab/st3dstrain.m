function [strain]=st3dstrain(dudx)
% strain tensor
%
% strain = (dudx + dudx.')/2;

nsource = size(dudx,3);
strain = zeros(3,3,nsource);

strain(1,1,:) = dudx(1,1,:) + dudx(1,1,:);
strain(2,2,:) = dudx(2,2,:) + dudx(2,2,:);
strain(3,3,:) = dudx(3,3,:) + dudx(3,3,:);

strain(1,2,:) = dudx(1,2,:) + dudx(2,1,:);
strain(2,1,:) = dudx(1,2,:) + dudx(2,1,:);

strain(1,3,:) = dudx(1,3,:) + dudx(3,1,:);
strain(3,1,:) = dudx(1,3,:) + dudx(3,1,:);

strain(2,3,:) = dudx(2,3,:) + dudx(3,2,:);
strain(3,2,:) = dudx(2,3,:) + dudx(3,2,:);

strain = strain/2;
