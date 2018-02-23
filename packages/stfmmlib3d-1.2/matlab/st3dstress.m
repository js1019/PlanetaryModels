function [stress]=st3dstress(p,dudx)
% stress tensor, for mu = 1 
%
% stress = -p + 2*mu (dudx + dudx.')/2;

nsource = size(dudx,3);
stress = zeros(3,3,nsource);

stress(1,1,:) = dudx(1,1,:) + dudx(1,1,:);
stress(2,2,:) = dudx(2,2,:) + dudx(2,2,:);
stress(3,3,:) = dudx(3,3,:) + dudx(3,3,:);

stress(1,2,:) = dudx(1,2,:) + dudx(2,1,:);
stress(2,1,:) = dudx(1,2,:) + dudx(2,1,:);

stress(1,3,:) = dudx(1,3,:) + dudx(3,1,:);
stress(3,1,:) = dudx(1,3,:) + dudx(3,1,:);

stress(2,3,:) = dudx(2,3,:) + dudx(3,2,:);
stress(3,2,:) = dudx(2,3,:) + dudx(3,2,:);

%%%stress = stress;

stress(1,1,:) = squeeze(stress(1,1,:)) - p(:);
stress(2,2,:) = squeeze(stress(2,2,:)) - p(:);
stress(3,3,:) = squeeze(stress(3,3,:)) - p(:);

