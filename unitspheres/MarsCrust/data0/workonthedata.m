% work on the data

%% thickness
clear all; clc;
fname = '../data/crustal_thickness_gmm3_rm1_constructed_density.grd';
finfo = ncinfo(fname);

depth0 = ncread(fname,'z');
% downsample change from meter to km
depth = depth0(1:4:end,end:-4:1)'; 

%% 
fname = '../data/megtcb.txt';
topo0 = textread(fname);
topo1 = reshape(topo0,1440,720); 
% extend the data
topo = zeros(721,1441); 
topo(1:720,1:1440) = topo1'*1e-3;
topo(721,1:1440) = topo(1,1:1440); 
topo(1:720,1441) = topo(1:720,1); 
topo(721,1441) = topo(1,1);

%figure; subplot(2,1,1); imagesc(topo); subplot(2,1,2); imagesc(depth)

save('mars_crust','depth','topo');

