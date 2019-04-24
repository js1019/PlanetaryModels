% build basic uniform mesh on Unit Sphere
clear all; clc
%parpool('local');
addpath('../../packages/distmesh/')

Ellp = 589e-5;
%fd=@(p) dsphere(p,0,0,0,1);
fd=@(p) p(:,1).^2+p(:,2).^2+p(:,3).^2/(1.-Ellp)^2-1;

surf_d = 0.04; 
tic
[p,t]=distmeshsurface(fd,@huniform,surf_d,1.1*[-1,-1,-1;1,1,1]);
toc

fname = ['MarsEllp',num2str(size(t,1)),'.mat'];
save(fname,'Ellp','fd','p','surf_d','t');

% for reproducing
% surf_d = 0.36;  p 132;   t 260
% surf_d = 0.3;   p 198;   t 392
% surf_d = 0.2;   p 480;   t 956
% surf_d = 0.1;   p 1806;  t 3608
% surf_d = 0.08;  p 3k;    t 6k
% surf_d = 0.05;  p 7k;    t 15k
% surf_d = 0.03;  p 21k;   t 42k
% surf_d = 0.02;  p 47k;   t 94k
% surf_d = 0.015; p 84k;   t 167k
% surf_d = 0.01;  p 188k;  t 377k
% surf_d = 0.008; p 294k;  t 589k
% surf_d = 0.006; p 523k;  t 1047k?
% surf_d = 0.004; p 1177k; t 2353k time:361294.354594
% surf_d = 0.0035; p 1538k; t 3077k time: 771746.288223


