% visualize perturbed gravity 
clear all;clc;
addpath('../../modelbuilder/'); 
addpath('../../packages/fmmlib3d-1.2/matlab/')
Globals3D;

% input parameters ---------------
fmesh  = '../../Demos/PREM/output/PREM3k/';
fout   = '../demos/PREM3k/';
fbase  = 'prem_3L_3k.1';
fdtail = '0.1000000_1.000000';

JOB = 2; pOrder = 2; nproc = 48; nth = 4; 
Radial = 6.371E3;

% memory control 
Gpsiz = 1e7; 
% end input parameters ------------------

% read eigenvectors
tic;
readEigvct; 

% calculate density changes 
StartUp3D;
sources;
toc; 
% compute the PG 
computePG; 
toc; 
% visualize the results
visualPGall;
toc;
