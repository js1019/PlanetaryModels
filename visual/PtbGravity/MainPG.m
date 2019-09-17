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

% for the density file
fname0 = '../../Demos/PREM/output/PREM3k/';
fname1  = [fname0,fbase]; 
% memory control 
Gpsiz = 1e7; 
% scaling scale the total size
scaling = 6371; 
% end input parameters ------------------

% read eigenvectors
readEigvct; 

% calculate density changes 
StartUp3D;
sources;

% compute the PG 
computePG; 

% visualize the results
visualPGall;