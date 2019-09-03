% plot reference interfaces
clear all; clc;
addpath('../../../modelbuilder/');
load ../../../radialmodels/moon/Moon6L.mat;

R1 = RD(2,1); R2 = RD(3,1); R3 = RD(4,1);
R4 = RD(5,1); R5 = RD(6,1); R6 = RD(7,1);

load ../../../discontinuities/data/Sph956.mat
p2 = p*R1; %/R1; 

vtk_write_general('Moon_ICB_org_face.vtk','test',p2,t);




