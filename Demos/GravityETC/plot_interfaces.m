% plot reference interfaces
clear all; clc;
addpath('../../modelbuilder/');
load ../../radialmodels/PREM/prem3L_noocean.mat

R1 = RD(1,1); R2 = RD(2,1); R3 = RD(3,1);

load ../../discontinuities/data/Sph94k.mat
p2 = p*R1; %/R1; 

vtk_write_general(['./Surf94k_org_face.vtk'],'test',p2,t);




