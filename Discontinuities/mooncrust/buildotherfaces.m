% draw inner/outer core boundary
clear all; clc;

R = 1737.15;
r0 = 330.0;

load ../data/Sph94k.mat

p = p*r0/R;

fsvtk = './vtk/Mweberoutercore';

vtktrisurf(t, p(:,1), p(:,2), p(:,3),'cmb', fsvtk);