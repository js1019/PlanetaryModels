% Purpose : Setup script, building operators, grid, metric,
%           and connectivity tables for 3D meshes of tetrahedra.
Globals3D;

% Just use N = 1;
N = 1;

% Definition of constants
Np = (N+1)*(N+2)*(N+3)/6; Nfp = (N+1)*(N+2)/2; Nfaces=4; NODETOL = 1e-5;

% Compute nodal set
[x,y,z] = Nodes3D(N); 
%nodel_plot_3d(x,y,z,0.08);
[r,s,t] = xyztorst(x,y,z);
%%%%% x,y,z: alpha-opt nodes' coordiantes in equiv-tetrahegen
%%%%% r,s,t: nodes' coordinates inside reference tetrahedron, i.e. xi,eta,zeta

% Build reference element matrices
V = Vandermonde3D(N,r,s,t); invV = inv(V); %%%%% Vandermonde Matrix & inverse
MassMatrix = invV'*invV;
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);
%%% {\partial \over \partial x_i} = {\partial r \over \partial x_i}Dr +
%%% {\partial s \over \partial x_i}Ds + {\partial r \over \partial x_i}Dt;\\
%%% Dr V = Vr; 
%%%%% Dr,Ds,Dt: Stiffmatrix*Massmatrix, i.e. M * K_x/y/z


% Build connectivity matrix
%[EToE, EToF] = tiConnect3D(EToV); 
%%%%% EToE: neighbouring element's ID, correspond to neigh
%%%%% EToF: neighbouting element's face, correspond to nmap

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)'; vd = EToV(:,4)';
x = 0.5*(-(1+r+s+t)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc)+(1+t)*VX(vd));
y = 0.5*(-(1+r+s+t)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc)+(1+t)*VY(vd));
z = 0.5*(-(1+r+s+t)*VZ(va)+(1+r)*VZ(vb)+(1+s)*VZ(vc)+(1+t)*VZ(vd));
%%%%% Global coordinates of Np nodes in each element
%%%%% x(in_ele_nodeid(1:Np), ele_id(1:Nele))

% find all the nodes that lie on each edge
                                                         %     t
%fmask1   = find( abs(1+t) < NODETOL)';                   %     |  4
%fmask2   = find( abs(1+s) < NODETOL)';                   %   2 |____s
%fmask3   = find( abs(1+r+s+t) < NODETOL)';               %     /
%fmask4   = find( abs(1+r) < NODETOL)';                   %   r/  1

fmask1   = find( abs(1+r+s+t) < NODETOL)'; 
fmask2   = find( abs(1+r) < NODETOL)'; 
fmask3   = find( abs(1+s) < NODETOL)';
fmask4   = find( abs(1+t) < NODETOL)';

Fmask  = [fmask1;fmask2;fmask3;fmask4]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :); Fz = z(Fmask(:), :);
%%%%% fmaski(1:Nfp): Local IDs' of nodes on i^th face

% Create surface integral terms
% LIFT = Lift3D(r, s, t);

% calculate geometric factors
[rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = GeometricFactors3D(x,y,z,Dr,Ds,Dt);
%%%%% J:Jacobian
%%%%% J=[xr,xs,xt;yr,ys,yt;zr,zs,zt]
%%%%% J^{-1}=[rx,sx,tx;ry,sy,ty;rz,sz,tz]

% calculate geometric factors
[nx, ny, nz, sJ] = Normals3D();
Fscale = sJ./(J(Fmask,:));
%%%%% nx,ny,nz: coordinates of outer normal vector
%%%%% n(in_ele_nodeid(1:4Nfp), ele_id(1:Nele))
%%%%% sJ: scale of each face to reference triangle (rectangular)

% Build connectivity maps
%[vmapM, vmapP, vmapB, mapB] = BuildMaps3D();


% Compute weak operators (could be done in preprocessing to save time)
%[Vr, Vs, Vt] = GradVandermonde3D(N, r, s, t);
%Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V'); Dtw = (V*Vt')/(V*V');
%%%%% Drw,Dsw,Dtw: Stiffness Matrix

% for 2D surfaces
Fm = Fmask(:,1); faceR = s(Fm); faceS = t(Fm);
V2D = Vandermonde2D(N, faceR, faceS);  vin = inv(V2D);
Mass2D = vin'*vin;

