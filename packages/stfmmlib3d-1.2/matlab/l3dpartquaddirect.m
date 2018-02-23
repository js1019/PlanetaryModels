function [U]=l3dpartquaddirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifpot,iffld,ifhess,ntarget,target,ifpottarg,iffldtarg,ifhesstarg)
%LFMM3DPARHESSTARG Laplace interactions in R^3, direct evaluation.
%
% Laplace FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interactions) and interactions with targets.
%
% [U]=L3DPARTQUADDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFQUAD,QUADSTR,QUADVEC);
%
% [U]=L3DPARTQUADDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFQUAD,QUADSTR,QUADVEC,...
%         IFPOT,IFFLD,IFHESS);
%
% [U]=L3DPARTQUADDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFQUAD,QUADSTR,QUADVEC,...
%         IFPOT,IFFLD,IFHESS,...
%         NTARGET,TARGET);
%
% [U]=L3DPARTQUADDIRECT(NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFQUAD,QUADSTR,QUADVEC,...
%         IFPOT,IFFLD,IFHESS,...
%         NTARGET,TARGET,IFPOTTARG,IFFLDTARG,IFHESSTARG);
%
%
% This subroutine evaluates the Laplace potential, field, and hessian due
% to a collection of charges, dipoles, and quadrupoles. We use
%
%     pot = charge / r + 
%           dipstr*  (dipvec(1)*U_x + dipvec(2)*U_y + dipvec(3)*U_z) +
%           quadstr* (quadvec(1)*V_xx + quadvec(2)*V_yy + quadvec(3)*V_zz+
%                     quadvec(4)*V_xy + quadvec(5)*V_xz + quadvec(6)*V_yz)
%
%     fld = -grad(pot)
%     hess = (potxx,potyy,potzz,potxy,potxz,potyz)
%
%     U_x = dx/r^3, U_y = dy/r^3, U_z = dz/r^3
%
%     V_xx = (-1/r^3 + 3*dx**2/r^5)
%     V_xy = 3*dx*dy/r^5
%     V_xz = 3*dx*dz/r^5
%     V_yy = (-1/r^3 + 3*dy**2/r^5)
%     V_yz = 3*dy*dz/r^5
%     V_zz = (-1/r^3 + 3*dz**2/r^5)
%
% for the Green's function, without the (1/4 pi) scaling. 
% Self-interactions are not-included.
%
%
% Input parameters:
% 
% nsource - number of sources
% source - real (3,nsource): source locations
% ifcharge - charge computation flag
%
%         0 => do not compute
%         1 => include charge contribution
% 
% charge - complex (nsource): charge strengths 
% ifdipole - dipole computation flag
%
%         0 => do not compute
%         1 => include dipole contributions
% 
% dipole - complex (nsource): dipole strengths
% dipvec - real (3,source): dipole orientation vectors
% ifquad - quadrupole computation flag
%
%         0 => do not compute
%         1 => include quadrupole contributions
% 
% quadstr - complex (nsource): quadrupole strengths
% quadvec - real (6,source): quadrupole orientation vectors
%
% ifpot - potential computation flag, 1 => compute the potential, otherwise no
% iffld - field computation flag, 1 => compute the field, otherwise no
% ifhess - hessian computation flag, 1 => compute the hessian, otherwise no
%
% ntarget - number of targets
% target - real (3,ntarget): target locations
%
% ifpottarg - target potential computation flag, 
%      1 => compute the target potential, otherwise no
% iffldtarg - target field computation flag, 
%      1 => compute the target field, otherwise no
% ifhesstarg - target hessian computation flag, 
%      1 => compute the target hessian, otherwise no
%
% Output parameters: 
%
% U.pot - complex (nsource) - potential at source locations
% U.fld - complex (3,nsource) - field (i.e. -gradient) at source locations
% U.hess - complex (6,nsource) - hessian at source locations
% U.pottarg - complex (ntarget) - potential at target locations
% U.fldtarg - complex (3,ntarget) - field (i.e. -gradient) at target locations
% U.hesstarg - complex (6,ntarget) - hessian at target locations
%
% U.ier - error return code
%
%             ier=0     =>  normal execution

if( nargin == 10 ) 
  ifpot = 1;
  iffld = 1;
  ifhess = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
  ifhesstarg = 0;
end

if( nargin == 13 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
  ifhesstarg = 0;
end

if( nargin == 15 ) 
  ifpottarg = 1;
  iffldtarg = 1;
  ifhesstarg = 1;
end

ifcharge = double(ifcharge); ifdipole = double(ifdipole);
ifquad = double(ifquad);
ifpot = double(ifpot); iffld = double(iffld); ifhess = double(ifhess);
ifpottarg = double(ifpottarg); iffldtarg = double(iffldtarg); 
ifhesstarg = double(ifhesstarg);

pot=0;
fld=zeros(3,1);
hess=zeros(6,1);
pottarg=0;
fldtarg=zeros(3,1);
hesstarg=zeros(6,1);

if( ifpot == 1 ), pot=zeros(1,nsource)+1i*zeros(1,nsource); end;
if( iffld == 1 ), fld=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifhess == 1 ), hess=zeros(6,nsource)+1i*zeros(6,nsource); end;
if( ifpottarg == 1 ), pottarg=zeros(1,ntarget)+1i*zeros(1,ntarget); end;
if( iffldtarg == 1 ), fldtarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;
if( ifhesstarg == 1 ), hesstarg=zeros(6,ntarget)+1i*zeros(6,ntarget); end;

ier=0;

mex_id_ = 'l3dpartquaddirect(i int[x], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[pot, fld, hess, pottarg, fldtarg, hesstarg] = stfmm3d_r2012a(mex_id_, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifquad, quadstr, quadvec, ifpot, pot, iffld, fld, ifhess, hess, ntarget, target, ifpottarg, pottarg, iffldtarg, fldtarg, ifhesstarg, hesstarg, 1, 3, nsource, 1, 1, 3, nsource, 1, 6, nsource, 1, 1, 1, 1, 1, 1, 1);


if( ifpot == 1 ), U.pot=pot; end
if( iffld == 1 ), U.fld=fld; end
if( ifhess == 1 ), U.hess=hess; end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
if( ifhesstarg == 1 ), U.hesstarg=hesstarg; end
U.ier=ier;



