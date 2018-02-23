function [U]=sthfmm3dpart(iprec,itype,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg)
%STHFMM3DPART Stokes half space particle target FMM in R^3.
%
% Stokes half space FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interaction) and interactions with targets.
%
% No slip (zero-velocity) boundary condition at z=0
%
% [U]=STHFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,SIGMA_DV,IFPOT,IFGRAD);
%
% [U]=STHFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,SIGMA_DV,IFPOT,IFGRAD,...
%         NTARGET,TARGET,IFPOTTARG,IFGRADTARG);
%
%
% This subroutine evaluates the Stokes potential (velocity/pressure) 
% and velocity gradient due
% to a collection of Stokes single and double forces. We use
%

%
%       ifsingle=1, stokeslet, f = sigma_sl
%       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
%       p = [r_j / r^3] f_j
%
%       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
%       u_i = [3 r_i r_j r_k / r^5] n_k g_j
%       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
%
%       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
%       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
%       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
%
%       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
%       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
%       p = 0
%
%       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
%                   g = sigma_dl, n = sigma_dv
%       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
%             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
%       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
%
% for the Green's function, without the (1/4 pi) scaling.  
%
% Half-space Green's function is the combination of direct arrival, 
% image contribution and Papkovich-Neuber correction.
%
% Self-interactions due to *direct-arrival* are not-included.
%
%
% Input parameters:
% 
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
%
% itype - half space Green's function evaluation flag
%         1 => include both direct arrival and image contribution
%         2 => include image contribution only
%
% nsource - number of sources
% source - double (3,nsource): source locations
% ifsingle - single force computation flag
%
%         0 => do not compute
%         1 => include Stokes single force contribution
% 
% sigma_sl - double (3,nsource): single force strengths
% ifdouble - double force computation flag
%
%         0 => do not compute
%         1 => include Stokes double force contribution
%         2 => include Stokes stresslet contribution
%         3 => include Stokes rotlet contribution
%         4 => include Stokes doublet contribution
% 
% sigma_dl - double (3,nsource): double force strengths
% sigma_dv - double (3,nsource): double force orientation vectors 
%
% ifpot - velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgrad - velocity gradient computation flag, 
%         1 => compute the velocity gradient, otherwise no
%
% ntarget - number of targets
% target - double (3,ntarget): target locations
%
% ifpottarg - target velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgradtarg - target velocity gradient computation flag, 
%         1 => compute the velocity gradient, otherwise no
%
%
% Output parameters: 
%
% U.pot - double (3,nsource) - velocity field at source locations
% U.pre - double (nsource) - pressure at source locations
% U.grad - double (3,3,nsource) - velocity gradient at source locations
% U.pottarg - double (3,ntarget) - velocity field at targets
% U.pretarg - double (ntarget) - pressure at targets
% U.gradtarg - double (3,3,ntarget) - velocity gradient at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%

if( nargin == 9 ) 
  ifpot = 1;
  ifgrad = 1;
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 11 ) 
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 13 ) 
  ifpottarg = 1;
  ifgradtarg = 1;
end

ifsingle = double(ifsingle); ifdouble = double(ifdouble);
ifpot = double(ifpot); ifgrad = double(ifgrad); 
ifpottarg = double(ifpottarg); ifgradtarg = double(ifgradtarg); 

pot=zeros(3,nsource);
pre=zeros(1,nsource);
grad=zeros(3,3,nsource);
if( ntarget > 0 ),
pottarg=zeros(3,ntarget);
pretarg=zeros(1,ntarget);
gradtarg=zeros(3,3,ntarget);
else
pottarg=zeros(3,1);
pretarg=zeros(1,1);
gradtarg=zeros(3,3,1);
end
ier=0;

if( ntarget == 0 ) 
mex_id_ = 'sthfmm3dpartself(io int[x], i int[x], i int[x], i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i double[xx], i int[x], io double[], io double[], i int[x], io double[])';
[ier, pot, pre, grad] = stfmm3d_r2012a(mex_id_, ier, iprec, itype, nsource, source, ifsingle, sigma_sl, ifdouble, sigma_dl, sigma_dv, ifpot, pot, pre, ifgrad, grad, 1, 1, 1, 1, 3, nsource, 1, 1, 3, nsource, 1, 1);
else
mex_id_ = 'sthfmm3dparttarg(io int[x], i int[x], i int[x], i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i double[xx], i int[x], io double[], io double[], i int[x], io double[], i int[x], i double[], i int[x], io double[], io double[], i int[x], io double[])';
[ier, pot, pre, grad, pottarg, pretarg, gradtarg] = stfmm3d_r2012a(mex_id_, ier, iprec, itype, nsource, source, ifsingle, sigma_sl, ifdouble, sigma_dl, sigma_dv, ifpot, pot, pre, ifgrad, grad, ntarget, target, ifpottarg, pottarg, pretarg, ifgradtarg, gradtarg, 1, 1, 1, 1, 3, nsource, 1, 1, 3, nsource, 1, 1, 1, 1, 1);
end

if( ifpot == 1 ), U.pot=pot; end
if( ifpot == 1 ), U.pre=pre; end
if( ifgrad == 1 ), U.grad=reshape(grad,3,3,nsource); end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( ifpottarg == 1 ), U.pretarg=pretarg; end
if( ifgradtarg == 1 ), U.gradtarg=reshape(gradtarg,3,3,ntarget); end
U.ier=ier;


