function [U]=sthfmm3dpart_matlab(iprec,itype,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg)
%STHFMM3DPART Stokes half space particle target FMM in R^3.
%
% Stokes half space FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interaction) and interactions with targets.
%
% No slip (zero-velocity) boundary condition at z=0
%
% [U]=STHFMM3DPART(IPREC,ITYPE,NSOURCE,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,SIGMA_DV,IFPOT,IFGRAD);
%
% [U]=STHFMM3DPART(IPREC,ITYPE,NSOURCE,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,SIGMA_DV,IFPOT,IFGRAD,...
%         NTARGET,TARGET,IFPOTTARG,IFGRADTARG);
%
%
% This subroutine evaluates the Stokes potential and gradient due
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
% for the free-space Green's function, without the (1/4 pi) scaling.  
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
% ifgrad - gradient computation flag, 
%         1 => compute the gradient, otherwise no
%
% ntarget - number of targets
% target - double (3,ntarget): target locations
%
% ifpottarg - target velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgradtarg - target gradient computation flag, 
%         1 => compute the gradient, otherwise no
%
%
% Output parameters: 
%
% U.pot - double (3,nsource) - velocity field at source locations
% U.pre - double (nsource) - pressure at source locations
% U.grad - double (3,3,nsource) - gradient at source locations
% U.pottarg - double (3,ntarget) - velocity field at targets
% U.pretarg - double (ntarget) - pressure at targets
% U.gradtarg - double (3,3,ntarget) - gradient at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%


if( nargin == 8 ) 
  ifpot = 1;
  ifgrad = 1;
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 10 ) 
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 12 ) 
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

if_use_fmm = 1;

%
%  Direct arrival
%
if( itype == 1 ),

if( if_use_fmm == 1 ),
U_direct = stfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl,...
        ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,...
        ntarget,target,ifpottarg,ifgradtarg);
else
U_direct = st3dpartdirect(nsource,source,ifsingle,sigma_sl,...
        ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,...
        ntarget,target,ifpottarg,ifgradtarg);
end

end

%
%  Image
%

source_image = source; source_image(3,:) = -source_image(3,:);

sigma_sl_image = sigma_sl; sigma_sl_image(3,:) = -sigma_sl(3,:); 

sigma_dl_image = sigma_dl; sigma_dl_image(3,:) = -sigma_dl(3,:); 
sigma_dv_image = sigma_dv; sigma_dv_image(3,:) = -sigma_dv(3,:); 

%
%  Join target and source lists for image processing
%
ntarget_list = ntarget + nsource;
target_list = [target source];

ifpot0 = 0; ifgrad0 = 0;
ifpottarg0 = 0; ifgradtarg0 = 0;
if( ifpot == 1 || ifpottarg == 1 ), ifpottarg0 = 1; end
if( ifgrad == 1 || ifgradtarg == 1 ), ifgradtarg0 = 1; end

if( if_use_fmm == 1 ),
U_image = stfmm3dpart(iprec,nsource,source_image,ifsingle,sigma_sl_image,...
        ifdouble,sigma_dl_image,sigma_dv_image,ifpot0,ifgrad0,...
        ntarget_list,target_list,ifpottarg0,ifgradtarg0);
else
U_image = st3dpartdirect(nsource,source_image,ifsingle,sigma_sl_image,...
        ifdouble,sigma_dl_image,sigma_dv_image,ifpot0,ifgrad0,...
        ntarget_list,target_list,ifpottarg0,ifgradtarg0);
end

%
%  Papkovich-Neuber potential
%

if( ifsingle == 1 ),

  ifcharge = 1;
  charge = sigma_sl_image(3,:);

  ifdipole = 1;
  dipstr = ones(1,nsource);
  dipvec = repmat(source(3,:),3,1) .* sigma_sl_image;

else

  ifcharge = 0;
  charge = zeros(1,nsource);
  ifdipole = 0;

  dipstr = zeros(1,nsource);
  dipvec = zeros(3,nsource);

end

ifquad  = 0;
quadstr = zeros(1,nsource);
quadvec = zeros(6,nsource);



if( ifdouble == 1 ),
  ifdipole = 1;
  dipstr = ones(1,nsource);
  dipvec = dipvec + [zeros(2,nsource); ... 
         2*(sigma_dl_image(1,:).*sigma_dv_image(1,:)+...
            sigma_dl_image(2,:).*sigma_dv_image(2,:)+...
            sigma_dl_image(3,:).*sigma_dv_image(3,:))];
elseif( ifdouble == 2 ),
elseif( ifdouble == 3 || ifdouble == 4 ),
  ifdipole = 1;
  dipstr = ones(1,nsource);
  dipvec = dipvec + (-repmat(sigma_dv_image(3,:),3,1).*sigma_dl_image+...
            +repmat(sigma_dl_image(3,:),3,1).*sigma_dv_image)*2;
else
end

if( ifdouble == 1 || ifdouble == 2 || ifdouble == 4 ),
  ifquad  = 1;
  quadstr = source(3,:)*2;
  quadvec(1,:) = sigma_dl_image(1,:).*sigma_dv_image(1,:); 
  quadvec(2,:) = sigma_dl_image(2,:).*sigma_dv_image(2,:); 
  quadvec(3,:) = sigma_dl_image(3,:).*sigma_dv_image(3,:); 
  quadvec(4,:) = sigma_dl_image(1,:).*sigma_dv_image(2,:)+...
                 sigma_dl_image(2,:).*sigma_dv_image(1,:);
  quadvec(5,:) = sigma_dl_image(1,:).*sigma_dv_image(3,:)+...
                 sigma_dl_image(3,:).*sigma_dv_image(1,:);
  quadvec(6,:) = sigma_dl_image(2,:).*sigma_dv_image(3,:)+...
                 sigma_dl_image(3,:).*sigma_dv_image(2,:);
end

ifpot0 = 0; iffld0 = 0; ifhess0 = 0;
ifpottarg0 = 1; iffldtarg0 = 1; ifhesstarg0 = 0;
if( ifgrad == 1 || ifgradtarg == 1 ), ifhesstarg0 = 1; end;

if( if_use_fmm == 1 ),
H=lfmm3dpartquad(iprec,nsource,source_image,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifpot0,iffld0,ifhess0,...
        ntarget_list,target_list,ifpottarg0,iffldtarg0,ifhesstarg0);
else
H=l3dpartquaddirect(nsource,source_image,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifpot0,iffld0,ifhess0,...
        ntarget_list,target_list,ifpottarg0,iffldtarg0,ifhesstarg0);
end
if( ifpottarg0 == 1 ), H.pottarg = real(H.pottarg); end;
if( iffldtarg0 == 1 ), H.fldtarg = real(H.fldtarg); end;
if( ifhesstarg0 == 1 ), H.hesstarg = real(H.hesstarg); end;


if( ifpot == 1 ),

s = (ntarget+1) : (ntarget+nsource);

F.pot = -repmat(source(3,:),3,1) .* H.fldtarg(:,s) ...
            -[zeros(2,nsource); H.pottarg(:,s)];
F.pre = -2 .* H.fldtarg(3,s);

if( itype == 1 ),
U.pot = U_direct.pot - U_image.pottarg(:,s) - F.pot;
U.pre = U_direct.pre - U_image.pretarg(:,s) - F.pre;
else
U.pot = - U_image.pottarg(:,s) - F.pot;
U.pre = - U_image.pretarg(:,s) - F.pre;
end

end


if( ifpottarg == 1 ),

t = 1:ntarget;

F.pottarg = -repmat(target(3,:),3,1) .* H.fldtarg(:,t) ...
            -[zeros(2,ntarget); H.pottarg(:,t)];
F.pretarg = -2 .* H.fldtarg(3,t);

if( itype == 1 ),
U.pottarg = U_direct.pottarg - U_image.pottarg(:,t) - F.pottarg;
U.pretarg = U_direct.pretarg - U_image.pretarg(:,t) - F.pretarg;
else
U.pottarg = - U_image.pottarg(:,t) - F.pottarg;
U.pretarg = - U_image.pretarg(:,t) - F.pretarg;
end

end


if( 2 == 2 ),
% to do
if( ifgrad == 1 ),

s = (ntarget+1) : (ntarget+nsource);

F.grad = zeros(3,3,nsource);
F.grad(1,1,:) = H.hesstarg(1,s).*source(3,:);
F.grad(2,2,:) = H.hesstarg(2,s).*source(3,:);
F.grad(3,3,:) = H.hesstarg(3,s).*source(3,:);
F.grad(1,2,:) = H.hesstarg(4,s).*source(3,:);
F.grad(1,3,:) = H.hesstarg(5,s).*source(3,:);
F.grad(2,3,:) = H.hesstarg(6,s).*source(3,:);
F.grad(2,1,:) = H.hesstarg(4,s).*source(3,:);
F.grad(3,1,:) = H.hesstarg(5,s).*source(3,:);
F.grad(3,2,:) = H.hesstarg(6,s).*source(3,:);

F.grad(1:3,3,:) = reshape(F.grad(1:3,3,:),3,nsource) - H.fldtarg(:,s);
F.grad(3,1:3,:) = reshape(F.grad(3,1:3,:),3,nsource) + H.fldtarg(:,s);

if( itype == 1 ),
U.grad = U_direct.grad - U_image.gradtarg(:,:,s) - F.grad;
else
U.grad = - U_image.gradtarg(:,:,s) - F.grad;
end

end


if( ifgradtarg == 1 ),

t = 1:ntarget;

F.gradtarg = zeros(3,3,ntarget);
F.gradtarg(1,1,:) = H.hesstarg(1,t).*target(3,:);
F.gradtarg(2,2,:) = H.hesstarg(2,t).*target(3,:);
F.gradtarg(3,3,:) = H.hesstarg(3,t).*target(3,:);
F.gradtarg(1,2,:) = H.hesstarg(4,t).*target(3,:);
F.gradtarg(1,3,:) = H.hesstarg(5,t).*target(3,:);
F.gradtarg(2,3,:) = H.hesstarg(6,t).*target(3,:);
F.gradtarg(2,1,:) = H.hesstarg(4,t).*target(3,:);
F.gradtarg(3,1,:) = H.hesstarg(5,t).*target(3,:);
F.gradtarg(3,2,:) = H.hesstarg(6,t).*target(3,:);

F.gradtarg(1:3,3,:) = reshape(F.gradtarg(1:3,3,:),3,ntarget) - H.fldtarg(:,t);
F.gradtarg(3,1:3,:) = reshape(F.gradtarg(3,1:3,:),3,ntarget) + H.fldtarg(:,t);

if( itype == 1 ),
U.gradtarg = U_direct.gradtarg - U_image.gradtarg(:,:,t) - F.gradtarg;
else
U.gradtarg = - U_image.gradtarg(:,:,t) - F.gradtarg;
end

end

end
