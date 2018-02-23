%
%  Test Stokes particle FMMs in R^3
%

% itype - half space Green's function evaluation flag
%         1 => include both direct arrival and image contribution
%         2 => include image contribution only
itype = 1

nsource = 2000

source = zeros(3,nsource);

idist=1;

if( idist == 1 ),
theta=rand(1,nsource)*pi;
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi).*sin(theta);
source(2,:)=.5*sin(phi).*sin(theta);
source(3,:)=.5*cos(theta);
source(3,:)=source(3,:) + .5;
end

if( idist == 2 ),
source(1,:)=rand(1,nsource);
source(2,:)=rand(1,nsource);
source(3,:)=rand(1,nsource);
end

%
%  timings
%

ifsingle=1;
sigma_sl = rand(3,nsource);
ifdouble=4;
sigma_dl = rand(3,nsource);
sigma_dv = rand(3,nsource);


ifsingle
ifdouble
ifpot = 1
ifgrad = 1


ntarget = min(10,nsource);
target = source(:,1:nsource);
target(1,:) = target(1,:) + 10;
[ndim,ntarget] = size(target);

target(3,:) = 0;

%%%ntarget = 0;

ntarget
ifpottarg = 1
ifgradtarg = 1

if( ntarget == 0 ),
ifpottarg = 0
ifgradtarg = 0
end


disp('')
'Stokes half-space particle target FMM in R^3 (Fortran)'

tic
iprec=1
[U]=sthfmm3dpart(iprec,itype,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);
total_time=toc

'Stokes half-space particle target FMM in R^3 (Matlab)'

tic
iprec=1
[F]=sthfmm3dpart_matlab(iprec,itype,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);
total_time=toc


if( ifpot ), U.pot=U.pot/(4*pi); end
if( ifpot ), U.pre=U.pre/(4*pi); end
if( ifgrad ), U.grad=U.grad/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( ifpot ), F.pre=F.pre/(4*pi); end
if( ifgrad ), F.grad=F.grad/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
rms_error_pre = norm((U.pre - F.pre),2)/sqrt(nsource)
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2)
rel_error_pre = norm((U.pre - F.pre),2)/norm((F.pre),2)
end

if( ifgrad ),
%rms_grad = norm(reshape(F.grad,9,nsource),2)/sqrt(nsource)
rms_error_grad = norm(reshape(U.grad - F.grad,9,nsource),2)/sqrt(nsource)
rel_error_grad = norm(reshape(U.grad - F.grad,9,nsource),2)/norm(reshape(F.grad,9,nsource),2)
end
%%%break;

if( ifpottarg ), U.pottarg=U.pottarg/(4*pi); end
if( ifpottarg ), U.pretarg=U.pretarg/(4*pi); end
if( ifgradtarg ), U.gradtarg=U.gradtarg/(4*pi); end

if( ifpottarg ), F.pottarg=F.pottarg/(4*pi); end
if( ifpottarg ), F.pretarg=F.pretarg/(4*pi); end
if( ifgradtarg ), F.gradtarg=F.gradtarg/(4*pi); end

if( ifpottarg ),
%rms_pottarg = norm((F.pottarg),2)/sqrt(nsource)
rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(ntarget)
rms_error_pretarg = norm((U.pretarg - F.pretarg),2)/sqrt(ntarget)
norm_pottarg = norm((F.pottarg),2)
rel_error_pottarg = norm((U.pottarg - F.pottarg),2)/norm((F.pottarg),2)
rel_error_pretarg = norm((U.pretarg - F.pretarg),2)/norm((F.pretarg),2)
end

if( ifgradtarg ),
rms_gradtarg = norm(reshape(F.gradtarg,9,ntarget),2)/sqrt(ntarget)
rms_error_gradtarg = ...
    norm(reshape(U.gradtarg - F.gradtarg,9,ntarget),2)/sqrt(ntarget)
rel_error_gradtarg = ...
    norm(reshape(U.gradtarg - F.gradtarg,9,ntarget),2)/ ...
    norm(reshape(F.gradtarg,9,ntarget),2)
end
%%%break;

