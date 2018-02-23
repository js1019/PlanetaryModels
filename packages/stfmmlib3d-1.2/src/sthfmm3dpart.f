cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains the FMM routines for Stokes particle
c        potentials in half space in R^3. 
c
c        No slip (zero-velocity) boundary condition at z=0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      
c      sthfmm3dpartself - Stokes half space FMM in R^3: 
c         evaluate all pairwise particle interactions 
c         (ignoring self-interaction)
c
c      sthfmm3dparttarg - Stokes half space FMM in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c      sth3dpartdirect - Stokes half space interactions in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine sthfmm3dpartself(ier,iprec,itype,
     $     nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad)
c
c     FMM calculation subroutine for Stokes N-body problem.
c
c     Half space Green's function.
c     No slip (zero-velocity) boundary condition at z=0
c
c     This subroutine evaluates the Stokes potential and gradient due
c     to a collection of Stokes single and double forces. We use
c
c       ifsingle=1, stokeslet, f = sigma_sl
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
c       u_i = [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
c       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 0
c
c       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
c                   g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
c             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c     for the free-space Green's function, without the (1/4 pi) scaling.  
c
c     Half-space Green's function is the combination of direct arrival, 
c     image contribution and Papkovich-Neuber correction.
c
c
c     input:
c
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c
c     iprec:  FMM precision flag
c     itype: half space Green's function evaluation flag
c         1 => include both direct arrival and image contribution
c         2 => include image contribution only
c
C     OUTPUT:
C
c     pot(3,nparts) = velocity at source locations
c     pre(nparts) = pressure at source locations
c     grad(3,3,nparts) = grad at source locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 pot(3,nparts),pre(nparts),train(3,3,nparts)
        integer nparts,ntargs

        ntargs=0
        ifpottarg=0
        ifgradtarg=0
        call sthfmm3dparttarg
     $     (ier,iprec,itype,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)

        return
        end
c
c
c
c
c
        subroutine sthfmm3dparttarg(ier,iprec,itype,
     $     nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,ntarget,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
c
c     Stokes interactions in R^3: evaluate all pairwise particle
c     interactions (excluding self interactions) and interactions with
c     targets using the direct O(N^2) algorithm.
c
c     Half space Green's function.
c     No slip (zero-velocity) boundary condition at z=0
c
c     This subroutine evaluates the Stokes potential and gradient due
c     to a collection of Stokes single and double forces. We use
c
c       ifsingle=1, stokeslet, f = sigma_sl
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
c       u_i = [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
c       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 0
c
c       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
c                   g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
c             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c     for the free-space Green's function, without the (1/4 pi) scaling.  
c
c     Half-space Green's function is the combination of direct arrival, 
c     image contribution and Papkovich-Neuber correction.
c
c
c       INPUT:
c
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       sigma_dv(3,nsource) - orientation of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifpot - velocity/pressure computation flag
c       ifgrad - grad computation flag
c       ifpottarg - target velocity/pressure computation flag
c       ifgradtarg - target grad computation flag
c
c       iprec: FMM precision flag
c
c       itype: half space Green's function evaluation flag
c         1 => include both direct arrival and image contribution
c         2 => include image contribution only
c
c       OUTPUT:
c
c       pot(3,nsource) - velocity at source locations
c       pre(nsource) - pressure at source locations
c       grad(3,3,nsource) - grad at source locations
c       pottarg(3,ntarget) - velocity at target locations
c       pretarg(ntarget) - pressure at target locations
c       gradtarg(3,3,ntarget) - grad at target locations
c
        real *8 source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8, allocatable :: source_image(:,:)
        real *8, allocatable :: sigma_sl_image(:,:)
        real *8, allocatable :: sigma_dl_image(:,:)
        real *8, allocatable :: sigma_dv_image(:,:)
c
        real *8, allocatable :: pot0(:,:)
        real *8, allocatable :: pre0(:)
        real *8, allocatable :: grad0(:,:,:)
c
        real *8, allocatable :: target0(:,:)
        real *8, allocatable :: pottarg0(:,:)
        real *8, allocatable :: pretarg0(:)
        real *8, allocatable :: gradtarg0(:,:,:)
c
        complex *16, allocatable :: charge(:)
        complex *16, allocatable :: dipstr(:)
        complex *16, allocatable :: quadstr(:)
        real *8, allocatable :: dipvec(:,:)
        real *8, allocatable :: quadvec(:,:)
c
        complex *16, allocatable :: cpot0(:)
        complex *16, allocatable :: cfld0(:,:)
        complex *16, allocatable :: chess0(:,:)
c
        complex *16, allocatable :: cpottarg0(:)
        complex *16, allocatable :: cfldtarg0(:,:)
        complex *16, allocatable :: chesstarg0(:,:)
c
c
        do i=1,nsource
        if( ifpot .eq. 1) then
           pot(1,i)=0
           pot(2,i)=0
           pot(3,i)=0
           pre(i)=0
        endif
        if( ifgrad .eq. 1) then
           grad(1,1,i)=0
           grad(2,1,i)=0
           grad(3,1,i)=0
           grad(1,2,i)=0
           grad(2,2,i)=0
           grad(3,2,i)=0
           grad(1,3,i)=0
           grad(2,3,i)=0
           grad(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) then
           pottarg(1,i)=0
           pottarg(2,i)=0
           pottarg(3,i)=0
           pretarg(i)=0
        endif
        if( ifgradtarg .eq. 1) then
           gradtarg(1,1,i)=0
           gradtarg(2,1,i)=0
           gradtarg(3,1,i)=0
           gradtarg(1,2,i)=0
           gradtarg(2,2,i)=0
           gradtarg(3,2,i)=0
           gradtarg(1,3,i)=0
           gradtarg(2,3,i)=0
           gradtarg(3,3,i)=0
        endif
        enddo
c
c
c   Direct arrival (free space Green's function)
c
        if( itype .eq. 1 ) then
c
        call stfmm3dparttarg(ier,iprec,
     $     nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,ntarget,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
        endif
c
c   Image
c       
        allocate( source_image(3,nsource) )
        do i=1,nsource
        source_image(1,i)=+source(1,i)
        source_image(2,i)=+source(2,i)
        source_image(3,i)=-source(3,i)
        enddo
c
        if( ifsingle .eq. 1 ) then
        allocate( sigma_sl_image(3,nsource) )
        do i=1,nsource
        sigma_sl_image(1,i)=+sigma_sl(1,i)
        sigma_sl_image(2,i)=+sigma_sl(2,i)
        sigma_sl_image(3,i)=-sigma_sl(3,i)
        enddo
        else
        allocate( sigma_sl_image(3,1) )
        endif

        if( ifdouble .ge. 1 ) then
        allocate( sigma_dl_image(3,nsource) )
        allocate( sigma_dv_image(3,nsource) )
        do i=1,nsource
        sigma_dl_image(1,i)=+sigma_dl(1,i)
        sigma_dl_image(2,i)=+sigma_dl(2,i)
        sigma_dl_image(3,i)=-sigma_dl(3,i)
        sigma_dv_image(1,i)=+sigma_dv(1,i)
        sigma_dv_image(2,i)=+sigma_dv(2,i)
        sigma_dv_image(3,i)=-sigma_dv(3,i)
        enddo
        else
        allocate( sigma_dl_image(3,1) )
        allocate( sigma_dv_image(3,1) )
        endif

c
c   Join target and source lists for image processing
c
        allocate( pot0(3, 1) )
        allocate( pre0(1) )
        allocate( grad0(3,3, 1) )
c
        ntarget0 = ntarget+nsource
        allocate( target0(3, ntarget+nsource) )
c
        do j=1,ntarget
        target0(1,j)=target(1,j)
        target0(2,j)=target(2,j)
        target0(3,j)=target(3,j)
        enddo
        do j=1,nsource
        target0(1,j+ntarget)=source(1,j)
        target0(2,j+ntarget)=source(2,j)
        target0(3,j+ntarget)=source(3,j)
        enddo
c
        allocate( pottarg0(3, ntarget+nsource) )
        allocate( pretarg0(ntarget+nsource) )
        allocate( gradtarg0(3,3, ntarget+nsource) )

        ifpot0 = 0
        ifgrad0 = 0
        ifpottarg0 = 0
        ifgradtarg0 = 0
        if( ifpot .eq. 1 .or. ifpottarg .eq. 1 ) ifpottarg0 = 1
        if( ifgrad .eq. 1 .or. ifgradtarg .eq. 1 ) ifgradtarg0 = 1

        call stfmm3dparttarg(ier,iprec,
     $     nsource,source_image,
     $     ifsingle,sigma_sl_image,
     $     ifdouble,sigma_dl_image,sigma_dv_image,
     $     ifpot0,pot0,pre0,ifgrad0,grad0,ntarget0,
     $     target0,ifpottarg0,pottarg0,pretarg0,
     $     ifgradtarg0,gradtarg0)

c
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
        do j=1,nsource
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)-pottarg0(1,j+ntarget)
        pot(2,j)=pot(2,j)-pottarg0(2,j+ntarget)
        pot(3,j)=pot(3,j)-pottarg0(3,j+ntarget)
        pre(j)=pre(j)-pretarg0(j+ntarget)
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)-gradtarg0(1,1,j+ntarget)
        grad(2,1,j)=grad(2,1,j)-gradtarg0(2,1,j+ntarget)
        grad(3,1,j)=grad(3,1,j)-gradtarg0(3,1,j+ntarget)
        grad(1,2,j)=grad(1,2,j)-gradtarg0(1,2,j+ntarget)
        grad(2,2,j)=grad(2,2,j)-gradtarg0(2,2,j+ntarget)
        grad(3,2,j)=grad(3,2,j)-gradtarg0(3,2,j+ntarget)
        grad(1,3,j)=grad(1,3,j)-gradtarg0(1,3,j+ntarget)
        grad(2,3,j)=grad(2,3,j)-gradtarg0(2,3,j+ntarget)
        grad(3,3,j)=grad(3,3,j)-gradtarg0(3,3,j+ntarget)
        endif
        enddo
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c       
        do j=1,ntarget
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)-pottarg0(1,j)
        pottarg(2,j)=pottarg(2,j)-pottarg0(2,j)
        pottarg(3,j)=pottarg(3,j)-pottarg0(3,j)
        pretarg(j)=pretarg(j)-pretarg0(j)
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)-gradtarg0(1,1,j)
        gradtarg(2,1,j)=gradtarg(2,1,j)-gradtarg0(2,1,j)
        gradtarg(3,1,j)=gradtarg(3,1,j)-gradtarg0(3,1,j)
        gradtarg(1,2,j)=gradtarg(1,2,j)-gradtarg0(1,2,j)
        gradtarg(2,2,j)=gradtarg(2,2,j)-gradtarg0(2,2,j)
        gradtarg(3,2,j)=gradtarg(3,2,j)-gradtarg0(3,2,j)
        gradtarg(1,3,j)=gradtarg(1,3,j)-gradtarg0(1,3,j)
        gradtarg(2,3,j)=gradtarg(2,3,j)-gradtarg0(2,3,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)-gradtarg0(3,3,j)
        endif
        enddo

        endif

c
c   Papkovich-Neuber potential
c

        if( ifsingle .eq. 1 ) then

        ifcharge = 1
        allocate( charge(nsource) )
        do i=1,nsource
        charge(i) = sigma_sl_image(3,i)
        enddo

        ifdipole = 1
        allocate( dipstr(nsource) )
        allocate( dipvec(3,nsource) )
        do i=1,nsource
        dipstr(i) = 1
        dipvec(1,i) = sigma_sl_image(1,i) * source(3,i)
        dipvec(2,i) = sigma_sl_image(2,i) * source(3,i)
        dipvec(3,i) = sigma_sl_image(3,i) * source(3,i)
        enddo

        else

        ifcharge = 0
        allocate( charge(nsource) )
        do i=1,nsource
        charge(i) = 0
        enddo

        ifdipole = 0
        allocate( dipstr(nsource) )
        allocate( dipvec(3,nsource) )
        do i=1,nsource
        dipstr(i) = 0
        dipvec(1,i) = 0
        dipvec(2,i) = 0
        dipvec(3,i) = 0
        enddo

        endif


        ifquad = 0
        allocate( quadstr(nsource) )
        allocate( quadvec(6,nsource) )
        do i=1,nsource
        quadstr(i) = 0
        quadvec(1,i) = 0
        quadvec(2,i) = 0
        quadvec(3,i) = 0
        quadvec(4,i) = 0
        quadvec(5,i) = 0
        quadvec(6,i) = 0
        enddo



        if( ifdouble .eq. 1 ) then

        ifdipole = 1
        do i=1,nsource
        dipstr(i) = 1
        dipvec(3,i) = dipvec(3,i) + 2*(
     $     sigma_dl_image(1,i) * sigma_dv_image(1,i)+
     $     sigma_dl_image(2,i) * sigma_dv_image(2,i)+
     $     sigma_dl_image(3,i) * sigma_dv_image(3,i))
        enddo

        endif

        if( ifdouble .eq. 2 .or. ifdouble .eq. 3 
     $     .or. ifdouble .eq. 4 ) then

        ifdipole = 1
        do i=1,nsource
        dipstr(i) = 1
        dipvec(1,i) = dipvec(1,i) + 2*(
     $     -sigma_dl_image(1,i) * sigma_dv_image(3,i)
     $     +sigma_dl_image(3,i) * sigma_dv_image(1,i))
        dipvec(2,i) = dipvec(2,i) + 2*(
     $     -sigma_dl_image(2,i) * sigma_dv_image(3,i)
     $     +sigma_dl_image(3,i) * sigma_dv_image(2,i))
        dipvec(3,i) = dipvec(3,i) + 2*(
     $     -sigma_dl_image(3,i) * sigma_dv_image(3,i)
     $     +sigma_dl_image(3,i) * sigma_dv_image(3,i))
        enddo

        endif

        if( ifdouble .eq. 1 .or. ifdouble .eq. 2 
     $     .or. ifdouble .eq. 4 ) then

        ifquad = 1
        do i=1,nsource
        quadstr(i) = source(3,i)*2
        quadvec(1,i) = sigma_dl_image(1,i)*sigma_dv_image(1,i) 
        quadvec(2,i) = sigma_dl_image(2,i)*sigma_dv_image(2,i) 
        quadvec(3,i) = sigma_dl_image(3,i)*sigma_dv_image(3,i) 
        quadvec(4,i) = sigma_dl_image(1,i)*sigma_dv_image(2,i)+
     $     sigma_dl_image(2,i)*sigma_dv_image(1,i)
        quadvec(5,i) = sigma_dl_image(1,i)*sigma_dv_image(3,i)+
     $     sigma_dl_image(3,i)*sigma_dv_image(1,i)
        quadvec(6,i) = sigma_dl_image(2,i)*sigma_dv_image(3,i)+
     $     sigma_dl_image(3,i)*sigma_dv_image(2,i)
        enddo

        endif



c
c   Join target and source lists for image processing
c
        allocate( cpot0(1) )
        allocate( cfld0(3,1) )
        allocate( chess0(6,1) )
c
        allocate( cpottarg0(ntarget+nsource) )
        allocate( cfldtarg0(3,ntarget+nsource) )
        allocate( chesstarg0(6,ntarget+nsource) )

        ifpot0 = 0
        iffld0 = 0
        ifhess0 = 0
        ifpottarg0 = 1
        iffldtarg0 = 1
        ifhesstarg0 = 0
        if( ifgrad .eq. 1 .or. ifgradtarg .eq. 1 ) ifhesstarg0 = 1
        
        call lfmm3dpartquadtarg(ier,iprec,nsource,
     $     source_image,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifquad,quadstr,quadvec,
     $     ifpot0,cpot0,iffld0,cfld0,ifhess0,chess0,ntarget0,
     $     target0,ifpottarg0,cpottarg0,iffldtarg0,cfldtarg0,
     $     ifhesstarg0,chesstarg0)
c
c
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
        do j=1,nsource
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+cfldtarg0(1,j+ntarget)*source(3,j)
        pot(2,j)=pot(2,j)+cfldtarg0(2,j+ntarget)*source(3,j)
        pot(3,j)=pot(3,j)+cfldtarg0(3,j+ntarget)*source(3,j)
        pot(3,j)=pot(3,j)+cpottarg0(j+ntarget)
        pre(j)=pre(j)+2*cfldtarg0(3,j+ntarget)
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)-chesstarg0(1,j+ntarget)*source(3,j)
        grad(2,1,j)=grad(2,1,j)-chesstarg0(4,j+ntarget)*source(3,j)
        grad(3,1,j)=grad(3,1,j)-chesstarg0(5,j+ntarget)*source(3,j)
        grad(1,2,j)=grad(1,2,j)-chesstarg0(4,j+ntarget)*source(3,j)
        grad(2,2,j)=grad(2,2,j)-chesstarg0(2,j+ntarget)*source(3,j)
        grad(3,2,j)=grad(3,2,j)-chesstarg0(6,j+ntarget)*source(3,j)
        grad(1,3,j)=grad(1,3,j)-chesstarg0(5,j+ntarget)*source(3,j)
        grad(2,3,j)=grad(2,3,j)-chesstarg0(6,j+ntarget)*source(3,j)
        grad(3,3,j)=grad(3,3,j)-chesstarg0(3,j+ntarget)*source(3,j)
        grad(1,3,j)=grad(1,3,j)+cfldtarg0(1,j+ntarget)
        grad(2,3,j)=grad(2,3,j)+cfldtarg0(2,j+ntarget)
        grad(3,3,j)=grad(3,3,j)+cfldtarg0(3,j+ntarget)
        grad(3,1,j)=grad(3,1,j)-cfldtarg0(1,j+ntarget)
        grad(3,2,j)=grad(3,2,j)-cfldtarg0(2,j+ntarget)
        grad(3,3,j)=grad(3,3,j)-cfldtarg0(3,j+ntarget)
        endif
        enddo
c
        endif

c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c
        do j=1,ntarget
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+cfldtarg0(1,j)*target(3,j)
        pottarg(2,j)=pottarg(2,j)+cfldtarg0(2,j)*target(3,j)
        pottarg(3,j)=pottarg(3,j)+cfldtarg0(3,j)*target(3,j)
        pottarg(3,j)=pottarg(3,j)+cpottarg0(j)
        pretarg(j)=pretarg(j)+2*cfldtarg0(3,j)
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)-chesstarg0(1,j)*target(3,j)
        gradtarg(2,1,j)=gradtarg(2,1,j)-chesstarg0(4,j)*target(3,j)
        gradtarg(3,1,j)=gradtarg(3,1,j)-chesstarg0(5,j)*target(3,j)
        gradtarg(1,2,j)=gradtarg(1,2,j)-chesstarg0(4,j)*target(3,j)
        gradtarg(2,2,j)=gradtarg(2,2,j)-chesstarg0(2,j)*target(3,j)
        gradtarg(3,2,j)=gradtarg(3,2,j)-chesstarg0(6,j)*target(3,j)
        gradtarg(1,3,j)=gradtarg(1,3,j)-chesstarg0(5,j)*target(3,j)
        gradtarg(2,3,j)=gradtarg(2,3,j)-chesstarg0(6,j)*target(3,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)-chesstarg0(3,j)*target(3,j)
        gradtarg(1,3,j)=gradtarg(1,3,j)+cfldtarg0(1,j)
        gradtarg(2,3,j)=gradtarg(2,3,j)+cfldtarg0(2,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)+cfldtarg0(3,j)
        gradtarg(3,1,j)=gradtarg(3,1,j)-cfldtarg0(1,j)
        gradtarg(3,2,j)=gradtarg(3,2,j)-cfldtarg0(2,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)-cfldtarg0(3,j)
        endif
        enddo
c
        endif

        return
        end
c
c
c
c
c
        subroutine sth3dpartdirect(itype,
     $     nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,ntarget,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
c
c     Stokes interactions in R^3: evaluate all pairwise particle
c     interactions (excluding self interactions) and interactions with
c     targets using the direct O(N^2) algorithm.
c
c     Half space Green's function.
c     No slip (zero-velocity) boundary condition at z=0
c
c     This subroutine evaluates the Stokes potential and gradient due
c     to a collection of Stokes single and double forces. We use
c
c       ifsingle=1, stokeslet, f = sigma_sl
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
c       u_i = [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
c       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 0
c
c       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
c                   g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
c             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c     for the free-space Green's function, without the (1/4 pi) scaling.  
c
c     Half-space Green's function is the combination of direct arrival, 
c     image contribution and Papkovich-Neuber correction.
c
c
c       INPUT:
c
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       sigma_dv(3,nsource) - orientation of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifpot - velocity/pressure computation flag
c       ifgrad - grad computation flag
c       ifpottarg - target velocity/pressure computation flag
c       ifgradtarg - target grad computation flag
c
c       itype: half space Green's function evaluation flag
c         1 => include both direct arrival and image contribution
c         2 => include image contribution only
c
c       OUTPUT:
c
c       pot(3,nsource) - velocity at source locations
c       pre(nsource) - pressure at source locations
c       grad(3,3,nsource) - grad at source locations
c       pottarg(3,ntarget) - velocity at target locations
c       pretarg(ntarget) - pressure at target locations
c       gradtarg(3,3,ntarget) - grad at target locations
c
        real *8 source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8, allocatable :: source_image(:,:)
        real *8, allocatable :: sigma_sl_image(:,:)
        real *8, allocatable :: sigma_dl_image(:,:)
        real *8, allocatable :: sigma_dv_image(:,:)
c
        real *8, allocatable :: pot0(:,:)
        real *8, allocatable :: pre0(:)
        real *8, allocatable :: grad0(:,:,:)
c
        real *8, allocatable :: target0(:,:)
        real *8, allocatable :: pottarg0(:,:)
        real *8, allocatable :: pretarg0(:)
        real *8, allocatable :: gradtarg0(:,:,:)
c
        complex *16, allocatable :: charge(:)
        complex *16, allocatable :: dipstr(:)
        complex *16, allocatable :: quadstr(:)
        real *8, allocatable :: dipvec(:,:)
        real *8, allocatable :: quadvec(:,:)
c
        complex *16, allocatable :: cpot0(:)
        complex *16, allocatable :: cfld0(:,:)
        complex *16, allocatable :: chess0(:,:)
c
        complex *16, allocatable :: cpottarg0(:)
        complex *16, allocatable :: cfldtarg0(:,:)
        complex *16, allocatable :: chesstarg0(:,:)
c
c
        do i=1,nsource
        if( ifpot .eq. 1) then
           pot(1,i)=0
           pot(2,i)=0
           pot(3,i)=0
           pre(i)=0
        endif
        if( ifgrad .eq. 1) then
           grad(1,1,i)=0
           grad(2,1,i)=0
           grad(3,1,i)=0
           grad(1,2,i)=0
           grad(2,2,i)=0
           grad(3,2,i)=0
           grad(1,3,i)=0
           grad(2,3,i)=0
           grad(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) then
           pottarg(1,i)=0
           pottarg(2,i)=0
           pottarg(3,i)=0
           pretarg(i)=0
        endif
        if( ifgradtarg .eq. 1) then
           gradtarg(1,1,i)=0
           gradtarg(2,1,i)=0
           gradtarg(3,1,i)=0
           gradtarg(1,2,i)=0
           gradtarg(2,2,i)=0
           gradtarg(3,2,i)=0
           gradtarg(1,3,i)=0
           gradtarg(2,3,i)=0
           gradtarg(3,3,i)=0
        endif
        enddo
c
c
c   Direct arrival (free space Green's function)
c
        if( itype .eq. 1 ) then
c
        call st3dpartdirect(
     $     nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,ntarget,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
        endif
c
c   Image
c       
        allocate( source_image(3,nsource) )
        do i=1,nsource
        source_image(1,i)=+source(1,i)
        source_image(2,i)=+source(2,i)
        source_image(3,i)=-source(3,i)
        enddo
c
        if( ifsingle .eq. 1 ) then
        allocate( sigma_sl_image(3,nsource) )
        do i=1,nsource
        sigma_sl_image(1,i)=+sigma_sl(1,i)
        sigma_sl_image(2,i)=+sigma_sl(2,i)
        sigma_sl_image(3,i)=-sigma_sl(3,i)
        enddo
        else
        allocate( sigma_sl_image(3,1) )
        endif

        if( ifdouble .ge. 1 ) then
        allocate( sigma_dl_image(3,nsource) )
        allocate( sigma_dv_image(3,nsource) )
        do i=1,nsource
        sigma_dl_image(1,i)=+sigma_dl(1,i)
        sigma_dl_image(2,i)=+sigma_dl(2,i)
        sigma_dl_image(3,i)=-sigma_dl(3,i)
        sigma_dv_image(1,i)=+sigma_dv(1,i)
        sigma_dv_image(2,i)=+sigma_dv(2,i)
        sigma_dv_image(3,i)=-sigma_dv(3,i)
        enddo
        else
        allocate( sigma_dl_image(3,1) )
        allocate( sigma_dv_image(3,1) )
        endif

c
c   Join target and source lists for image processing
c
        allocate( pot0(3, 1) )
        allocate( pre0(1) )
        allocate( grad0(3,3, 1) )
c
        ntarget0 = ntarget+nsource
        allocate( target0(3, ntarget+nsource) )
c
        do j=1,ntarget
        target0(1,j)=target(1,j)
        target0(2,j)=target(2,j)
        target0(3,j)=target(3,j)
        enddo
        do j=1,nsource
        target0(1,j+ntarget)=source(1,j)
        target0(2,j+ntarget)=source(2,j)
        target0(3,j+ntarget)=source(3,j)
        enddo
c
        allocate( pottarg0(3, ntarget+nsource) )
        allocate( pretarg0(ntarget+nsource) )
        allocate( gradtarg0(3,3, ntarget+nsource) )

        ifpot0 = 0
        ifgrad0 = 0
        ifpottarg0 = 0
        ifgradtarg0 = 0
        if( ifpot .eq. 1 .or. ifpottarg .eq. 1 ) ifpottarg0 = 1
        if( ifgrad .eq. 1 .or. ifgradtarg .eq. 1 ) ifgradtarg0 = 1

        call st3dpartdirect(
     $     nsource,source_image,
     $     ifsingle,sigma_sl_image,
     $     ifdouble,sigma_dl_image,sigma_dv_image,
     $     ifpot0,pot0,pre0,ifgrad0,grad0,ntarget0,
     $     target0,ifpottarg0,pottarg0,pretarg0,
     $     ifgradtarg0,gradtarg0)

c
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
        do j=1,nsource
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)-pottarg0(1,j+ntarget)
        pot(2,j)=pot(2,j)-pottarg0(2,j+ntarget)
        pot(3,j)=pot(3,j)-pottarg0(3,j+ntarget)
        pre(j)=pre(j)-pretarg0(j+ntarget)
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)-gradtarg0(1,1,j+ntarget)
        grad(2,1,j)=grad(2,1,j)-gradtarg0(2,1,j+ntarget)
        grad(3,1,j)=grad(3,1,j)-gradtarg0(3,1,j+ntarget)
        grad(1,2,j)=grad(1,2,j)-gradtarg0(1,2,j+ntarget)
        grad(2,2,j)=grad(2,2,j)-gradtarg0(2,2,j+ntarget)
        grad(3,2,j)=grad(3,2,j)-gradtarg0(3,2,j+ntarget)
        grad(1,3,j)=grad(1,3,j)-gradtarg0(1,3,j+ntarget)
        grad(2,3,j)=grad(2,3,j)-gradtarg0(2,3,j+ntarget)
        grad(3,3,j)=grad(3,3,j)-gradtarg0(3,3,j+ntarget)
        endif
        enddo
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c       
        do j=1,ntarget
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)-pottarg0(1,j)
        pottarg(2,j)=pottarg(2,j)-pottarg0(2,j)
        pottarg(3,j)=pottarg(3,j)-pottarg0(3,j)
        pretarg(j)=pretarg(j)-pretarg0(j)
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)-gradtarg0(1,1,j)
        gradtarg(2,1,j)=gradtarg(2,1,j)-gradtarg0(2,1,j)
        gradtarg(3,1,j)=gradtarg(3,1,j)-gradtarg0(3,1,j)
        gradtarg(1,2,j)=gradtarg(1,2,j)-gradtarg0(1,2,j)
        gradtarg(2,2,j)=gradtarg(2,2,j)-gradtarg0(2,2,j)
        gradtarg(3,2,j)=gradtarg(3,2,j)-gradtarg0(3,2,j)
        gradtarg(1,3,j)=gradtarg(1,3,j)-gradtarg0(1,3,j)
        gradtarg(2,3,j)=gradtarg(2,3,j)-gradtarg0(2,3,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)-gradtarg0(3,3,j)
        endif
        enddo

        endif

c
c   Papkovich-Neuber potential
c

        if( ifsingle .eq. 1 ) then

        ifcharge = 1
        allocate( charge(nsource) )
        do i=1,nsource
        charge(i) = sigma_sl_image(3,i)
        enddo

        ifdipole = 1
        allocate( dipstr(nsource) )
        allocate( dipvec(3,nsource) )
        do i=1,nsource
        dipstr(i) = 1
        dipvec(1,i) = sigma_sl_image(1,i) * source(3,i)
        dipvec(2,i) = sigma_sl_image(2,i) * source(3,i)
        dipvec(3,i) = sigma_sl_image(3,i) * source(3,i)
        enddo

        else

        ifcharge = 0
        allocate( charge(nsource) )
        do i=1,nsource
        charge(i) = 0
        enddo

        ifdipole = 0
        allocate( dipstr(nsource) )
        allocate( dipvec(3,nsource) )
        do i=1,nsource
        dipstr(i) = 0
        dipvec(1,i) = 0
        dipvec(2,i) = 0
        dipvec(3,i) = 0
        enddo

        endif


        ifquad = 0
        allocate( quadstr(nsource) )
        allocate( quadvec(6,nsource) )
        do i=1,nsource
        quadstr(i) = 0
        quadvec(1,i) = 0
        quadvec(2,i) = 0
        quadvec(3,i) = 0
        quadvec(4,i) = 0
        quadvec(5,i) = 0
        quadvec(6,i) = 0
        enddo



        if( ifdouble .eq. 1 ) then

        ifdipole = 1
        do i=1,nsource
        dipstr(i) = 1
        dipvec(3,i) = dipvec(3,i) + 2*(
     $     sigma_dl_image(1,i) * sigma_dv_image(1,i)+
     $     sigma_dl_image(2,i) * sigma_dv_image(2,i)+
     $     sigma_dl_image(3,i) * sigma_dv_image(3,i))
        enddo

        endif

        if( ifdouble .eq. 2 .or. ifdouble .eq. 3 
     $     .or. ifdouble .eq. 4 ) then

        ifdipole = 1
        do i=1,nsource
        dipstr(i) = 1
        dipvec(1,i) = dipvec(1,i) + 2*(
     $     -sigma_dl_image(1,i) * sigma_dv_image(3,i)
     $     +sigma_dl_image(3,i) * sigma_dv_image(1,i))
        dipvec(2,i) = dipvec(2,i) + 2*(
     $     -sigma_dl_image(2,i) * sigma_dv_image(3,i)
     $     +sigma_dl_image(3,i) * sigma_dv_image(2,i))
        dipvec(3,i) = dipvec(3,i) + 2*(
     $     -sigma_dl_image(3,i) * sigma_dv_image(3,i)
     $     +sigma_dl_image(3,i) * sigma_dv_image(3,i))
        enddo

        endif

        if( ifdouble .eq. 1 .or. ifdouble .eq. 2 
     $     .or. ifdouble .eq. 4 ) then

        ifquad = 1
        do i=1,nsource
        quadstr(i) = source(3,i)*2
        quadvec(1,i) = sigma_dl_image(1,i)*sigma_dv_image(1,i) 
        quadvec(2,i) = sigma_dl_image(2,i)*sigma_dv_image(2,i) 
        quadvec(3,i) = sigma_dl_image(3,i)*sigma_dv_image(3,i) 
        quadvec(4,i) = sigma_dl_image(1,i)*sigma_dv_image(2,i)+
     $     sigma_dl_image(2,i)*sigma_dv_image(1,i)
        quadvec(5,i) = sigma_dl_image(1,i)*sigma_dv_image(3,i)+
     $     sigma_dl_image(3,i)*sigma_dv_image(1,i)
        quadvec(6,i) = sigma_dl_image(2,i)*sigma_dv_image(3,i)+
     $     sigma_dl_image(3,i)*sigma_dv_image(2,i)
        enddo

        endif



c
c   Join target and source lists for image processing
c
        allocate( cpot0(1) )
        allocate( cfld0(3,1) )
        allocate( chess0(6,1) )
c
        allocate( cpottarg0(ntarget+nsource) )
        allocate( cfldtarg0(3,ntarget+nsource) )
        allocate( chesstarg0(6,ntarget+nsource) )

        ifpot0 = 0
        iffld0 = 0
        ifhess0 = 0
        ifpottarg0 = 1
        iffldtarg0 = 1
        ifhesstarg0 = 0
        if( ifgrad .eq. 1 .or. ifgradtarg .eq. 1 ) ifhesstarg0 = 1
        
        call l3dpartquaddirect(nsource,
     $     source_image,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifquad,quadstr,quadvec,
     $     ifpot0,cpot0,iffld0,cfld0,ifhess0,chess0,ntarget0,
     $     target0,ifpottarg0,cpottarg0,iffldtarg0,cfldtarg0,
     $     ifhesstarg0,chesstarg0)
c
c
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
        do j=1,nsource
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+cfldtarg0(1,j+ntarget)*source(3,j)
        pot(2,j)=pot(2,j)+cfldtarg0(2,j+ntarget)*source(3,j)
        pot(3,j)=pot(3,j)+cfldtarg0(3,j+ntarget)*source(3,j)
        pot(3,j)=pot(3,j)+cpottarg0(j+ntarget)
        pre(j)=pre(j)+2*cfldtarg0(3,j+ntarget)
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)-chesstarg0(1,j+ntarget)*source(3,j)
        grad(2,1,j)=grad(2,1,j)-chesstarg0(4,j+ntarget)*source(3,j)
        grad(3,1,j)=grad(3,1,j)-chesstarg0(5,j+ntarget)*source(3,j)
        grad(1,2,j)=grad(1,2,j)-chesstarg0(4,j+ntarget)*source(3,j)
        grad(2,2,j)=grad(2,2,j)-chesstarg0(2,j+ntarget)*source(3,j)
        grad(3,2,j)=grad(3,2,j)-chesstarg0(6,j+ntarget)*source(3,j)
        grad(1,3,j)=grad(1,3,j)-chesstarg0(5,j+ntarget)*source(3,j)
        grad(2,3,j)=grad(2,3,j)-chesstarg0(6,j+ntarget)*source(3,j)
        grad(3,3,j)=grad(3,3,j)-chesstarg0(3,j+ntarget)*source(3,j)
        grad(1,3,j)=grad(1,3,j)+cfldtarg0(1,j+ntarget)
        grad(2,3,j)=grad(2,3,j)+cfldtarg0(2,j+ntarget)
        grad(3,3,j)=grad(3,3,j)+cfldtarg0(3,j+ntarget)
        grad(3,1,j)=grad(3,1,j)-cfldtarg0(1,j+ntarget)
        grad(3,2,j)=grad(3,2,j)-cfldtarg0(2,j+ntarget)
        grad(3,3,j)=grad(3,3,j)-cfldtarg0(3,j+ntarget)
        endif
        enddo
c
        endif

c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c
        do j=1,ntarget
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+cfldtarg0(1,j)*target(3,j)
        pottarg(2,j)=pottarg(2,j)+cfldtarg0(2,j)*target(3,j)
        pottarg(3,j)=pottarg(3,j)+cfldtarg0(3,j)*target(3,j)
        pottarg(3,j)=pottarg(3,j)+cpottarg0(j)
        pretarg(j)=pretarg(j)+2*cfldtarg0(3,j)
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)-chesstarg0(1,j)*target(3,j)
        gradtarg(2,1,j)=gradtarg(2,1,j)-chesstarg0(4,j)*target(3,j)
        gradtarg(3,1,j)=gradtarg(3,1,j)-chesstarg0(5,j)*target(3,j)
        gradtarg(1,2,j)=gradtarg(1,2,j)-chesstarg0(4,j)*target(3,j)
        gradtarg(2,2,j)=gradtarg(2,2,j)-chesstarg0(2,j)*target(3,j)
        gradtarg(3,2,j)=gradtarg(3,2,j)-chesstarg0(6,j)*target(3,j)
        gradtarg(1,3,j)=gradtarg(1,3,j)-chesstarg0(5,j)*target(3,j)
        gradtarg(2,3,j)=gradtarg(2,3,j)-chesstarg0(6,j)*target(3,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)-chesstarg0(3,j)*target(3,j)
        gradtarg(1,3,j)=gradtarg(1,3,j)+cfldtarg0(1,j)
        gradtarg(2,3,j)=gradtarg(2,3,j)+cfldtarg0(2,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)+cfldtarg0(3,j)
        gradtarg(3,1,j)=gradtarg(3,1,j)-cfldtarg0(1,j)
        gradtarg(3,2,j)=gradtarg(3,2,j)-cfldtarg0(2,j)
        gradtarg(3,3,j)=gradtarg(3,3,j)-cfldtarg0(3,j)
        endif
        enddo
c
        endif

        return
        end
c
c
c
c
c
