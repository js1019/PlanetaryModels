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
c    $Date: 2012-03-31 11:17:40 -0400 (Sat, 31 Mar 2012) $
c    $Revision: 2875 $
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains the FMM routines for Stokes particle
c        potentials in free space in R^3. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      stfmm3dpartself - Stokes FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c      stfmm3dparttarg - Stokes FMM in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c      st3dpartdirect - Stokes interactions in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine stfmm3dpartself
     $     (ier,iprec,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad)
c
c     FMM calculation subroutine for Stokes N-body problem.
c
c     \delta u = \grad p, div u = 0, mu = 1.
c
c     Free space Stokes Green's functions:
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
c     without the (1/4 pi) scaling.
c
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
        call stfmm3dparttarg
     $     (ier,iprec,nparts,source,
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
c*********************************
      subroutine stfmm3dparttarg
     $     (ier,iprec,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
c     FMM calculation subroutine for Stokes N-body problem.
c
c     Free space Stokes Green's function.
c
c     INPUT:
c
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     pre(nparts) = pressure at source locations
c     grad(3,3,nparts) = grad at source locations
c     pottarg(3,ntargs) = velocity at target locations
c     pretarg(ntargs) = pressure at target locations
c     gradtarg(3,3,ntargs) = grad at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 pot(3,nparts),pre(nparts),grad(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 pottarg(3,ntargs),pretarg(ntargs),
     $     gradtarg(3,3,ntargs)
        integer nparts,ntargs
        real *8, allocatable :: w(:)
c       
        ier=0
        lused=0
c
c       ... allocate work arrays
c
        icharge=lused+1
        lcharge=2*nparts
        lused=lused+lcharge

        idipstr=lused+1
        ldipstr=2*nparts
        lused=lused+ldipstr

        idipvec=lused+1
        ldipvec=3*nparts
        lused=lused+ldipvec

        icpot=lused+1
        lcpot=2*nparts
        lused=lused+lcpot

        icfld=lused+1
        lcfld=2*3*nparts
        lused=lused+lcfld

        ichess=lused+1
        lchess=2*6*nparts
        lused=lused+lchess

        ihessmatr=lused+1
        lhessmatr=3*3*nparts
        lused=lused+lhessmatr

        icpottarg=lused+1
        lcpottarg=2*ntargs
        lused=lused+lcpottarg

        icfldtarg=lused+1
        lcfldtarg=2*3*ntargs
        lused=lused+lcfldtarg

        ichesstarg=lused+1
        lchesstarg=2*6*ntargs
        lused=lused+lchesstarg

        ihessmatrtarg=lused+1
        lhessmatrtarg=3*3*ntargs
        lused=lused+lhessmatrtarg
c
        allocate(w(lused+10),stat=ier)
        if( ier .ne. 0 ) return
c
        call stfmm3dparttargmain_fast
     $     (ier,iprec,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg,
     $     w(icharge),w(idipstr),w(idipvec),
     $     w(icpot),w(icfld),w(ichess),w(ihessmatr),
     $     w(icpottarg),w(icfldtarg),w(ichesstarg),w(ihessmatrtarg))
c
c     reconstruct FMM data structure and account for all local
c     interactions using quadrature routines - no interactions are saved
c     in the present version.
c     
        call stfmm3dparttarg0
     $     (ier,iprec,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
        if( ier .ne. 0 ) return
c
        return
        end
c
c
c
c
c
c*********************************
      subroutine stfmm3dparttargmain_fast
     $     (ier,iprec,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg,
     $     charge,dipstr,dipvec,cpot,cfld,chess,hessmatr,
     $     cpottarg,cfldtarg,chesstarg,hessmatrtarg)
c
c     FMM calculation subroutine for Stokes N-body problem
c
c     4 Laplace FMM calls
c
c     3 Laplace FMM calls for rotlet
c     7 Laplace FMM calls for doublet
c
c
c     INPUT:
c
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     pre(nparts) = pressure at source locations
c     grad(3,3,nparts) = grad at source locations
c     pottarg(3,ntargs) = velocity at target locations
c     pretarg(ntargs) = pressure at target locations
c     gradtarg(3,3,ntargs) = grad at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 pot(3,nparts),pre(nparts),grad(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 pottarg(3,ntargs),pretarg(ntargs),
     $     gradtarg(3,3,ntargs)
        integer nparts,ntargs
c       
        complex *16 charge(1)
        complex *16 dipstr(1)
        real *8 dipvec(3,1)
c
        complex *16 cpot(1)
        complex *16 cfld(3,1)
        complex *16 chess(6,1) 
        real *8 hessmatr(3,3,1)

        complex *16 cpottarg(1)
        complex *16 cfldtarg(3,1)
        complex *16 chesstarg(6,1) 
        real *8 hessmatrtarg(3,3,1)
c
c
        do k=1,nparts
c
        if( ifpot .eq. 1 ) then
        pot(1,k) = 0.0d0
        pot(2,k) = 0.0d0
        pot(3,k) = 0.0d0
        pre(k) = 0.0d0
        endif
c
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        hessmatr(i,j,k) = 0.0d0
        enddo
        enddo
        endif        
c
        enddo
c
c
        do k=1,ntargs
c
        if( ifpottarg .eq. 1 ) then
        pottarg(1,k) = 0.0d0
        pottarg(2,k) = 0.0d0
        pottarg(3,k) = 0.0d0
        pretarg(k) = 0.0d0
        endif
c
        if( ifgradtarg .eq. 1 ) then
        do i=1,3
        do j=1,3
        hessmatrtarg(i,j,k) = 0.0d0
        enddo
        enddo
        endif
c
        enddo
c
c
c
        ifpot0=0
        iffld0=0
        ifhess0=0
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
        ifpot0=1
        iffld0=1
        endif
        if( ifgrad .eq. 1 ) ifhess0=1

        ifpottarg0=0
        iffldtarg0=0
        ifhesstarg0=0
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
        ifpottarg0=1
        iffldtarg0=1
        endif
        if( ifgradtarg .eq. 1 ) ifhesstarg0=1
c
c
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
        do j = 1,3

        ifcharge=0
        ifdipole=0

        do k = 1,nparts
            charge(k) = 0
            dipstr(k) = 0
            dipvec(1,k) = 0
            dipvec(2,k) = 0
            dipvec(3,k) = 0
            if( ifsingle .eq. 1 ) then
            charge(k) = sigma_sl(j,k)/2
            ifcharge=1
            endif
            if( ifdouble .eq. 1 .or. ifdouble .eq. 2 
     $         .or. ifdouble .eq. 4) then
            dipstr(k) = 1
            dipvec(1,k) = sigma_dv(1,k)*sigma_dl(j,k)
            dipvec(2,k) = sigma_dv(2,k)*sigma_dl(j,k)
            dipvec(3,k) = sigma_dv(3,k)*sigma_dl(j,k)
            dipvec(1,k) = dipvec(1,k)+sigma_dl(1,k)*sigma_dv(j,k)
            dipvec(2,k) = dipvec(2,k)+sigma_dl(2,k)*sigma_dv(j,k)
            dipvec(3,k) = dipvec(3,k)+sigma_dl(3,k)*sigma_dv(j,k)
            dipvec(1,k) = dipvec(1,k)/2
            dipvec(2,k) = dipvec(2,k)/2
            dipvec(3,k) = dipvec(3,k)/2
            ifdipole=1
            endif
        enddo

        if( ifcharge .ne. 0 .or. ifdipole .ne. 0 ) then

        call lfmm3dparthesstftarg(ier,iprec,
     $     nparts,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot0,cpot,iffld0,cfld,ifhess0,chess,
     $     ntargs,target,ifpottarg0,cpottarg,iffldtarg0,cfldtarg,
     $     ifhesstarg0,chesstarg)

        call stfmm3dlap1(nparts,j,cpot,cfld,chess,
     $     source,ifpot,pot,pre,ifgrad,grad)
        call stfmm3dlap1(ntargs,j,cpottarg,cfldtarg,chesstarg,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)

        endif

        enddo

c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
        ifcharge=0
        ifdipole=0
c        
        do k = 1,nparts
          charge(k) = 0
          dipstr(k) = 0
          dipvec(1,k) = 0
          dipvec(2,k) = 0
          dipvec(3,k) = 0
          if( ifsingle .eq. 1 ) then
          charge(k) = 
     $      (sigma_sl(1,k)*source(1,k)+
     $       sigma_sl(2,k)*source(2,k)+
     $       sigma_sl(3,k)*source(3,k))/2
          ifcharge = 1
          endif
          if( ifdouble .eq. 2 .or. ifdouble .eq. 4 ) then
          charge(k) = charge(k) + 
     $        (sigma_dl(1,k)*sigma_dv(1,k) + 
     1         sigma_dl(2,k)*sigma_dv(2,k) + 
     2         sigma_dl(3,k)*sigma_dv(3,k))
          ifcharge = 1
          endif
          if( ifdouble .eq. 1 .or. ifdouble .eq. 2 
     $       .or. ifdouble .eq. 4 ) then
          dipstr(k) = 1
          dipvec(1,k) = sigma_dv(1,k)*
     $        (sigma_dl(1,k)*source(1,k) + 
     1         sigma_dl(2,k)*source(2,k) + 
     2         sigma_dl(3,k)*source(3,k) )
          dipvec(2,k) = sigma_dv(2,k)*
     $        (sigma_dl(1,k)*source(1,k) + 
     1         sigma_dl(2,k)*source(2,k) + 
     2         sigma_dl(3,k)*source(3,k) )
          dipvec(3,k) = sigma_dv(3,k)*
     $        (sigma_dl(1,k)*source(1,k) + 
     1         sigma_dl(2,k)*source(2,k) + 
     2         sigma_dl(3,k)*source(3,k) )
          dipvec(1,k) = dipvec(1,k) + sigma_dl(1,k)*
     $        (sigma_dv(1,k)*source(1,k) + 
     1         sigma_dv(2,k)*source(2,k) + 
     2         sigma_dv(3,k)*source(3,k))
          dipvec(2,k) = dipvec(2,k) + sigma_dl(2,k)*
     $        (sigma_dv(1,k)*source(1,k) + 
     1         sigma_dv(2,k)*source(2,k) + 
     2         sigma_dv(3,k)*source(3,k))
          dipvec(3,k) = dipvec(3,k) + sigma_dl(3,k)*
     $        (sigma_dv(1,k)*source(1,k) + 
     1         sigma_dv(2,k)*source(2,k) + 
     2         sigma_dv(3,k)*source(3,k))
          dipvec(1,k) = dipvec(1,k)/2
          dipvec(2,k) = dipvec(2,k)/2
          dipvec(3,k) = dipvec(3,k)/2
          ifdipole = 1
          endif
        enddo

        if( ifcharge .ne. 0 .or. ifdipole .ne. 0 ) then

        call lfmm3dparthesstftarg(ier,iprec,
     $     nparts,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot0,cpot,iffld0,cfld,ifhess0,chess,
     $     ntargs,target,ifpottarg0,cpottarg,iffldtarg0,cfldtarg,
     $     ifhesstarg0,chesstarg)

        call stfmm3dlap2(nparts,cpot,cfld,chess,
     $     ifpot,pot,ifgrad,grad)
        call stfmm3dlap2(ntargs,cpottarg,cfldtarg,chesstarg,
     $     ifpottarg,pottarg,ifgradtarg,gradtarg)

        endif


        if( ifdouble .eq. 3 .or. ifdouble .eq. 4 ) then
c
c       ... rotlet part
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
        do j = 1,3
c
        ifcharge=0
        ifdipole=0
c        
        do k = 1,nparts
          charge(k) = 0
          dipstr(k) = 0
          dipvec(1,k) = 0
          dipvec(2,k) = 0
          dipvec(3,k) = 0
          if( ifdouble .eq. 3 .or. ifdouble .eq. 4 ) then
          dipstr(k) = 1
          dipvec(1,k) = sigma_dv(1,k)*sigma_dl(j,k)
          dipvec(2,k) = sigma_dv(2,k)*sigma_dl(j,k)
          dipvec(3,k) = sigma_dv(3,k)*sigma_dl(j,k)
          dipvec(1,k) = dipvec(1,k)-sigma_dl(1,k)*sigma_dv(j,k)
          dipvec(2,k) = dipvec(2,k)-sigma_dl(2,k)*sigma_dv(j,k)
          dipvec(3,k) = dipvec(3,k)-sigma_dl(3,k)*sigma_dv(j,k)
          ifdipole = 1
          endif
        enddo

        if( ifcharge .ne. 0 .or. ifdipole .ne. 0 ) then

        ifhess0=0
        ifhesstarg0=0
        call lfmm3dparthesstftarg(ier,iprec,
     $     nparts,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot0,cpot,iffld0,cfld,ifhess0,chess,
     $     ntargs,target,ifpottarg0,cpottarg,iffldtarg0,cfldtarg,
     $     ifhesstarg0,chesstarg)

        call stfmm3dlap3(nparts,j,cpot,cfld,
     $     ifpot,pot,ifgrad,grad)
        call stfmm3dlap3(ntargs,j,cpottarg,cfldtarg,
     $     ifpottarg,pottarg,ifgradtarg,gradtarg)

        endif

        enddo

        endif

        return
        end
c
c
c
c
c
        subroutine stfmm3dparttarg0(ier,iprec,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntarget,target,
     $     ifpottarg,pottarg,pretarg,ifgradtarg,gradtarg)
c
c     FMM calculation subroutine for Stokes N-body problem
c
c     Direct evaluation routine, for local interactions only
c
c     INPUT:
c
c     nsource = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     grad(3,3,nparts) = grad at source locations
c     pottarg(3,ntargs) = velocity at target locations
c     gradtarg(3,3,ntargs) = grad at target locations
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,1)
        real *8 sigma_sl(3,1)
        real *8 sigma_dl(3,1)
        real *8 sigma_dv(3,1)
        real *8 pot(3,1)
        real *8 pre(1)
        real *8 grad(3,3,1)
        real *8 target(3,1)
        real *8 pottarg(3,1)
        real *8 pretarg(1)
        real *8 gradtarg(3,3,1)
c
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
c
        real *8 timeinfo(10)
c       
        real *8 center(3)
c       
        integer laddr(2,200)
        real *8 scale(0:200)
        real *8 bsize(0:200)
        integer nterms(0:200)
c       
        integer box(20)
        real *8 center0(3),corners0(3,8)
c       
        integer box1(20)
        real *8 center1(3),corners1(3,8)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
c
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c       ... build the oct-tree
c       
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c       
        call prin2('epsfmm=*',epsfmm,1)
c
        if( iprec .eq. -2 ) nbox=40
        if( iprec .eq. -1 ) nbox=50
        if( iprec .eq. 0 ) nbox=80
        if( iprec .eq. 1 ) nbox=160
        if( iprec .eq. 2 ) nbox=400
        if( iprec .eq. 3 ) nbox=800
        if( iprec .eq. 4 ) nbox=1200
        if( iprec .eq. 5 ) nbox=1400
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        call prinf('nbox=*',nbox,1)
c
c
c     create oct-tree data structure
c
        ntot = 100*(nsource+ntarget)+10000
        do ii = 1,10
           allocate (wlists(ntot))
           call lfmm3dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
              goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return          
        endif
c
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
c       
        call prin2('scale=*',scale,nlev+1)
c       
c       
c       carve up workspace further
c
c     isourcesort is pointer for sorted source locations
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
c
        isourcesort = lused7
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget

        isigma_slsort = itargetsort+ltargetsort
        if (ifsingle.eq.1) then
          lsigma_slsort = 3*nsource
        else
          lsigma_slsort = 3
        endif
        isigma_dlsort = isigma_slsort+lsigma_slsort
        if (ifdouble.ge.1) then
          lsigma_dlsort = 3*nsource
        else
          lsigma_dlsort = 3
        endif
        isigma_dvsort = isigma_dlsort+lsigma_dlsort
        if (ifdouble.ge.1) then
          lsigma_dvsort = 3*nsource
        else
          lsigma_dvsort = 3
        endif

        lused7 = isigma_dvsort+lsigma_dvsort
c
c
c       ... allocate the potential and field arrays
c
c
        ipot = lused7 
        if( ifpot .eq. 1) then
        lpot = 2*(3*nsource)
        else
        lpot=6
        endif
        lused7=lused7+lpot
c      
        ipre = lused7 
        if( ifpot .eq. 1) then
        lpre = 2*(nsource)
        else
        lpre=2
        endif
        lused7=lused7+lpre
c      
        igrad = lused7
        if( ifgrad .eq. 1) then
        lgrad = 2*(3*3*nsource)
        else
        lgrad= 2*3*3
        endif
        lused7=lused7+lgrad
c      
        ipottarg = lused7
        if( ifpottarg .eq. 1) then
        lpottarg = 2*(3*ntarget)
        else
        lpottarg=6
        endif
        lused7=lused7+lpottarg
c      
        ipretarg = lused7
        if( ifpottarg .eq. 1) then
        lpretarg = 2*(ntarget)
        else
        lpretarg=2
        endif
        lused7=lused7+lpretarg
c      
        igradtarg = lused7
        if( ifgradtarg .eq. 1) then
        lgradtarg = 2*(3*3*ntarget)
        else
        lgradtarg= 2*3*3
        endif
        lused7=lused7+lgradtarg
c      
c
        call prinf(' lused7 is *',lused7,1)
c
c   
c       ... allocate temporary arrays
c
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                   lused7 is *',lused7,1)
           ier = 8
           return          
        endif
c

        call l3dreordertarg
     $     (nsource,source,wlists(iisource),w(isourcesort))
        if( ifsingle .eq. 1 ) then
        call l3dreordertarg
     $     (nsource,sigma_sl,wlists(iisource),w(isigma_slsort))
        endif
        if( ifdouble .ge. 1 ) then
        call l3dreordertarg
     $     (nsource,sigma_dl,wlists(iisource),w(isigma_dlsort))
        call l3dreordertarg
     $     (nsource,sigma_dv,wlists(iisource),w(isigma_dvsort))
        endif
        
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
c
c
c
        call stfmm3dparttarg0_evalloc(ier,iprec,
     $     nsource,w(isourcesort),
     $     ifsingle,w(isigma_slsort),
     $     ifdouble,w(isigma_dlsort),w(isigma_dvsort),
     $     ifpot,w(ipot),w(ipre),ifgrad,w(igrad),
     $     ntarget,w(itargetsort),
     $     ifpottarg,w(ipottarg),w(ipretarg),
     $     ifgradtarg,w(igradtarg),
     $     nboxes,laddr,nlev,wlists(iwlists),lwlists)
c
        call prinf('lwlists=*',lwlists,1)
        call prinf('lused total =*',lused7,1)
c       
        call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call st3dptsort(nsource,wlists(iisource),w(ipot),pot)
        if(ifpot .eq. 1) 
     $     call st3dftsort(nsource,wlists(iisource),w(ipre),pre)
        if(ifgrad .eq. 1) 
     $     call st3dstsort(nsource,wlists(iisource),w(igrad),grad)
c
        if(ifpottarg .eq. 1 )
     $     call st3dptsort(ntarget,wlists(iitarget),
     $     w(ipottarg),pottarg)
        if(ifpottarg .eq. 1 )
     $     call st3dftsort(ntarget,wlists(iitarget),
     $     w(ipretarg),pretarg)
        if(ifgradtarg .eq. 1) 
     $     call st3dstsort(ntarget,wlists(iitarget),
     $     w(igradtarg),gradtarg)
c       
        return
        end
c
c
c
c
c
        subroutine st3dptsort(n,isource,potsort,pot)
        implicit real *8 (a-h,o-z)
        integer isource(1)
        real *8 pot(3,1),potsort(3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        do m=1,3
ccc        pot(m,isource(i))=potsort(m,i)
        pot(m,isource(i))=pot(m,isource(i))+potsort(m,i)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine st3dftsort(n,isource,presort,pre)
        implicit real *8 (a-h,o-z)
        integer isource(1)
        real *8 pre(1),presort(1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
ccc        pre(isource(i))=presort(i)
        pre(isource(i))=pre(isource(i))+presort(i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine st3dstsort(n,isource,gradsort,grad)
        implicit real *8 (a-h,o-z)
        integer isource(1)
        real *8 grad(3,3,1),gradsort(3,3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        do j=1,3
        do m=1,3
ccc        grad(m,j,isource(i))=gradsort(m,j,i)
        grad(m,j,isource(i))=
     $     grad(m,j,isource(i))+gradsort(m,j,i)
        enddo
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine stfmm3dparttarg0_evalloc(ier,iprec,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntarget,target,
     $     ifpottarg,pottarg,pretarg,ifgradtarg,gradtarg,
     $     nboxes,laddr,nlev,wlists,lwlists)
        implicit real *8 (a-h,o-z)
        real *8 source(3,1)
        real *8 sigma_sl(3,1)
        real *8 sigma_dl(3,1)
        real *8 sigma_dv(3,1)
        real *8 pot(3,1)
        real *8 pre(1)
        real *8 grad(3,3,1)
        real *8 target(3,1)
        real *8 pottarg(3,1)
        real *8 pretarg(1)
        real *8 gradtarg(3,3,1)
c
        integer laddr(2,200)
        integer list(10 000)
c
        real *8 timeinfo(10)
c
        real *8 wlists(1)
c
        integer box(20)
        real *8 center0(3),corners0(3,8)
        integer box1(20)
        real *8 center1(3),corners1(3,8)
c
ccc        save
c
c     
c       ... set the velocity and grad to zero
c
        do i=1,nsource
        if (ifpot .eq. 1) then
        pot(1,i)=0
        pot(2,i)=0
        pot(3,i)=0
        pre(i)=0
        endif
        if (ifgrad .eq. 1) then
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
        if (ifpottarg .eq. 1) then
        pottarg(1,i)=0
        pottarg(2,i)=0
        pottarg(3,i)=0
        pretarg(i)=0
        endif
        if (ifgradtarg .eq. 1) then
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
        do i=1,10
        timeinfo(i)=0
        enddo
c
        call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        ifprint=0
        if (ifprint .eq. 1) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .eq. 1) then
               call prinf('npts=*',npts,1)
            endif
        endif
c
c
        if (nkids .eq. 0 ) then
c
c       ... evaluate self interactions
c
        call stfmm3dpart_direct_self(box,
     $     source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .eq. 1) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c    
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c
               call stfmm3dpart_direct(box1,box,
     $            source,
     $            ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $            ifpot,pot,pre,ifgrad,grad,
     $            target,ifpottarg,pottarg,pretarg,
     $            ifgradtarg,gradtarg)
c
 6203           continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
c
ccc        call prin2('inside fmm, pot=*',pot,3*nsource)
ccc        call prin2('inside fmm, pottarg=*',pottarg,3*ntarget)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        call prin2('timeinfo=*',timeinfo,8)
c       
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
c       
        return
        end
c
c
c
c
c
        subroutine stfmm3dpart_direct(box,box1,
     $     source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        real *8 source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8 pot0(3),grad0(3,3)
c
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(14),box1(14)+box1(15)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifsingle .eq. 1 ) then
        call green3sup_eval(source(1,i),
     $     sigma_sl(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call green3stp_arb_eval(ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(16),box1(16)+box1(17)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifsingle .eq. 1 ) then
        call green3sup_eval(source(1,i),
     $     sigma_sl(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call green3stp_arb_eval(ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine stfmm3dpart_direct_self(box,
     $     source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        real *8 source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8 pot0(3),grad0(3,3)
c
c       ... sources
c
        ione=1
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(14),box(14)+box(15)-1
        do i=box(14),box(14)+box(15)-1
c
        if( i .eq. j ) cycle
c
        if (ifsingle .eq. 1 ) then
        call green3sup_eval(source(1,i),
     $     sigma_sl(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call green3stp_arb_eval(ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c       
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(16),box(16)+box(17)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifsingle .eq. 1 ) then
        call green3sup_eval(source(1,i),
     $     sigma_sl(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call green3stp_arb_eval(ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine st3dpartdirect(
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
c     \delta u = \grad p, div u = 0, mu = 1.
c
c     Free space Stokes Green's functions:
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
c     without the (1/4 pi) scaling.
c
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
        real *8 pot0(3),grad0(3,3)
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
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,pot0,pre0,grad0)
        do 6550 j=1,nsource
        do 6540 i=1,nsource
        if( i .eq. j ) goto 6540

        if (ifsingle .eq. 1 ) then
        call green3sup_eval(source(1,i),
     $     sigma_sl(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif

        if (ifdouble .ge. 1) then
        call green3stp_arb_eval(ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
c
 6540   continue
 6550   continue
c
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,pot0,pre0,grad0)
        do j=1,ntarget
        do i=1,nsource
        if (ifsingle .eq. 1 ) then
        call green3sup_eval(source(1,i),
     $     sigma_sl(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
c
        if (ifdouble .ge. 1) then
        call green3stp_arb_eval(ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
C$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
