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
c    $Date: 2012-02-19 15:58:32 -0500 (Sun, 19 Feb 2012) $
c    $Revision: 2721 $
c
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Laplace potentials, fields, and
c     hessians due to point charges and dipoles.  (FORTRAN 90 VERSION)
c
c     FAR FIELD CONTRIBUTION ONLY.
c
c     lfmm3dparthesstf - Laplace FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm3dparthesstftarg - Laplace FMM in R^3: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the Laplace particle FMM in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine lfmm3dparthesstf(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
c       Laplace FMM in R^3: this subroutine evaluates the harmonic
c       potential, field, and hessian due to a collection due to a
c       collection of point charges and dipoles (all pairwise particle
c       interactions, ignoring self-interaction).  We use (1/r) for the
c       Green's function, without the (1/4 pi)
c       scaling. Self-interactions are not included.
c   
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
c
        dimension source(3,1)
        complex *16 charge(1)
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1)
        dimension w(1)
c
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
        complex *16 hesstarg(6,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        iffldtarg=0
        ifhesstarg=0
c
        call lfmm3dparthesstftarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dparthesstftarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dparttargmain.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of sources
c       source: real *8 (3,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
c       dipvec: real *8 (3,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pot: complex *16 (nsource): potential at source locations
c       fld: complex *16 (3,nsource): field (-gradient) at source locations
c       hess: complex *16 (6,nsource): hessian at triangle centroids 
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c       hesstarg: complex *16 (6,ntarget): hessian at target locations
c
cf2py   intent(out) ier
cf2py   intent(in) iprec
cf2py   intent(in) nsource
cf2py   intent(in) source
cf2py   intent(in) ifcharge,charge
cf2py   check(!ifcharge || (shape(charge,0) == nsource))  charge
cf2py   depend(nsource)  charge
cf2py   intent(in) ifdipole,dipvec,dipstr
cf2py   check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
cf2py   depend(nsource)  dipstr
cf2py   intent(in) ifpot,iffld,ifhess
cf2py   intent(out) pot,fld,hess
cf2py   intent(in) ifpottarg, iffldtarg, ifhesstarg
cf2py   intent(in) target
cf2py   intent(in) ntarget
cf2py   check((!ifpottarg && !iffldtarg && !ifhesstarg) || (shape(target,0)==3 && shape(target,1) == ntarget))  target
cf2py   check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
cf2py   check((!iffldtarg) || (shape(fldtarg,0)==3 && shape(fldtarg,1) == ntarget))  fldtarg
cf2py   check((!ifhesstarg) || (shape(hesstarg,0)==6 && shape(hesstarg,1) == ntarget))  hesstarg
c
c       (F2PY workaround: *targ must be input because f2py
c       refuses to allocate zero-size output arrays.)
c
cf2py   intent(in,out) pottarg,fldtarg, hesstarg
c
        dimension source(3,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        dimension dipvec(3,nsource)
        complex *16 ima
        complex *16 pot(nsource)
        complex *16 fld(3,nsource)
        complex *16 hess(6,nsource)
        dimension target(3,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 fldtarg(3,ntarget)
        complex *16 hesstarg(6,ntarget)
c
        dimension timeinfo(10)
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 evels of refinment, which is 
c     more than enough for any non-pathological case.
c
 
        dimension laddr(2,200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
        real *8, allocatable :: wrmlexp(:)
        complex *16 ptemp,ftemp(3)
c       
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
c       
        done=1
        pi=4*atan(done)
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
        ifprint=1
c
c     set fmm tolerance based on iprec flag.
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
        if (ifprint .ge. 1) call prin2('epsfmm=*',epsfmm,1)
c
c
c     set criterion for box subdivision (number of sources per box)
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
        if (ifprint .ge. 1) call prinf('nbox=*',nbox,1)
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
        if (ifprint .ge. 1) call prin2('scale=*',scale,nlev+1)
c       
c
c       carve up workspace further
c
c     isourcesort is pointer for sorted source coordinates
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c       
        isourcesort = lused7 + 5
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 3*nsource
          ldipstr = 2*nsource
        else
          ldipvec = 3
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr
c
c
c       ... allocate the potential and field arrays
c
        ipot = lused7
        lpot = 2*nsource
        lused7=lused7+lpot
c       
        ifld = lused7
        if( iffld .eq. 1) then
        lfld = 2*(3*nsource)
        else
        lfld=6
        endif
        lused7=lused7+lfld
c      
        ihess = lused7
        if( ifhess .eq. 1) then
        lhess = 2*(6*nsource)
        else
        lhess=12
        endif
        lused7=lused7+lhess
c      
        ipottarg = lused7
        lpottarg = 2*ntarget
        lused7=lused7+lpottarg
c       
        ifldtarg = lused7
        if( iffldtarg .eq. 1) then
        lfldtarg = 2*(3*ntarget)
        else
        lfldtarg=6
        endif
        lused7=lused7+lfldtarg
c      
        ihesstarg = lused7
        if( ifhesstarg .eq. 1) then
        lhesstarg = 2*(6*ntarget)
        else
        lhesstarg=12
        endif
        lused7=lused7+lhesstarg
c      
        if(ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
c
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c      
        nmax = 0
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call l3dterms(epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        nquad=2*nmax        
c       
        ixnodes = lused7 
        iwts = ixnodes + nquad
        lused7 = iwts + nquad
c
c
        if(ifprint.eq.1) call prinf('nterms=*',nterms,nlev+1)
        if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c   
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7 
        imptemp = iiaddr + 2*nboxes
        lmptemp = (nmax+1)*(2*nmax+1)*2
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                  lused7 is *',lused7,1)
           ier = 8
           return
        endif
c
c     reorder sources, targets so that each box holds
c     contiguous list of source/target numbers.
c
        call l3dreorder(nsource,source,ifcharge,charge,wlists(iisource),
     $     ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 
c       
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        if(ifprint.eq.1) then
          call prinf('finished reordering=*',ier,1)
          call prinf('ier=*',ier,1)
          call prinf('nboxes=*',nboxes,1)
          call prinf('nlev=*',nlev,1)
          call prinf('nboxes=*',nboxes,1)
          call prinf('lused7=*',lused7,1)
        endif
c
        ifinit=1
        call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
c
ccc        call prin2('xnodes=*',xnodes,nquad)
ccc        call prin2('wts=*',wts,nquad)
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call l3dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if (ifprint .ge. 1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                  lused7 is *',lused7,1)
           ier = 16
           return
        endif
c
c       
ccc        do i=lused7+1,lused7+1+100
ccc        w(i)=777
ccc        enddo
c
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far field and local fields (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
        ifevalfar=1
        ifevalloc=0
c
        call lfmm3dparthesstftargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,w(isourcesort),w(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),iffld,w(ifld),ifhess,w(ihess),
     $     ntarget,w(itargetsort),w(iitarget),
     $     ifpottarg,w(ipottarg),iffldtarg,w(ifldtarg),
     $     ifhesstarg,w(ihesstarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     w(ixnodes),w(iwts),nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if(ifprint.eq.1) then
        call prinf('lwlists=*',lwlists,1)
        call prinf('lused total =*',lused7,1)
        call prin2('memory / point = *',(lused7)/dble(nsource),1)
        endif
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
        if(iffld .eq. 1) 
     $     call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
        if(ifhess .eq. 1) 
     $     call l3dhsort(nsource,wlists(iisource),w(ihess),hess)
c
        if(ifpottarg .eq. 1 )
     $     call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
        if(iffldtarg .eq. 1) 
     $     call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
        if(ifhesstarg .eq. 1) 
     $     call l3dhsort(ntarget,wlists(iitarget),w(ihesstarg),hesstarg)
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dparthesstftargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ifhess,hess,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,xnodes,wts,nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension sourcesort(3,1), isource(1)
        complex *16 chargesort(1)
        complex *16 dipstrsort(1)
        dimension dipvecsort(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1)
        dimension targetsort(3,1), itarget(1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        complex *16 hesstarg(6,1)
        dimension wlists(1)
        dimension iaddr(2,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension xnodes(nquad),wts(nquad)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3),htemp(6)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(40 000)
        dimension nterms_eval(4,0:200)
c
        real *8, allocatable :: scarray_local(:)
        real *8, allocatable :: scarray_mpole(:)
c
        data ima/(0.0d0,1.0d0)/

c
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
c     
c       ... set the potential, field, and hessian to zero
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(1,i)=0
           hess(2,i)=0
           hess(3,i)=0
           hess(4,i)=0
           hess(5,i)=0
           hess(6,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(1,i)=0
           hesstarg(2,i)=0
           hesstarg(3,i)=0
           hesstarg(4,i)=0
           hesstarg(5,i)=0
           hesstarg(6,i)=0
        endif
        enddo
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c
c       ... initialize Legendre function evaluation routines
c
        nlege=100
        lw7=40 000
        call ylgndrfwini(nlege,wlege,lw7,lused7)
c
        do i=0,nlev
        do itype=1,4
        call l3dterms_eval(itype,epsfmm,
     1       nterms_eval(itype,i),ier)
        enddo
        enddo
c
        if(ifprint .ge. 2)
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
        do ibox = 1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        enddo
c
c
        if (ifprint .ge. 1) call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
        do 1300 ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        level=box(1)
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(14)
c        npts=box(15)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(15)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(14)),box(15))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(15) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	    radius = (corners0(1,1) - center0(1))**2
	    radius = radius + (corners0(2,1) - center0(2))**2
	    radius = radius + (corners0(3,1) - center0(3))**2
	    radius = sqrt(radius)
c
            call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))

            if( ifcharge .eq. 1 ) then
c
            call l3dformmp_add_trunc(ier,scale(level),
     1         sourcesort(1,box(14)),chargesort(box(14)),
     $         npts,center0,
     $         nterms(level),nterms_eval(1,level),
     2         rmlexp(iaddr(1,ibox)),wlege,nlege)        
c
            endif
c 
            if (ifdipole .eq. 1 ) then

            call l3dformmp_dp_add_trunc(ier,scale(level),
     $         sourcesort(1,box(14)),
     1         dipstrsort(box(14)),dipvecsort(1,box(14)),
     $         npts,center0,nterms(level),nterms_eval(1,level),
     2         rmlexp(iaddr(1,ibox)),wlege,nlege)
            
            endif
         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
        if (ifprint .ge. 1) call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and fields directly
c 
         do 3251 ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(level,npts,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            level1=box1(1)
c
            ifdirect3 = 0
            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
               npts=box(15)
c
               if( ifcharge .eq. 1 ) then
c
               call l3dformta_add_trunc(ier,scale(level1),
     1            sourcesort(1,box(14)),chargesort(box(14)),
     $            npts,center1,
     $            nterms(level1),nterms_eval(1,level1),
     2            rmlexp(iaddr(2,jbox)),wlege,nlege)
c
               endif
c
               if( ifdipole .eq. 1 ) then

               call l3dformta_dp_add_trunc(ier,scale(level1),
     1            sourcesort(1,box(14)),dipstrsort(box(14)),
     2            dipvecsort(1,box(14)),npts,center1,
     3            nterms(level1),nterms_eval(1,level1),
     $            rmlexp(iaddr(2,jbox)),wlege,nlege)

               endif
c
            else

            call lfmm3dparthess_direct(box,box1,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

            endif
         enddo
C$OMP END PARALLEL DO
c
 3251    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c
        if(ifprint .ge. 1)
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0
        if (ifhess.eq.1) ifprune_list2 = 0
        call lfmm3d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,ifprune_list2)

c
        allocate( scarray_mpole(0:100000) )
        call l3dmpevalhessdini(nterms(0),scarray_mpole)
c       
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and fields directly
c
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,level,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box(1)
c
            ifdirect4 = 0
c
            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then
               do j=box1(14),box1(14)+box1(15)-1
                  if( ifhess .eq. 1 ) then
c                  call l3dmpevalhess(scale(level),center0,
c     $               rmlexp(iaddr(1,ibox)),nterms(level),
c     $               sourcesort(1,j),
c     $               ptemp,iffld,ftemp,ifhess,htemp,
c     $               ier)
                  call l3dmpevalhessd_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),
     $               sourcesort(1,j),
     $               ptemp,iffld,ftemp,ifhess,htemp,
     $               scarray_mpole,wlege,nlege)
                  else
                  call l3dmpeval_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),nterms(level),
     $               sourcesort(1,j),
     $               ptemp,iffld,ftemp,
     $               wlege,nlege,ier)
                  endif
                  if( ifpot .eq. 1 ) pot(j)=pot(j)+ptemp
                  if( iffld .eq. 1 ) then
                     fld(1,j)=fld(1,j)+ftemp(1)
                     fld(2,j)=fld(2,j)+ftemp(2)
                     fld(3,j)=fld(3,j)+ftemp(3)
                  endif
                  if (ifhess .eq. 1) then
                  hess(1,j)=hess(1,j)+htemp(1)
                  hess(2,j)=hess(2,j)+htemp(2)
                  hess(3,j)=hess(3,j)+htemp(3)
                  hess(4,j)=hess(4,j)+htemp(4)
                  hess(5,j)=hess(5,j)+htemp(5)
                  hess(6,j)=hess(6,j)+htemp(6)
                  endif
               enddo
               do j=box1(16),box1(16)+box1(17)-1
                  if( ifhesstarg .eq. 1 ) then
c                  call l3dmpevalhess(scale(level),center0,
c     $               rmlexp(iaddr(1,ibox)),nterms(level),
c     $               targetsort(1,j),
c     $               ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
c     $               ier)
                  call l3dmpevalhessd_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),
     $               targetsort(1,j),
     $            ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
     $            scarray_mpole,wlege,nlege)
                  else
                  call l3dmpeval_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),nterms(level),
     $               targetsort(1,j),
     $               ptemp,iffldtarg,ftemp,
     $               wlege,nlege,ier)
                  endif
                  if( ifpottarg .eq. 1 ) pottarg(j)=pottarg(j)+ptemp
                  if( iffldtarg .eq. 1 ) then
                     fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                     fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                     fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
                  endif
                  if (ifhesstarg .eq. 1) then
                  hesstarg(1,j)=hesstarg(1,j)+htemp(1)
                  hesstarg(2,j)=hesstarg(2,j)+htemp(2)
                  hesstarg(3,j)=hesstarg(3,j)+htemp(3)
                  hesstarg(4,j)=hesstarg(4,j)+htemp(4)
                  hesstarg(5,j)=hesstarg(5,j)+htemp(5)
                  hesstarg(6,j)=hesstarg(6,j)+htemp(6)
                  endif
               enddo
            else
            
            call lfmm3dparthess_direct(box,box1,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)

            endif
        enddo
C$OMP END PARALLEL DO
 3252   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c
        allocate( scarray_local(0:100000) )
        call l3dtaevalhessdini(nterms(0),scarray_local)

        if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(15)
c       
        if (level .ge. 2) then
            do j=box(14),box(14)+box(15)-1
               if( ifhess .eq. 1 ) then
c                  call l3dtaevalhess(scale(level),center0,
c     $               rmlexp(iaddr(2,ibox)),nterms(level),
c     $               sourcesort(1,j),
c     $               ptemp,iffld,ftemp,ifhess,htemp,
c     $               ier)
               call l3dtaevalhessd_trunc(scale(level),center0,
     $              rmlexp(iaddr(2,ibox)),nterms(level),sourcesort(1,j),
     $              ptemp,iffld,ftemp,ifhess,htemp,
     $              scarray_local,wlege,nlege)
               else
               call l3dtaeval_trunc(scale(level),center0,
     $            rmlexp(iaddr(2,ibox)),nterms(level),nterms(level),
     $            sourcesort(1,j),
     $            ptemp,iffld,ftemp,
     $            wlege,nlege,ier)
               endif
               if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
               if (iffld .eq. 1) then
                  fld(1,j)=fld(1,j)+ftemp(1)
                  fld(2,j)=fld(2,j)+ftemp(2)
                  fld(3,j)=fld(3,j)+ftemp(3)
               endif
               if (ifhess .eq. 1) then
               hess(1,j)=hess(1,j)+htemp(1)
               hess(2,j)=hess(2,j)+htemp(2)
               hess(3,j)=hess(3,j)+htemp(3)
               hess(4,j)=hess(4,j)+htemp(4)
               hess(5,j)=hess(5,j)+htemp(5)
               hess(6,j)=hess(6,j)+htemp(6)
               endif
            enddo
c
            do j=box(16),box(16)+box(17)-1
               if( ifhesstarg .eq. 1 ) then
c                  call l3dtaevalhess(scale(level),center0,
c     $               rmlexp(iaddr(2,ibox)),nterms(level),
c     $               targetsort(1,j),
c     $               ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
c     $               ier)
               call l3dtaevalhessd_trunc(scale(level),center0,
     $              rmlexp(iaddr(2,ibox)),nterms(level),
     $              targetsort(1,j),
     $              ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
     $              scarray_local,wlege,nlege)
               else
               call l3dtaeval_trunc(scale(level),center0,
     $            rmlexp(iaddr(2,ibox)),nterms(level),nterms(level),
     $            targetsort(1,j),
     $            ptemp,iffldtarg,ftemp,
     $            wlege,nlege,ier)
               endif
               if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
               if (iffldtarg .eq. 1) then
                  fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                  fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                  fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
               endif
               if (ifhesstarg .eq. 1) then
               hesstarg(1,j)=hesstarg(1,j)+htemp(1)
               hesstarg(2,j)=hesstarg(2,j)+htemp(2)
               hesstarg(3,j)=hesstarg(3,j)+htemp(3)
               hesstarg(4,j)=hesstarg(4,j)+htemp(4)
               hesstarg(5,j)=hesstarg(5,j)+htemp(5)
               hesstarg(6,j)=hesstarg(6,j)+htemp(6)
               endif
            enddo
        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
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
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        call lfmm3dparthess_direct_self(box,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
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
            call lfmm3dparthess_direct(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
c
 6203        continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if(ifprint .ge. 1) then
        call prin2('timeinfo=*',timeinfo,8)       
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        endif
c       
        return
        end
c
c
c
c
c
