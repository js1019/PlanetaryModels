c
c     testing code for FMM - tests charges and dipoles against
c     O(N^2) direct method 
c
c
        implicit real *8 (a-h,o-z)
        parameter(lw=20 000 000)
        dimension w(lw)
c       
        call test_part(w,lw)
c
        stop
        end
c
c
c
c
c
        subroutine test_part(w,lw)
        implicit real *8 (a-h,o-z)
        dimension source(3,1 000 000)
        complex *16 charge(1 000 000)
        complex *16 dipstr(1 000 000)
        dimension dipvec(3,1 000 000)
        complex *16 quadstr(1 000 000)
        dimension quadvec(6,1 000 000)
        complex *16 pot(1 000 000)
        complex *16 fld(3,1 000 000)
        complex *16 hess(6,1 000 000)
c       
        complex *16 pot2(1 000 000)
        complex *16 fld2(3,1 000 000)
        complex *16 hess2(6,1 000 000)
c       
        dimension target(3,2 000 000)
        complex *16 pottarg(2 000 000)
        complex *16 fldtarg(3,2 000 000)
        complex *16 hesstarg(6,2 000 000)
c
        complex *16 ptemp,ftemp(3),htemp(6)
c       
ccc        parameter(lw=120 000 000)
        dimension w(1)
        save
c       
        complex *16 ima
        complex *16 zk
        data ima/(0.0d0,1.0d0)/
c
c
        done=1
        pi=4*atan(done)
c
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
ccc        call prini(0,13)
c
        print *, 'ENTER n'
        read *, nsource
c
c
        call prinf('nsource=*',nsource,1)
c
        zk = 1.2d0 + ima*0.1d0
c
c
        do i=1,lw
        w(i)=0
        enddo
c
        idist=3
c
        if( idist .eq. 1 ) then
c
c       ... construct randomly located charge distribution on a unit cube
c       
        do i=1,nsource
        source(1,i)=hkrand(0)
        source(2,i)=hkrand(0)
        source(3,i)=hkrand(0)
        source(1,i)=source(1,i)-0.5
        source(2,i)=source(2,i)-0.5
        source(3,i)=source(3,i)-0.5
        enddo
c
        endif
c
        if( idist .eq. 2 ) then
c
c       ... construct charge distribution on a curve in R^3
c       
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        source(1,i)=sin(1.1*a)
        source(2,i)=cos(2.1*a)
        source(3,i)=cos(3.1*a)
        enddo
c
        endif
c       
        if( idist .eq. 3 ) then
c
c       ... construct randomly located charge distribution on a unit sphere
c 
        done=1
        pi=4*atan(done)
c
        d=hkrand(0)
        do i=1,nsource
c
c        source(1,i)=hkrand(0)
c        source(2,i)=hkrand(0)
c        source(3,i)=hkrand(0)
c        source(1,i)=source(1,i)-0.5
c        source(2,i)=source(2,i)-0.5
c        source(3,i)=source(3,i)-0.5
c        rr=source(1,i)**2+source(2,i)**2+source(3,i)**2
c        rr=sqrt(rr)
c        source(1,i)=source(1,i)/rr
c        source(2,i)=source(2,i)/rr
c        source(3,i)=source(3,i)/rr
c
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        source(1,i)=.5d0*cos(phi)*sin(theta)
        source(2,i)=.5d0*sin(phi)*sin(theta)
        source(3,i)=.5d0*cos(theta)
c
        d=hkrand(0)
        d=hkrand(0)
        enddo
c
        endif
c
        if( idist .eq. 4 ) then
c
c       ... construct the grid of charges
c
        ngrid=nsource-1
        kk=0
        do i=1,ngrid+1
        do j=1,ngrid+1
        do k=1,ngrid+1
        kk=kk+1
        source(1,kk)=(i-1.0)/ngrid
        source(2,kk)=(j-1.0)/ngrid
        source(3,kk)=(k-1.0)/ngrid
        source(1,kk)=source(1,kk)+hkrand(0)*.0001
        source(2,kk)=source(2,kk)+hkrand(0)*.0001
        source(3,kk)=source(3,kk)+hkrand(0)*.0001
ccc        call prin2('source=*',source(1,kk),3)
        enddo
        enddo
        enddo
        nsource=kk
c       
        call prinf('after grid build, nsource=*',nsource,1)
c
        endif
c
c        do i=1,nsource
c        source(1,i)=source(1,i)*10
c        source(2,i)=source(2,i)*10
c        source(3,i)=source(3,i)*10
c        enddo


c       ... set up the targets
c
        if( idist .eq. 1 .or. idist .eq. 4 ) then
        do i=1,nsource
        target(1,i)=source(1,i) + 0.1
        target(2,i)=source(2,i)
        target(3,i)=source(3,i)
        enddo
        ntarget=nsource
        do i=1,nsource
        target(1,i+nsource)=source(1,i) 
        target(2,i+nsource)=source(2,i) - 0.2
        target(3,i+nsource)=source(3,i)
        enddo
        ntarget=nsource*2
        endif
c
c
        if( idist .eq. 2 ) then
c
c       ... construct target distribution on a curve in R^3
c       
        ntarget=nsource*4
        do i=1,ntarget
        a=2*pi*dble(i)/ntarget
        target(1,i)=sin(1.1*a)/2
        target(2,i)=cos(2.1*a)/2
        target(3,i)=cos(3.1*a)/2
        enddo
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct target distribution on a unit sphere 
c       highly oversampled 
c
        ntarget=nsource*1
        do i=1,ntarget
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        target(1,i)=.5d0*cos(phi)*sin(theta) + 1
        target(2,i)=.5d0*sin(phi)*sin(theta)
        target(3,i)=.5d0*cos(theta)
        d=hkrand(0)
        d=hkrand(0)
        enddo
        endif
c
        call prinf('ntarget=*',ntarget,1)
c       
ccc        iw=71
ccc        read(iw,*) nsource
ccc        do i=1,nsource
ccc        read(iw,*) source(1,i),source(2,i),source(3,i)
ccc        enddo
c
ccc        call prinf('read data from file iw=*',iw,1)
ccc        call prinf('nsource=*',n,1)
c
c
        iprec=1
c       
        call prinf('iprec=*',iprec,1)
c       
        ifpot=1
        iffld=1
        ifhess=1
c
        ifcharge=1
        ifdipole=0
        ifquad=0
c
        ifpottarg=1
        iffldtarg=1
        ifhesstarg=1
c
        itest=1
c       
        if (itest .eq. 1) then
c       
        if (ifcharge .eq. 1 ) then
c
        do i=1,nsource
        charge(i)=1+ima
        enddo
c
        endif
c       
        if (ifdipole .eq. 1) then
c       
        do i=1,nsource
           dipstr(i)=1+ima
           dipvec(1,i)=1
           dipvec(2,i)=2
           dipvec(3,i)=3
        enddo
c
        endif
c
        if (ifquad .eq. 1) then
c       
        do i=1,nsource
           quadstr(i)=1+ima
           quadvec(1,i)=1
           quadvec(2,i)=2
           quadvec(3,i)=3
           quadvec(4,i)=4
           quadvec(5,i)=5
           quadvec(6,i)=6
        enddo
c
        endif
c
        endif
c
        if (itest .eq. 2) then
c
        if (ifcharge .eq. 1 ) then
c
        do i=1,nsource
        charge(i)=0
        enddo
c
        endif
c
        if (ifdipole .eq. 1) then
c
        do i=1,nsource
           dipstr(i)=0
           dipvec(1,i)=0
           dipvec(2,i)=0
           dipvec(3,i)=0
        enddo
c
        endif
c
        if (ifquad .eq. 1 ) then
c
        do i=1,nsource
           quadstr(i)=0
           quadvec(1,i)=0
           quadvec(2,i)=0
           quadvec(3,i)=0
           quadvec(4,i)=0
           quadvec(5,i)=0
           quadvec(6,i)=0
        enddo
c
        endif
c
        if (ifcharge .eq. 1 ) then
c
        do i=1,1
        charge(i)=1
        enddo
c
        endif
c
        if (ifdipole .eq. 1) then
c
        do i=1,1
           dipstr(i)=1
           dipvec(1,i)=1
           dipvec(2,i)=2
           dipvec(3,i)=3
        enddo
c
        if (ifquad .eq. 1) then
c       
        do i=1,1
           quadstr(i)=1
           quadvec(1,i)=1
           quadvec(2,i)=2
           quadvec(3,i)=3
           quadvec(4,i)=4
           quadvec(5,i)=5
           quadvec(6,i)=6
        enddo
c
        endif
        endif
c
        endif
c
        t1=second()
C$        t1=omp_get_wtime()
c       
        call lfmm3dpartquadtarg(ier,iprec,
     $     nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifquad,quadstr,quadvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c       
        t2=second()
C$        t2=omp_get_wtime()
c       
c       
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        call prin2('after fmm, time (sec)=*',t2-t1,1)
ccc        call prin2('after fmm, speed (points/sec)=*',nsource/(t2-t1),1)
        call prin2('after fmm, speed (points+targets/sec)=*',
     $     (nsource+ntarget)/(t2-t1),1)
c       
c
ccc        m=nsource
        m=min(nsource,10)
c
c
        ifprint=0
        if (ifprint .eq. 1) then
        call prin2('source=*',source,3*nsource)
        endif

        ifprint=0
        if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('after fmm, pot=*',pot,2*m)
        if( iffld.eq.1 ) call prin2('after fmm, fld=*',fld,3*2*m)
        if( ifhess.eq.1 ) call prin2('after fmm, hess=*',hess,6*2*m)
        endif
c
c
c
        do i=1,nsource
        if (ifpot .eq. 1) pot2(i)=0
        if (iffld .eq. 1) then
           fld2(1,i)=0
           fld2(2,i)=0
           fld2(3,i)=0
        endif
        if (ifhess .eq. 1) then
           hess2(1,i)=0
           hess2(2,i)=0
           hess2(3,i)=0
           hess2(4,i)=0
           hess2(5,i)=0
           hess2(6,i)=0
        endif
        enddo
c        
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 7160 j=1,m
        do 7150 i=1,nsource       
        if( i .eq. j ) goto 7150
        if( ifcharge .eq. 1 ) then
        call lpotfld3dhess(iffld,ifhess,source(1,i),charge(i),
     $     source(1,j),
     1     ptemp,ftemp,htemp)
        if (ifpot .eq. 1) pot2(j)=pot2(j)+ptemp
        if (iffld .eq. 1) then
           fld2(1,j)=fld2(1,j)+ftemp(1)
           fld2(2,j)=fld2(2,j)+ftemp(2)
           fld2(3,j)=fld2(3,j)+ftemp(3)
        endif
        if (ifhess .eq. 1) then
           hess2(1,j)=hess2(1,j)+htemp(1)
           hess2(2,j)=hess2(2,j)+htemp(2)
           hess2(3,j)=hess2(3,j)+htemp(3)
           hess2(4,j)=hess2(4,j)+htemp(4)
           hess2(5,j)=hess2(5,j)+htemp(5)
           hess2(6,j)=hess2(6,j)+htemp(6)
        endif
        endif
        if (ifdipole .eq. 1) then
           call lpotfld3dhess_dp(iffld,ifhess,source(1,i),
     $     dipstr(i),dipvec(1,i),
     $     source(1,j),ptemp,ftemp,htemp)
           if (ifpot .eq. 1) pot2(j)=pot2(j)+ptemp
           if (iffld .eq. 1) then
              fld2(1,j)=fld2(1,j)+ftemp(1)
              fld2(2,j)=fld2(2,j)+ftemp(2)
              fld2(3,j)=fld2(3,j)+ftemp(3)
           endif
           if (ifhess .eq. 1) then
           hess2(1,j)=hess2(1,j)+htemp(1)
           hess2(2,j)=hess2(2,j)+htemp(2)
           hess2(3,j)=hess2(3,j)+htemp(3)
           hess2(4,j)=hess2(4,j)+htemp(4)
           hess2(5,j)=hess2(5,j)+htemp(5)
           hess2(6,j)=hess2(6,j)+htemp(6)
           endif
        endif
        if (ifquad .eq. 1) then
           call lpotfld3dhess_qp(iffld,ifhess,source(1,i),
     $     quadstr(i),quadvec(1,i),
     $     source(1,j),ptemp,ftemp,htemp)
           if (ifpot .eq. 1) pot2(j)=pot2(j)+ptemp
           if (iffld .eq. 1) then
              fld2(1,j)=fld2(1,j)+ftemp(1)
              fld2(2,j)=fld2(2,j)+ftemp(2)
              fld2(3,j)=fld2(3,j)+ftemp(3)
           endif
           if (ifhess .eq. 1) then
           hess2(1,j)=hess2(1,j)+htemp(1)
           hess2(2,j)=hess2(2,j)+htemp(2)
           hess2(3,j)=hess2(3,j)+htemp(3)
           hess2(4,j)=hess2(4,j)+htemp(4)
           hess2(5,j)=hess2(5,j)+htemp(5)
           hess2(6,j)=hess2(6,j)+htemp(6)
           endif
        endif
c
 7150   continue
 7160   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('directly, pot=*',pot2,2*m)
        if( iffld.eq.1 ) call prin2('directly, fld=*',fld2,3*2*m)
        if( ifhess.eq.1 ) call prin2('directly, hess=*',hess2,6*2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpot .eq. 1)  then
        call h3derror(pot,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in potential=*',rerr,1)
        endif
c
        if (iffld .eq. 1) then
ccc        call prin2('after fmm, fld=*',fld2,2*m*3)
ccc        call prin2('directly, fld=*',fld2,2*m*3)
        call h3derror(fld,fld2,3*m,aerr,rerr)
ccc         call prin2('absolute L2 error in field=*',aerr,1)
        call prin2('relative L2 error in field=*',rerr,1)
        endif
c       
        if (ifhess .eq. 1) then
ccc        call prin2('after fmm, hess=*',hess2,2*m*6)
ccc        call prin2('directly, hess=*',hess2,2*m*6)
        call h3derror(hess,hess2,6*m,aerr,rerr)
ccc         call prin2('absolute L2 error in hessian=*',aerr,1)
        call prin2('relative L2 error in hessian=*',rerr,1)
        endif
c       
c
        do i=1,ntarget
        if (ifpottarg .eq. 1) pot2(i)=0
        if (iffldtarg .eq. 1) then
           fld2(1,i)=0
           fld2(2,i)=0
           fld2(3,i)=0
        endif
        if (ifhesstarg .eq. 1) then
           hess2(1,i)=0
           hess2(2,i)=0
           hess2(3,i)=0
           hess2(4,i)=0
           hess2(5,i)=0
           hess2(6,i)=0
        endif
        enddo
c        
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 8160 j=1,m
        do 8150 i=1,nsource        
        if( ifcharge .eq. 1 ) then
        call lpotfld3dhess(iffldtarg,ifhesstarg,
     $     source(1,i),charge(i),target(1,j),
     1     ptemp,ftemp,htemp)
        if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
        if (iffldtarg .eq. 1) then
           fld2(1,j)=fld2(1,j)+ftemp(1)
           fld2(2,j)=fld2(2,j)+ftemp(2)
           fld2(3,j)=fld2(3,j)+ftemp(3)
        endif
        if (ifhesstarg .eq. 1) then
           hess2(1,j)=hess2(1,j)+htemp(1)
           hess2(2,j)=hess2(2,j)+htemp(2)
           hess2(3,j)=hess2(3,j)+htemp(3)
           hess2(4,j)=hess2(4,j)+htemp(4)
           hess2(5,j)=hess2(5,j)+htemp(5)
           hess2(6,j)=hess2(6,j)+htemp(6)
        endif
        endif
        if (ifdipole .eq. 1) then
           call lpotfld3dhess_dp(iffldtarg,ifhesstarg,
     $     source(1,i),dipstr(i),dipvec(1,i),
     $     target(1,j),ptemp,ftemp,htemp)
           if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
           if (iffldtarg .eq. 1) then
              fld2(1,j)=fld2(1,j)+ftemp(1)
              fld2(2,j)=fld2(2,j)+ftemp(2)
              fld2(3,j)=fld2(3,j)+ftemp(3)
           endif
           if (ifhesstarg .eq. 1) then
           hess2(1,j)=hess2(1,j)+htemp(1)
           hess2(2,j)=hess2(2,j)+htemp(2)
           hess2(3,j)=hess2(3,j)+htemp(3)
           hess2(4,j)=hess2(4,j)+htemp(4)
           hess2(5,j)=hess2(5,j)+htemp(5)
           hess2(6,j)=hess2(6,j)+htemp(6)
           endif
        endif
        if (ifquad .eq. 1) then
           call lpotfld3dhess_qp(iffldtarg,ifhesstarg,
     $     source(1,i),quadstr(i),quadvec(1,i),
     $     target(1,j),ptemp,ftemp,htemp)
           if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
           if (iffldtarg .eq. 1) then
              fld2(1,j)=fld2(1,j)+ftemp(1)
              fld2(2,j)=fld2(2,j)+ftemp(2)
              fld2(3,j)=fld2(3,j)+ftemp(3)
           endif
           if (ifhesstarg .eq. 1) then
           hess2(1,j)=hess2(1,j)+htemp(1)
           hess2(2,j)=hess2(2,j)+htemp(2)
           hess2(3,j)=hess2(3,j)+htemp(3)
           hess2(4,j)=hess2(4,j)+htemp(4)
           hess2(5,j)=hess2(5,j)+htemp(5)
           hess2(6,j)=hess2(6,j)+htemp(6)
           endif
        endif
c
 8150   continue
 8160   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        if (ifprint .eq. 1) then
        if( ifpottarg.eq.1 ) 
     $     call prin2('after fmm, pottarg=*',pottarg,2*m)
        if( iffldtarg.eq.1 ) 
     $     call prin2('after fmm, fldtarg=*',fldtarg,3*2*m)
        if( ifhesstarg.eq.1 ) 
     $     call prin2('after fmm, hesstarg=*',hesstarg,6*2*m)
        endif

        if (ifprint .eq. 1) then
        if (ifpottarg .eq. 1) 
     $     call prin2('directly, pottarg=*',pot2,2*m)
        if( iffldtarg.eq.1 ) 
     $     call prin2('directly, fldtarg=*',fld2,3*2*m)
        if( ifhesstarg.eq.1 ) 
     $     call prin2('directly, hesstarg=*',hess2,6*2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpottarg .eq. 1) then
        call h3derror(pottarg,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in target potential=*',rerr,1)
        endif
c
        if (iffldtarg .eq. 1) then
        call h3derror(fldtarg,fld2,3*m,aerr,rerr)
ccc         call prin2('absolute L2 error in field=*',aerr,1)
        call prin2('relative L2 error in target field=*',rerr,1)
        endif
c       
        if (ifhesstarg .eq. 1) then
        call h3derror(hesstarg,hess2,6*m,aerr,rerr)
ccc         call prin2('absolute L2 error in hessian=*',aerr,1)
        call prin2('relative L2 error in target hessian=*',rerr,1)
        endif
        
c       
        return
        end
c
c
c
c
c
        subroutine h3dmperr(mpole1,mpole2,nterms,d)
        implicit real *8 (a-h,o-z)
c       
        complex *16 mpole1(0:nterms,-nterms:nterms)
        complex *16 mpole2(0:nterms,-nterms:nterms)
c       
        d=0
c       
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole1(n,m)-mpole2(n,m))**2
        enddo
        enddo
c       
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
c       
        return
        end
c
c
c
c
c
        subroutine h3dmpnorm(mpole,nterms,d)
        implicit real *8 (a-h,o-z)
c
        complex *16 mpole(0:nterms,-nterms:nterms)
c
        d=0
c
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole(n,m))**2
        enddo
        enddo
c
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
c
        return
        end
c
c
c
c
c
        subroutine h3derror(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        complex *16 pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
        d=d+abs(pot1(i)-pot2(i))**2
        a=a+abs(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
c
