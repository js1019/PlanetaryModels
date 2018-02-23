c
c       Test the fundamental solutions (a.k.a Green's functions) 
c       for free space Stokes 3d problem 
c
c       The governing Stokes equations are checked by numerical differentiation
c       (please run in extended precision to get double precision accuracy)
c
c       Test parallel plates
c       
        implicit real *8 (a-h,o-z)
        dimension df(3),du(3),rnorm(3)
        dimension xyz(3),source(3),target(3),uvec(3)
        dimension umatr1(3,3),umatr2(3,3,3),uout(3)
        dimension pmatr1(3),pmatr2(3,3)
        dimension strain(3,3),stress(3,3)
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        du(1)=1
        du(2)=2
        du(3)=3
c
        df(1)=-4
        df(2)=2
        df(3)=-5
c
        dx=3
        dy=2
        dz=1
        d=sqrt(dx*dx+dy*dy+dz*dz)
c        
        rnorm(1)=dx/d
        rnorm(2)=dy/d
        rnorm(3)=dz/d
c        
        call prin2('df=*',df,3)
        call prin2('du=*',du,3)
        call prin2('rnorm=*',rnorm,3)
c
c
        source(1)=1
        source(2)=2
        source(3)=3
c
        target(1)=-1
        target(2)=3
        target(3)=0
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)

c
        call prin2('xyz=*',xyz,3)
c
        call green3sup(xyz,df,uvec,p)
        call prin2('uvec (U)=*',uvec,3)
        call prin2('p (U)=*',p,1)
c
        call green3stp(xyz,du,rnorm,uvec,p)
        call prin2('uvec (T)=*',uvec,3)
        call prin2('p (T)=*',p,1)
c
c       
c       ... simple image
c       
        source(1)=1
        source(2)=2
        source(3)=-3
c
        target(1)=-1
        target(2)=3
        target(3)=0
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        df(1)=-df(1)
        df(2)=-df(2)
        df(3)=+df(3)
c
c
cccc        call prin2('xyz=*',xyz,3)
c
        call green3sup(xyz,df,uvec,p)
        call prin2('uvec image(U)=*',uvec,3)
        call prin2('p image(U)=*',p,1)
c
        du(1)=-du(1)
        du(2)=-du(2)
        du(3)=+du(3)
c
        rnorm(1)=-rnorm(1)
        rnorm(2)=-rnorm(2)
        rnorm(3)=+rnorm(3)
c
        call green3stp(xyz,du,rnorm,uvec,p)
        call prin2('uvec image(T)=*',uvec,3)
        call prin2('p image(T)=*',p,1)
c
c       
        source(1)=1
        source(2)=2
        source(3)=-3
c
        target(1)=-1
        target(2)=3
        target(3)=0
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        df(1)=-4
        df(2)=2
        df(3)=-5
c
        call prinf('check the first fundamental solution U*',i,0)
c
c       ... check the first fundamental solution U
c
c       ... first derivatives of velocity
c
        call test1u(xyz,df,umatr1,pmatr1)
        call prin2('umatr=*',umatr1,3*3)
        call prin2('pmatr=*',pmatr1,3)
c
c       ... second derivatives of velocity
c
        call test2u(xyz,df,umatr2,pmatr2)
        call prin2('umatr=*',umatr2,3*3*3)
        call prin2('pmatr=*',pmatr2,3*3)
c
c       ... check Stokes equations numerically
c
        call test3(umatr1,umatr2,pmatr1,uout,pout)
c
c       ... this should be approximately zero
c       
        call prin2('and the value of Stokes operator,lap=*',uout,3)
        call prin2('and the value of Stokes operator,div=*',pout,1)
c
c
        source(1)=1
        source(2)=2
        source(3)=-3
c
        target(1)=-1
        target(2)=3
        target(3)=0
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        du(1)=-4
        du(2)=2
        du(3)=-5
c
        dx=3
        dy=-2
        dz=-1
        d=sqrt(dx*dx+dy*dy+dz*dz)
c        
        rnorm(1)=dx/d
        rnorm(2)=dy/d
        rnorm(3)=dz/d
c        
        call prinf('check the second fundamental solution T*',i,0)
c
c       ... check the first fundamental solution T
c
c       ... first derivatives of velocity
c
        call test1t(xyz,du,rnorm,umatr1,pmatr1)
        call prin2('umatr=*',umatr1,3*3)
        call prin2('pmatr=*',pmatr1,3)
c
c       ... second derivatives of velocity
c
        call test2t(xyz,du,rnorm,umatr2,pmatr2)
        call prin2('umatr=*',umatr2,3*3*3)
        call prin2('pmatr=*',pmatr2,3*3)
c
c       ... check Stokes equations numerically
c
        call test3(umatr1,umatr2,pmatr1,uout,pout)
c
c       ... this should be approximately zero
c       
        call prin2('and the value of Stokes operator,lap=*',uout,3)
        call prin2('and the value of Stokes operator,div=*',pout,1)
c
ccc     stop
c
c
c       ... layered media, single interface at z=0
c
c       half space green's function
c
c       ... evaluate the velocity
c       
        call prin2('====================*',i,0)
        call prin2('u at the interface*',i,0)
c
        source(1)=1.2345d0
        source(2)=2
        source(3)=3
c
        target(1)=-1.137d0
        target(2)=3
        target(3)=0
c
        df(1)=-4
        df(2)=2
        df(3)=-5
c
        call green3suph_brute(source,target,df,uvec,p)
        call prin2('uvec (U)=*',uvec,3)
        call prin2('p (U)=*',p,1)
c
ccc        stop
c
        call prinf('check the first fundamental solution U*',i,0)
c
c       ... check the first fundamental solution U
c
c       ... first derivatives of velocity
c
        call test1uh(source,target,df,umatr1,pmatr1)
        call prin2('umatr1=*',umatr1,3*3)
        call prin2('pmatr1=*',pmatr1,3)
c
c       ... second derivatives of velocity
c
        call test2uh(source,target,df,umatr2,pmatr2)
        call prin2('umatr2=*',umatr2,3*3*3)
        call prin2('pmatr2=*',pmatr2,3*3)
c
c       ... check Stokes equations numerically
c
        call test3(umatr1,umatr2,pmatr1,uout,pout)
c
c       ... this should be approximately zero
c       
        call prin2('and the value of Stokes operator,lap=*',uout,3)
        call prin2('and the value of Stokes operator,div=*',pout,1)
c
c
        source(1)=1.2345d0
        source(2)=2
        source(3)=3
c
        target(1)=-1.137d0
        target(2)=3
        target(3)=0
c
        du(1)=-4
        du(2)=2
        du(3)=-5
c
        dx=3
        dy=-2
        dz=-1
        d=sqrt(dx*dx+dy*dy+dz*dz)
c        
        rnorm(1)=dx/d
        rnorm(2)=dy/d
        rnorm(3)=dz/d
c
        do 5000 ifdouble = 1,4
        call prinf('ifdouble=*',ifdouble,1)

        call green3stph_arb_brute
     $     (ifdouble,source,target,du,rnorm,uvec,p)
        call prin2('uvec (T)=*',uvec,3)
        call prin2('p (T)=*',p,1)
c
        call prinf('check the second fundamental solution T*',i,0)
c
c       ... check the second fundamental solution T
c
c       ... first derivatives of velocity
c
        call test1th(ifdouble,source,target,du,rnorm,umatr1,pmatr1)
        call prin2('umatr1=*',umatr1,3*3)
        call prin2('pmatr1=*',pmatr1,3)
c
c       ... second derivatives of velocity
c
        call test2th(ifdouble,source,target,du,rnorm,umatr2,pmatr2)
        call prin2('umatr2=*',umatr2,3*3*3)
        call prin2('pmatr2=*',pmatr2,3*3)
c
c       ... check Stokes equations numerically
c
        call test3(umatr1,umatr2,pmatr1,uout,pout)
c
c       ... this should be approximately zero
c       
        call prin2('and the value of Stokes operator,lap=*',uout,3)
        call prin2('and the value of Stokes operator,div=*',pout,1)
c
 5000   continue
c
        stop
        end
c
c
c
c
c
        subroutine test1u(xyz,df,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension df(3),xyz(3)
        dimension xyz1(3),xyz2(3),fvec1(3),fvec2(3)
        dimension umatr(3,3),pmatr(3)
c
c       ... first derivative of free-space SLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=xyz(1)+h
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)-h
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)
c
        call green3sup(xyz1,df,fvec1,q1)
        call green3sup(xyz2,df,fvec2,q2)
c
        umatr(1,1)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,1)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,1)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(1) = (q1-q2)/(2*h)
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)+h
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)-h
        xyz2(3)=xyz(3)
c
        call green3sup(xyz1,df,fvec1,q1)
        call green3sup(xyz2,df,fvec2,q2)
c
        umatr(1,2)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,2)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,2)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(2) = (q1-q2)/(2*h)
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)+h
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)-h
c
        call green3sup(xyz1,df,fvec1,q1)
        call green3sup(xyz2,df,fvec2,q2)
c
        umatr(1,3)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,3)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,3)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(3) = (q1-q2)/(2*h)
c
        return
        end
c
c
c
        subroutine test2u(xyz,df,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension df(3),xyz(3)
        dimension xyz1(3),xyz2(3)
        dimension fvec1(3,3),fvec2(3,3),pvec1(3),pvec2(3)
        dimension umatr(3,3,3),pmatr(3,3)
c
c       ... second derivative of free-space SLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=xyz(1)+h
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)-h
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)
c
        call test1u(xyz1,df,fvec1,pvec1)
        call test1u(xyz2,df,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,1)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,1)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)+h
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)-h
        xyz2(3)=xyz(3)
c
        call test1u(xyz1,df,fvec1,pvec1)
        call test1u(xyz2,df,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,2)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,2)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)+h
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)-h
c
        call test1u(xyz1,df,fvec1,pvec1)
        call test1u(xyz2,df,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,3)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,3)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        return
        end
c
c
c
        subroutine test1uh(source,target,df,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension df(3),source(3),target(3)
        dimension xyz1(3),xyz2(3),fvec1(3),fvec2(3)
        dimension umatr(3,3),pmatr(3)
c
c       ... first derivative of half-space SLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=target(1)+h
        xyz1(2)=target(2)
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)-h
        xyz2(2)=target(2)
        xyz2(3)=target(3)
c
        call green3suph_brute(source,xyz1,df,fvec1,q1)
        call green3suph_brute(source,xyz2,df,fvec2,q2)
c
        umatr(1,1)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,1)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,1)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(1) = (q1-q2)/(2*h)
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)+h
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)-h
        xyz2(3)=target(3)
c
        call green3suph_brute(source,xyz1,df,fvec1,q1)
        call green3suph_brute(source,xyz2,df,fvec2,q2)
c
        umatr(1,2)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,2)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,2)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(2) = (q1-q2)/(2*h)
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)
        xyz1(3)=target(3)+h
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)
        xyz2(3)=target(3)-h
c
        call green3suph_brute(source,xyz1,df,fvec1,q1)
        call green3suph_brute(source,xyz2,df,fvec2,q2)
c
        umatr(1,3)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,3)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,3)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(3) = (q1-q2)/(2*h)
c
        return
        end
c
c
c
        subroutine test2uh(source,target,df,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension df(3),xyz(3)
        dimension source(3),target(3)
        dimension xyz1(3),xyz2(3)
        dimension fvec1(3,3),fvec2(3,3),pvec1(3),pvec2(3)
        dimension umatr(3,3,3),pmatr(3,3)
c
c       ... second derivative of half-space SLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=target(1)+h
        xyz1(2)=target(2)
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)-h
        xyz2(2)=target(2)
        xyz2(3)=target(3)
c
        call test1uh(source,xyz1,df,fvec1,pvec1)
        call test1uh(source,xyz2,df,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,1)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,1)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)+h
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)-h
        xyz2(3)=target(3)
c
        call test1uh(source,xyz1,df,fvec1,pvec1)
        call test1uh(source,xyz2,df,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,2)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,2)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)
        xyz1(3)=target(3)+h
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)
        xyz2(3)=target(3)-h
c
        call test1uh(source,xyz1,df,fvec1,pvec1)
        call test1uh(source,xyz2,df,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,3)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,3)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        return
        end
c
c
c
        subroutine test1th
     $     (ifdouble,source,target,du,rnorm,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension du(3),rnorm(3),source(3),target(3)
        dimension xyz1(3),xyz2(3),fvec1(3),fvec2(3)
        dimension umatr(3,3),pmatr(3)
c
c       ... first derivative of half-space SLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=target(1)+h
        xyz1(2)=target(2)
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)-h
        xyz2(2)=target(2)
        xyz2(3)=target(3)
c
        call green3stph_arb_brute
     $     (ifdouble,source,xyz1,du,rnorm,fvec1,q1)
        call green3stph_arb_brute
     $     (ifdouble,source,xyz2,du,rnorm,fvec2,q2)
c
        umatr(1,1)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,1)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,1)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(1) = (q1-q2)/(2*h)
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)+h
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)-h
        xyz2(3)=target(3)
c
        call green3stph_arb_brute
     $     (ifdouble,source,xyz1,du,rnorm,fvec1,q1)
        call green3stph_arb_brute
     $     (ifdouble,source,xyz2,du,rnorm,fvec2,q2)
c
        umatr(1,2)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,2)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,2)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(2) = (q1-q2)/(2*h)
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)
        xyz1(3)=target(3)+h
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)
        xyz2(3)=target(3)-h
c
        call green3stph_arb_brute
     $     (ifdouble,source,xyz1,du,rnorm,fvec1,q1)
        call green3stph_arb_brute
     $     (ifdouble,source,xyz2,du,rnorm,fvec2,q2)
c
        umatr(1,3)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,3)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,3)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(3) = (q1-q2)/(2*h)
c
        return
        end
c
c
c
        subroutine test2th
     $     (ifdouble,source,target,du,rnorm,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension du(3),rnorm(3),xyz(3)
        dimension source(3),target(3)
        dimension xyz1(3),xyz2(3)
        dimension fvec1(3,3),fvec2(3,3),pvec1(3),pvec2(3)
        dimension umatr(3,3,3),pmatr(3,3)
c
c       ... second derivative of half-space SLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=target(1)+h
        xyz1(2)=target(2)
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)-h
        xyz2(2)=target(2)
        xyz2(3)=target(3)
c
        call test1th(ifdouble,source,xyz1,du,rnorm,fvec1,pvec1)
        call test1th(ifdouble,source,xyz2,du,rnorm,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,1)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,1)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)+h
        xyz1(3)=target(3)
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)-h
        xyz2(3)=target(3)
c
        call test1th(ifdouble,source,xyz1,du,rnorm,fvec1,pvec1)
        call test1th(ifdouble,source,xyz2,du,rnorm,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,2)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,2)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=target(1)
        xyz1(2)=target(2)
        xyz1(3)=target(3)+h
c
        xyz2(1)=target(1)
        xyz2(2)=target(2)
        xyz2(3)=target(3)-h
c
        call test1th(ifdouble,source,xyz1,du,rnorm,fvec1,pvec1)
        call test1th(ifdouble,source,xyz2,du,rnorm,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,3)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,3)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        return
        end
c
c
c
        subroutine test3(umatr1,umatr2,pmatr1,uout,qout)
        implicit real *8 (a-h,o-z)
        dimension du(3),rnorm(3),xyz(3)
        dimension umatr1(3,3),umatr2(3,3,3),pmatr1(3)
        dimension uout(3),delta(3,3)
c
c       ... check the Stokes equations
c
        do i=1,3
        uout(i)=0
        enddo
c
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
        do i=1,3
        delta(i,i)=1
        enddo
c
c
        do i=1,3
c
        do j=1,3
        uout(i)=uout(i)+umatr2(i,j,j)
        enddo
        uout(i)=uout(i)-pmatr1(i)
c
        enddo
c
c
        qout=umatr1(1,1)+umatr1(2,2)+umatr1(3,3)
c
        return
        end
c
c        
c
        subroutine test1t(xyz,du,rnorm,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension du(3),rnorm(3),xyz(3)
        dimension xyz1(3),xyz2(3),fvec1(3),fvec2(3)
        dimension umatr(3,3),pmatr(3)
c
c       ... first derivative of free-space DLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=xyz(1)+h
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)-h
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)
c
        call green3stp(xyz1,du,rnorm,fvec1,q1)
        call green3stp(xyz2,du,rnorm,fvec2,q2)
c
        umatr(1,1)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,1)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,1)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(1) = (q1-q2)/(2*h)
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)+h
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)-h
        xyz2(3)=xyz(3)
c
        call green3stp(xyz1,du,rnorm,fvec1,q1)
        call green3stp(xyz2,du,rnorm,fvec2,q2)
c
        umatr(1,2)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,2)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,2)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(2) = (q1-q2)/(2*h)
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)+h
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)-h
c
        call green3stp(xyz1,du,rnorm,fvec1,q1)
        call green3stp(xyz2,du,rnorm,fvec2,q2)
c
        umatr(1,3)=(fvec1(1)-fvec2(1))/(2*h)
        umatr(2,3)=(fvec1(2)-fvec2(2))/(2*h)
        umatr(3,3)=(fvec1(3)-fvec2(3))/(2*h)
        pmatr(3) = (q1-q2)/(2*h)
c
        return
        end
c
c
c
        subroutine test2t(xyz,du,rnorm,umatr,pmatr)
        implicit real *8 (a-h,o-z)
        dimension du(3),rnorm(3),xyz(3)
        dimension xyz1(3),xyz2(3)
        dimension fvec1(3,3),fvec2(3,3),pvec1(3),pvec2(3)
        dimension umatr(3,3,3),pmatr(3,3)
c
c       ... second derivative of free-space DLP (numerical approximation)
c
        h=1.0e-5
ccc        h=1.0e-10
c
        xyz1(1)=xyz(1)+h
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)-h
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)
c
        call test1t(xyz1,du,rnorm,fvec1,pvec1)
        call test1t(xyz2,du,rnorm,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,1)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,1)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)+h
        xyz1(3)=xyz(3)
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)-h
        xyz2(3)=xyz(3)
c
        call test1t(xyz1,du,rnorm,fvec1,pvec1)
        call test1t(xyz2,du,rnorm,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,2)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,2)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        xyz1(1)=xyz(1)
        xyz1(2)=xyz(2)
        xyz1(3)=xyz(3)+h
c
        xyz2(1)=xyz(1)
        xyz2(2)=xyz(2)
        xyz2(3)=xyz(3)-h
c
        call test1t(xyz1,du,rnorm,fvec1,pvec1)
        call test1t(xyz2,du,rnorm,fvec2,pvec2)
c
        do i=1,3
        do j=1,3
        umatr(i,j,3)=(fvec1(i,j)-fvec2(i,j))/(2*h)
        enddo
        pmatr(i,3)=(pvec1(i)-pvec2(i))/(2*h)
        enddo
c
        return
        end
c
c
c
