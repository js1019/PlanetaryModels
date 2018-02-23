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
c    $Date: 2012-04-10 19:48:49 -0400 (Tue, 10 Apr 2012) $
c    $Revision: 2891 $
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of half-space Stokes Green's functions in R^3,
c       that satisfy u=0 at the half space interface z=0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c**********************************************************************C
        subroutine green3suph_brute_v1(source,target,df,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Naive implementation of Stokes SLP with zero boundary condition
c       on lower half-space.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),uvec(3),df(3)
        dimension source(3),target(3),delta(3,3),rr(3),df_image(3)
c
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
c
c       === PART 1 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df(j)
        enddo
        enddo
c
        uout(1)=uvec(1) * (+0.5d0)
        uout(2)=uvec(2) * (+0.5d0)
        uout(3)=uvec(3) * (+0.5d0)
c
        pout=(dx*df(1)+dy*df(2)+dz*df(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df_image(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df_image(j)
        enddo
        enddo
c
        uout(1)=uout(1)-uvec(1) * (+0.5d0)
        uout(2)=uout(2)-uvec(2) * (+0.5d0)
        uout(3)=uout(3)-uvec(3) * (+0.5d0)
c
        pout=pout-(dx*df_image(1)+dy*df_image(2)+dz*df_image(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c
c
c       === PART 3 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
c
        uvec(i)=uvec(i)+(-delta(i,j)/cd**3) *df_image(j)
        uvec(i)=uvec(i)+(3*rr(i)*rr(j)/cd**5) *df_image(j)
c
        enddo
        enddo
c
        uout(1)=uout(1)+source(3)*target(3)*uvec(1)
        uout(2)=uout(2)+source(3)*target(3)*uvec(2)
        uout(3)=uout(3)+source(3)*target(3)*uvec(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(3)=uout(3)+source(3)*(rr(1)/cd**3) *df_image(1)
        uout(3)=uout(3)+source(3)*(rr(2)/cd**3) *df_image(2)
        uout(3)=uout(3)+source(3)*(rr(3)/cd**3) *df_image(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(1)=uout(1)-target(3)*(rr(1)/cd**3) *df(3)
        uout(2)=uout(2)-target(3)*(rr(2)/cd**3) *df(3)
        uout(3)=uout(3)-target(3)*(rr(3)/cd**3) *df(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(3)=uout(3)-1/cd *df(3)
c
c
        pout = pout + source(3)*uvec(3) *2
        pout = pout - rr(3)/cd**3*df(3) *2
c
ccc        call prin2('pout=*',pout,1)
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
c
        return
        end
c
c
c
c
c
c 
c**********************************************************************C
        subroutine green3suph_brute_v2(source,target,df,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Naive implementation of Stokes SLP with zero boundary condition
c       on lower half-space, yet another decomposition.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),uvec(3),df(3)
        dimension source(3),target(3),delta(3,3),rr(3),df_image(3)
c
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
        uout(1)=0
        uout(2)=0
        uout(3)=0
        pout=0

        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        
c
c       === PART 1 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df(j)
        enddo
        enddo
c
        uout(1)=uvec(1) * (+0.5d0)
        uout(2)=uvec(2) * (+0.5d0)
        uout(3)=uvec(3) * (+0.5d0)
c
        pout=+(dx*df(1)+dy*df(2)+dz*df(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
 5200   continue
c
c       === PART 2 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df(j)
        enddo
        enddo
c
        uout(1)=uout(1)-uvec(1) * (+0.5d0)
        uout(2)=uout(2)-uvec(2) * (+0.5d0)
        uout(3)=uout(3)-uvec(3) * (+0.5d0)
c
        pout=pout-(dx*df(1)+dy*df(2)+dz*df(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c
 5300   continue
c
c       === PART 3 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
c
        uvec(i)=uvec(i)+(-delta(i,j)/cd**3) *df_image(j)
        uvec(i)=uvec(i)+(3*rr(i)*rr(j)/cd**5) *df_image(j)
c
        enddo
        enddo
c
        uout(1)=uout(1)+source(3)*rr(3)*uvec(1)
        uout(2)=uout(2)+source(3)*rr(3)*uvec(2)
        uout(3)=uout(3)+source(3)*rr(3)*uvec(3)
c
        uout(1)=uout(1)-source(3)*source(3)*uvec(1)
        uout(2)=uout(2)-source(3)*source(3)*uvec(2)
        uout(3)=uout(3)-source(3)*source(3)*uvec(3)
c
        uout(3)=uout(3)+source(3)*(rr(1)/cd**3) *df_image(1)
        uout(3)=uout(3)+source(3)*(rr(2)/cd**3) *df_image(2)
        uout(3)=uout(3)+source(3)*(rr(3)/cd**3) *df_image(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(1)=uout(1)+source(3)*(rr(1)/cd**3) *df(3)
        uout(2)=uout(2)+source(3)*(rr(2)/cd**3) *df(3)
        uout(3)=uout(3)+source(3)*(rr(3)/cd**3) *df(3)
c
ccc        call prin2('uout=*',uout,3)
c
        pout = pout + source(3)*uvec(3) *2
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
c
        return
        end
c
c
c**********************************************************************C
        subroutine green3suph_brute(source,target,df,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes SLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),df(3)
        dimension source(3),target(3),df_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,pot,fld(3)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3sup(xyz,df,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... stokeslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        call green3sup(xyz,df_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... dipole image
c
        dipstr = source(3)
        call lpotfld3d_dp(iffld,sourceim,dipstr,df_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... charge image
c
        charge = df_image(3)
        call lpotfld3d(iffld,sourceim,charge,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... dlp stresslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
c       ... dipole image
c
        dipstr = rnorm_image(1)*du_image(1)+
     $           rnorm_image(2)*du_image(2)+
     $           rnorm_image(3)*du_image(3)
        dipstr = dipstr*2
        dipvec(1)=0
        dipvec(2)=0
        dipvec(3)=1
        call lpotfld3d_dp(iffld,sourceim,dipstr,dipvec,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_arb_brute(
     $     ifdouble,source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
        if( ifdouble .eq. 1 ) then
        call green3stph_brute1
     $     (source,target,du,rnorm,uout,pout)
        endif

        if( ifdouble .eq. 2 ) then
        call green3stph_brute2
     $     (source,target,du,rnorm,uout,pout)
        endif

        if( ifdouble .eq. 3 ) then
        call green3stph_brute3
     $     (source,target,du,rnorm,uout,pout)
        endif

        if( ifdouble .eq. 4 ) then
        call green3stph_brute4
     $     (source,target,du,rnorm,uout,pout)
        endif

c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute1(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... dlp stresslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
c       ... dipole image
c
        dipstr = rnorm_image(1)*du_image(1)+
     $           rnorm_image(2)*du_image(2)+
     $           rnorm_image(3)*du_image(3)
        dipstr = dipstr*2
        dipvec(1)=0
        dipvec(2)=0
        dipvec(3)=1
        call lpotfld3d_dp(iffld,sourceim,dipstr,dipvec,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute2(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6),strain(3,3)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp_stresslet_dlp(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... stresslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_stresslet_dlp(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute3(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp_rotlet(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... rotlet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_rotlet(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... dipole image
c
        dipstr = -rnorm_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,du_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... dipole image
c
        dipstr = du_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,rnorm_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute4(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp_doublet(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... doublet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_doublet(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... dipole image
c
        dipstr = -rnorm_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,du_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... dipole image
c
        dipstr = du_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,rnorm_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
        return
        end
c
c
c
c
