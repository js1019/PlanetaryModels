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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of various (free space) Stokes
c       Green's functions in R^3.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine green3sup(xyz,df,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes single layer Green's function
c
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static single force df at the origin.
c
c       stokeslet
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       df (real *8) - the strength of the single force source 
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),df(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
	cdinv=1/cd
        uvec(1)=df(1)*cdinv
        uvec(2)=df(2)*cdinv
        uvec(3)=df(3)*cdinv
c
        cdinv3=cdinv**3
c
c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
	dd=dd*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
        uvec(1)=uvec(1) * (+0.5d0)
        uvec(2)=uvec(2) * (+0.5d0)
        uvec(3)=uvec(3) * (+0.5d0)
c
        p=+dd
c
        return
        end
c
c
c
c
c
        subroutine green3sup_eval(source,df,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes single layer Green's function
c
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static single force df at the origin.
c
c       stokeslet
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       df (real *8) - the strength of the single force source 
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),df(3),dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
        cdinv=1/cd
        uvec(1)=df(1)*cdinv
        uvec(2)=df(2)*cdinv
        uvec(3)=df(3)*cdinv
c
        cdinv3=cdinv**3
c
        if( ifgrad .eq. 1 ) then

        px=dx*cdinv3
        py=dy*cdinv3
        pz=dz*cdinv3

        dudx(1,1)=-df(1)*px
        dudx(1,2)=-df(1)*py
        dudx(1,3)=-df(1)*pz
        dudx(2,1)=-df(2)*px
        dudx(2,2)=-df(2)*py
        dudx(2,3)=-df(2)*pz
        dudx(3,1)=-df(3)*px
        dudx(3,2)=-df(3)*py
        dudx(3,3)=-df(3)*pz

        endif

c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
        dd=dd*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
        if( ifgrad .eq. 1 ) then

        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv
c
        dudx(1,1)=dudx(1,1)+df(1)*px+dd*(1-3*cx*cx)
        dudx(1,2)=dudx(1,2)+df(2)*px+dd*( -3*cx*cy)
        dudx(1,3)=dudx(1,3)+df(3)*px+dd*( -3*cx*cz)
        dudx(2,1)=dudx(2,1)+df(1)*py+dd*( -3*cy*cx)
        dudx(2,2)=dudx(2,2)+df(2)*py+dd*(1-3*cy*cy)
        dudx(2,3)=dudx(2,3)+df(3)*py+dd*( -3*cy*cz)
        dudx(3,1)=dudx(3,1)+df(1)*pz+dd*( -3*cz*cx)
        dudx(3,2)=dudx(3,2)+df(2)*pz+dd*( -3*cz*cy)
        dudx(3,3)=dudx(3,3)+df(3)*pz+dd*(1-3*cz*cz)
c
        endif

        uvec(1)=uvec(1) * (+0.5d0)
        uvec(2)=uvec(2) * (+0.5d0)
        uvec(3)=uvec(3) * (+0.5d0)
c
        if( ifgrad .eq. 1 ) then

        dudx(1,1)=dudx(1,1)*(+0.5d0)
        dudx(1,2)=dudx(1,2)*(+0.5d0)
        dudx(1,3)=dudx(1,3)*(+0.5d0)
        dudx(2,1)=dudx(2,1)*(+0.5d0)
        dudx(2,2)=dudx(2,2)*(+0.5d0)
        dudx(2,3)=dudx(2,3)*(+0.5d0)
        dudx(3,1)=dudx(3,1)*(+0.5d0)
        dudx(3,2)=dudx(3,2)*(+0.5d0)
        dudx(3,3)=dudx(3,3)*(+0.5d0)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        p=dd
c
        return
        end
c
c
c
c
c
        subroutine green3stp(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part 
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0
c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0

        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        endif

c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        cdinv5=1/cd**5
        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
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
        subroutine green3stp_arb_eval(itype,
     $     source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3)
        dimension rnorm(3),du(3),dudx(3,3),grad(3,3)

        if( itype .eq. 1 ) then
        call green3stp_stresslet_dlp_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        if( itype .eq. 2 ) then
        call green3stp_stresslet_sym_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        if( itype .eq. 3 ) then
        call green3stp_rotlet_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        if( itype .eq. 4 ) then
        call green3stp_doublet_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        return
        end
c
c
c
c
c       
        subroutine green3stp_stresslet_dlp(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part 
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0
c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_stresslet_dlp_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0

        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        endif

c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        cdinv5=1/cd**5
        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
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
        subroutine green3stp_stresslet_sym(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes stresslet (type 2), forms symmetric part of Stokes doublet
c
c       u_i = [-r_i /r^3] n_j f_j + [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3),utmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dd=dd/cd**3
c
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
        uvec(1)=utmp(1)
        uvec(2)=utmp(2)
        uvec(3)=utmp(3)
c
        p=0
c
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_stresslet_sym_eval(source,du,rnorm,
     $     target,uvec,p,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes stresslet (type 2), forms symmetric part of Stokes doublet
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),utmp(3)
        dimension dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dd=dd/cd**3
c
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
        uvec(1)=utmp(1)
        uvec(2)=utmp(2)
        uvec(3)=utmp(3)
c
        p=0
c
        if( ifgrad .eq. 1 ) then

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        dudx(1,1)=-dd
        dudx(2,2)=-dd
        dudx(3,3)=-dd        

        dudx(1,1)=dudx(1,1)+3*dd*cx*cx
        dudx(1,2)=dudx(1,2)+3*dd*cy*cx
        dudx(1,3)=dudx(1,3)+3*dd*cz*cx
        dudx(2,1)=dudx(2,1)+3*dd*cx*cy
        dudx(2,2)=dudx(2,2)+3*dd*cy*cy
        dudx(2,3)=dudx(2,3)+3*dd*cz*cy
        dudx(3,1)=dudx(3,1)+3*dd*cx*cz
        dudx(3,2)=dudx(3,2)+3*dd*cy*cz
        dudx(3,3)=dudx(3,3)+3*dd*cz*cz

        endif

c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then

        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        cdinv5=1/cd**5
        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
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
        subroutine green3stp_rotlet(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes rotlet
c       u_i = [r_j n_j /r^3] f_i - [r_j f_j/ r^3] n_i
c       p = 0
c
c       (forms part of the Stokes doublet)
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3),utmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=+dd*du(1)
        utmp(2)=+dd*du(2)
        utmp(3)=+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
        return
        end
c
c
c
c
c
        subroutine green3stp_rotlet_eval(source,du,rnorm,
     $     target,uvec,p,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes rotlet
c       u_i = [r_j n_j /r^3] f_i - [r_j f_j/ r^3] n_i
c       p = 0
c
c       (forms part of the Stokes doublet)
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),utmp(3)
        dimension dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=+dd*du(1)
        utmp(2)=+dd*du(2)
        utmp(3)=+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        dudx(1,1)=dudx(1,1)+rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)+rnorm(2)*du(1)
        dudx(1,3)=dudx(1,3)+rnorm(3)*du(1)
        dudx(2,1)=dudx(2,1)+rnorm(1)*du(2)
        dudx(2,2)=dudx(2,2)+rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)+rnorm(3)*du(2)
        dudx(3,1)=dudx(3,1)+rnorm(1)*du(3)
        dudx(3,2)=dudx(3,2)+rnorm(2)*du(3)
        dudx(3,3)=dudx(3,3)+rnorm(3)*du(3)

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv
        
        dd=3*(cx*rnorm(1)+cy*rnorm(2)+cz*rnorm(3))
        dudx(1,1)=dudx(1,1)-dd*cx*du(1)
        dudx(1,2)=dudx(1,2)-dd*cy*du(1)
        dudx(1,3)=dudx(1,3)-dd*cz*du(1)
        dudx(2,1)=dudx(2,1)-dd*cx*du(2)
        dudx(2,2)=dudx(2,2)-dd*cy*du(2)
        dudx(2,3)=dudx(2,3)-dd*cz*du(2)
        dudx(3,1)=dudx(3,1)-dd*cx*du(3)
        dudx(3,2)=dudx(3,2)-dd*cy*du(3)
        dudx(3,3)=dudx(3,3)-dd*cz*du(3)

        dudx(1,1)=dudx(1,1)-rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)-rnorm(1)*du(2)
        dudx(1,3)=dudx(1,3)-rnorm(1)*du(3)
        dudx(2,1)=dudx(2,1)-rnorm(2)*du(1)
        dudx(2,2)=dudx(2,2)-rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)-rnorm(2)*du(3)
        dudx(3,1)=dudx(3,1)-rnorm(3)*du(1)
        dudx(3,2)=dudx(3,2)-rnorm(3)*du(2)
        dudx(3,3)=dudx(3,3)-rnorm(3)*du(3)

        dd=3*(cx*du(1)+cy*du(2)+cz*du(3))
        dudx(1,1)=dudx(1,1)+dd*cx*rnorm(1)
        dudx(1,2)=dudx(1,2)+dd*cy*rnorm(1)
        dudx(1,3)=dudx(1,3)+dd*cz*rnorm(1)
        dudx(2,1)=dudx(2,1)+dd*cx*rnorm(2)
        dudx(2,2)=dudx(2,2)+dd*cy*rnorm(2)
        dudx(2,3)=dudx(2,3)+dd*cz*rnorm(2)
        dudx(3,1)=dudx(3,1)+dd*cx*rnorm(3)
        dudx(3,2)=dudx(3,2)+dd*cy*rnorm(3)
        dudx(3,3)=dudx(3,3)+dd*cz*rnorm(3)

        dudx(1,1)=dudx(1,1)*cdinv3
        dudx(1,2)=dudx(1,2)*cdinv3
        dudx(1,3)=dudx(1,3)*cdinv3
        dudx(2,1)=dudx(2,1)*cdinv3
        dudx(2,2)=dudx(2,2)*cdinv3
        dudx(2,3)=dudx(2,3)*cdinv3
        dudx(3,1)=dudx(3,1)*cdinv3
        dudx(3,2)=dudx(3,2)*cdinv3
        dudx(3,3)=dudx(3,3)*cdinv3
        
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
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
        subroutine green3stp_doublet(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function
c
c       Stokes doublet = symmetric stresslet (type 2) + rotlet
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3),utmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=utmp(1)+dd*du(1)
        utmp(2)=utmp(2)+dd*du(2)
        utmp(3)=utmp(3)+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        cdinv5=1/cd**5
        dd=3*dd*cdinv5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_doublet_eval(source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function
c
c       Stokes doublet = stresslet (type 2) + rotlet
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),utmp(3)
        dimension dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=utmp(1)+dd*du(1)
        utmp(2)=utmp(2)+dd*du(2)
        utmp(3)=utmp(3)+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

c       ... symmetric part
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dudx(1,1)=-dd
        dudx(2,2)=-dd
        dudx(3,3)=-dd

        dudx(1,1)=dudx(1,1)+3*dd*cx*cx
        dudx(1,2)=dudx(1,2)+3*dd*cy*cx
        dudx(1,3)=dudx(1,3)+3*dd*cz*cx
        dudx(2,1)=dudx(2,1)+3*dd*cx*cy
        dudx(2,2)=dudx(2,2)+3*dd*cy*cy
        dudx(2,3)=dudx(2,3)+3*dd*cz*cy
        dudx(3,1)=dudx(3,1)+3*dd*cx*cz
        dudx(3,2)=dudx(3,2)+3*dd*cy*cz
        dudx(3,3)=dudx(3,3)+3*dd*cz*cz

c       ... antisymmetric part
        dudx(1,1)=dudx(1,1)+rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)+rnorm(2)*du(1)
        dudx(1,3)=dudx(1,3)+rnorm(3)*du(1)
        dudx(2,1)=dudx(2,1)+rnorm(1)*du(2)
        dudx(2,2)=dudx(2,2)+rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)+rnorm(3)*du(2)
        dudx(3,1)=dudx(3,1)+rnorm(1)*du(3)
        dudx(3,2)=dudx(3,2)+rnorm(2)*du(3)
        dudx(3,3)=dudx(3,3)+rnorm(3)*du(3)

        dd=3*(cx*rnorm(1)+cy*rnorm(2)+cz*rnorm(3))
        dudx(1,1)=dudx(1,1)-dd*cx*du(1)
        dudx(1,2)=dudx(1,2)-dd*cy*du(1)
        dudx(1,3)=dudx(1,3)-dd*cz*du(1)
        dudx(2,1)=dudx(2,1)-dd*cx*du(2)
        dudx(2,2)=dudx(2,2)-dd*cy*du(2)
        dudx(2,3)=dudx(2,3)-dd*cz*du(2)
        dudx(3,1)=dudx(3,1)-dd*cx*du(3)
        dudx(3,2)=dudx(3,2)-dd*cy*du(3)
        dudx(3,3)=dudx(3,3)-dd*cz*du(3)

        dudx(1,1)=dudx(1,1)-rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)-rnorm(1)*du(2)
        dudx(1,3)=dudx(1,3)-rnorm(1)*du(3)
        dudx(2,1)=dudx(2,1)-rnorm(2)*du(1)
        dudx(2,2)=dudx(2,2)-rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)-rnorm(2)*du(3)
        dudx(3,1)=dudx(3,1)-rnorm(3)*du(1)
        dudx(3,2)=dudx(3,2)-rnorm(3)*du(2)
        dudx(3,3)=dudx(3,3)-rnorm(3)*du(3)

        dd=3*(cx*du(1)+cy*du(2)+cz*du(3))
        dudx(1,1)=dudx(1,1)+dd*cx*rnorm(1)
        dudx(1,2)=dudx(1,2)+dd*cy*rnorm(1)
        dudx(1,3)=dudx(1,3)+dd*cz*rnorm(1)
        dudx(2,1)=dudx(2,1)+dd*cx*rnorm(2)
        dudx(2,2)=dudx(2,2)+dd*cy*rnorm(2)
        dudx(2,3)=dudx(2,3)+dd*cz*rnorm(2)
        dudx(3,1)=dudx(3,1)+dd*cx*rnorm(3)
        dudx(3,2)=dudx(3,2)+dd*cy*rnorm(3)
        dudx(3,3)=dudx(3,3)+dd*cz*rnorm(3)

        dudx(1,1)=dudx(1,1)*cdinv3
        dudx(1,2)=dudx(1,2)*cdinv3
        dudx(1,3)=dudx(1,3)*cdinv3
        dudx(2,1)=dudx(2,1)*cdinv3
        dudx(2,2)=dudx(2,2)*cdinv3
        dudx(2,3)=dudx(2,3)*cdinv3
        dudx(3,1)=dudx(3,1)*cdinv3
        dudx(3,2)=dudx(3,2)*cdinv3
        dudx(3,3)=dudx(3,3)*cdinv3

        endif

c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        cdinv5=1/cd**5
        dd=3*dd*cdinv5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then
c
        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
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
