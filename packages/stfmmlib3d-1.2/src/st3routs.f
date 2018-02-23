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
c     Stokes potential evaluation and decomposition routines     
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Auxiliary routines for Stokes FMM
c
c
        subroutine stfmm3dlap1(npts,j,pot,fld,hess,
     $     point,ifptfrc,ptfrc,ptpre,ifstrain,hessmatr)
        implicit real *8 (a-h,o-z)
        complex *16 pot(1),fld(3,1),hess(6,1)
        real *8 point(3,1),ptfrc(3,1),ptpre(1),hessmatr(3,3,1)
        real *8 rpot,temp1(3),temp2(3,3)
c
        if( ifptfrc .eq. 0 .and. ifstrain .eq. 0 ) return
c
        do k=1,npts
c
        rpot=pot(k)
        temp1(1)=fld(1,k)
        temp1(2)=fld(2,k)
        temp1(3)=fld(3,k)
        if( ifstrain .eq. 1 ) then
        temp2(1,1)=hess(1,k)
        temp2(2,2)=hess(2,k)
        temp2(3,3)=hess(3,k)
        temp2(1,2)=hess(4,k)
        temp2(1,3)=hess(5,k)
        temp2(2,3)=hess(6,k)
        temp2(2,1)=hess(4,k)
        temp2(3,1)=hess(5,k)
        temp2(3,2)=hess(6,k)
        endif
c
        if( ifptfrc .eq. 1 ) then
        ptfrc(j,k)=ptfrc(j,k)+rpot
        ptfrc(1,k)=ptfrc(1,k)+point(j,k)*temp1(1)
        ptfrc(2,k)=ptfrc(2,k)+point(j,k)*temp1(2)
        ptfrc(3,k)=ptfrc(3,k)+point(j,k)*temp1(3)
        ptpre(k)=ptpre(k)+temp1(j)*2
        endif
        if( ifstrain .eq. 1 ) then
        hessmatr(j,1,k)=hessmatr(j,1,k)-temp1(1)
        hessmatr(j,2,k)=hessmatr(j,2,k)-temp1(2)
        hessmatr(j,3,k)=hessmatr(j,3,k)-temp1(3)
        hessmatr(1,j,k)=hessmatr(1,j,k)+temp1(1)
        hessmatr(2,j,k)=hessmatr(2,j,k)+temp1(2)
        hessmatr(3,j,k)=hessmatr(3,j,k)+temp1(3)
        hessmatr(1,1,k)=hessmatr(1,1,k)-point(j,k)*temp2(1,1)
        hessmatr(1,2,k)=hessmatr(1,2,k)-point(j,k)*temp2(1,2)
        hessmatr(1,3,k)=hessmatr(1,3,k)-point(j,k)*temp2(1,3)
        hessmatr(2,1,k)=hessmatr(2,1,k)-point(j,k)*temp2(2,1)
        hessmatr(2,2,k)=hessmatr(2,2,k)-point(j,k)*temp2(2,2)
        hessmatr(2,3,k)=hessmatr(2,3,k)-point(j,k)*temp2(2,3)
        hessmatr(3,1,k)=hessmatr(3,1,k)-point(j,k)*temp2(3,1)
        hessmatr(3,2,k)=hessmatr(3,2,k)-point(j,k)*temp2(3,2)
        hessmatr(3,3,k)=hessmatr(3,3,k)-point(j,k)*temp2(3,3)
        endif
c
        enddo
c
        return
        end
c
c
c
c
c
        subroutine stfmm3dlap2(npts,pot,fld,hess,
     $     ifptfrc,ptfrc,ifstrain,hessmatr)
        implicit real *8 (a-h,o-z)
        complex *16 pot(1),fld(3,1),hess(6,1)
        real *8 ptfrc(3,1),hessmatr(3,3,1)
        real *8 rpot,temp1(3),temp2(3,3)
c
        if( ifptfrc .eq. 0 .and. ifstrain .eq. 0 ) return
c
        do k=1,npts
c
        rpot=pot(k)
        temp1(1)=fld(1,k)
        temp1(2)=fld(2,k)
        temp1(3)=fld(3,k)
        if( ifstrain .eq. 1 ) then
        temp2(1,1)=hess(1,k)
        temp2(2,2)=hess(2,k)
        temp2(3,3)=hess(3,k)
        temp2(1,2)=hess(4,k)
        temp2(1,3)=hess(5,k)
        temp2(2,3)=hess(6,k)
        temp2(2,1)=hess(4,k)
        temp2(3,1)=hess(5,k)
        temp2(3,2)=hess(6,k)
        endif
c
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,k)=ptfrc(1,k)-temp1(1)
        ptfrc(2,k)=ptfrc(2,k)-temp1(2)
        ptfrc(3,k)=ptfrc(3,k)-temp1(3)
        endif
c
        if( ifstrain .eq. 1 ) then
        hessmatr(1,1,k)=hessmatr(1,1,k)+temp2(1,1) 
        hessmatr(1,2,k)=hessmatr(1,2,k)+temp2(1,2) 
        hessmatr(1,3,k)=hessmatr(1,3,k)+temp2(1,3)
        hessmatr(2,1,k)=hessmatr(2,1,k)+temp2(2,1)
        hessmatr(2,2,k)=hessmatr(2,2,k)+temp2(2,2)
        hessmatr(2,3,k)=hessmatr(2,3,k)+temp2(2,3)
        hessmatr(3,1,k)=hessmatr(3,1,k)+temp2(3,1)
        hessmatr(3,2,k)=hessmatr(3,2,k)+temp2(3,2)
        hessmatr(3,3,k)=hessmatr(3,3,k)+temp2(3,3)
        endif
c
        enddo
c
        return
        end
c
c
c
c
c
        subroutine stfmm3dlap3(npts,j,pot,fld,
     $     ifptfrc,ptfrc,ifstrain,hessmatr)
        implicit real *8 (a-h,o-z)
        complex *16 pot(1),fld(3,1)
        real *8 ptfrc(3,1),hessmatr(3,3,1)
        real *8 rpot,temp1(3),temp2(3,3)
c
        if( ifptfrc .eq. 0 .and. ifstrain .eq. 0 ) return
c
        do k=1,npts
c
        rpot=pot(k)
        if( ifstrain .eq. 1 ) then
        temp1(1)=fld(1,k)
        temp1(2)=fld(2,k)
        temp1(3)=fld(3,k)
        endif
c
        if( ifptfrc .eq. 1 ) then
        ptfrc(j,k)=ptfrc(j,k)+rpot
        endif
        if( ifstrain .eq. 1 ) then
        hessmatr(j,1,k)=hessmatr(j,1,k)-temp1(1)
        hessmatr(j,2,k)=hessmatr(j,2,k)-temp1(2)
        hessmatr(j,3,k)=hessmatr(j,3,k)-temp1(3)
        endif
c
        enddo
c
        return
        end
c
c
c
c
c
