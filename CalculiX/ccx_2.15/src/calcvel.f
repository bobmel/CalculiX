!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine calcvel(ne,nactdoh,vel,bv,neq,nef)
!
!     stores the velocities into field vel
!
      implicit none
!
      integer ne,neq,nactdoh(*),i,j,nef
!
      real*8 vel(nef,0:7),bv(neq,3)
!
      do i=1,ne
         do j=1,3
            vel(i,j)=bv(i,j)
         enddo
c         write(*,*) 'calcvel ',i,(vel(i,j),j=1,3)
      enddo
!     
      return
      end
