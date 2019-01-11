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
      subroutine calcrhoel(nef,vel,rhcon,nrhcon,ielmat,ntmat_,
     &  ithermal,mi)
!
!     calculation of rho in the element centers (incompressible
!     fluids)
!
      implicit none
!
      integer nef,i,nrhcon(*),imat,ithermal,ntmat_,mi(*),
     &  ielmat(mi(3),*)
!
      real*8 t1l,vel(nef,0:7),rhcon(0:1,ntmat_,*)
!     
c$omp parallel default(none)
c$omp& shared(nef,vel,ielmat,rhcon,nrhcon,ntmat_,ithermal)
c$omp& private(i,t1l,imat)
c$omp do
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmat(1,i)
         call materialdata_rho(rhcon,nrhcon,imat,vel(i,5),t1l,ntmat_,
     &            ithermal)
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
