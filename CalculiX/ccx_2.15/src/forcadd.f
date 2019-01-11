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
      subroutine forcadd(node,i,val,nodeforc,ndirforc,xforc,
     &  nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,co,
     &  ikforc,ilforc,isector,add,user,idefforc,ipompc,nodempc,
     &  nmpc,ikmpc,ilmpc,labmpc)
!
!     adds a cload condition to the data base
!
      implicit none
!
      logical add,user
!
      character*20 labmpc(*)
!
      integer nodeforc(2,*),ndirforc(*),node,i,nforc,nforc_,j,
     &  iamforc(*),iamplitude,nam,ntrans,inotr(2,*),itr,idf(3),
     &  ikforc(*),ilforc(*),idof,id,k,isector,idefforc(*),ipompc(*),
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*)
!
      real*8 xforc(*),val,trab(7,*),a(3,3),co(3,*)
!
      if(ntrans.eq.0) then
         itr=0
      else
         itr=inotr(1,node)
      endif
!
!     checking for boundary conditions on rotational dofs of
!     distributing couplings 
!
      if((i.ge.4).and.(i.le.6)) then
!
!        rotational dof
!
         idof=8*(node-1)+i
         call nident(ikmpc,idof,nmpc,id)
         if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
               if(labmpc(ilmpc(id))(1:14).eq.'ROTTRACOUPLING') then
                  node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
                  i=nodempc(2,nodempc(3,ipompc(ilmpc(id))))
                  itr=0
               endif
            endif
         endif
      endif
!
!     change: transformations on rotations are taken into account
!     by the normal of the mean rotation MPC, not by expanding the
!     MPC in Carthesian coordinates
!
c      if((itr.eq.0).or.(i.eq.0).or.(i.gt.3)) then
!
!        no transformation applies to the node
!
         idof=8*(node-1)+i
         call nident(ikforc,idof,nforc,id)
         if(id.gt.0) then
            do
               if(ikforc(id).eq.idof) then
                  k=ilforc(id)
                  if(nodeforc(2,k).eq.isector) then
                     if(add.or.(idefforc(k).eq.1)) then
                        if(nam.gt.0) then
                           if(iamforc(k).ne.iamplitude) then
                              write(*,*) '*ERROR in forcadd:'
                              write(*,*) '       it is not allowed to '
                              write(*,*)'       define two concentrated'
                              write(*,*) '       loads/fluxes'
                              write(*,*) '       different amplitudes '
                              write(*,*) '       in one step'
                              write(*,*) 'node: ',node,' dof:',i
                              call exit(201)
                           endif
                        endif
                        xforc(k)=xforc(k)+val
                     else
                        xforc(k)=val
                        if(.not.user) idefforc(k)=1
                     endif
                     if(nam.gt.0) iamforc(k)=iamplitude
                     return
                  endif
                  id=id-1
                  if(id.eq.0) exit
               else
                  exit
               endif
            enddo
         endif
c
         nforc=nforc+1
         if(nforc.gt.nforc_) then
            write(*,*) '*ERROR in forcadd: increase nforc_'
            call exit(201)
         endif
         nodeforc(1,nforc)=node
         nodeforc(2,nforc)=isector
         ndirforc(nforc)=i
         xforc(nforc)=val
         if(.not.user) idefforc(nforc)=1
         if(nam.gt.0) iamforc(nforc)=iamplitude
!
!        updating ikforc and ilforc
!            
         do j=nforc,id+2,-1
            ikforc(j)=ikforc(j-1)
            ilforc(j)=ilforc(j-1)
         enddo
         ikforc(id+1)=idof
         ilforc(id+1)=nforc
c      else
c!
c!        a transformation applies
c!
c         call transformatrix(trab(1,itr),co(1,node),a)
c!
c         do j=1,3
c            idf(j)=0
c            idof=8*(node-1)+j
c            call nident(ikforc,idof,nforc,id)
c            if(id.gt.0) then
c               do
c                  if(ikforc(id).eq.idof) then
c                     k=ilforc(id)
c                     if(nodeforc(2,k).eq.isector) then
c                        idf(j)=ilforc(id)
c                        exit
c                     endif
c                     id=id-1
c                     if(id.eq.0) exit
c                  else
c                     exit
c                  endif
c               enddo
c            endif
c         enddo
c!
c         if((idf(1).ne.0).and.(.not.user)) then
c!
c!        a force was previously applied to this node. The component
c!        in direction i is filtered out (F.n) and replaced by the new
c!        value
c!
c!        if an amplitude is selected, it applies to all components
c!        of the force in the node. No separate amplitudes are allowed.
c!
c!        val=val-F.n (substract previous contribution of the force
c!        in direction n from the actual value specified by the user)
c!
c            if((.not.add).and.(idefforc(idf(i)).ne.1))
c     &        val=val-xforc(idf(1))*a(1,i)-xforc(idf(2))*a(2,i)
c     &             -xforc(idf(3))*a(3,i)
c!
c            xforc(idf(1))=xforc(idf(1))+val*a(1,i)
c            xforc(idf(2))=xforc(idf(2))+val*a(2,i)
c            xforc(idf(3))=xforc(idf(3))+val*a(3,i)
c!
c!           only first entry is tagged
c!
c            if(idefforc(idf(i)).eq.0) then
c!
c!              no force was defined before for this node and direction
c!              within the actual step
c!
c               idefforc(idf(i))=1
c            else
c!
c!              a force was defined before within the actual step for
c!              this node and direction: check whether amplitude is
c!              the same, if any
c!
c               if(nam.gt.0) then
c                  if((iamforc(idf(1)).ne.iamplitude).or.
c     &                 (iamforc(idf(2)).ne.iamplitude).or.
c     &                 (iamforc(idf(3)).ne.iamplitude)) then
c                     write(*,*) '*ERROR in forcadd:'
c                     write(*,*) '       it is not allowed to '
c                     write(*,*) '       define two concentrated'
c                     write(*,*) '       loads/fluxes with'
c                     write(*,*) '       different amplitudes '
c                     write(*,*) '       in one step'
c                     write(*,*) 'node: ',node,' dof:',i
c                     call exit(201)
c                  endif
c               endif
c            endif
c            if(nam.gt.0) then
c               do j=1,3
c                  iamforc(idf(j))=iamplitude
c               enddo
c            endif
c         else
c            do j=1,3
c               nforc=nforc+1
c               if(nforc.gt.nforc_) then
c                  write(*,*) '*ERROR in forcadd: increase nforc_'
c                  call exit(201)
c               endif
c               nodeforc(1,nforc)=node
c               nodeforc(2,nforc)=isector
c               ndirforc(nforc)=j
c               if(user) then
c                  xforc(nforc)=val
c               else
c                  xforc(nforc)=val*a(j,i)
c               endif
c               idefforc(nforc)=1
c               if(nam.gt.0) iamforc(nforc)=iamplitude
c!
c!              updating ikforc and ilforc
c! 
c               idof=8*(node-1)+j
c               call nident(ikforc,idof,nforc-1,id)
c               do k=nforc,id+2,-1
c                  ikforc(k)=ikforc(k-1)
c                  ilforc(k)=ilforc(k-1)
c               enddo
c               ikforc(id+1)=idof
c               ilforc(id+1)=nforc
c            enddo
c         endif
c      endif
!
      return
      end

