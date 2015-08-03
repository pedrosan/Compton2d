c     Parallelized version.  Each process reads in particles
c     from its own census file for the previous time step
c     and tracks them.
c     J. Finke, 5 May 2005
      subroutine imcfield2d 
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer l, sender
c
c
c
c     master node
      if(myid.eq.master) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'myid=', myid, ' calling field_send_synch'
      call field_send_synch
c      write(*,*) 'myid=', myid, ' calling field_calc'
c      call field_calc
c      write(*,*) 'myid=', myid, ' called field_calc'
c
      do 10 l = 1, (numprocs-1)
c         write(*,*) 'l=', l, ' calling field_recv, numprocs=', numprocs
         call field_recv(sender)
c         write(*,*) 'l=', l, ' sender=', sender, ' called field_recv'
 10   continue
c
c
c     slave node
      else
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) 'myid=', myid, ' calling field_recv_synch'
      call field_recv_synch
c
c      write(*,*) 'myid=', myid, ' calling field_calc'
      call field_calc
c      write(*,*) 'myid=', myid, ' calling field_send'
      call field_send
c      write(*,*) 'myid=', myid, ' called field_sendq'
c
c
c     end of imcfield2d
      endif
c      write(*,*) 'myid=', myid, ' in imcfield ecens=', ecens(1,1)
      return
      end
c
c
c
c
c     Read in photons from census file (previous time step) and track them.
c     J. Finke, 2 May 2005
      subroutine field_calc
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c       integer ucens, iucens, ducens,
c       parameter(ucens = 1000)
c       parameter(iucens = 5000)  ! Xuhui why 5000?
c       parameter(ducens = 6000)
c       parameter(ucens = 50000000)
c       parameter(iucens = 300000000)
c       parameter(ducens = 300000000) ! Xuhui ucens
c
       integer ncrem, nread, lwad, lwai, ndxin
c
       double precision cdt
c
c
c
c      
       ncrem = ndxout
       ndxout = 0 ! The re-initialization of ndxout. used to be in imcgen Xuhui 08/15/11
c       
c
       if (ncrem.eq.0) goto 900
c
       if (ncrem.ge.ucens) then
          write(*,*)'too many photons needed to read'
          stop ! Xuhui cens
          nread = ucens
          ncrem = ncrem - ucens
       else
          nread = ncrem
          ncrem = 0
       endif
c
     
c      call read_cens(nread) ! Xuhui cens
       dbufin(:)=dbufout(:)
       ibufin(:)=ibufout(:)
       ndxin = 0
c  
 110   continue
c      
         lwad = 6*ndxin
         lwai = 6*ndxin
         rpre = dbufin(lwad+1)
         zpre = dbufin(lwad+2)
         wmu = dbufin(lwad+3)
         phi = dbufin(lwad+4)
         ew = dbufin(lwad+5)
         xnu = dbufin(lwad+6)
         jgpsp = ibufin(lwai+1)
         jgplc = ibufin(lwai+2)
         jgpmu = ibufin(lwai+3)
         jph = ibufin(lwai+4)
         kph = ibufin(lwai+5)
         seeds(jph,kph) = ibufin(lwai+6)
         call initialize_rand(jph,kph) 
         dcen = c_light*dt(1)
c
         if (wmu.gt.0.99999999d0) wmu = 0.99999999d0
         if (wmu.lt.-0.99999999d0) wmu = -0.99999999d0
c
         ndxin = ndxin+1
c
         call imctrk2d(-1)
c         
         if(ndxin .ge. nread) then
            if (ncrem.le.0) goto 900
c            if (ncrem.ge.ucens) then
c               nread = ucens
c               ncrem = ncrem - ucens
c            else
               nread = ncrem
               ncrem = 0
c            endif
c           call read_cens(nread) ! Xuhui cens
            ndxin = 0
         endif
      
         goto 110
c
 900    continue
c
        return
        end
c
c
c
c
c     This subroutine is here for synchronization purposes.  In the
c     future, if something is needed for sending in imcfield,
c     it may be added here.
c     J. Finke, 9 May 2005
      subroutine field_send_synch
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer node
c
      do 10 node=1, (numprocs-1)
         call MPI_SEND(MPI_BOTTOM, 0, MPI_INTEGER, node, myid,
     1        MPI_COMM_WORLD, ierr)
 10   continue
c
      return
      end
c
c
c
c     This subroutine is here for synchronization purposes.  In the
c     future, if something is needed for recieving in imcfield,
c     it may be added here.
c     J. Finke, 9 May 2005
      subroutine field_recv_synch
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer dummy
      integer status(MPI_STATUS_SIZE)
c
c
      call MPI_RECV(dummy, 1, MPI_INTEGER, master, MPI_ANY_TAG,
     1     MPI_COMM_WORLD, status, ierr)
c
c
      return
      end
c
c
c
c
c     Sends number of the node that completed its field 
c     calculation.  Paired with field_recv.
c     J. Finke, 2 May 2005
      subroutine field_send
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c
      call MPI_SEND(master, 1, MPI_DOUBLE_PRECISION,
     1     master, myid, MPI_COMM_WORLD, ierr)
c
      return
      end
c
c
c
c     Recieves number of the node that completed its field 
c     calculation.  Paired with field_send.
c     J. Finke, 2 May 2005
      subroutine field_recv(sender)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer sender
c
      integer status(MPI_STATUS_SIZE)
      double precision dummy
c
c
      call MPI_RECV(dummy, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      sender = status(MPI_TAG)
c
      return
      end
c
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Fri Jun 16 12:26:36 EDT 2006
c version: 2
c Name: J. Finke
c fibran seed is now stored in census files. 
c imcfield is now determinable.     
