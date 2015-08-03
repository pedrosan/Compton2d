c     Broadcasts to the slaves data needed for the volume calculation.
c     J. Finke, 18 May 2005
      subroutine vol_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer status(MPI_STATUS_SIZE)
      double precision esurf
c
c
c     broadcast integers
c
c     common block timeindex
c      call MPI_BCAST(t, 1, MPI_INTEGER, master, 
c     1      MPI_COMM_WORLD, ierr)
c     common block izones
c      call MPI_BCAST(nz, 1, MPI_INTEGER, master, 
c     1      MPI_COMM_WORLD, ierr)
c      call MPI_BCAST(nr, 1, MPI_INTEGER, master, 
c     1      MPI_COMM_WORLD, ierr)
c
c     broadcast doubles
c
c     common block zones
c      call MPI_BCAST(z, jmax, MPI_DOUBLE_PRECISION, master, 
c     1      MPI_COMM_WORLD, ierr)
c      call MPI_BCAST(r, kmax, MPI_DOUBLE_PRECISION, master, 
c     1      MPI_COMM_WORLD, ierr)
c     common block times
c      call MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION, master, 
c     1      MPI_COMM_WORLD, ierr)
c      call MPI_BCAST(dt, 2, MPI_DOUBLE_PRECISION, master, 
c     1      MPI_COMM_WORLD, ierr)
c
c     common block eweights
c      write(*,*) 'myid=',myid,' in vol_bcast'
      call MPI_BCAST(ewsv, jmax*kmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)   
c     common block vol_em
      call MPI_BCAST(E_ph, n_vol, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(eps_Th, n_vol*jmax*kmax, MPI_DOUBLE_PRECISION, 
     1      master, MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c     Creates struct type for sending jobs to the slaves.
c     J. Finke, 18 May 2005
      subroutine vol_create_job_type
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c

      integer lengths(2), displacements(2), types(2), address(2)

c

c
      lengths(1) = int_ans_size
      lengths(2) = double_ans_size
      call MPI_ADDRESS( integer_ans(1), address(1), ierr)
      call MPI_ADDRESS( double_ans(1), address(2), ierr )
      displacements(1) = 0
      displacements(2) = address(2) - address(1)
      types(1) = MPI_INTEGER
      types(2) = MPI_DOUBLE_PRECISION

      call MPI_TYPE_STRUCT(2, lengths, displacements, types, 
     1     vol_job_type, ierr)
      call MPI_TYPE_COMMIT(vol_job_type, ierr)
c
      return
      end
c
c
c
c     Sends volume job for zone to available_proc.  
c     Paired with vol_recv_job.
c     J. Finke, 18 May 2005
      subroutine vol_send_job(available_proc, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer available_proc, zone
c
c
      integer j, k, i
c
      integer status(MPI_STATUS_SIZE)
c
c
      call get_j_k(zone, j, k, nr)
c
c     put items to send in arrays
c
c     common block random
      integer_ans(1) = seeds(j,k)
c     common block ph_numbers
      integer_ans(2) = nsv(j,k)
c     common block zones
      double_ans(1) = zsurf(j,k)
c     common block vol_em
      double_ans(2) = Eloss_th(j,k)
      double_ans(3) = Eloss_tot(j,k)
      do 10 i = 1, n_vol
         double_ans(3+i) = eps_tot(i,j,k)
 10   continue
c
      call MPI_SEND(integer_ans, 1, vol_job_type, available_proc, 
     1     zone, MPI_COMM_WORLD, ierr)
c
c
      return
      end
c

      subroutine vol_send_end_signal(available_proc,zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer node
c
c
      integer end_signal, zone, available_proc
      integer status(MPI_STATUS_SIZE)
c
c
      end_signal = jmax*kmax+1
c      write(*,*) 'myid=', myid, ' in send_end_signal, end_signal=', 
c     1     end_signal
      write(*,*)'sender=',available_proc,
     1 'send end. end_signal=',end_signal
      call MPI_SEND(integer_ans, 1, vol_job_type, available_proc, 
     1     end_signal, MPI_COMM_WORLD, ierr)
            write(*,*)'sender=',available_proc,'end successfully sent'
c
      return
      end
c
c     Recieves volume job for zone from master node.  
c     Paired with vol_send_job.
c     J. Finke, 18 May 2005
      subroutine vol_recv_job(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer zone
c
c
      integer j, k, i
c
      integer status(MPI_STATUS_SIZE)
      integer end_signal
c
      integer lengths(2), displacements(2), types(2), address(2)
c
c
      integer v
      logical flag
c
c
      end_signal = jmax*kmax+1
c
c
c      write(*,*)'myid=',myid,'will recieve'
      call MPI_RECV(integer_ans, 1, vol_job_type, master,MPI_ANY_TAG,
     1     MPI_COMM_WORLD, status, ierr)
c      write(*,*)'myid=',myid,'recieved'
      zone = status(MPI_TAG)
c      write(*,*)'zone tag=',zone
      if(zone.eq.end_signal) return
      call get_j_k(zone, j, k, nr)
c
c     common block random
      seeds(j,k) = integer_ans(1)
c     common block ph_numbers
      nsv(j,k) = integer_ans(2)
c     common block zones
      zsurf(j,k) = double_ans(1)
c     common block vol_em
      Eloss_th(j,k) = double_ans(2)
      Eloss_tot(j,k) = double_ans(3)

      do 10 i = 1, n_vol
         eps_tot(i,j,k) = double_ans(3+i)
 10   continue
c
c
c
      return
      end
c
c
c
c     Sends result of volume calculation for zone to master node.  
c     Paired with vol_recv_job.
c     J. Finke, 18 May 2005
      subroutine vol_send_result(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer zone
c
      integer junk
c
c
c      write(*,*)'myid=',myid,'will send result'
      call MPI_SEND(junk, 0, MPI_INTEGER, master, zone,
     1     MPI_COMM_WORLD, ierr)
c      write(*,*)'myid=',myid,'result sent'
c
c
      return
      end
c
c
c
c     Recieves result of volume calculation for zone from sender.
c     paired with vol_send_result.
c     J. Finke, 18 May 2005
      subroutine vol_recv_result(sender, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer sender, zone
c
      integer junk
      integer status(MPI_STATUS_SIZE)
c
c
c      write(*,*)'myid=',myid,'will recieve result'
      call MPI_RECV(junk, 0, MPI_INTEGER, MPI_ANY_SOURCE, 
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
c      write(*,*)'myid=',myid,'result recieved'
c
      sender = status(MPI_SOURCE)
      zone = status(MPI_TAG)
c
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:36:26 EDT 2006
c version: 2
c Name: J. Finke
c Changed common block 'random'.     
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Thu Oct  5 13:22:38 EDT 2006
c version: 3
c Name: J. Finke
c Changed routines vol_create_job_type, vol_send_job, and vol_recv_job 
c so that 
c eps_tot is sent to the slave nodes.  
