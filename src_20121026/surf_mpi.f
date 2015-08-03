c     This subroutine will broadcast to all processes the data needed
c     for the inner and outer surfaces part of imcsurf2d.
c     J. Finke,  6 April 2005
      subroutine z_surf_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c
      integer status(MPI_STATUS_SIZE)

      double precision esurf
c
c     broadcast integers
c
c     common block timeindex
      call MPI_BCAST(ti, 1, MPI_INTEGER, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block ps
      call MPI_BCAST(pair_switch, 1, MPI_INTEGER, master, 
     1     MPI_COMM_WORLD, ierr)
c     common block ienergy
      call MPI_BCAST(nphtotal, 1, MPI_INTEGER, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nph_lc, 1, MPI_INTEGER, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block imubins
      call MPI_BCAST(nmu, 1, MPI_INTEGER, master, 
     1      MPI_COMM_WORLD, ierr)
c
c     broadcast doubles
c
c     common block photonfield
      call MPI_BCAST(E_field, nphfield, MPI_DOUBLE_PRECISION, master,
     1     MPI_COMM_WORLD, ierr)
c      call MPI_BCAST(n_field, jmax*kmax*nphfield, MPI_DOUBLE_PRECISION, 
c     1     master, MPI_COMM_WORLD, ierr) 
c     This broadcast made all the n_field in slave nodes zero. ! Xuhui 10/15/10

c     common block zone_quantities
      call MPI_BCAST(tea, jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
c     common block nontherm
      call MPI_BCAST(amxwl, jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gmin, jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gmax, jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(p_nth, jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gbar_nth, jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
c     common block bdtemp
      call MPI_BCAST(tbbi, jmax*ntmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tbbo, jmax*ntmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block zones
      call MPI_BCAST(z, jmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(r, kmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block times
      call MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dt, 2, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block eweights
      call MPI_BCAST(ewsurfi, jmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ewsurfo, jmax*ntmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block energy
      call MPI_BCAST(hu, nphomax+1, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(Elcmin, nphlcmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(Elcmax, nphlcmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block Pnth
      call MPI_BCAST(Pnt, num_nt*jmax*kmax, MPI_DOUBLE_PRECISION, 
     1      master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gnt, num_nt, MPI_DOUBLE_PRECISION, 
     1      master, MPI_COMM_WORLD, ierr)
c     common block ddh
      call MPI_BCAST(Ed_abs, kmax, MPI_DOUBLE_PRECISION, 
     1      master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(Ed_in, kmax, MPI_DOUBLE_PRECISION, 
     1      master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(Ed_ref, kmax, MPI_DOUBLE_PRECISION, 
     1      master, MPI_COMM_WORLD, ierr)
c     common block vol_em
      call MPI_BCAST(E_ph, n_vol, MPI_DOUBLE_PRECISION, master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(kappa_tot, n_vol*jmax*kmax, MPI_DOUBLE_PRECISION, 
     1     master, MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c
c     z_surf_send_job sends the job of the vertical surface number
c     (js) to the available process.  It must be paired with 
c     z_surf_recv_job.
c     J. Finke, 6 April 2005
      subroutine z_surf_send_job(available_proc, js)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer available_proc, js
c
      integer ans_size
      parameter(ans_size = 63)
c
      integer status(MPI_STATUS_SIZE)
      integer ans(ans_size), i
c
c
c     common block random
      ans(1) = zseeds(js)
c     common block s_fnames
      do 10 i=1, 30
         ans(i+1)    = ichar( i_fname(js,ti)(i:i) )
         ans(i+30+1) = ichar( o_fname(js,ti)(i:i) )
 10   continue
c     common block ph_numbers
      ans(62) = nsurfi(js)
      ans(63) = nsurfo(js)
c
      call MPI_SEND(ans, ans_size, MPI_INTEGER, 
     1     available_proc, js, MPI_COMM_WORLD, ierr)
c
c
c
      return
      end
c
c
c
c     surf_send_end_signal sends the end signal to a node when
c     there are no more calculations to perform.
c     J. Finke, 6 April 2005
      subroutine surf_send_end_signal(node)
      implicit none      
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer node
c
      integer end_signal
      integer status(MPI_STATUS_SIZE)
c
c
      end_signal = jmax*kmax+1
      call MPI_SEND(0, 1, MPI_INTEGER, node, 
     1     end_signal, MPI_COMM_WORLD, ierr)
c
      return
      end
c
c
c
c     z_surf_recv recieves the job of the vertical zone, js, from
c     the master node.  It must be paired with z_surf_send.
c     J. Finke, 6 April 2005
      subroutine z_surf_recv_job(js)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer js
c
      integer ans_size
      parameter(ans_size = 63)
c
      integer end_signal
      integer status(MPI_STATUS_SIZE)
      integer ans(ans_size), i
c
      end_signal = jmax*kmax+1
c
      call MPI_RECV(ans, ans_size, MPI_INTEGER, master, 
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      js = status(MPI_TAG)
      if(js.eq.end_signal) then
         js = end_signal
         return
      endif
c
c
c     common block random
      zseeds(js) = ans(1)
c     common block s_fnames
      do 10 i=1, 30
         i_fname(js,ti)(i:i)  = char( ans(i+1) )
         o_fname(js,ti)(i:i)  = char( ans(i+30+1) )
 10   continue
c     common block ph_numbers
      nsurfi(js) = ans(62)
      nsurfo(js) = ans(63)
c
c
      return
      end
c
c
c     Sends the result of a vertical zone (js) calculation to the master
c     node.  Must be paired with z_surf_recv_result.
c     J. Finke, 27 April 2005
      subroutine z_surf_send_result(js)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer js
c
      integer ans_size
      parameter (ans_size = 2)
c
      integer ans(ans_size)
      integer status(MPI_STATUS_SIZE)
c
c
      ans(1) = npsurfz(js)
      call MPI_SEND(ans, ans_size, MPI_INTEGER, master, js,
     1     MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c     Recieves the result of a vertical zone (js) from the
c     slave node (sender).  Must be paired with z_surf_send_result.
c     J. Finke, 27 April 2005
      subroutine z_surf_recv_result(sender, js)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer js, sender
c
      integer ans_size
      parameter (ans_size = 2)
c
      integer ans(ans_size)
      integer status(MPI_STATUS_SIZE)
c
c
      call MPI_RECV(ans, ans_size, MPI_INTEGER, MPI_ANY_SOURCE,
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
c
      sender = status(MPI_SOURCE)
      js = status(MPI_TAG)
c
      npsurfz(js) = ans(1)
c
c
      return
      end
c
c
c
c
c     This subroutine will broadcast to all processes the data needed
c     for the upper and lower surfaces part of imcsurf2d.
c     J. Finke,  27 April 2005
      subroutine r_surf_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer status(MPI_STATUS_SIZE)
      double precision esurf
c
c
c
c
c     broadcast doubles
c
c     common block bdtemp
      call MPI_BCAST(tbbu, jmax*ntmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tbbl, jmax*ntmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c     common block eweights
      call MPI_BCAST(ewsurfu, jmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ewsurfl, jmax*ntmax, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c
c     r_surf_send_job sends the job of the radial surface number
c     (ks) to the available process.  It must be paired with 
c     r_surf_recv_job.
c     J. Finke, 11 April 2005
      subroutine r_surf_send_job(available_proc, ks)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer available_proc, ks
c
      integer ans_size
      parameter(ans_size = 63)
c
      integer status(MPI_STATUS_SIZE)
      integer ans(ans_size), i
c
c     common block random
      ans(1) = rseeds(ks)
c     common block s_fnames
      do 10 i=1, 30
         ans(i+1)    = ichar( u_fname(ks,ti)(i:i) )
         ans(i+30+1) = ichar( l_fname(ks,ti)(i:i) )
 10   continue
c     common block ph_numbers
      ans(62) = nsurfu(ks)
      ans(63) = nsurfl(ks)
c
      call MPI_SEND(ans, ans_size, MPI_INTEGER, 
     1     available_proc, ks, MPI_COMM_WORLD, ierr)
c
c
c
      return
      end
c
c
c     r_surf_recv recieves the job of the radial zone, ks, from
c     the master node.  It must be paired with r_surf_send.
c     J. Finke, 11 April 2005
      subroutine r_surf_recv_job(ks)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer ks
c
      integer ans_size
      parameter(ans_size = 63)
c
      integer end_signal
      integer status(MPI_STATUS_SIZE)
      integer ans(ans_size), i
c
      end_signal = jmax*kmax+1
c
      call MPI_RECV(ans, ans_size, MPI_INTEGER, master, 
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ks = status(MPI_TAG)
      if(ks.eq.end_signal) then
         ks = end_signal
         return
      endif
c
c
c     common block random
      rseeds(ks) = ans(1)
c     common block s_fnames
      do 10 i=1, 30
         u_fname(ks,ti)(i:i)  = char( ans(i+1) )
         l_fname(ks,ti)(i:i)  = char( ans(i+30+1) )
 10   continue
c     common block ph_numbers
      nsurfu(ks) = ans(62)
      nsurfl(ks) = ans(63)
c
c
      return
      end
c
c
c
c     r_surf_recv recieves the job of the radial zone, ks, from
c     the master node.  It must be paired with r_surf_recv result.
c     J. Finke, 6 April 2005
      subroutine r_surf_send_result(ks)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer ks
c
      integer ans_size
      parameter (ans_size = 2)
c
      integer ans(ans_size)
      integer status(MPI_STATUS_SIZE)
c
      ans(1) = npsurfr(ks)
      call MPI_SEND(ans, ans_size, MPI_INTEGER, master, ks,
     1     MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c     Recieves the result of a radial zone (ks) from the
c     slave node (sender).  Must be paired with r_surf_send_result.
c     J. Finke, 27 April 2005
      subroutine r_surf_recv_result(sender, ks)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer ks, sender
c
      integer ans_size
      parameter (ans_size = 2)
c
      integer ans(ans_size)
      integer status(MPI_STATUS_SIZE)
c

c
      call MPI_RECV(ans, ans_size, MPI_INTEGER, MPI_ANY_SOURCE,
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
c
      sender = status(MPI_SOURCE)
      ks = status(MPI_TAG)
c
      npsurfr(ks) = ans(1)
c
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:33:56 EDT 2006
c version: 2
c Name: J. Finke
c changed common block 'random' in many subroutines.  
