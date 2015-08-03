c      
c
c
c
c
c     This broadcasts the names of the event file (eventfile)
c     to the slaves, as well as a few other things.
c     J. Finke, 19 July 2005
c     updated again, name_bcast (in xec2d) and setup_bcast combined to make
c     new setup_bcast.
c     J. Finke, 7 February 2006.
      subroutine setup_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer i
      double precision junk, ran1
      integer start_address, address
      integer displacements(64), types(64), lengths(64)
      integer setup_bcast_type
c
c
c
c
c     broadcast filenames names
      call MPI_BCAST(eventfile, 30, MPI_CHARACTER, master,
     1     MPI_COMM_WORLD, ierr)
c
c     broadcast tc
      call MPI_BCAST(T_const, 1, MPI_INTEGER, master, 
     1      MPI_COMM_WORLD, ierr)
c
c     broadcast times
      call MPI_BCAST(tstop, 1, MPI_DOUBLE_PRECISION, master, 
     1      MPI_COMM_WORLD, ierr)
c
c     imubins
      call MPI_BCAST(nmu, 1, MPI_INTEGER, master,
     1     MPI_COMM_WORLD, ierr)
c
c     zone_quantities
      call MPI_BCAST(n_e, jmax*kmax, MPI_DOUBLE_PRECISION, master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(B_field,jmax*kmax,MPI_DOUBLE_PRECISION,
     1     master,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea,jmax*kmax,MPI_DOUBLE_PRECISION,
     1     master,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna,jmax*kmax,MPI_DOUBLE_PRECISION,
     1     master,MPI_COMM_WORLD, ierr)
c
c     icr
      call MPI_BCAST(cr_sent, 1, MPI_INTEGER, master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(upper_sent, 1, MPI_INTEGER, master,
     1     MPI_COMM_WORLD, ierr)
c
c     star_switch
      call MPI_BCAST(star_switch, 1, MPI_INTEGER, master, 
     1     MPI_COMM_WORLD, ierr)
c
c     common block random
      call MPI_BCAST(rand_switch, 1, MPI_INTEGER, master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(rseed, 1, MPI_INTEGER, master,
     1     MPI_COMM_WORLD, ierr)
      rseed = rseed + myid*84725
      if(rand_switch.eq.2) then
         do 10 i=1,myid+10
            junk = ran1(rseed)
c            write(*,*) i,' myid=',myid,' junk=',junk,' rseed=',rseed
 10      continue
      endif
c
c     common block photonfield
      call MPI_BCAST(E_field, nphfield, MPI_DOUBLE_PRECISION, master, 
     1     MPI_COMM_WORLD, ierr)
c
c     common block Pnth and npos
      call MPI_BCAST(Pnt,jmax*kmax*num_nt,MPI_DOUBLE_PRECISION,master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gnt, num_nt, MPI_DOUBLE_PRECISION, master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(f_nt,jmax*kmax*num_nt,MPI_DOUBLE_PRECISION,master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(n_pos,jmax*kmax*num_nt,MPI_DOUBLE_PRECISION,master,
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dg, 1, MPI_DOUBLE_PRECISION, master,
     1     MPI_COMM_WORLD, ierr)
c
c     integers
c
c     common block i_flare
      call MPI_ADDRESS(cf_sentinel, start_address, ierr)
      displacements(1) = 0
      lengths(1) = 1
      types(1) = MPI_INTEGER

c     common block izones
      call MPI_ADDRESS(nr, address, ierr)
      displacements(2) = address - start_address
      lengths(2) = 1
      types(2) = MPI_INTEGER

      call MPI_ADDRESS(nz, address, ierr)
      displacements(3) = address - start_address
      lengths(3) = 1
      types(3) = MPI_INTEGER

c     common block ps
      call MPI_ADDRESS(pair_switch, address, ierr)
      displacements(4) = address - start_address
      lengths(4) = 1
      types(4) = MPI_INTEGER
c
c     doubles
c
c     common block zones.
      call MPI_ADDRESS(z, address, ierr)
      displacements(5) = address - start_address
      lengths(5) = jmax
      types(5) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(r, address, ierr)
      displacements(6) = address - start_address
      lengths(6) = kmax
      types(6) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(rmin, address, ierr)
      displacements(7) = address - start_address
      lengths(7) = 1
      types(7) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(zmin, address, ierr)
      displacements(8) = address - start_address
      lengths(8) = 1
      types(8) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(vol, address, ierr)
      displacements(9) = address - start_address
      lengths(9) = kmax*jmax
      types(9) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(zsurf, address, ierr)
      displacements(10) = address - start_address
      lengths(10) = kmax*jmax
      types(10) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfu, address, ierr)
      displacements(11) = address - start_address
      lengths(11) = kmax
      types(11) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfl, address, ierr)
      displacements(12) = address - start_address
      lengths(12) = kmax
      types(12) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfo, address, ierr)
      displacements(13) = address - start_address
      lengths(13) = jmax
      types(13) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfi, address, ierr)
      displacements(14) = address - start_address
      lengths(14) = kmax
      types(14) = MPI_DOUBLE_PRECISION

c     common block d_flare
      call MPI_ADDRESS(r_flare, address, ierr)
      displacements(15) = address - start_address
      lengths(15) = 1
      types(15) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(z_flare, address, ierr)
      displacements(16) = address - start_address
      lengths(16) = 1
      types(16) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(t_flare, address, ierr)
      displacements(17) = address - start_address
      lengths(17) = 1
      types(17) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sigma_r, address, ierr)
      displacements(18) = address - start_address
      lengths(18) = 1
      types(18) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sigma_z, address, ierr)
      displacements(19) = address - start_address
      lengths(19) = 1
      types(19) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sigma_t, address, ierr)
      displacements(20) = address - start_address
      lengths(20) = 1
      types(20) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(flare_amp, address, ierr)
      displacements(21) = address - start_address
      lengths(21) = 1
      types(21) = MPI_DOUBLE_PRECISION
c
c     common block inj
      call MPI_ADDRESS(inj_g1, address, ierr)
      displacements(22) = address - start_address
      lengths(22) = 1
      types(22) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_g2, address, ierr)
      displacements(23) = address - start_address
      lengths(23) = 1
      types(23) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_p, address, ierr)
      displacements(24) = address - start_address
      lengths(24) = 1
      types(24) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_t, address, ierr)
      displacements(25) = address - start_address
      lengths(25) = 1
      types(25) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_v, address, ierr)
      displacements(26) = address - start_address
      lengths(26) = 1
      types(26) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_gg, address, ierr)
      displacements(27) = address - start_address
      lengths(27) = 1
      types(27) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_sigma, address, ierr)
      displacements(28) = address - start_address
      lengths(28) = 1
      types(28) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_switch, address, ierr)
      displacements(29) = address - start_address
      lengths(29) = 1
      types(29) = MPI_INTEGER
      call MPI_ADDRESS(inj_dis, address, ierr)
      displacements(30) = address - start_address
      lengths(30) = 1
      types(30) = MPI_INTEGER
      call MPI_ADDRESS(g2var_switch, address, ierr)
      displacements(31) = address - start_address
      lengths(31) = 1
      types(31) = MPI_INTEGER
      call MPI_ADDRESS(r_esc, address, ierr)
      displacements(32) = address - start_address
      lengths(32) = 1
      types(32) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(r_acc, address, ierr)
      displacements(33) = address - start_address
      lengths(33) = 1
      types(33) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_L, address, ierr)
      displacements(34) = address - start_address
      lengths(34) = 1
      types(34) = MPI_DOUBLE_PRECISION  
      call MPI_ADDRESS(pick_sw, address, ierr)
      displacements(35) = address - start_address
      lengths(35) = 1
      types(35) = MPI_INTEGER
      call MPI_ADDRESS(pick_rate, address, ierr)
      displacements(36) = address - start_address
      lengths(36) = 1
      types(36) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(g_bulk, address, ierr)
      displacements(37) = address - start_address
      lengths(37) = 1
      types(37) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(R_blr, address, ierr)
      displacements(38) = address - start_address
      lengths(38) = 1
      types(38) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(fr_blr, address, ierr)
      displacements(39) = address - start_address
      lengths(39) = 1
      types(39) = MPI_DOUBLE_PRECISION  
      call MPI_ADDRESS(R_ir, address, ierr)
      displacements(40) = address - start_address
      lengths(40) = 1
      types(40) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(fr_ir, address, ierr)
      displacements(41) = address - start_address
      lengths(41) = 1
      types(41) = MPI_DOUBLE_PRECISION    
      call MPI_ADDRESS(R_disk, address, ierr)
      displacements(42) = address - start_address
      lengths(42) = 1
      types(42) = MPI_DOUBLE_PRECISION 
      call MPI_ADDRESS(d_jet, address, ierr)
      displacements(43) = address - start_address
      lengths(43) = 1
      types(43) = MPI_DOUBLE_PRECISION     
      call MPI_ADDRESS(split1, address, ierr)
      displacements(44) = address - start_address
      lengths(44) = 1
      types(44) = MPI_INTEGER
      call MPI_ADDRESS(split2, address, ierr)
      displacements(45) = address - start_address
      lengths(45) = 1
      types(45) = MPI_INTEGER
      call MPI_ADDRESS(split3, address, ierr)
      displacements(46) = address - start_address
      lengths(46) = 1
      types(46) = MPI_INTEGER
      call MPI_ADDRESS(spl3_trg, address, ierr)
      displacements(47) = address - start_address
      lengths(47) = 1
      types(47) = MPI_INTEGER



c
      call MPI_TYPE_STRUCT(47, lengths, displacements, types, 
     1     setup_bcast_type, ierr)
      call MPI_TYPE_COMMIT(setup_bcast_type, ierr)
      call MPI_BCAST(cf_sentinel, 1, setup_bcast_type,
     1      master, MPI_COMM_WORLD, ierr)
c
c
c
      return
      end
c
c
c
c
c
c
c     This subroutine calculates the j and k (radial and vertical
c     zone numbers) from the zone number, zone.
      subroutine get_j_k(zone, j, k, r)
      implicit none
      integer zone, j, k, r
c
      j = (zone-1)/r+1
      k = zone - r*(j-1)
      return 
      end
c
c
c     
c     This calculates the zone number from j and k, the radial
c     and vertical zone numbers.
      integer function get_zone(j, k, r)
      implicit none
      integer j, k, r
c
      get_zone = k+r*(j-1)
      return 
      end
c
c
c
c
c
c     This routine makes the types used in FP_bcast to efficiently
c     broadcast to the processes.
c     J. Finke, 7 February 2006
      subroutine make_FP_bcast_type
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c
      integer status(MPI_STATUS_SIZE)

c     the above Asurfs are not used in fp_calc but are in the common
c     block with z, r etc.


c     Below are libraries from the subroutine coulomb.
c     variables for building data type FP_bcast_type
      integer start_address, address
      integer displacements(28), types(28), lengths(28)
      integer FP_Bcast_type3
c

c
c
c     variables common from outside update.
c
c      write(*,*) 'myid=', myid, ' in make_fp_bcast_type'
c

c     doubles from outside update.
c
c     common block times
c      call MPI_ADDRESS(time, start_address, ierr)
c      displacements(1) = 0
c      lengths(1) = 1
c      types(1) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(dt, address, ierr)
c      displacements(2) = address - start_address
c      lengths(2) = 2
c      types(2) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(t0, start_address, ierr)
      displacements(1) = 0
      lengths(1) = ntmax
      types(1) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(t1, address, ierr)
      displacements(2) = address - start_address
      lengths(2) = ntmax
      types(2) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(tstop, address, ierr)
c      displacements(5) = address - start_address
c      lengths(5) = 1
c      types(5) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(mcdt, address, ierr)
      displacements(3) = address - start_address
      lengths(3) = 1
      types(3) = MPI_DOUBLE_PRECISION

c     common block fic.
      call MPI_ADDRESS(F_IC, address, ierr)
      displacements(4) = address - start_address
c      lengths(7) = 100
      lengths(4) = num_nt*nphfield
      types(4) = MPI_DOUBLE_PRECISION

c      call MPI_TYPE_INDEXED(2, lengths, displacements, 
c     1     MPI_DOUBLE_PRECISION, FP_bcast_type, ierr)
      call MPI_TYPE_STRUCT(4, lengths, displacements, 
     1     types, FP_bcast_type, ierr)
      call MPI_TYPE_COMMIT(FP_bcast_type, ierr)

c     common block photonfield.
c      call MPI_ADDRESS(n_field, start_address, ierr)
c      displacements(1) = 0
c      lengths(1) = nphfield*jmax*kmax
c      types(1) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(E_field, start_address, ierr)
      displacements(1) = 0
      lengths(1) = nphfield
      types(1) = MPI_DOUBLE_PRECISION

c      call MPI_TYPE_STRUCT(2, lengths, displacements, 
c     1     types, FP_bcast_type2, ierr)
c      call MPI_TYPE_COMMIT(FP_bcast_type2, ierr)

c     common block n_pair
c      call MPI_ADDRESS(n_pos, start_address, ierr)
c      displacements(1) = 0
c      lengths(1) =  kmax*jmax*num_nt
c      types(1) = MPI_DOUBLE_PRECISION

c     common block Pnth
c      call MPI_ADDRESS(Pnt, address, ierr)
c      displacements(2) = address - start_address
c      lengths(2) = kmax*jmax*num_nt
c      types(2) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(gnt, address, ierr)
c      displacements(3) = address - start_address
c      lengths(3) = num_nt
c      types(3) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(f_nt, address, ierr)
c      displacements(4) = address - start_address
c      lengths(6) =  kmax*jmax*num_nt
c      types(6) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(dg, address, ierr)
      displacements(2) = address - start_address
      lengths(2) = 1
      types(2) = MPI_DOUBLE_PRECISION

c     common block dnpp
c      call MPI_ADDRESS(dn_pp, address, ierr)
c      displacements(6) = address - start_address
c      lengths(8) =  kmax*jmax*num_nt
c      lengths(8) = 1
c      types(8) = MPI_DOUBLE_PRECISION

c     common block annihil
c      call MPI_ADDRESS(dne_pa, address, ierr)
c      displacements(7) = address - start_address
c      lengths(9) =  kmax*jmax*num_nt
c      types(9) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(dnp_pa, address, ierr)
c      displacements(8) = address - start_address
c      lengths(10) = kmax*jmax*num_nt
c      types(10) = MPI_DOUBLE_PRECISION
c
c     variables common only between update, FP_calc, and photon_fill
c
c     common block d_update
c      call MPI_ADDRESS(Te_new, address, ierr)
c      displacements(9) = address - start_address
c      lengths(11) = kmax*jmax
c      types(11) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(dt_new, address, ierr)
      displacements(3) = address - start_address
      lengths(3) = 1
      types(3) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(lnL, address, ierr)
c      displacements(11) = address - start_address
c      lengths(13) = 1
c      types(13) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(dT_max, address, ierr)
      displacements(4) = address - start_address
      lengths(4) = 1
      types(4) = MPI_DOUBLE_PRECISION

c     common block to_fp_calc
      call MPI_ADDRESS(E_tot_old, address, ierr)
      displacements(5) = address - start_address
      lengths(5) = 1
      types(5) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(E_tot_new, address, ierr)
      displacements(6) = address - start_address
      lengths(6) = 1
      types(6) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(hr_st_total, address, ierr)
      displacements(7) = address - start_address
      lengths(7) = 1
      types(7) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(hr_total, address, ierr)
      displacements(8) = address - start_address
      lengths(8) = 1
      types(8) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(f_t_implicit, address, ierr)
      displacements(9) = address - start_address
      lengths(9) = 1
      types(9) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(dr, address, ierr)
c      displacements(20) = address - start_address
c      lengths(20) = 1
c      types(20) = MPI_DOUBLE_PRECISION

c      call MPI_ADDRESS(dz, address, ierr)
c      displacements(21) = address - start_address
c      lengths(21) = 1
c      types(21) = MPI_DOUBLE_PRECISION

      call MPI_TYPE_STRUCT(9, lengths, displacements, 
     1     types, FP_bcast_type2, ierr)
      call MPI_TYPE_COMMIT(FP_bcast_type2, ierr)
c
c
c

c      call MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION,
c     1      master, MPI_COMM_WORLD, ierr)
c      call MPI_BCAST(time, 1, FP_bcast_type, 
c     1      master, MPI_COMM_WORLD, ierr)
c      call MPI_BCAST(n_field, 1, FP_bcast_type2, 
c     1      master, MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c
c
c
c     This subroutine will broadcast to all of the processes
c     data needed for the Fokker-Planck routine.
c     J. Finke, 7 March 2005
      subroutine FP_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c
      integer status(MPI_STATUS_SIZE)
c     the above Asurfs are not used in fp_calc but are in the common
c     block with z, r etc.
c
c
c      write(*,*) 'myid=', myid, ' in fp_bcast'
c
c
c      call MPI_BCAST(time, 1, FP_bcast_type, 
c     1      master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(t0, 1, FP_bcast_type, 
     1      master, MPI_COMM_WORLD, ierr) ! Xuhui 5/28/09
      call MPI_BCAST(E_field, 1, FP_bcast_type2, 
     1      master, MPI_COMM_WORLD, ierr)
c
c
c     
      return
      end
c
c
c
c     FP_send_job sends the zone number 
c     from the master node to the available slave.  It
c     must be paired with FP_recv_job by the slave node.
c     J. Finke,  7 March 2005
      subroutine FP_send_job(available_proc, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer available_proc, zone
c
      integer ans_size
      parameter (ans_size = num_nt*7+nphfield+30)
c
      integer code, j, k, i
      integer status(MPI_STATUS_SIZE)
c
      double precision ans(ans_size)
c
c      write(*,*) 'send job myid=', myid, ' zone=', zone
      call get_j_k(zone, j, k, nr)
      code = j*100 + k
c
c     common block zone_quantities
      ans(1) = tea(j,k)
      ans(2) = tna(j,k)
      ans(3) = n_e(j,k)
      ans(4) = B_field(j,k)
c     common block vol_em
      ans(5) = Eloss_sy(j,k)
      ans(6) = Eloss_cy(j,k)
      ans(7) = Eloss_br(j,k)
      ans(8) = Eloss_th(j,k)
      ans(9) = Eloss_tot(j,k)
c     common block deposition
      ans(10) = edep(j,k)
      ans(11) = prdep(j,k)
c     common block cens_energy
      ans(12) = ecens(j,k)
c     common block ecold
      ans(13) = ec_old(j,k)
c     common block fpair
      ans(14) = f_pair(j,k)
c     common block turbulence.
      ans(15) = q_turb(j,k)
      ans(16) = turb_lev(j,k)
c     common block nontherm.  These variables are modified by
c     FP_calc and update and returned to xec2d.
      ans(17) = amxwl(j,k)
      ans(18) = gmin(j,k)
      ans(19) = gmax(j,k)
      ans(20) = gbar_nth(j,k)
      ans(21) = N_nth(j,k)
      ans(22) = p_nth(j,k)
c     common block Pnth
      do 10 i=1, num_nt
         ans(22+i) = Pnt(j, k, i)
         ans(num_nt+i+22) = gnt(i)
         ans(num_nt+num_nt+i+22) = f_nt(j,k,i)
 10   continue
c     common block photonfield
      do 20 i=1, nphfield
         ans(num_nt*3+22+i) = n_field(i,j,k)
 20   continue
c      write(*,*) 'myid=',myid,' n_field=',n_field(50,1,1)
c      stop
c     common block n_pair
      do 30 i=1, num_nt
         ans(num_nt*3+nphfield+22+i) = n_pos(j,k,i)
 30   continue
c     common block dnpp
      do 40 i=1, num_nt
         ans(num_nt*4+nphfield+22+i) = dn_pp(j,k,i)
 40   continue
c     common block annihil
      do 50 i=1, num_nt
         ans(num_nt*5+nphfield+22+i) = dne_pa(j,k,i)
         ans(num_nt*6+nphfield+22+i) = dnp_pa(j,k,i)
 50   continue
c      ans(num_nt*3+23) = dg
c     common block to_fp_calc
c      ans(num_nt*3+24) = E_tot_old
c      ans(num_nt*3+25) = E_tot_new
c      ans(num_nt*3+26) = hr_st_total
c      ans(num_nt*3+27) = hr_total
c      ans(num_nt*3+28) = f_t_implicit
c      ans(num_nt*3+29) = dr
c      ans(num_nt*3+30) = dz
c
      call MPI_SEND(ans, ans_size, MPI_DOUBLE_PRECISION, 
     1     available_proc, code, MPI_COMM_WORLD, ierr)
c      write(*,*) 'myid=', myid, ' end of fp_send_job'
c
c
      return
      end
c
c
c     
c
c     This subroutine send the signal to a node that there are
c     no more calculations for it to perform.
c     J. Finke 18 March 2005
      subroutine FP_send_end_signal(node)
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
c      write(*,*) 'myid=', myid, ' in send_end_signal, end_signal=', 
c     1     end_signal
      call MPI_SEND(0.d0, 1, MPI_DOUBLE_PRECISION, node, 
     1     end_signal, MPI_COMM_WORLD, ierr)
c
      return
      end
c
c
c
c
c
c     FP_recv_job recieves the zone number 
c     from the master node by the slave node.  It must 
c     be paired with FP_send_job or FP_send_end_signal 
c     in the master node.
c     J. Finke, 7 March 2005
      subroutine FP_recv_job(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer zone, get_zone
c
      integer ans_size
      parameter (ans_size = num_nt*7+nphfield+30)
c
      integer j, k, i, code, end_signal
      integer status(MPI_STATUS_SIZE)
      double precision ans(ans_size)
c
c      write(*,*) 'recv job myid=', myid, ' zone=', zone
      end_signal = jmax*kmax+1
c      write(*,*) 'myid=', myid, ' in fp_recv_job'
c
      call MPI_RECV(ans, ans_size, MPI_DOUBLE_PRECISION, master, 
     1     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      code = status(MPI_TAG)
      if(code.eq.end_signal) then
         zone = end_signal
         return
      else
         j=code/100
         k=code-j*100
         zone = get_zone(j, k, nr)
      endif
c
c     common block zone_quantities
      tea(j,k) = ans(1)
      tna(j,k) = ans(2)
      n_e(j,k) = ans(3)
      B_field(j,k) = ans(4)
c      write(*,*) 'myid=', myid, ' in fp_recv_job, j=', j      
c      write(*,*) 'myid=', myid, ' in fp_recv_job, k=', k
c      write(*,*) 'myid=', myid, ' in fp_recv_job, zone=', zone
c      tea(j,k) = ans
c      write(*,*) 'myid=', myid, ' tea(j,k) =', tea(j,k)

c     common block vol_em
      Eloss_sy(j,k) = ans(5)
      Eloss_cy(j,k) = ans(6)
      Eloss_br(j,k) = ans(7)
      Eloss_th(j,k) = ans(8)
      Eloss_tot(j,k) = ans(9)
c     common block deposition
      edep(j,k) = ans(10)
      prdep(j,k) = ans(11)
c     common block cens_energy
      ecens(j,k) = ans(12)
c     common block ecold
      ec_old(j,k) = ans(13)
c     common block fpair
      f_pair(j,k) = ans(14)
c     common block turbulence.
      q_turb(j,k) = ans(15)
      turb_lev(j,k) = ans(16)
c     common block nontherm.  These variables are modified by
c     FP_calc and update and returned to xec2d.
      amxwl(j,k) = ans(17)
      gmin(j,k) = ans(18)
      gmax(j,k) = ans(19)
      gbar_nth(j,k) = ans(20)
      N_nth(j,k) = ans(21)
      p_nth(j,k) = ans(22)
c     common block Pnth
      do 10 i=1, num_nt
         Pnt(j,k,i) = ans(22+i)
         gnt(i) = ans(num_nt+i+22)
         f_nt(j,k,i) = ans(num_nt+num_nt+i+22)
 10   continue
c     common block photonfield
      do 20 i=1, nphfield
         n_field(i,j,k) = ans(num_nt*3+22+i) 
 20   continue
c     common block n_pair
      do 30 i=1, num_nt
         n_pos(j,k,i) = ans(num_nt*3+nphfield+22+i) 
 30   continue
c     common block dnpp
      do 40 i=1, num_nt
         dn_pp(j,k,i) = ans(num_nt*4+nphfield+22+i)
 40   continue
c     common block annihil
      do 50 i=1, num_nt
         dne_pa(j,k,i) = ans(num_nt*5+nphfield+22+i) 
         dnp_pa(j,k,i) = ans(num_nt*6+nphfield+22+i) 
 50   continue
c      write(*,*) 'recv j,k=',j,k,' dne_pa=',dne_pa(j,k,num_nt),
c     1     ' ans=', ans(num_nt*5+nphfield+22+num_nt) 
c      write(*,*) 'recv j,k=',j,k,' dnp_pa=',dnp_pa(j,k,num_nt),
c     1     ' ans=', ans(num_nt*6+nphfield+22+num_nt) 
c      write(*,*) 'myid=',myid,' j,k=',j,k,' dn_pp=',dn_pp(j,k,num_nt),
c     1     ' dne_pa=',dne_pa(j,k,num_nt),' dnp_pa=',dnp_pa(j,k,num_nt)
c      dg = ans(num_nt*3+23)
c     common block to_fp_calc
c      E_tot_old = ans(num_nt*3+24)
c      E_tot_new = ans(num_nt*3+25)
c      hr_st_total = ans(num_nt*3+26)
c      hr_total = ans(num_nt*3+27)
c      f_t_implicit = ans(num_nt*3+28)
c      dr = ans(num_nt*3+29)
c      dz = ans(num_nt*3+30)
c
c
c
c      write(*,*) 'myid=', myid, ' end of recv_job'
      return
      end
c
c
c
c     FP_send_result sends the results from FP_calc from the
c     slave node to the master node.  It must be paired with
c     FP_recv_result.
c     J. Finke, 7 March 2005
      subroutine FP_send_result(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer zone
c
      integer ans_size
      parameter (ans_size = num_nt*4+30)
c
      integer code, i, j, k
c
      integer status(MPI_STATUS_SIZE)
      double precision ans(ans_size)
c
c     nontherm variables, except for gbar_nth and N_nth, are modified
c     by FP_calc and send back to the master node by this subroutine.
c     tea is calculated in FP_calc and returned to the master node.
c     to_fp_calc has variables modified by FP_calc and returned
c     to the master node with this subroutine.
c      write(*,*) 'send result myid=', myid, ' zone=', zone
      call get_j_k(zone, j, k, nr)
c
c     set up the code for the zone numbers.  They are sent
c     to the master node as the tag.
      code = j*100 + k
c
c      common block nontherm
      ans(1) = amxwl(j,k)
      ans(2) = gmax(j,k)
      ans(3) = gmin(j,k)
      ans(4) = p_nth(j,k)
c     electron temp.
      ans(5) = tea(j,k)
c
c     common block Pnth
      do 10 i=1, num_nt
         ans(i+5) = Pnt(j, k, i)
         ans(i+5+num_nt) = gnt(i)
         ans(i+5+num_nt+num_nt) = f_nt(j,k,i)
         ans(i+5+num_nt*3) = n_pos(j,k,i) 
 10   continue
      ans(6+num_nt*4) = dg
c     common block to_fp_calc
c      ans(7+num_nt*3) = E_tot_old
c      ans(8+num_nt*3) = E_tot_new
c      ans(9+num_nt*3) = hr_total
c      ans(10+num_nt*3) = hr_st_total
      ans(11+num_nt*4) = f_t_implicit
      ans(12+num_nt*4) = dr
      ans(13+num_nt*4) = dz
c     common block d_update.
      ans(14+num_nt*4) = Te_new(j,k)
c     common block zone quantities and npos.
      ans(15+num_nt*4) = n_e(j,k)
      ans(16+num_nt*4) = B_field(j,k)
      ans(17+num_nt*4) = tea(j,k)
      ans(18+num_nt*4) = tna(j,k) ! Xuhui
c      write(*,*) 'myid=', myid, ' ans=', ans
c      write(*,*) 'myid=', myid, ' Te_new(j,k)=', Te_new(j,k)
c
c      write(*,*) 'Pnt(1,1,50) =', Pnt(1,1,50)
c      write(*,*) 'Pnt(1,1,70) =', Pnt(1,1,70)
c      write(*,*) 'f_nt(1,1,50) =', f_nt(1,1,50)
c      write(*,*) 'f_nt(1,1,70) =', f_nt(1,1,70)
c      write(*,*) 'Te_new=', Te_new(j,k)
      call MPI_SEND(ans, ans_size, MPI_DOUBLE_PRECISION, 
     1     master, code, MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
c
c
c     FP_recv_result recieves the results of FP_calc from the
c     slave node.  It must be paired with FP_send_job.
c     J. Finke, 7 March 2005
      subroutine FP_recv_result(sender, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
      integer zone, get_zone, sender
c
      integer ans_size
      parameter (ans_size = num_nt*4+30)
c
      integer code, i, j, k
c
      integer status(MPI_STATUS_SIZE)
      double precision ans(ans_size)
c
c     nontherm variables, except for gbar_nth and N_nth, are modified
c     by FP_calc and recieved by the master node by this subroutine.
c     tea is calculated in FP_calc and recieved by the master node.
c     to_fp_calc has variables modified by FP_calc and recieved
c     by the master node with this subroutine.
c
c      write(*,*) 'recv result myid=', myid, ' zone=', zone
c    recieving common block nontherm
      call MPI_RECV(ans, ans_size, MPI_DOUBLE_PRECISION, 
     1             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
     1             status, ierr)
c     j and k, the zone numbers, are stored as code in the tag.
      sender = status(MPI_SOURCE)
      code = status(MPI_TAG)

      j = code/100
      k = code-j*100
      zone = get_zone(j, k, nr)

      amxwl(j,k) = ans(1)
      gmax(j,k) = ans(2)
      gmin(j,k) = ans(3)
      p_nth(j,k) = ans(4)
c     electron temp.
      tea(j,k) = ans(5)
c     common block Pnth
      do 10 i=1, num_nt
         Pnt(j,k,i) = ans(i+5)
         gnt(i) = ans(i+num_nt+5)
         f_nt(j,k,i) = ans(i+num_nt+num_nt+5)
         n_pos(j,k,i) = ans(i+3*num_nt +5) 
 10   continue
      dg = ans(4*num_nt+6)
c     common block to_fp_calc
c      E_tot_old = ans(3*num_nt+7)
c      E_tot_new = ans(3*num_nt+8)
c      hr_total = ans(3*num_nt+9)
c      write(*,*) 'myid=', myid, ' recv hr_total=', hr_total
c      hr_st_total = ans(3*num_nt+10)
      f_t_implicit = ans(4*num_nt+11)
      dr = ans(4*num_nt+12)
      dz = ans(4*num_nt+13)
c      write(*,*) 'myid=', myid, ' recv dz=', dz
c     common block d_update.
c      write(*,*) 'myid=', myid, ' recv ans=', ans
      Te_new(j,k) = ans(4*num_nt+14)
c     common block zone quantities
      n_e(j,k) = ans(4*num_nt+15)
      B_field(j,k) = ans(4*num_nt+16)
      tea(j,k) = ans(4*num_nt+17)
      tna(j,k) = ans(4*num_nt+18) ! Xuhui 4/1/09
c      write(*,*) 'myid=', myid, ' recv Te_new=', Te_new(j,k)
c
c
c      write(*,*) 'Pnt(1,1,50) =', Pnt(1,1,50)
c      write(*,*) 'Pnt(1,1,70) =', Pnt(1,1,70)
c      write(*,*) 'f_nt(1,1,50) =', f_nt(1,1,50)
c      write(*,*) 'f_nt(1,1,70) =', f_nt(1,1,70)
c      write(*,*) 'Te_new=', Te_new(j,k)
      return
      end

      
c     This routine broadcast the results after update to every slave node.
      subroutine FP_end_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'

      integer zone
c
c
      integer code, i, j, k
c
      integer status(MPI_STATUS_SIZE)
c
c     nontherm variables, except for gbar_nth and N_nth, are modified
c     by FP_calc and send back to the master node by this subroutine.
c     tea is calculated in FP_calc and returned to the master node.
c     to_fp_calc has variables modified by FP_calc and returned
c     to the master node with this subroutine.
cc
c     common block nontherm
      write(*,*)'myid, in FP_end',myid
      call MPI_BCAST(amxwl, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gmax, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gmin, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(P_nth, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
c     electron temp.
      call MPI_BCAST(amxwl, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
c     common block Pnth and n_pair
      call MPI_BCAST(Pnt, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gnt, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(f_nt, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dg, 1, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(n_pos, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
c     common block to_fp_calc
c      call MPI_BCAST(f_t_implicit, 1, MPI_DOUBLE_PRECISION,
c     1               master, MPI_COMM_WORLD, ierr)
      write(*,*)'myid, in the mid of FP_end',myid
      call MPI_BCAST(dr, 1, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dz, 1, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
c     common block d_update.
      call MPI_BCAST(Te_new, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
c     common block zone quantities and npos
      call MPI_BCAST(n_e, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(B_field, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna, jmax*kmax, MPI_DOUBLE_PRECISION,
     1               master, MPI_COMM_WORLD, ierr)
! Xuhui
c
c
      return
      end

c
c
c
c
c
c
      subroutine coulomb_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c
      integer status(MPI_STATUS_SIZE)
c
c
c
c     broadcast number of zones to all processes
      call MPI_BCAST(nz, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nr, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
c
c     broadcast coulomb scattering libraries to all processes
      call MPI_BCAST(lib_dg_ce, num_nt+num_temp_max, 
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lib_dg_cp, num_nt+num_temp_max, 
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lib_disp_ce, num_nt+num_temp_max, 
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lib_disp_cp, num_nt+num_temp_max, 
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(rate_gm1, num_nt, 
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea_min, 1, 
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea_max, 1,
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna_min, 1,
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna_max, 1,
     1     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
c
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:22:41 EDT 2006
c version: 2
c Name: J. Finke
c Commented out 'write' statements.     
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Wed Jul  5 12:10:18 EDT 2006
c version: 3
c Name: J. Finke
c Get rid of sending and recieving 'hr_total' and 
c 'hr_st_total' in 'FP_send_result' and 'FP_recv_result'. 
c This facilitates the 
c adding up of 'hr_total' and 'hr_st_total' in routine 
c 'E_add_up' in update2d.f      
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Sun Aug  6 18:55:16 EDT 2006
c version: 4
c Name: J. Finke
c In setup_bcast, it now broadcasts star_switch.   
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Wed Nov  1 15:18:29 EST 2006
c version: 5
c Name: J. Finke
c Changed setup_bcast so that it now broadcasts rand_switch. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Wed Nov  1 16:49:51 EST 2006
c version: 6
c Name: J. Finke
c Changed setup_bcast so now it broadcasts rseed. 
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c Mon Apr 06 14:13 CDT 2009
c Name: Xuhui Chen
c setup_bcast now broadcasts Pnt, f_nt, n_pos, dg, B_field
c and the injection parameters.
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Mon Apr 06 14:36 CDT 2009
c Name: Xuhui Chen
c Established the FP_end_bcast subroutine to broadcast the
c updated electron informations to every slave node.
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
