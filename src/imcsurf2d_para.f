      subroutine imcsurf2d
c     This subroutine convert the input energy from external source into photonsc     and write them into census file. !comment by Xuhui Chen 7/8/08
c     

  
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
c
      integer i, n, js, ks, nxs
      integer npsurf, isurfmu
c     MPI variables
      integer num_sent, end_signal, sender
      integer zone
c
      double precision esurf
      double precision psi, wv1, wv2, tpl
c
      double precision fibran, t_average, x1, x2
c
c

c
c
c     MPI Initialization 
c      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      num_sent = 0
      end_signal = jmax*kmax+1
c      write(*,*) 'myid=', myid, 'end_signal=', end_signal
c
c
c
c     master node part
      if(myid.eq.master) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c      write(*,*) 'beginning of imcsurf2d'
      isurfmu = 1
      npsurf = 0
      esurf = 0.d0
c
c      do 120 n = 1, nphtotal
c 120  F(n) = 0.d0
c
c
c    Scan through time bins to find 
c      the correct time index t 
c
       if (ncycle.eq.0) then
          ti = 1
          goto 135
       endif
       t_average = time + 5.d-1*dt(1)
c
      do 130 ti = 1, ntime
         if (t1(ti).gt.t_average) goto 135
 130  continue
 135  continue
c
c
c
c     Vertical boundaries
c
      call z_surf_bcast
c
c     This sends the first round of zones to the slaves for processing.
         do 300 js = 1, min(nz, numprocs-1)
            zone = js
            call z_surf_send_job(js, zone)
c            write(*,*) 'myid=', myid, ' sent zone=', zone
            num_sent = num_sent + 1
 300     continue
c
c     As slaves complete processing a zone, this recieves the results
c     and sends the next zone to the slaves.
         do 302 js=1, nz
c            write(*,*) 'js=', js
            call z_surf_recv_result(sender, zone)
            npsurf = npsurf + npsurfz(js)
            if(num_sent.lt.nz) then
c           if there are still more zones to send, send next zone
c           to the available node.
               zone = num_sent + 1
               call z_surf_send_job(sender, zone)
               num_sent = num_sent + 1
            else
               call surf_send_end_signal(sender)
            endif
 302     continue
c
c
c
c
c     Radial boundaries
c
         num_sent = 0
         call r_surf_bcast
c
c     This sends the first round of zones to the slaves for processing.
         do 400 ks = 1, min(nr, numprocs-1)
            zone = ks
c            write(*,*) 'myid=', myid, ' before send zone=', zone,
c     1           ' nsurfl(zone)=', nsurfl(zone)
            call r_surf_send_job(ks, zone)
c            write(*,*) 'myid=', myid, ' sent zone=', zone
            num_sent = num_sent + 1
 400     continue
c
c     As slaves complete processing a zone, this recieves the results
c     and sends the next zone to the slaves.
         do 402 ks=1, nr
            call r_surf_recv_result(sender, zone)
c            write(*,*) 'master recieved sender=', sender, ' zone=', zone
c     1           , ' ks=', ks
            npsurf = npsurf + npsurfr(ks)
            if(num_sent.lt.nr) then
c           if there are still more zones to send, send next zone
c           to the available node.
               zone = num_sent + 1
c            write(*,*) 'myid=', myid, ' before send zone=', zone,
c     1           ' nsurfl(zone)=', nsurfl(zone)
               call r_surf_send_job(sender, zone)
c            write(*,*) 'myid=', myid, ' sent zone=', zone
               num_sent = num_sent + 1
            else
               call surf_send_end_signal(sender)
            endif
 402     continue
c
c
c     End of master part.
c
c
c      open(unit=10, file='planck_test.dat', status='unknown')
c
c
c      close(10)             
c
 800  format(' rpre = ',e14.7,', zpre = ',e14.7)
 801  format(' xnu = ',e14.7)
 802  format(' phi = ',e14.7,',  wmu = ',e14.7)
 803  format(' jgpsp = ',i3,',  jgplc = ',i2)
 804  format(' jgpmu = ',i2)
 805  format(' jgpsp = ',i3,', ew = ',e14.7)
 806  format('     F = ',e14.7,', xnu = ',e14.7)
 807  format('     F = ',e14.7,', Delta (F) = ',e14.7)
c
c
c
c
c     slave nodes part.
      else if(myid.ne.master) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c     vertical zones part
c
      call z_surf_bcast
c
c     if there are more nodes than work skip this node.
      if(myid.gt.nz) goto 990
c
c     as long as the node doesn't recieve the end_signal, it
c     will keep tracking photons for the zones.
c      write(*,*) 'myid=', myid, ' before recv zone=', zone
 991  call z_surf_recv_job(zone)
      if(zone.eq.end_signal) goto 990
      call z_surf_calc(zone)
c      write(*,*) 'myid=', myid, ' js=', zone, ' after z_surf_calc'
      call z_surf_send_result(zone)
c      write(*,*) 'myid=', myid, ' zone=', zone, ' result was sent'
      goto 991
 990  continue
c      write(*,*) 'myid=', myid, ' z recieved end signal'
c      stop
c
c
c     radial zones part.
c
      call r_surf_bcast
c
c     if there are more nodes than work skip this node.
      if(myid.gt.nr) goto 995
c
c     as long as the node doesn't recieve the end_signal, it
c     will keep tracking photons for the zones.
c      write(*,*) 'myid=', myid, ' before recv zone=', zone
 996  call r_surf_recv_job(zone)
c      write(*,*) 'myid=', myid, ' zone=', zone, ' recved job', 
c     1     ' nsurfl(zone)=', nsurfl(zone)
      if(zone.eq.end_signal) goto 995
c      write(*,*) 'myid=', myid, ' before r_surf_calc, ncycle=', ncycle
      call r_surf_calc(zone)
c      if( (myid.eq.2).or.(myid.eq.3) ) write(*,*) 'myid=', myid, 
c     1     ' after r_surf_calc, zone=', zone
c      write(*,*) 'myid=', myid, ' ks=', zone, ' after r_surf_calc'
      call r_surf_send_result(zone)
c      if(myid.eq.2) write(*,*) 'myid=', myid, ' zone=', 
c     1     zone, ' result was sent'
      goto 996
 995  continue
c      write(*,*) 'myid=', myid, ' r recieved end signal'
c
c
c
c     end of slave part
      endif
c
c
c
c     end of imcsurf2d
c      write(*,*) 'myid=', myid, ' end imsurf2d ecens=', ecens(1,1)
 900  return
      end
c
c
c
c     This subroutine calculates and tracks the photons produced 
c     at the inner and outer surface boundaries of the zones.
c     J. Finke 4 April 2005
      subroutine z_surf_calc(js)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer js
c
      integer nxs
      integer npsurf
c
      double precision fibran, t_average, x1, x2
      double precision esurf
      double precision tpl, F(nphomax)
c
c
      call initialize_zrand(js)
c      write(*,*) 'js=',js,' rseed=',rseed
c      write(*,*) 'js=',js,' rand=',fibran()
c      write(*,*) 'js=',js,' rand_switch=',rand_switch
c
c      rseed = 8535
c         added time greater than t0 to make sure the external radiation begins after a certain time ! Xuhui 1/28/10
              if (tbbi(js, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti))
     1            call file_sp(i_fname(js, ti))

              if (nsurfi(js).le.0) goto 201
              npsurfz(js) = 0
              do 200 nxs = 1, nsurfi(js)
c                 write(*,*) 'js=', js, ' nxs=', nxs
                 jph = js
                 npsurfz(js) = npsurfz(js) + 1
c
c                 wmu = 2.d0*drand(0) - 1.d0
                 wmu = 2.d0*fibran() - 1.d0
                 if (wmu.gt.0.9999999999d0) wmu = 0.9999999999d0
                 if (wmu.lt.-0.9999999999d0) wmu = -0.9999999999d0
c
c                 phi = -1.1d1/7.d0 + 2.2d1/7.d0*drand(0)
                 phi = -1.1d1/7.d0 + 2.2d1/7.d0*fibran()
                 if (phi.lt.-1.5707963d0) phi = -1.57079063d0
                 if (phi.gt.1.5707963d0) phi = 1.5707963d0
                 rpre = rmin
                 if (jph.eq.1) then          
                      zpre = z(1)*fibran()
                 else
                      zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                 endif
                 ew = ewsurfi(jph)
                 esurf = esurf+ew
c                 dcen = drand(0)*c_light*dt(1)
                 dcen = fibran()*c_light*dt(1)
c
                 if (tbbi(jph,ti).gt.0.d0) then
                    tpl = tbbi(jph,ti)
                    call planck(tpl)
                 else
                    call file_sample()
                 endif
c
c                 F(jgpsp) = F(jgpsp) + ew
c                
                 kph = 1
c
                 call imctrk2d(-1)
c
 200          continue
 201          continue
c
              if (tbbo(js, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti))
     1           call file_sp(o_fname(js, ti))
              do 250 nxs = 1, nsurfo(js)
                  jph = js
                  npsurfz(js) = npsurfz(js) + 1
c                  wmu = 2.d0*drand(0) - 1.d0
                  wmu = 2.d0*fibran() - 1.d0
                  if (wmu.gt.0.9999999999d0) wmu = 0.9999999999d0
                  if (wmu.lt.-0.9999999999d0) wmu = -0.9999999999d0
c
                  rpre = r(nr)
                  if (jph.eq.1) then          
                      zpre = z(1)*fibran()
                  else
                      zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                  endif
c                  x1 = drand(0)
c                  x2 = drand(0)
                  x1 = fibran()
                  x2 = fibran()
                  if (x1.lt.0.5) then
                     phi = 1.1d1/7.d0 + (1.1d1/7.d0)*x2
                     if (phi.lt.1.57079638d0) phi = 1.57079638d0
                  else 
                     phi = -1.1d1/7.d0 - (1.1d1/7.d0)*x2
                     if (phi.gt.-1.57079638d0) phi = -1.57079638d0
                  endif 
                  ew = ewsurfo(jph)  
                  esurf = esurf+ew
c                  dcen = drand(0)*c_light*dt(1)
                  dcen = fibran()*c_light*dt(1)
c
                  if (tbbo(jph,ti).gt.0.d0) then
                     tpl = tbbo(jph,ti)
                     call planck(tpl)
                  else
                     call file_sample()
                  endif
c
c                  F(jgpsp) = F(jgpsp) + ew
c                 
                  kph = nr
c
                  call imctrk2d(-1)
c
 250           continue
c
c
c
c             write(*,*) 'end js=',js,' rand=',fibran()
               return
               end
c
c
c
c     This subroutine calculates and tracks the photons produced 
c     at the upper and lower surface boundaries of the zones.
c     J. Finke 4 April 2005
      subroutine r_surf_calc(ks)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer ks
c

      integer npsurf, nxs
c
      double precision psi
      double precision fibran, t_average, x1, x2
      double precision esurf
      double precision tpl, F(nphomax)

c
c

c
c
      call initialize_rrand(ks)

             if (tbbu(ks, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti) )
     1             call file_sp(u_fname(ks, ti))

             npsurfr(ks) = 0
             do 330 nxs = 1, nsurfu(ks)
                 kph = ks
                 npsurfr(ks) = npsurfr(ks) + 1
c                 
c                 wmu = -drand(0)
c                 phi = 4.4d1/7.d0*drand(0)
c                 
                 wmu = -fibran()
                 if (wmu.gt.0.9999999999d0) wmu = 0.9999999999d0
                 if (wmu.lt.-0.9999999999d0) wmu = -0.9999999999d0
c
                 phi = 2.d0*pi*fibran()
c                 psi = drand(0)
                 psi = fibran()
                 Zpre= z(nz)    
                 if (kph.eq.1) then
                       rpre = dsqrt((rmin**2) + psi*(r(kph)**2
     1                                             - rmin**2))
                 else
                       rpre = dsqrt((r(kph-1)**2) + psi*(r(kph)**2 
     1                                               - r(kph-1)**2))
                 endif            
                 ew = ewsurfu(kph)
                 esurf = esurf+ew
c
c                 dcen = drand(0)*c_light*dt(1)
                 dcen = fibran()*c_light*dt(1)
c
                 if (tbbu(kph,ti).gt.0.d0) then
                    tpl = tbbu(kph,ti)
                    call planck(tpl)
                 else
                    call file_sample()
                 endif
c
c                 F(jgpsp) = F(jgpsp) + ew
c
                 jph = nz
c
c                 write(*,*)
c                 write(*,*) 'Upper surface; particle no. ', npsurf
c                 write(*,*) 'xnu = ',xnu
c                 write(*,*) 'wmu = ',wmu
c                 write(*,*) 'phi = ',phi
c                 write(*,*) 'rpre = ',rpre
c                 write(*,*) 'zpre = ',zpre
c                 write(*,*) 'dcen = ',dcen
c                 write(*,*) 'ew = ',ew
c                 write(*,*) 'jph = ',jph
c                 write(*,*) 'kph = ',kph
c                 write(*,*) 'jgpsp = ',jgpsp
c                 write(*,*) 'jgplc = ',jgplc
c                 write(*,*) 'jgpmu = ',jgpmu
c
                 call imctrk2d(-1)
c
 330          continue
c
c              write(4,*) 'Lower boundary, zone ',ks,': ',nsurfl(ks),
c     1                   ' photons'
c              write(*,*) 'Lower boundary, zone ',ks,': ',nsurfl(ks),
c     1                   ' photons'
c
              if (tbbl(ks, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti))
     1                call file_sp(l_fname(ks, ti))

              do 360 nxs = 1, nsurfl(ks)
c                 write(*,*) 'nxs=', nxs
                 kph = ks
                 npsurfr(ks) = npsurfr(ks) + 1
c
c                  wmu = drand(0)
                 if(star_switch.eq.1) then
                    wmu = 9.9999999d-1
                 else
c                   all the beamed external radiation are in the up direction
c                    wmu = fibran()
                     wmu = 9.9999999d-1 ! Xuhui 2/2/10
                    if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                    if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                 endif
                  if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                  if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
c
c
c                 wv1 = drand(0)
c                 wv2 = drand(0)
c                 wmu = dmax1(wv1, wv2)
c
c                 phi = (4.4d1/7.d0)*drand(0)
                 phi = 2.d0*pi*fibran()  
                 psi = fibran()
                 Zpre = zmin
                 if (kph.eq.1) then
                       rpre = dsqrt((rmin**2)+ psi*(r(kph)**2
     1                                            - rmin**2))
                 else
                       rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2
     1                                              - r(kph-1)**2))
                 endif             
                 ew = ewsurfl(kph)
                 esurf = esurf+ew
c
                 if (tbbl(kph,ti).gt.0.d-20) then
                    tpl = tbbl(kph,ti)
                    call planck(tpl)
                 else
                    call file_sample()
                 endif
c                 dcen = drand(0)*c_light*dt(1)
                 dcen = fibran()*c_light*dt(1)
c
c                 F(jgpsp) = F(jgpsp) + ew
c
                 jph = 1
c
                 if(ncycle.eq.2) then
                 if((nxs.gt.1542).and.(nxs.lt.1545)) then
c     problem btwn 1543 and 1544
c                 if(nxs.gt.nsurfl(ks)-10) then
c                 if((ks.eq.2).and.(jph.eq.1))  then
c                 if(nxs.le.10) then
c                write(*,*)
c                write(*,*) 'Lower surface; particle no. ', npsurfr(ks)
c                write(*,*) 'nxs=',nxs
c                write(*,*) 'Lower surface; particle no. ', nxs
c                write(*,*) 'xnu = ',xnu
c                write(*,*) 'wmu = ',wmu
c                write(*,*) 'phi = ',phi
c                write(*,*) 'rpre = ',rpre
c                write(*,*) 'zpre = ',zpre
c                write(*,*) 'dcen = ',dcen
c                write(*,*) 'ew = ',ew
c                write(*,*) 'jph = ',jph
c                write(*,*) 'kph = ',kph
c                write(*,*) 'ks=', ks
c                write(*,*) 'jgpsp = ',jgpsp
c                write(*,*) 'jgplc = ',jgplc
c                write(*,*) 'jgpmu = ',jgpmu
c                write(*,*) 'rand=',fibran()
c                stop
c                write(*,*) 'Ed_ref=', Ed_ref(ks)
c                write(*,*) 'rseed=', rseed
             endif
             endif
c             endif
c
c
                 call imctrk2d(-1)
c              
 360       continue
c
c
c
c
      return
      end
c
c
c
c
c     Read external spectrum file;
c     calculate flux normalizations
c     and probability distributions
c
c
      subroutine file_sp(fname)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      character *30 fname
c
      integer,parameter :: ndisk = 10
       
c
      integer i, n_input, j, ijump
      double precision mu_disk(ndisk), dmu(ndisk)
      double precision r_d(ndisk),
     1  Ltot_disk, Ftot_blr, Ftot_ir, Ftot_blr_norm, Ftot_ir_norm
      double precision Id_file(nfmax,ndisk), F_como(nfmax,ndisk),
     1  L_disk(nfmax), F_blr(nfmax), F_ir(nfmax)
c
      double precision Isum
      double precision beta, dopp(ndisk)
c
c

c
      n_input = 25

      open(n_input, file=fname, status='unknown') ! Xuhui Chen 05/08
      i = 0
 100  i = i + 1
      read(unit=n_input, err=120, fmt=*) 
     1 E_file(i), L_disk(i), F_blr(i), F_ir(i)
c      L_disk(i)=L_disk(i)*max(1.d0,7.3d0*exp(-abs(time-8.d6)/2.d6)) ! Xuhui 1/13/10
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculation of radiative intensity from the accretion disk.
c     followed the calculation in Ghisellini & Tavecchio 2009.
c
c     dopp is the doppler factor of the photons travel from the
c     accretion disk to the jet.
      do j=1,ndisk
        r_d(j)=R_disk*(j-0.5)/ndisk
      enddo
      beta = sqrt(1.d0-1.d0/g_bulk**2)
      mu_disk(:) = d_jet/sqrt(d_jet**2+r_d(:)**2)
      do j=1,ndisk
         dmu(j) = d_jet/sqrt(d_jet**2+(R_disk*(j-1)/ndisk)**2)-
     1   d_jet/sqrt(d_jet**2+(R_disk*j/ndisk)**2)
      enddo
      dopp(:) = g_bulk*(1.d0-beta*mu_disk(:))
      Id_file(i,:)=L_disk(i)/(4*pi**2*R_disk**2*mu_disk(:))
      F_como(i,:)=2*pi*dopp(:)*Id_file(i,:)*dmu(:)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((E_file(i).gt.0.d0).and.(i.le.(nfmax-1))) goto 100
 120  close(n_input)
      nfile = i - 1

c     The following procedure has assumed E_file to be evenly spaced
c     in logarithm. It moves the new spectrum to energy grids of its own
c      F_file(:)=1.d-29
c      do j=1,ndisk
c       if(dopp(j).ge.1.d0)then
c         do ijump=0,nfmax-1
c            if(E_file(ijump+1).ge.E_file(1)*dopp(j))exit
c         enddo
c         do i=1,nfmax-ijump
c            F_file(i+ijump)=F_file(i+ijump)+F_como(i,j)
c         enddo
c       else
c         do ijump=0,nfmax-1
c            if(E_file(nfmax-ijump).le.E_file(nfmax)*dopp(j))exit
c         enddo
c         do i=1,nfmax-ijump
c            F_file(i)=F_file(i)+F_como(i+ijump,j)
c         enddo
c       endif
c      enddo

c     Ltot_disk is the total power of the disk radiation.
c     E_file is the median of every energy grid, so sqrt() is divided for
c     the actual power.
      Ltot_disk = 0.d0
      Ftot_blr = 0.d0
      Ftot_ir = 0.d0
      do i=1,nfmax-1
        Ltot_disk = Ltot_disk+L_disk(i)*(E_file(i+1)-E_file(i))
        Ftot_blr = Ftot_blr+F_blr(i)*(E_file(i+1)-E_file(i))
        Ftot_ir = Ftot_ir+F_ir(i)*(E_file(i+1)-E_file(i))
      enddo
      Ltot_disk = Ltot_disk/sqrt(E_file(2)/E_file(1))
c     Ftot_blr is the total power of the BLR radiation of the input file.
c     Ftot_ir is the total power of the torus radiation of the input file.
      Ftot_blr = Ftot_blr/sqrt(E_file(2)/E_file(1))
      Ftot_ir = Ftot_ir/sqrt(E_file(2)/E_file(1))
c     Ftot_blr_norm and Ftot_ir_norm are the total power the BLR and torus 
c     radiation should have according
c     to Ghisellini & Madau 1996 (timed c for flux)
      Ftot_blr_norm = 17.d0/48.d0/pi*g_bulk**2*fr_blr*Ltot_disk/R_blr**2
      Ftot_ir_norm =  1.d0/4.d0/pi*g_bulk**2*fr_ir*Ltot_disk/R_ir**2
c     F_blr in the end is the normalized BLR and torus flux.
      F_blr(:) = F_blr(:)/Ftot_blr*Ftot_blr_norm
      F_ir(:) = F_ir(:)/Ftot_ir*Ftot_ir_norm


c************** Here it decides which of the above calculation counts*******
c      F_file(:)=F_file(:)+F_blr(:)
      F_file(:)=F_blr(:)+F_ir(:)

      if (nfile.lt.2) then
         write(*,*) 'Error in input spectrum file:'
         write(*,*) 'Less than 2 lines of input read!'
         write(*,*)
         stop
      endif
c
c      do 130 i = 1, nfile
c         write(*,*) i, E_file(i), F_file(i)
c 130  continue
c
      Isum = 0.d0
      do 150 i = 1,nfile-1
         alpha(i) = log(F_file(i+1)/F_file(i))
     1             /log(E_file(i+1)/E_file(i))
         a1(i) = alpha(i) + 1.d0
         if (a1(i).gt.2.d1) a1(i) = 2.d1
         if (a1(i).lt.-2.d1) a1(i) = -2.d1
         if (dabs(a1(i)).lt.1.d-3) then
            I_file(i) = F_file(i)*E_file(i)
     1                 *log(E_file(i+1)/E_file(i))
         else
            I_file(i) = F_file(i)*E_file(i)
     1                *((E_file(i+1)/E_file(i))**a1(i)
     2                  - 1.d0)/a1(i)
         endif
         Isum = Isum + I_file(i)
         P_file(i) = Isum
 150  continue
c
      do 200 i = 1, nfile-1
         P_file(i) = P_file(i)/Isum
 200  continue

      int_file = Isum
c
      return
      end
c
c
c
c
c     Subroutine to sample a photon energy
c     from a user-specified input spectrum
c
c
      subroutine file_sample
      implicit none
      include 'general.pa'
      include 'commonblock.f'
c
c
      integer i
      integer jbot, jtop, jmid, m, n
c
      double precision hubot, hutop, Isum
      double precision fibran, x1, x2
c
c      write(*,*) 'in file_sample'
c      write(*,*) 'xnu=', xnu
c      write(*,*) 'nfile=', nfile
c
c      x1 = drand(0)
      x1 = fibran()
      do 210 i = 1, nfile-1
 210  if (P_file(i).gt.x1) goto 220
 220  continue
c
c      x2 = drand(0)
      x2 = fibran()
      xnu = E_file(i)*((a1(i)*I_file(i)*x2
     1     /(F_file(i)*E_file(i)) + 1.d0)**(1./a1(i)))
c
c
c
c     Determine photon group of new photon energy
c     in photon energy bin structure for spectrum
c
      jbot = 1
      jtop = nphtotal + 1
      hubot = 1.000001d0*hu(1)
      hutop =.999999d0*hu(jtop)
      if( xnu .ge. hutop ) then
c         xnu = hutop
c         jgp = numgps
         jgpsp = 0
         go to 300
      endif
      if( xnu .le. hubot ) then
c         xnu = hubot
c         jgp = 1
         jgpsp = 0
         go to 300
      endif
c
  140 continue
      jmid = (jbot+jtop)/2
      if( jmid .eq. jbot ) go to 160
      if(xnu .eq. hu(jmid)) go to 160
      if(xnu .lt. hu(jmid)) then
         jtop = jmid
      else
         jbot = jmid
      endif
      go to 140
c
  160 continue
      jgpsp = jmid
c
  300 continue
c
c
c     Determine photon group of new photon energy
c   in photon energy bin structure for light curves
c
      jgplc = 0
      do 400 m = 1, nph_lc
         if ((xnu.gt.Elcmin(m)).and.
     1       (xnu.le.Elcmax(m))) then
            jgplc = m
            goto 600
         endif
  400 continue
c
  600 continue
c
c
c    Determine angular bin of photon
c
      do 700 n = 1, nmu
         if (wmu.le.mu(n)) then
            jgpmu = n
            goto 900
         endif
  700 continue
      jgpmu = nmu
c
  900 continue      
c
      return
      end

c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:25:40 EDT 2006
c version: 2
c Name: J. Finke
c Changed call to 'initialize_rand' to 'initialize_zrand' and 'initialize_rrand'. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Sun Aug  6 18:51:03 EDT 2006
c version: 3
c Name: J. Finke
c Added possibility of using upper surface to represent 
c a star. Upper surface photons can come in 
c parallel to each other, representing the star, instead 
c of in random directions.     
