c     everything is done in master node. Then some of them are broadcasted to
c     slave nodes 
      subroutine setup
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      integer i, j, k, n, m, m1
      integer status(MPI_STATUS_SIZE)
c
      double precision delj, delk
      double precision dmu
      double precision dE, priE, secE
      double precision dist_max
c____________________________________________________________________
      double precision sum_E
      double precision inj_sum, inj_E, inject_ne(num_nt), inj_y,
     1                 gamma0(num_nt), vol_tot !local varible
c_____________________________________________________________________

c
c___________________________________________________________
c
c
c     most of setup is done by the master node only.
      if(myid.eq.master) then
c
c
c
c     Call imcdate to load Compton cross-section averages
c       and average energy exchange data into memory
c
c
      call imcdate
c
c
c      Calculate Compton reflection matrix
c
      call Pref_calc
      call Wref_calc
c

c
c        Initialize time and time step
c
      ncycle = 0 ! Xuhui 5/15/11
      time = 0.d0
c      dist_max = dmax1((2.d0*(r(nr) - rmin)), z(nz))
      dist_max = dmin1(r(nr)/nr,z(nz)/nz) ! Xuhui 5/11/09
      dt(1) = mcdt*dist_max/inj_v
      dt(2) = dt(1)
      write(4,*)
      write(4, 100) dt(1)
 100  format('Initial time step: ',e14.7,' s')
c
c
c      Calculate spatial zoning (linear)
c
      j = 1
      k = 1
      delj = z(nz)/dble(nz)
      delk= (r(nr) - rmin)/dble(nr)
      zmin = 0.d0
      z(1) = delj
      r(1) = rmin + delk
c
      do 120 j = 1, nz-1
 120  z(j+1) = z(j) + delj
c
      do 130 k = 1, nr-1
 130  r(k+1) = r(k) + delk
c
      write(4,*)
      do 133 j = 1, nz
 133  write(4,810) j,z(j)
      do 136 k = 1, nr
 136  write(4,815) k,r(k)
c
      vol(1,1) = pi*(r(1)**2 - rmin**2)*z(1)
      zsurf(1,1) = 2.d0*pi*((r(1) + rmin)*z(1) 
     1                     + (r(1)**2 - rmin**2))
      do 138 k = 2, nr
         vol(1,k) = pi*(r(k)**2 - r(k-1)**2)*z(1)
         zsurf(1,k) = 2.d0*pi*((r(k) + r(k-1))*z(1) 
     1                        + (r(k)**2 - r(k-1)**2))
 138  continue
      do 143 j = 2, nz
         vol(j, 1) = pi*(r(1)**2 - rmin**2)
     1                        *(z(j) - z(j-1))
         zsurf(j, 1) = 2.d0*pi*((r(1) + rmin)*(z(j) - z(j-1))
     1                         + (r(1)**2 - rmin**2))
         do 142 k = 2, nr
            vol(j,k) = pi*(r(k)**2 - r(k-1)**2)
     1                          *(z(j) - z(j-1))
            zsurf(j, k) = 2.d0*pi*((r(k) + r(k-1))
     1                             *(z(j) - z(j-1))
     2                            + (r(k)**2 - r(k-1)**2))
 142     continue
 143  continue
c
      Asurfu(1) = pi*(r(1)**2 - rmin**2)
      Asurfl(1) = Asurfu(1)
      Asurfi(1) = 2.d0*pi*rmin*z(1)
      Asurfo(1) = 2.d0*pi*r(nr)*z(1)
c
      do 1100 j = 2, nz
         Asurfi(j) = 2.d0*pi*rmin*(z(j) - z(j-1))
         Asurfo(j) = 2.d0*pi*r(nr)*(z(j) - z(j-1))
 1100 continue
      do 1110 k = 2, nr
         Asurfu(k) = pi*(r(k)**2 - r(k-1)**2)
         Asurfl(k) = Asurfu(k)
 1110 continue
c
c      write(*,*) 'Partial surfaces calculated.'
c
c
c       Calculate parameters of nonthermal distribution:
c
      sum_E = 0.d0
      do 1150 j = 1, nz
         do 1140 k = 1, nr
            f_pair(j,k) = 0.d0
            do 1130 i = 1, num_nt
               dne_pa(j, k, i) = 0.d0
               dnp_pa(j, k, i) = 0.d0
 1130       continue
c            write(*,*) 'j = ',j,'; k = ',k,'; calling gam_min'
            call gam_min(j, k, tea(j,k))
c            write(*,*) 'j = ',j,'; k = ',k,'; calling P_nontherm'
            call P_nontherm(j, k)
c            write(*,*) 'Returned from P_nontherm.'
             do i = 1, num_nt-1
                sum_E = sum_E + 
     1        f_nt(j,k,i)*(gnt(i+1)-gnt(i))*(gnt(i)+1)*n_e(j,k)*vol(j,k)
             enddo
 1140    continue
 1150 continue
      write(*,*) 'Initial Energy amount:', sum_E*8.186d-7
c
c      write(*,*) 'Nonthermal parameters calculated.'
c
c
c       Calculating the set-up angle for the problem 
c                     (linear)
c
        dmu = 2.d0/dble(nmu)
        mu(1) = dmu - 1.d0
c        
        do 145 n = 2, nmu
          mu(n) = mu(n-1)+dmu
 145    continue
c
        write(4,*)
        do 150 n = 1, nmu
 150    write(4,805) n, mu(n)
c
c
c  Calculating the Energy set-up grid for the problem 
c              (logarithmic)
c        
        i = 1
        do 160 m = 1, nphreg
           priE  = log(Ephmax(m)/Ephmin(m))
           secE = priE/dble(nphbins(m))
           dE    = exp(secE)
           hu(i) = Ephmin(m)
           do 170  m1=1, nphbins(m)
              i = i+1
              hu(i) = hu(i-1)*dE
  170      continue
  160   continue
c
c         write(4,*)
c         do 180 i = 1, nphtotal+1
c  180    write(4,800) i, hu(i)
c
c
c          Array initializations
c
         do 210 j = 1, nz
            do 200 k = 1, nr
               ecens(j,k) = 0.d0
               npcen(j,k) = 0
  200       continue
  210    continue
c
         do 220 k = 1, nr
            Ed_in(k) = 0.d0
            Ed_abs(k) = 0.d0
            Ed_ref(k) = 0.d0
  220    continue
c
c
c      Initialize arrays for internal photon fields
c        for gamma-gamma absorption calculations
c
         dE = dexp(log(1.d2)/dble(n_gg))
         E_gg(1) = 5.d1
         do 300 i = 1, n_gg
            if (i.gt.1) E_gg(i) = E_gg(i-1)*dE
            do 280 j = 1, nz
               do 260 k = 1, nr
                  n_ph(i, j, k) = 0.d0
                  k_gg(i, j, k) = 0.d0
 260           continue
 280        continue     
 300     continue
c
c
c
c      Initialize arrays for internal photon
c       fields for Compton loss calculations
c
c      dE = dexp(log(5.d5)/dble(nphfield))
      dE = dexp(log(1.d20)/dble(nphfield)) ! Xuhui
c      E_field(1) = 1.d-3
      E_field(1) = 1.d-10 ! Xuhui
      do i = 2, nphfield
        E_field(i) = E_field(i-1)*dE
      enddo
c
c
c     Calculate table of IC loss kernel values
c
      call IC_loss
c
c
c     Initialize spectral and light curve output arrays
c
         do 400 n = 1, nmu
            do 380 i = 1, nphtotal
 380          fout(n, i) = 0.d0
 400     continue
c
c
c     Initialize sums for temperature distribution output
c
         time_sum = 0.d0
         do 450 j = 1, nz
            do 440 k = 1, nr
               T_sum(j,k) = 0.d0
 440        continue
 450     continue
c   
c
c     end master-only part; master and slaves call setup_bcast.
      endif

      call setup_bcast
c
c   
c
c
c
c
c          Output formats for log file
c
  800    format('hu(',i3,') = ',e14.7,' keV')
  805    format('mu(',i3,') = ',e14.7)
  810    format('Vertical zone boundary z(',i2,') = ',e14.7,' cm')
  815    format('Radial zone boundary r(',i2,') = ',e14.7,' cm')
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:20:38 EDT 2006
c version: 2
c Name: J. Finke
c Fixed bug with T_const.     
