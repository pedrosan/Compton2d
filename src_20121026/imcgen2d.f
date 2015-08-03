      subroutine imcgen2d
c     this subroutine calculate some general stuff, like the radiative energy   
c     loss in volume, the total energy of the system, and photon numbers to be
c     produced this time step. 
c     All works are done in the master node
c     !comment by Xuhui Chen 7/8/08

      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f'
c
      double precision sigma
      parameter(sigma = 1.0267d24)
c
      integer i, j, k, npwant, nst_1st, nt, t, n
      integer n_new
      integer n_input
      integer syc_switch
c
      double precision sum_g_1, t_average
      double precision fac, bingo, fcens, gamma_R, fbias, bias
      double precision fas(jmax, kmax)
      double precision p_1, l_min, isy, Th, uB
      double precision g_av, Th_e, N_e_nt, f_rel_br
      double precision kgg_calc, int_sy, McDonald
      double precision gamma_bar
      double precision Eloss_pa

c
c
c
c        Formats
c
  10  format('Zone no. ',i2,',',i2)
  15  format('Eloss_sy = ',e14.7,' erg')
  20  format('Eloss_cy = ',e14.7,' erg')
  25  format('Eloss_th = ',e14.7,' erg')
  30  format('Eloss_br = ',e14.7,' erg (f_rel = ',e14.7,')')
  35  format('Eloss_pa = ',e14.7,' erg')
  40  format('Ingoing energy at inner boundary, zone ',i2,': ',e14.7,
     1       ' erg')
  41  format('Ingoing energy at outer boundary, zone ',i2,': ',e14.7,
     1       ' erg')
  42  format('Ingoing energy at upper boundary, zone ',i2,': ',e14.7,
     1       ' erg')
  43  format('Ingoing energy at lower boundary, zone ',i2,': ',e14.7, 
     1       ' erg')
  50  format('Eloss_total = ',e14.7,' erg')
  55  format('Total ingoing energy: bingo = ',e14.7,' erg')
  60  format('Photons produced at inner boundary, zone ',i2,': ',i7)
  61  format('Photons produced at outer boundary, zone ',i2,': ',i7)
  62  format('Photons produced at upper boundary, zone ',i2,': ',i7)
  63  format('Photons produced at lower boundary, zone ',i2,': ',i7)
  64  format('Photons produced in zone (',i2,',',i2,'): ',i7)
  70  format(e14.7,1x,e14.7)
  71  format(e14.7,1x,e14.7,1x,e14.7)
c 
  90  format('Boundary temp. at upper boundary, zone ',i2,': ',
     1        e14.7,' keV')
  91  format('Boundary temp. at lower boundary, zone ',i2,': ',
     1        e14.7,' keV')
  92  format('Boundary temp. at inner boundary, zone ',i2,': ',
     1        e14.7,' keV')
  93  format('Boundary temp. at outer boundary, zone ',i2,': ',
     1        e14.7,' keV')
  94  format('Ed_abs = ',e14.7,' ergs; Ed_ref = ',e14.7,' ergs.')
c
c
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
c
c     master node
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(myid.eq.master) then
c
c
c
c
c
      n_input = 19
c
c
c       Initilaizing arrays
c
c      ndxout = 0
      do 100 j = 1, nz
         do 110  k = 1, nr
            npcen(j,k) = 0
            fas(j,k) = 0.d0
            edep(j,k) = 0.d0
            if (ncycle.eq.0) then
               ec_old(j,k) = 0.d0
            else
               ec_old(j,k) = ecens(j,k)
            endif
            prdep(j,k) = 0.d0
 110     continue
 100  continue 
c
      do 120 n = 1, nmu
         do 115 i = 1, nph_lc
 115        edout(n, i) = 0.d0
 120  continue
c
      bingo = 0.d0 
c
c
c     Calculate (time-dependent) surface emission
c
      if (ncycle.eq.0) then
         t = 1
         goto 135
      endif
      t_average = time + 5.d-1*dt(1)
c
      do 130 t = 1, ntime
         if (t1(t).gt.t_average) goto 135
 130  continue
 135  continue
 
 
c
c
      do 140 j = 1, nz
c         added time greater than t0 to make sure the external radiation begins after a certain time ! Xuhui 1/28/10
          if (tbbi(j,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t))then
             call file_sp(i_fname(j,t))
             erini(j) = dt(1)*Asurfi(j)*int_file
             write(4, 40) j, erini(j)
          else
             erini(j) = dt(1)*Asurfi(j)*sigma*(tbbi(j, t)**4.d0)
             if (tbbi(j,t).gt.0.d0) then
                write(4, 92) j, tbbi(j, t)
                write(4, 40) j, erini(j)
             endif
          endif
c
          if (tbbo(j,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t)) then
             call file_sp(o_fname(j,t))
             erino(j) = dt(1)*Asurfo(j)*int_file
             write(4, 41) j, erino(j)
          else
             erino(j) = dt(1)*Asurfo(j)*sigma*(tbbo(j, t)**4.d0)
             if (tbbo(j,t).gt.0.d0) then
                write(4, 93) j, tbbo(j, t)
                write(4, 41) j, erino(j)
             endif
          endif
          erlko(j) = 0.d0
          erlki(j) = 0.d0
 140  continue    
c
      do 150 k = 1, nr          
          Ed_abs(k) = Ed_in(k) - Ed_ref(k)
          if (tbbu(k,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t)) then
             call file_sp(u_fname(k,t))
             erinu(k) = dt(1)*Asurfu(k)*int_file
             write(4, 42) k, erinu(k)
          else
             erinu(k) = dt(1)*Asurfu(k)*sigma*(tbbu(k, t)**4.d0)
             if (star_switch.eq.1) 
     1            erinu(k) = erinu(k)*(Rstar/dist_star)**2.d0
c            if weight the luminosity by the radius of the star
c            and the distance to the star.
c            J. Finke 20 July 2006.
             if (tbbu(k,t).gt.0.d0) then
                write(4, 90) k, tbbu(k,t)
                write(4, 42) k, erinu(k)
             endif

          endif
c
          if (tbbl(k,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t)) then
             call file_sp(l_fname(k,t))
             erinl(k) = dt(1)*Asurfl(k)*int_file
             write(4, 43) k, erinl(k)
          else
             erinl(k) = dt(1)*Asurfl(k)*sigma*(tbbl(k, t)**4.d0)
             if (dh_sentinel.eq.1) then
                erinl(k) = erinl(k) + Ed_abs(k)*dt(1)/dt(2)
                if (tbbl(k,t).lt.1.d-20) erinl(k) = 0.d0
                write(4, 94) Ed_abs(k), Ed_ref(k)
             endif
             if (tbbl(k,t).gt.0.d0) then
                write(4, 91) k, tbbl(k,t)
                write(4, 43) k, erinl(k)
             endif
          endif
          erlku(k) = 0.d0
          erlkl(k) = 0.d0
c
 150  continue
c
c      Calculate volume emissivities etc.
c
c
      open(16, file='n_ph1.dat', status='unknown')
      open(17, file='nph1_smooth.dat', status='unknown')
      open(18, file='n_ph2.dat', status='unknown')
      open(19, file='nph2_smooth.dat', status='unknown')
c
      do 400 j = 1, nz
         do 390 k = 1, nr
c
            Eloss_cy(j,k) = 0.d0
            Eloss_sy(j,k) = 0.d0
            Eloss_br(j,k) = 0.d0
            Eloss_th(j,k) = 0.d0
            Eloss_pa = 0.d0
            Eloss_tot(j,k) = 0.d0
c
            write(4, 10) j,k
c
c          Calculate the Magnetic Field 
            if (ep_switch(j,k).eq.1) then
               Th = 1.957d-3*tea(j,k)
               if (Th.lt.1.d-2) then
                  uB = 1.5d0*Th + 7.5d0*(Th**2.d0)
               else
                  uB = McDonald(3.d0, (1.d0/Th))
     1                /McDonald(2.d0, (1.d0/Th)) - Th - 1.d0
               endif
               uB = uB*n_e(j,k)*8.176d-7*(1. + 2.d0*f_pair(j,k))
               B_field(j,k) = dsqrt(2.513d1*uB)
            else if (ep_switch(j,k).eq.2) then
               Th = 1.066d-6*tna(j,k)
               if (Th.lt.1.d-2) then
                  uB = 1.5d0*Th + 7.5d0*(Th**2.d0)
               else
                  uB = McDonald(3.d0, (1.d0/Th))
     1                /McDonald(2.d0, (1.d0/Th)) - Th - 1.d0
               endif
               uB = uB*n_e(j,k)*1.5d-3
               B_field(j,k) = dsqrt(2.513d1*uB)
            endif            
c
            if ((j.eq.1).and.(k.eq.1)) then
               l_min = dmin1(z(1), (r(1) - rmin))
            else if (j.eq.1) then
               l_min = dmin1(z(1), (r(k) - r(k-1)))
            else if (k.eq.1) then
               l_min = dmin1((z(j) - z(j-1)), (r(1) - rmin))
            else
               l_min = dmin1((z(j) - z(j-1)), (r(k) - r(k-1)))
            endif
c            
            call volume_em(j, k, tea(j,k), n_e(j,k), B_field(j,k),l_min)
            Th_e = tea(j,k)/5.11d2
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        syc_switch = 0
c        syc_switch 1 means syc using power-law fit. 0 means without power-law fit.
        if(syc_switch.eq.1)then
c           use power law fit to calculate synchrotron loss
            if (amxwl(j,k).gt.9.9999d-1) then
               Eloss_sy(j,k) = 0.d0
               goto 200
            endif

         p_1 = 1. - p_nth(j,k)
         if (dabs(p_1).gt.1.d-3) then
            N_e_nt = (1.d0 - amxwl(j,k))*p_1*n_e(j,k)
     1              /(gmax(j,k)**p_1 - gmin(j,k)**p_1)
         else
            N_e_nt = (1.d0 - amxwl(j,k))*n_e(j,k)
     1              /log(gmax(j,k)/gmin(j,k))
         endif

         g_av = gamma_bar(Th_e)
         gamma_R = 2.1d-3*sqrt(n_e(j,k))/(B_field(j,k)*dsqrt(g_av))

         isy = int_sy(gmin(j,k), gmax(j,k), p_nth(j,k), gamma_R)
         Eloss_sy(j,k) = 1.05838d-15*dt(1)*N_e_nt*(B_field(j,k)**2.d0)
     1      *vol(j,k)*isy
         write(*,*) 'j,k=',j,k,' Eloss_sy=',Eloss_sy(j,k)
         write(*,*) 'dt1=',dt(1),' N_e_nt=',N_e_nt
       else
c      use actual electron distribution to calculate synchrotron loss
          sum_g_1 = 0.d0
            do 170 i=1, num_nt-1
  170       sum_g_1 = sum_g_1 + ((gnt(i)+1.d0)**2.d0 - 1.d0)*f_nt(j,k,i)
     1                         *(gnt(i+1) - gnt(i))
          Eloss_sy(j,k) = 1.058d-15*n_e(j,k)*dt(1)*B_field(j,k)**2.d0
     1                    *sum_g_1*vol(j,k)
       endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Xuhui syc
c
 200     continue
c       write(*,*)'gmax=',gmax(j,k),'gmin=',gmin(j,k),'slope=',p_nth(j,k)
         if (Eloss_sy(j,k).gt.1.d-20) write(4,15) Eloss_sy(j,k)
c
         Eloss_cy(j,k) = dt(1)*vol(j,k)*Eloss_cy(j,k)
c
         if (Eloss_cy(j,k).gt.1.d-20) write(4,20) Eloss_cy(j,k)
c
         Eloss_th(j,k) = dt(1)*zsurf(j,k)*Eloss_th(j,k)
c
         if (Eloss_th(j,k).gt.1.d-20) write(4,25) Eloss_th(j,k)
c
         if (Th_e.gt.0.1d0) then
            f_rel_br = 1.41d0*dsqrt(Th_e)*(dlog(2.d0*Th_e) + 9.228d-1)
     1               - 1.d0
            f_rel_br = 1.d0 + (Th_e**2.d0)*f_rel_br/(1.d0 + 
     1           (Th_e**2.d0))
            if (f_rel_br.lt.1.d0) f_rel_br = 1.d0
         else
            f_rel_br = 1.d0
         endif
c
         Eloss_br(j,k) = 5.34d-24*vol(j,k)*dt(1)*amxwl(j,k)
     1                *sqrt(tea(j,k))*f_rel_br*(2.828d0*f_pair(j,k) 
     2                + 1.d0)*(n_e(j,k)**2.d0)         
c
         if (Eloss_br(j,k).gt.1.d-20) write(4,30) Eloss_br(j,k), 
     1                                            f_rel_br
c
         if (f_pair(j,k).lt.1.d-20) then
            Eloss_pa = 1.223d-20*vol(j,k)*dt(1)*f_pair(j,k)
     1                *(1.d0 + f_pair(j,k))*(n_e(j,k)**2.d0)
     2                /(1.d0/(1.d0 + 6.d0*Th_e)
     3              + Th_e/(dlog(1.123d0*Th_e + 1.d0) + 2.5d-1))
         endif
c
         if (Eloss_pa.gt.1.d-20) write(4,35) Eloss_pa
c
c
c         Eloss_tot(j,k) = Eloss_sy(j,k) + Eloss_br(j,k) + Eloss_cy(j,k) 
c     1                  + Eloss_th(j,k) + Eloss_pa
         Eloss_tot(j,k) = Eloss_sy(j,k)
c        Xuhui 6/1/09 deactivated any emission except synchrotron.
c
         fas(j,k) = Eloss_tot(j,k)
         write(4,50) Eloss_tot(j,k)
c
c          Check if blob is inside current zone;
c              if so, add blob emission
c
c         if ((r_blob.lt.xq(j,k)).and.(r_blob.gt.xq(j-1,k))) then
c            kT_bl = bc(3)*((xq(jm)/r_blob)**pT_b)
c            R_bl = h_r*xq(jm)
c            L_bl = 1.29d25*(kT_bl**4)*(R_bl**2)
c            E_bl = L_bl*1.d-8*dt(1)
c            fas(j) = Eloss_tot(j) + E_bl
c            p_blob(j) = E_bl/fas(j)
c
c
c         else
c            p_blob(j) = 0.d0
c            fas(j) = Eloss_tot(j)
c         endif
c
c
         if (pair_switch.eq.0) goto 385
c
c         write(*,*) 'Normalizing photon density spectra ...'
c
         do 370 i = 1, n_gg-1
             n_ph(i, j, k) = n_ph(i, j, k)
     1                      /(vol(j, k)*(E_gg(i+1) - E_gg(i)))
             if (j.eq.1.and.k.eq.1) then
                write(16,70) E_gg(i), dmax1(1.d-10, n_ph(i, j, k))
             else if (j.eq.2.and.k.eq.1) then
                write(18,70) E_gg(i), dmax1(1.d-10, n_ph(i, j, k))
             endif
  370    continue
c
c         write(*,*) 'Smoothing photon density spectra ...'
c
         call nph_smooth(j, k, tea(j, k))
c
c
c       Calculate pair production opacities and rates
c
c         write(*,*) 'Calculating gamma-gamma opacities ...'
c
         do 380 i = 1, n_gg-1
             k_gg(i, j, k) = kgg_calc(E_gg(i), j, k)
             if (k_gg(i, j, k).lt.1.d-50) k_gg(i, j, k) = 0.d0
c
             if (j.eq.1.and.k.eq.1) then
                write(17, 71) E_gg(i), dmax1(1.d-10, n_ph(i, j, k)),
     1                        dmax1(1.d-10, k_gg(i, j, k))
             else if (j.eq.2.and.k.eq.1) then
                write(19, 71) E_gg(i), dmax1(1.d-10, n_ph(i, j, k)),
     1                        dmax1(1.d-10, k_gg(i, j, k))
             endif
 380     continue
c
c         write(*,*) 'Calculating pair production ...'
c
         call pairprod(j,k)        
c
c         write(*,*) 'Pair production rates done.'
c
 385     continue
 390     continue
 400  continue
c
      close(16)
      close(17)
      close(18)
      close(19)
c
c
c      Calculate total input energy
c
      Emiss_old = Emiss_tot
      if(ncycle.eq.0)Emiss_old = 0.d0
      Emiss_tot = 0.d0
      do 420 j = 1, nz
         do 410 k = 1, nr
            bias = 1.d0
            bingo = bingo + bias*(ecens(j,k) + fas(j,k))
            Emiss_tot = Emiss_tot + fas(j,k)
  410    continue
  420 continue
c
      do 430 j = 1, nz
  430 bingo = bingo + erini(j) + erino(j)
      do 440 k = 1, nr
  440 bingo = bingo + erinu(k) + erinl(k)
      write(4,55) bingo
c
c
      nst=nst
c
      n_new = 0
c
      do 700 j = 1, nz
          nsurfi(j) = 0
          nsurfo(j) = 0
          if(tbbi(j,t).lt..0)nsurfi(j) = nst/nz
          if(tbbo(j,t).lt..0)nsurfo(j) = nst/nz
          n_new = n_new + nsurfi(j) + nsurfo(j)
 700  continue
       
      do 710 k = 1, nr
          nsurfu(k) = 0
          nsurfl(k) = 0
          if(tbbu(k,t).lt..0)nsurfu(k)=nst*(r(k)**2-r(k-1)**2)/r(nr)**2
          if(tbbl(k,t).lt..0)nsurfl(k)=nst*(r(k)**2-r(k-1)**2)/r(nr)**2
          n_new = n_new + nsurfu(k) + nsurfl(k)
 710   continue      
       
       do 720 j = 1, nz
           do 730 k = 1, nr
             nsv(j,k) = 0.5d0*nst*fas(j,k)/Emiss_tot !(r(k)**2-r(k-1)**2)/r(nr)**2/nz
                         ! Xuhui new way to distribute MC particles 4/1/11 Xuhui
             n_new = n_new + nsv(j, k)
             if (nsv(j,k).gt.0) then
                ewsv(j,k) = fas(j,k)/dble(nsv(j,k))
             else
                ewsv(j,k) = 0.d0
             endif
 730      continue
 720  continue
c
c Calculating the surface energy weights
c
c
       do 800 j = 1,nz
          if (nsurfi(j).gt.0) then
             ewsurfi(j) = erini(j)/dble(nsurfi(j))
          else
             ewsurfi(j) = 0.d0
          endif
          if (nsurfo(j).gt.0) then
             ewsurfo(j) = erino(j)/dble(nsurfo(j)) 
          else
             ewsurfo(j) = 0.d0
          endif
  800  continue                    
                      
       do 810 k =1,nr                   
          if (nsurfu(k).gt.0) then
             ewsurfu(k) = erinu(k)/dble(nsurfu(k))
          else
             ewsurfu(k) = 0.d0
          endif
          if (nsurfl(k).gt.0) then
             ewsurfl(k) = erinl(k)/dble(nsurfl(k))
          else 
             ewsurfl(k) = 0.d0
          endif
  810  continue
c
c    Introduce "bias" correction if more than
c     10*nst particles would be produced
c
      if (n_new.gt.(10*nst)) then
         fbias = dble(10*nst)/dble(n_new)
c      if ((n_new).gt.nst) then
c         fbias = dble(nst)/dble(n_new) ! Xuhui
         do 815 j = 1, nz
            nsurfi(j) = int(dble(nsurfi(j))*fbias)
            nsurfo(j) = int(dble(nsurfo(j))*fbias)
            ewsurfi(j) = ewsurfi(j)/fbias
            ewsurfo(j) = ewsurfo(j)/fbias
 815     continue
         do 820 k = 1, nr
            nsurfu(k) = int(dble(nsurfu(k))*fbias)
            nsurfl(k) = int(dble(nsurfl(k))*fbias)
            ewsurfu(k) = ewsurfu(k)/fbias
            ewsurfl(k) = ewsurfl(k)/fbias
 820     continue
c
        do 830 j = 1, nz
           do 825 k = 1, nr
              write(*,*)'fbias=',fbias ! Xuhui
c              nsv(j, k) = int(dble(nsv(j, k)/1000)*fbias)*1000
               nsv(j,k) = int(nsv(j,k)*fbias)
              ewsv(j, k) = ewsv(j, k)/fbias
 825       continue
 830    continue
c         
      endif
      do 835 j = 1, nz
         if (nsurfi(j).gt.10) write(4,60) j, nsurfi(j)
         if (nsurfo(j).gt.10) write(4,61) j, nsurfo(j)
 835  continue
      do 836 k = 1, nr
         if (nsurfu(k).gt.10) write(4,62) k, nsurfu(k)
         if (nsurfl(k).gt.10) write(4,63) k, nsurfl(k)
 836  continue  
c
      do 838 j = 1, nz
         do 837 k = 1, nr
            if (nsv(j,k).gt.0) write(4,64) j, k, nsv(j,k)
 837     continue
 838  continue
c
c
      do 850 j = 1, nz
        do 840 k = 1, nr
           ecens(j,k) = 0.d0
           do 839 i = 1, nphfield
 839       n_field(i, j, k) = 0.d0
 840    continue
 850  continue
c
      do 860 k = 1, nr
         Ed_abs(k) = 0.d0
         Ed_in(k) = 0.d0
         Ed_ref(k) = 0.d0
 860  continue
      do i=1,num_nt
        E_IC(i)=0.d0
      enddo
c
      call seed_zone(rseed)


c
c     slave node
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
c
      do 900 j=1, nz
          erlko(j) = 0.d0
          erlki(j) = 0.d0
 900  continue
      do 901 k=1, nr
          erlku(k) = 0.d0
          erlkl(k) = 0.d0
 901  continue
c
c      ndxout = 0
      do 851 j = 1, nz
        do 841 k = 1, nr
            npcen(j,k) = 0
            fas(j,k) = 0.d0
            edep(j,k) = 0.d0
            ecens(j,k) = 0.d0
           do 842 i = 1, nphfield
 842          n_field(i, j, k) = 0.d0
 841       continue
 851    continue
c
      do 716 n = 1, nmu
         do 715 i = 1, nph_lc
 715        edout(n, i) = 0.d0
 716     continue
c
      do 861 k = 1, nr
         Ed_abs(k) = 0.d0
         Ed_in(k) = 0.d0
         Ed_ref(k) = 0.d0
 861  continue
c
      do 210 j=1, nz
          erlko(j) = 0.d0
          erlki(j) = 0.d0      
 210   continue
      do 220 k=1, nr
          erlku(k) = 0.d0
          erlkl(k) = 0.d0
 220   continue
      do i=1,num_nt
        E_IC(i) =0.d0
      enddo
c
c
c
c     end of imcgen2d
      endif
      return 
      end
c
c
c
c
c       Subroutine to integrate the synchrotron
c           emissivity over photon energy
c
c
      double precision function int_sy(gmin, gmax, p, g_R)
      double precision gmin, gmax, p, g_R
c
      double precision sum, g, gs, dg, d, s, sd, x, y
c
      dg = 1.d-2
      d = 1.d0 + dg
      s = 1.d0 + 5.d-1*d
      sum = 0.d0
      g = gmin
      y = 2.d0 - p
c
 100  gs = g*s
      x = -g_R/gs
      if (x.gt.-1.d2) then
         sd = (gs**y)*dexp(x)
         sum = sum + g*dg*sd
      endif
      g = g*d
      if (g.lt.gmax) goto 100
c
      int_sy = sum
      return
      end
c
c 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:21:19 EDT 2006
c version: 2
c Name: J. Finke
c Fixed bug with T_const.     
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Sun Aug  6 18:31:12 EDT 2006
c version: 3
c Name: J. Finke
c Added possibility of using upper surface to represent 
c a star. Weighs the calculation of erinu by 
c star parameters.       
