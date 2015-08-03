c
      subroutine imcleak
      implicit none
      include 'general.pa'
      include 'commonblock.f'
c
c
c     NAME OF ROUTINE...imcleak
c
c     AUTHOR OF ROUTINE...E. H. Canfield July 27, 1983
c
c     PURPOSE...
c
c        Called by imctrk when particle hits inner or outer surface, or
c        by imcfield if hydro motion sweeps either surface past particle.
c
c        default mode:
c        Particles at the outer surface escape the system, and all leakage
c        and pinex/pipe data are updated.  Particles at the inner surface
c        are kept for tracking across the void (except for slabs).
c
c        Options allow specular reflection at surfaces, and killing particles
c        entering an interior void.
c
c        For reflection at outer surface, no energy escapes, but leakage
c        and pinex/pipe arrays are updated to show internal spectral data.
c
c     GLOSSARY OF IMPORTANT VARIABLES...
c
c        izbnd = 1  particle is at xq(1)
c              = jm particle is at xq(jm).
c        ibrad = 1  xq(1 ) is a reflecting surface.
c              = 2  xq(jm) is a reflecting surface.
c              = 3  both are reflecting surfaces.
c        lkin  = 1  particle is to be killed on entering interior void.
c        idead = 1  set if particle leaked from system.
c              = 3  set if particle was reflected.
c              = 4  set if particle is to be tracked across void.
c
c     LIST CHANGES HERE...
c
cccccc
c 09/02/88 spectral "enhancement" code for dermer.
c      accumulate data for <xnu> and <xnu**2>
c 11/23/04 implicit none added at begining.  J.D. Finke
cccccc
c/c
c
c
      integer ef_switch, n_in, n_out, i_0, i_1,
     1        jbot, jtop, jmid
      integer n, m
c     To calculate the spectrum incident on the upper and lower
c     boundaries, such as in photon bubble simulations, set
c     spec_switch to 1.  Otherwise, to calculate the spectrum
c     that leaves the region, set it to 0.
c     J. Finke 7/22/06
c
      double precision rnum, ewref, ewnew
      double precision hubot, hutop
      double precision f
      double precision xnu_sv
      double precision fibran
c


c
c
      xnu_sv = xnu
c
      if (kph.eq.0) then
c
c        Particle is at inner r-boundary:
c           If rmin > 0, absorb photon,
c         else continue on the other side
c
         if (rmin.gt.1.d-10) then
            idead = 1
            erlki(jph) = erlki(jph) + ew
c            incounter = incounter + 1
            goto 900
         else
            idead = 0
            phi = 1.d-6
            kph = 1
            goto 900
         endif
      endif
c
      if(jph .gt. 0) go to 100
c
c     particle is at lower surface:
c         Test for reflection.
c
      if (tbbl(kph,ti).gt.0.d0) then
         Ed_in(kph) = Ed_in(kph) + ew
         erlkl(kph) = erlkl(kph) + ew
         if (kph.eq.2) then
c         write(*,*) 'ew=',ew,' k=',k,' Ed_in=',Ed_in(k)
c         stop
         endif
      endif
c
c       Compton reflection at lower boundary
c       (cr_sent = 1 or cr_sent = 3)
c
      if ((cr_sent.eq.1).or.(cr_sent.eq.3).or.(cr_sent.eq.4)) then
c         dE_in = dE_in + ew
c
c      write(*,*) 'Compton reflection at lower boundary.'
c      write(*,*) 'xnu = ',xnu
c
         if ((tbbl(kph,ti).le.0.d0).or.(cr_sent.eq.4)) then
c            write(*,*) 'b4 refl wmu=', wmu
            wmu = -wmu
            if ( (jgpsp.gt.0).and.(spec_switch.eq.1) )
     1         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
c            write(*,*) 'j=', j, ' k=', k, ' refl.'
c            write(*,*) ' wmu=', wmu
c            stop
            goto 60
         endif
c
         n_in = 1
 10      if ((E_ref(n_in).lt.xnu).and.(n_in.lt.n_ref)) then
            n_in = n_in + 1
            goto 10
         endif
c
c         write(*,*) 'n_in = ',n_in
c
         n_out = 1
c         rnum = drand(0)
         rnum = fibran()
 20      if ((P_ref(n_out, n_in).lt.rnum).and.(n_out.lt.n_ref)) then
            n_out = n_out + 1
            goto 20
         endif
c
c         write(*,*) 'n_out = ',n_out
c
         if (n_out.gt.1) then
c            xnu = E_ref(n_out-1) 
c     1          + drand(0)*(E_ref(n_out) - E_ref(n_out-1))
            xnu = E_ref(n_out-1) 
     1          + fibran()*(E_ref(n_out) - E_ref(n_out-1))
         else
            xnu = E_ref(n_out)
         endif
c
c         write(*,*) 'xnu_new = ',xnu
c         write(*,*) 'W_abs = ',W_abs(n_out, n_in)
c
         ewnew = ew*W_abs(n_out, n_in)*xnu/xnu_sv
         Ed_ref(kph) = Ed_ref(kph) + ewnew
         wmu = dabs(wmu)
         ew = ewnew
c
c
 60      continue
         call get_bin
c
         idead = 0
         jph = 1
         goto 900
c
      else
c     no reflection, photon exits.
c     Added by J. Finke, 7 August 2006.
         if (ncycle.gt.0) then
            write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu, phi
            if (jgplc.gt.0) 
     1           edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
            if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )
     1           fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
         endif
         idead = 1
         goto 900
      endif
c
 105  format(6(e14.7,1x),e14.7)
c
c
c         particle is at outer surface.
c    Update all spectral and light curve data
c
  100 continue
c
      if (jph.eq.nz+1) goto 500
c
c       Particle is at outer r-boundary:
c         Test for Compton reflection
c
      erlko(jph) = erlko(jph) + ew
c      outcounter = outcounter + 1
c
      if ((wmu.gt.0.).or.(cr_sent.eq.0).or.
     1    (cr_sent.eq.1).or.(cr_sent.eq.4)) then
c
c             Write photon 
c       to spectrum and event file
c
         t_bound = time + dt(1) - rad_cp*dcen
         if (ncycle.gt.0) then
            write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu, phi
            if (jgplc.gt.0) 
     1         edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
            if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )
     1         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
         endif
         idead = 1
         goto 900
c
c
c
      else
c
c      Compton reflection at outer disk 
c
c         dE_in = dE_in + ew*f_ref
c
         n_in = 1
 110     if ((E_ref(n_in).lt.xnu).and.(n_in.lt.n_ref)) then
            n_in = n_in + 1
            goto 110
         endif
         n_out = 1
c         rnum = drand(0)
         rnum = fibran()
 120     if ((P_ref(n_out, n_in).lt.rnum).and.(n_out.lt.n_ref)) then
            n_out = n_out + 1
            goto 120
         endif
c
         if (n_out.gt.1) then
c            xnu = E_ref(n_out-1) 
c     1          + drand(0)*(E_ref(n_out) - E_ref(n_out-1))
            xnu = E_ref(n_out-1) 
     1          + fibran()*(E_ref(n_out) - E_ref(n_out-1))
         else
            xnu = E_ref(n_out)
         endif
c
c
         ewref = ew*W_abs(n_out, n_in)*xnu/xnu_sv
         ew = ewref
         if (dabs(wmu).gt.1.d-6) then 
            t_bound = time + dt(1) 
     1              - rad_cp*(dcen + zpre/wmu)
            zpre = 0.d0
            f = zpre*dsqrt(1. - wmu**2)/dabs(wmu)
            rpre = dsqrt(rpre**2 + f**2 + 2.d0*rpre*f*cos(phi))
         else 
            t_bound = 1.d20
         endif
c         wmu = drand(0)
         wmu = fibran()
c
c
 180     continue
         call get_bin
c
         if (ncycle.gt.0) then
            write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu, phi
            if (jgplc.gt.0) 
     1         edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
            if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )
     1         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
         endif
c
         idead = 1
         goto 900
c
c
      endif
c
c
c      Photon is at upper z surface.
c       Write photon to event file
c
c
 500  continue

c     Reflection at upper boundary added
c     J. Finke, 28 July 2005
c      if(upper_sent.eq.1) then
c         wmu = -wmu
c         if ( (jgpsp.gt.0).and.(spec_switch.eq.1) )
c     1        fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
c         do 270 n = 1, nmu
c            if(wmu.le.mu(n)) then
c               jgpmu = n
c               goto 180
c            endif
c 270     continue
c         jgpmu = nmu
c      endif
c      Xuhui deleted this because the goto command is getting warning from the
c      compiler

      erlku(kph) = erlku(kph) + ew
      t_bound = time + dt(1) - rad_cp*dcen
      if (ncycle.gt.0.and.wmu.lt.0.98) then
          write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu, phi
         ! wmu<0.98 to remove upper boudary record of large amount of external
         ! radiation coming out there. Xuhui 5/28/11
          if (jgplc.gt.0) 
     1         edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
          if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )
     1         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
      endif
c
      idead = 1
c
c
c
  900 continue
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c
c     This subroutine determines the spectral, lightcurve, and
c     anglular bin a photon should be in.
c     J. Finke 2 January 2007
      subroutine get_bin
      implicit none
      include 'general.pa'
      include 'commonblock.f'
c
      integer ef_switch, n_in, n_out, i_0, i_1,
     1        jbot, jtop, jmid
      integer m, n
      double precision hubot, hutop
c
c

c
c     Determine photon group of new photon energy
c     in photon energy bin structure for spectrum
c
      jbot = 1
      jtop = nphtotal + 1
      hubot = 1.000001*hu(1)
      hutop =.999999*hu(jtop)
      if( xnu .ge. hutop ) then
         jgpsp = 0
         go to 30
      endif
      if( xnu .le. hubot ) then
         jgpsp = 0
         go to 30
      endif
c
 24   continue
      jmid = (jbot+jtop)/2
      if( jmid .eq. jbot ) go to 26
      if(xnu .eq. hu(jmid)) go to 26
      if(xnu .lt. hu(jmid)) then
         jtop = jmid
      else
         jbot = jmid
      endif
      go to 24
c     
 26   continue
      jgpsp = jmid
c     
 30   continue
c
c
c     Determine photon group of new photon energy
c   in photon energy bin structure for light curves
c
      jgplc = 0
      do 40 m = 1, nph_lc
         if ((xnu.gt.Elcmin(m)).and.
     1        (xnu.le.Elcmax(m))) then
            jgplc = m
            goto 60
         endif
 40   continue
c
 60   continue
c
c
c    Determine angular bin of photon
c
      do 70 n = 1, nmu
         if (wmu.le.mu(n)) then
            jgpmu = n
            goto 80
         endif
 70   continue
      jgpmu = nmu
c
 80   continue
c     
c     
c
      return
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Fri Jun  2 15:14:30 EDT 2006
c version: 1
c Name: J. Finke
c Replaced ran1 with fibran and ran1. Added reflection 
c at lower boundary.      
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Fri Jun  2 15:22:22 EDT 2006
c version: 2
c Name: J. Finke
c Added fibran as well and ran1 for random 
c number generation. Added total reflection at lower boundary. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Sun Jun  4 19:37:16 EDT 2006
c version: 3
c Name: J. Finke
c added zseeds and rseeds     
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Thu Jul 27 13:07:43 EDT 2006
c version: 4
c Name: J. Finke
c Added spec_switch, so quick looks spectra can be 
c spectra incident on top and bottom boundaries.  
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Aug  8 12:12:06 EDT 2006
c version: 5
c Name: J. Finke
c Fixed bug at end, where photons at upper 
c z surface are read to event and spectra 
c files twice. Also changed leak to lower z 
c surface, so that photons there are read to 
c event and spectra files.     
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jan  2 20:15:39 EST 2007
c version: 5
c Name: J. Finke
c Created subroutine get_bin for determining bins of a 
c photon.        
