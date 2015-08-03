      subroutine planck(tpl)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
c
c     called by genimc, imcsurf, imcvol
c
c    AUTHOR...Gene Canfield
c
c     samples normalized planckian (or wien) spectrum
c
c    GLOSSARY...
c
c     input: tpl(j) = temperature array
c             iwien = 0  select from planck spectrum
c                   = 1  select from   wien spectrum
c             hu(l) = photon group boundary array [keV]
c     output:
c             xnu = photon energy [keV]
c           jgpsp = photon group (spectrum binning)
c           jgplc = light curve photon group
c           jgpmu = angular bin no.
c
c
      double precision tpl
c
      integer iwien
      integer jbot, jtop, jmid, m, n
c
      double precision hubot, hutop
      double precision u4, ap0, ap1, ap2, ap3, rn1
      double precision fibran
c
c
c
c
      iwien = 0
c
c   99 u4 = drand(0)
c      u4 = u4*drand(0)
c      u4 = u4*drand(0)
c      u4 = u4*drand(0)
c
   99 u4 = fibran()
      u4 = u4*fibran()
      u4 = u4*fibran()
      u4 = u4*fibran()
      if (u4.le.1.d-200) goto 99
      ap0 = -log(u4)
      ap1 = 1.d0
      ap2 = 1.d0
      ap3 = 1.d0
      if (iwien .eq. 1) go to 120
c      rn1 = 1.08232d0*drand(0)
      rn1 = 1.08232d0*fibran()
c
  100 continue
      if( rn1 .le. ap1 ) go to 120
      ap2 = ap2 + 1.d0
      ap3 = 1.d0/ap2
      ap1 = ap1 + ap3**4
      go to 100
c
  120 continue
      xnu = ap0*ap3*tpl
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
c      write(10, 800) rpre, zpre
c      write(10, 801) xnu
c      write(10, 802) phi, wmu
c      write(10, 803) jgpsp, jgplc
c      write(10, 804) jgpmu
c
c
      return
c
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:32:12 EDT 2006
c version: 2
c Name: J. Finke
c Changed common block 'random'.     
