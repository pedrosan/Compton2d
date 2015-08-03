        subroutine compb2d(i_gam)
        implicit none
        include 'general.pa'
        include 'commonblock.f'
c
        
        integer i_gam, m, n, icoms
        integer isurfmu
        integer jbot, jtop, bdotr, ncoll, jmid
       
        double precision wa, wb, swa, cazes, omege,  fuzz, cazs, caz    
        double precision phis, phie, rhos, alphas, theta, omeg
        double precision znue, znue3, betz, gamz, xxx, sz, 
     1                   games, gams, phat, omegs, omeges
        double precision costheta, cosphis, cos_alphas, tl,
     1                   tr, znues, znus, xknot, sintheta, cosphi
        double precision thetas, cosrhos, cosalphas, ucaz, ucazs,
     1                   uca, coscaz, znu
        double precision xnus, wmus, hubot, hutop
        double precision tt
        double precision emasskev, gamm, betb
        double precision cosdphi, dphi
        double precision fibran
c
c
c
c
        data fuzz / 1.d-10 /
c
        icoms = 1
        emasskev = 5.11d2      
        znu = xnu / emasskev
        tt = tea(jph,kph) / emasskev
c
c        
 100  continue        
c
c
c     Draw new particle energy and velocity
c      from electron population according
c       to FP solution (MB, 13/May/1999)
c
c      write(*,*) 'xnu = ',xnu
c      write(*,*) 'wmu = ',wmu
c      write(*,*) 'phi = ',phi
c
      call nth2d(jph, kph, gamm, betb, i_gam)
c
c      write(*,*) 'gamma = ',gamm
c      write(*,*) 'beta = ', betb
c
c
c     betb = sample of speed from isotropic relativistic
c     Maxwellian + non-thermal power-law distribution.
c     This prospective target may end up rejected.  next select
c     omeg = lab cosine of angle between electron and incident photon
c
c         omeg = 2.d0*drand(0) - 1.d0
         omeg = 2.d0*fibran() - 1.d0
         if (omeg.gt.9.9999999d-1) omeg = 9.9999999d-1
         if (omeg.lt.-9.9999999d-1) omeg = -9.9999999d-1
c         tl = drand(0)
         tl = fibran()
         tr = 0.5d0*(1.0d0 - betb*omeg)
         if(tl .gt. tr) omeg = -omeg

         if (omeg.gt.9.9999999d-1) omeg = 9.9999999d-1
         if (omeg.lt.-9.9999999d-1) omeg = -9.9999999d-1
c
c     znue = doppler shifted photon frequency.
c     accept target electron with prob. = xknot = total K-N/Thomson.
c     for small znue, use series to prevent bad roundoff.
c     (for znue = .01, error .lt. (544/7)*znue**5 = e-08)
c
      znue = (1.-betb*omeg)*znu*gamm
      if  (znue.lt.1.d-10) goto 100
      if(znue .le. 1.d-2) then
        xknot = 1.0d0-znue
     1         *(2.0d0-znue*(5.2d0-znue*(13.3d0-1.144d3*znue/3.5d1)))
      else
        znue3 = znue*znue*znue
        betz   = 1.0d0 + 2.0d0*znue
        gamz   = znue*(znue-2.0d0) - 2.0d0
        xxx = 4.0d0*znue + 2.0d0*znue3*(1.0d0+znue)/betz**2 
     1      + gamz*dlog(betz)
        xknot = 3.75d-1*xxx/znue3
      endif
c
c      write(*,*) 'xknot = ',xknot
c      write(*,*) 'znue = ',znue
c        
c      if(drand(0) .gt. xknot) go to 100
       if(fibran() .gt. xknot) go to 100
c     ncoll = ncoll + 1
c       
      betz = 1.0d0 + 2.0d0*znue
c
 200  continue
 202  sz = (1.0d0 + 2.0d0*znue*fibran())/betz
c     sz = (1.0 + 2.0*znue*drand(0))/betz
      games = 1.0d0 + (1.0d0-1.0d0/sz)/znue
      if ((1.d0 - games**2).lt.0.) goto 202
      tr = games**2 - 1.0d0 + sz + 1.0d0/sz
      phat = betz + 1.d0/betz
c     if(drand(0)*phat .gt. tr) go to 200
      if(fibran()*phat .gt. tr) go to 200
      znues = znue*sz
      
c      write(*,*) 'znues = ',znues
c
  210 continue
c      wa = drand(0)
c      wb = 2.*drand(0) - 1.
      wa = fibran()
      wb = 2.d0*fibran() - 1.d0
      swa = wa*wa+wb*wb
      if((swa.ge.1.0d0) .or. (swa.le.1.d-20)) go to 210
c
      cazes = (wa*wa-wb*wb)/swa
      
c      write(*,*) 'cazes = ',cazes
c        
c     omege = cosine of angle between photon and electron in rest frame.
c     omeges = cosine between elec. direction and scattered photon.
c
      omege = (omeg-betb)/(1.-betb*omeg)
      if (omege.gt.9.9999999d-1) omege = 9.9999999d-1
      if (omege.lt.-9.9999999d-1) omege = -9.9999999d-1
      omeges = games*omege + 
     1         cazes*dsqrt((1.-omege**2+fuzz)*(1.-games**2))
      if (omeges.gt.9.9999999d-1) omeges = 9.9999999d-1
      if (omeges.lt.-9.9999999d-1) omeges = -9.9999999d-1
c
c      write(*,*) 'omeges = ',omeges
c      write(*,*) 'omege = ',omege
c
c     xform back to lab frame.
c
c     omegs = cosine between electron and scattered photon in lab frame.
c     znus = scattered photon frequency in lab.
c     gams = cosine of scattering angle in lab.
c
      omegs = (omeges+betb)/(1.0d0+omeges*betb)
      if (omegs.gt.9.9999999d-1) omegs = 9.9999999d-1
      if (omegs.lt.-9.9999999d-1) omegs = -9.9999999d-1
      znus = (1.d0 + betb*omeges)*gamm*znues
      gams = 1.d0 - (znue - znues)/(znu*znus)
      if (gams.gt.9.9999999d-1) gams = 9.9999999d-1
      if (gams.lt.-9.9999999d-1) gams = -9.9999999d-1
c
c      write(*,*) 'omegs = ',omegs
c      write(*,*) 'znus = ',znus
c      write(*,*) 'gams = ',gams
c
c     calculate lab scattered direction vector relative to radius.
c     cazs = cos. of azimuthal scattering angle (uniform from symmetry.)
c     wmus = lab cosine between scattered photon and radius.
c
 220  continue
c      wa = drand(0)
c      wb = 2.*drand(0) - 1.d0
      wa = fibran()
      wb = 2.d0*fibran() - 1.d0
      swa = wa*wa + wb*wb
      if((swa.ge.1.0d0) .or. (swa.le.1.d-20)) go to 220
      cazs = (wa*wa-wb*wb)/swa
      if (cazs.gt.9.9999999d-1) cazs = 9.9999999d-1
      if (cazs.lt.-9.9999999d-1) cazs = -9.9999999d-1
      wmus = wmu*gams + cazs*dsqrt((1.0d0-gams**2)
     1      *(1.0d0-wmu**2 + fuzz))
      if (wmus.gt.9.9999999d-1) wmus = 9.9999999d-1
      if (wmus.lt.-9.9999999d-1) wmus = -9.9999999d-1
      xnus = znus*emasskev
c
c      write(*,*) 'cazs = ',cazs
c      write(*,*) 'wmus = ',wmus
c      write(*,*) 'xnus = ',xnus
c
c  Introducing the two dimensional scenario in compton scattering
c  Now using two sets of axes 'z' and 'r'.
c  phie =  Azimuthal scattering angle of the electron - randomly distributed.
c  phis =  Azimuthal scattering angle of the scattered photon with the 'r' axis.
c  phi =  
c  theta =  Angle between the 'z' axis and the electron.
c  rhos = Angle between the scattered photon and 'r' axis.
c
c
c      phie = 4.4d1*drand(0)/7.d0
c
c      caz = (omegs- gams *omeg)/ dsqrt((1.-omegs**2)*(1.-gams**2 ))
c      if (caz.gt.0.99999999d0) caz = 0.99999999d0
c      if (caz.lt.-0.99999999d0) caz = -0.99999999d0
c      ucaz  = acos(caz)
c      ucazs = acos(cazs) 
c      uca   = ucaz + ucazs
c      coscaz = cos(uca)
c      
c      write(*,*) 'caz = ',caz
c      write(*,*) 'uca = ',uca
c      write(*,*) 'coscaz = ',coscaz
c      
c      costheta = omeg*wmu + (dsqrt((1.-omeg**2)*(1.-wmu**2))*coscaz)
c      if (costheta.gt.0.99999999d0) costheta = 0.99999999d0
c      if (costheta.lt.-0.99999999d0) costheta = -0.99999999d0
c      theta = acos(costheta)
c      sintheta = sin(theta)
c      if (dabs(sintheta).lt.1.d-10) sintheta = 1.d-10
c      cosphi = (omegs - wmus*costheta)/(sintheta*dsqrt(1.-wmus**2))
c      if (cosphi.gt.0.99999999d0) cosphi = 0.99999999d0
c      if (cosphi.lt.-0.99999999d0) cosphi = -0.99999999d0
c      
c      write(*,*) 'costheta = ',costheta
c      write(*,*) 'sintheta = ',sintheta
c      write(*,*) 'cosphi = ',cosphi
c
c
c      phi = acos(cosphi)
c      phis = phie-phi
c      thetas = acos(wmus)
c
c
c      alphas = ( 2.2d1/1.4d1 - thetas)
c      cosphis= cos(phis)
c      cos_alphas= cos(alphas)
c      cosrhos = cosphis*cos_alphas
c      if (cosrhos.gt.0.99999999d0) cosrhos = 0.99999999d0
c      if (cosrhos.lt.-0.99999999d0) cosrhos = -0.99999999d0
c      rhos = acos(cosrhos)
c
       cosdphi = (gams - wmu*wmus)
     1          /dsqrt((1.0d0 - wmu**2)*(1.0d0 - wmus**2))
       if (cosdphi.gt.9.9999999d-1) cosdphi = 9.9999999d-1
       if (cosdphi.lt.-9.9999999d-1) cosdphi = -9.9999999d-1
       dphi = acos(cosdphi)
       phis = phi + dphi
c
c       write(*,*) 'cosdphi = ',cosdphi
c       write(*,*) 'dphi = ',dphi
c       write(*,*) 'phis = ',phis
c
c
  250 continue
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
      if( xnus .ge. hutop ) then
         jgpsp = 0
         go to 300
      endif
      if( xnus .le. hubot ) then
         jgpsp = 0
         go to 300
      endif
c
  140 continue
      jmid = (jbot+jtop)/2
      if( jmid .eq. jbot ) go to 160
      if(xnus .eq. hu(jmid)) go to 160
      if(xnus .lt. hu(jmid)) then
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
         if ((xnus.gt.Elcmin(m)).and.
     1       (xnus.le.Elcmax(m))) then
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
         if (wmus.le.mu(n)) then
            jgpmu = n
            goto 900
         endif
  700 continue
      jgpmu = nmu
c
c
c
 900  continue
      ew = ew*xnus/xnu
      xnu = xnus
      wmu = wmus
      phi = phis

c  
c      write(*,*) 'xnus = ',xnu
c      write(*,*) 'wmus = ',wmu
c      write(*,*) 'phis = ',phis
c
      return
      end
     
      
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:24:02 EDT 2006
c version: 2
c Name: J. Finke
c Fixed bug in common block 'random'.   
