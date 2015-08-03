c     This subroutine is a copy of imctrk2d, and it is called and calls
c     imctrk2d. However, it differs imctrk2d in several ways:
c     1.It does not split the photon before it is scattered.
c     2.It has an option of whether to scatter the photon. When the photon
c       was decided not to scatter in imctrk2d, imctrk will call it with 
c       this option 0.
c     3.When the scattering is turned off, edep will not be recorded.
c     4.Writting of photons to census can only happen in this subroutine.

      subroutine imctrk2d_sub(scat_flag)
      implicit none
      include 'ucens.pa'
      include 'general.pa'
      include 'commonblock.f'
c
c, ucens, iucens, ducens
      integer scat_flag, split2,split3
      parameter(split2=1000)
      parameter(split3=1000)
c      parameter(ucens = 1000)
c      parameter(iucens = 6000)
c      parameter(ducens = 6000)
c      parameter(ucens = 50000000)
c      parameter(iucens = 300000000)
c      parameter(ducens = 300000000) ! Xuhui ucens 
c
      logical gone
c
      integer i_gam
      integer i
      integer iii,pfrl(30),ii,ii2
      integer jcsv,kcsv,jgpspcsv,
     1        jgplccsv,jgpmucsv
      double precision ewcsv, xnucsv, wmucsv,rprecsv,
     1                 zprecsv, phicsv, dcencsv
      
      integer idead, nptrks, knew, jnew, kbnd, jbnd, npkill
      integer iscat, lwad, lwai
      integer i_gg, inout, npbnd
      integer isurfmu, nphtotal, nph_lc
      integer ikind
      integer Eta_switch
c
      integer incounter, outcounter, cmcount1,cmcount2, lkcount
      double precision ewold
      double precision Zr, Rr, eoutr(jmax, kmax)
      double precision disp, psq, disbr, dpbsq, Zbnd, rbnd,
     1                 xqsqleft, trld, trldb, dcol       
      double precision wtmin, wkth, wmue, ewsv
      double precision Egg_min
      double precision kgg, kgg0, kgg1
      double precision pair_enhance
      double precision theta, Eta, comac, enexc
      double precision sigabs, velfact, facomp
      double precision colmfp, mb_ran, denom, f
      double precision sigcom, sigcomex, sigsc, deleabs, delecomp
      double precision ewnew, xabs, ekill, ekillt, ewpl, sstar
      double precision delpr, exan, wmusv, xnusv, wmustar
      double precision fibran, rnew, znew
c
      common /frl/ pfrl !Xuhui
      common / count / incounter, outcounter, cmcount1,cmcount2,lkcount ! Xuhui cm
      common / itrk / idead, iscat 
c      common / sigmas / comac, enexc
c
c        common / debug / nlost,gone
c        
       idead = 0
       sigabs = 1.d-40
       velfact = 1.d0
       facomp = 1.d0
c       wkth = 1.d-2
       wkth = 1.d-3
c
c      idead = -1 from imcvol if rwlk turned on.
c              =  0 particle returned from rwlk alive.  track it.
c              =  1 killed in rwlk or imcleak.
c              =  2 retire rwlk particle to census.
c              =  3 particle was reflected from surface in imcleak.
c              =  4 particle leaked into void for tracking.
c
c     
         wtmin = wkth*ew
c
c         if ( idead. eq. 1) go to 230
c         
c        Sum up squared photon energies
c        and squares for recoil  
c         
c          sum_e = sum_e + 1.957d-3*ew*xnu
c          sum_e2 = sum_e2+ew*((1.957d-3*xnu)**2)
c
  100  continue
       sigabs = 1.d-40 ! Xuhui 5/24/09
c         mb_ran = drand(0)
         if(scat_flag .eq. 0)then
           mb_ran = 1.d-10
         else
           mb_ran = fibran()
c Xuhui 5/1/09
          if(rand_switch.eq.2)
     1       mb_ran=int(mb_ran*1.d6)/1.d6 + 1.d-6*fibran()
c added to combine two 32 bits random number into one 64 bits random number.

         endif
         if (mb_ran.gt.0.) then
             colmfp = -dlog(mb_ran)
         else 
             goto 100
         endif
c
  110  continue
c
       if (ew.lt.1.d-40) goto 900
       if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
       if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
c          
         if (pair_switch.eq.1) then
             pair_enhance = 1.d0 + 2.d0*f_pair(jph, kph)
         else
             pair_enhance = 1.d0
         endif
c  
         if (scat_flag .ne. 0)then
             call comtot(jph, kph, xnu, 1,comac,enexc)
         else
             comac = 0.
             enexc = 0.
         endif
c
         sigcom = velfact*comac*pair_enhance
         sigcomex = velfact*enexc*pair_enhance
         sigsc = sigcom
c
         if (kph.eq.1) then
            xqsqleft = rmin**2
         else
            xqsqleft = r(kph-1)**2
         endif
c         
c         
        if(scat_flag .ne. 0)then
           dcol = colmfp/sigsc
        else
           dcol =100*dmax1(r(nr),z(nz)) ! Xuhui 3/4/09
        endif

c        iii = int(dcol/1e16)
c        if(iii.le.30 .and. iii.gt.0)pfrl(iii) = pfrl(iii) +1 !Xuhui

        if ( dcen .le. dcol ) then
           trld = dcen
           ikind = 2
        else
          trld = dcol
          ikind = 3
        endif
c
c        if (xnu.gt.4.3d2) write(*,*) 'trld(col/cen) = ',trld
c
       Eta = cos(phi)
c     Eta_switch is +1 if phi is in quadrants I or II, and -1
c     if it is in quadrants III or IV.  This way
c     the quadrant of phi is not lost when taking its cosine.  
c     J. Finke, 6 Sept. 2005
       if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
          Eta_switch = 1
       else
          Eta_switch = -1
       endif
c       
       if (Eta.gt.9.9999999d-1) Eta = 9.9999999d-1
       if (Eta.lt.-9.9999999d-1) Eta = -9.9999999d-1
c
       disp = Eta*rpre
       psq = rpre*rpre*(1. - Eta**2)
c
       if ((Eta.lt.0.d0) .and. (psq.lt.xqsqleft))  then
c
          incounter = incounter + 1
          kbnd = kph-1
          inout = -1
          if (kph.gt.1) then
             rbnd = r(kph-1)
          else
             rbnd = rmin
          endif
c
       else
          outcounter = outcounter + 1
          kbnd = kph
          inout = 1
          rbnd = r(kph)

       endif     
c
       dpbsq = rbnd**2 - psq
       if (dpbsq.lt.1.d-6) dpbsq = 1.d-6
c
       disbr = dble(inout)*sqrt(dpbsq) - disp
       trldb = disbr/dsqrt(1.d0 - (wmu**2))
       f     = disbr
       Zr = Zpre + wmu*trldb
c
        if (jph.eq.1) then
c        
           if  ( Zr.gt.Z(jph)) then
c
c         nearest boundary = upward z boundary
c          
              Zbnd = Z(jph)
              knew = kph
              jnew = jph+1
              f = (zbnd - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt (rpre**2 + f**2 + 2.d0*rpre*f*Eta)
              rbnd =  Rr     
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
                               
           else if (Zr.lt.zmin) then
c
c        nearest boundary = downward z boundary
c       
              Zbnd = zmin
              knew = kph
              jnew = jph - 1
              f = (zmin - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt (rpre**2 + f**2 + 2.d0* rpre*f*Eta)           
              rbnd = Rr
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
     
           else 
c
c        nearest boundary = r boundary
c
               knew = kph + inout
               jnew = jph
               if (kbnd.gt.0) then
                  Rr = r(kbnd)
               else
                  Rr = rmin
               endif
               rbnd = Rr
               Zbnd = Zr
              
           endif
c
        else
c       j is greater than 1   cccccccccccccccccccccc
c
           if  ( Zr.gt.Z(jph)) then
c
c         nearest boundary = upward z boundary
c          
              Zbnd = Z(jph)
              knew = kph
              jnew = jph+1
              f = (zbnd - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt(rpre**2 + f**2 + 2.d0*rpre*f*Eta)
              rbnd =  Rr     
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
                           
           else if (Zr.lt.Z(jph-1)) then
c
c        nearest boundary = downward z boundary
c       
              Zbnd = Z(jph-1)
              knew = kph
              jnew = jph - 1
              f = (zbnd - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt(rpre**2 + f**2 + 2.d0*rpre*f*Eta)            
              rbnd = Rr
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
c
           else
c
c        nearest boundary = r boundary
c
               knew = kph + inout
               jnew = jph
               if (kbnd.gt.0) then
                  Rr = r(kbnd)
               else
                  Rr = rmin
               endif
               rbnd = Rr
               Zbnd = Zr

           endif
         endif
c
       if (trldb.lt.trld) then
c  
           ikind = 1
           trld = trldb
           rnew = rbnd
           znew = Zbnd
c
       else 
           jnew = jph
           knew = kph
           f = trld*dsqrt(1.d0 - wmu**2)
           rnew = dsqrt(f**2 + rpre**2 + 2.d0*f*rpre*Eta)
           znew = zpre + trld*wmu
c
       endif
c
       do 122 i = 1, n_vol-1
 122      if (xnu.lt.E_ph(i+1)) goto 123
 123   continue
c
       if (pair_switch.eq.0) goto 128
       do 127 i_gg = 1, n_gg
 127   if (xnu.lt.E_gg(i_gg)) goto 128
 128   continue
c
         sigabs = sigabs + pair_enhance*kappa_tot(i, jph, kph)
c
         if (pair_switch.eq.1) then
            if (xnu.lt.E_gg(1)) then
               kgg = (xnu/E_gg(1))*k_gg(1, jph, kph)
            else if (i_gg.eq.n_gg) then
               kgg = k_gg(n_gg, jph, kph)
            else
               kgg0 = k_gg(i_gg, jph, kph)
               kgg1 = k_gg(i_gg+1, jph, kph)
               kgg = kgg0 + (xnu - E_gg(i_gg))*(kgg1 - kgg0)
     1                     /(E_gg(i_gg+1) - E_gg(i_gg))
            endif
            sigabs = sigabs + kgg
         endif
         if (sigabs.lt.1.d-40) sigabs = 1.d-40
c
         xabs = sigabs*trld
c
         if (xabs.lt.100.) then
            ewnew = ew*dexp(-xabs)
        else
            ewnew = 0.
        endif
c        if ( ewnew .le. wtmin ) then
c           npkill = npkill + 1
c           ekill  = ekill + ewnew
c           ekillt = ekillt + ewnew
c           ewnew  = 0.
c       endif
c
c
c       if ((xnu.gt.47.d0).and.(sigabs.gt.1.d-40)) then
       if((xnu.gt.47.d0).and.(sigabs.gt.1.d-40).and.
     1       (pair_switch.eq.1))then ! Xuhui 5/19/09
          deleabs = (ew - ewnew)*(sigabs - k_gg(i_gg, jph, kph))/sigabs
       else
          deleabs = ew - ewnew
       endif
       if (deleabs.lt.1.d-50) deleabs = 1.d-50
c
c
c  ewpl = energy weight path length
c
       if ( xabs .le. .00001d0 ) then
         ewpl = ew*trld*(1.d0-.5d0*xabs)
         wmustar = wmu
      else
         ewpl = deleabs/sigabs
c  119    mb_ran = drand(0)
  119    mb_ran = fibran()
         if (mb_ran.lt.(ew/deleabs)) then
            sstar = -dlog(1. - mb_ran*deleabs/ew)/sigabs
         else
            goto 119
         endif
         
         denom = dsqrt (rpre**2 + 2.*wmu*rpre*sstar + sstar**2)
         wmustar =  (wmu*rpre + sstar)/denom
      endif
c
      delpr = deleabs*wmustar*c_light
      delecomp = facomp*sigcomex*ewpl
      if(scat_flag .ne. 0)then
        edep(jph, kph) = edep(jph, kph) + deleabs ! + delecomp Xuhui 2/15/11
        prdep(jph, kph) = prdep(jph, kph) + delpr
      endif ! Xuhui 3/4/09
c      
c      if (ewnew .le. 1.d-40) go to 900
      if(ewnew.le.wtmin) goto 900 !Xuhui 5/26/09
c
      ew = ewnew
      dcen = dcen - trld
c      if (rnew.gt.1.d-10) then
c         Eta = (trld - Eta*rpre)/rnew
c      else      
         Eta = (trld + Eta*rpre)/rnew
c      endif
      
      if (Eta.gt.0.999999999d0) Eta = 0.999999999d0
      if (Eta.lt.-0.999999999d0) Eta = -0.999999999d0
      phi = acos(Eta)
      if (eta_switch.eq.-1) phi = 2.d0*pi - phi
       if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
          Eta_switch = 1
       else
          Eta_switch = -1
       endif 
      rpre = rnew
      Zpre = znew
c      call bin_add(phi)  
c       
       if ( ikind .eq. 3 ) go to 210
       if ( ikind .eq. 2 ) go to 150
c
c     particle has reached boundary, r(kbnd) (ikind = 1)
c
c      npbnd = npbnd + 1
c      eoutr(jbnd, kbnd) = eoutr(jbnd, kbnd) + inout*ew/dt(1)
c
      colmfp = colmfp - sigsc*trld
c
      
 140  continue
      if((jnew .eq. nz+1) .or. (jnew.eq.0) .or.
     1   (knew .eq. nr+1) .or. (knew.eq.0)) then
c
         jph = jnew
         kph = knew
c
c         if(knew.eq.0) incounter = incounter + 1
c         if(knew.eq.nr+1) outcounter = outcounter + 1
c
         call imcleak
         lkcount = lkcount + 1 ! Xuhui
c
         if (idead .eq. 1) go to 900
c
         if (idead .eq. 4) then
            kph = 0
            go to 100
         endif
c
c
c         go to 110
         go to 100 ! Xuhui 5/24/09
      endif
c     
      kph = knew
      jph = jnew
c      go to 110 
      go to 100 ! Xuhui 5/24/09
c
c       case ikind = 2: write particle to census
c
c
 150  continue
c      if(ew.lt.wtmin)goto 900 ! Xuhui 5/18/09
c      if the energy weight is below wtmin, kill the photon.
      npcen(jph,kph) = npcen(jph,kph) + 1
 160  continue
      ecens(jph, kph) = ecens(jph, kph) + ew
c
      if (pair_switch.eq.1) then
         do 162 i_gg = 1, n_gg
 162        if (xnu.lt.E_gg(i_gg)) goto 163
 163     continue
         Egg_min = (E_gg(1)**2)/E_gg(2)
         if (xnu.gt.Egg_min) then
             n_ph(i_gg, jph, kph) = n_ph(i_gg, jph, kph)+(ew*6.25d8)/xnu
         endif
      endif
c
      do 165 i_gg = 1, nphfield-1
  165 if (xnu.lt.E_field(i_gg+1)) goto 166
  166 continue
      Egg_min = (E_field(1)**2)/E_field(2)
c
      if (xnu.gt.Egg_min) then
         n_field(i_gg, jph, kph) = n_field(i_gg, jph, kph) 
     1                       + 6.25d8*ew/xnu
      endif
c
       lwad = 6*ndxout
       lwai = 6*ndxout
       dbufout(lwad+1) = rpre
       dbufout(lwad+2) = Zpre
       dbufout(lwad+3) = wmu
       dbufout(lwad+4) = phi
       dbufout(lwad+5) = ew
       dbufout(lwad+6) = xnu
       ibufout(lwai+1) = jgpsp
       ibufout(lwai+2) = jgplc
       ibufout(lwai+3) = jgpmu
       ibufout(lwai+4) = jph
       ibufout(lwai+5) = kph
       ibufout(lwai+6) = int( fibran()*1.d5 )
       ndxout = ndxout + 1
       if(ndxout .ge. ucens) then
c           call write_cens(ndxout) ! Xuhui cens
          write(*,*)'too many photons'
          stop ! Xuhui cens
          ndxout = 0
       endif
       go to 900
c
c     case ikind = 3: Compton scattering
c
 210   continue
       wmusv = wmu
       xnusv = xnu
       ewsv = ew
       iscat = 0
c       
c       if ( drand(0)*sigsc .lt. sigcom ) iscat = 1
c

*****************XUHUI*******************************************************
c       write(*,*)'second scattering'
       ewcsv = ew
       dcencsv = dcen
       xnucsv = xnu
       zprecsv = zpre
       rprecsv = rpre
       wmucsv = wmu
       phicsv = phi
       jcsv = jph
       kcsv = kph
       jgpspcsv = jgpsp
       jgplccsv = jgplc
       jgpmucsv = jgpmu
ccccccccccccccccccccccc this subsection is to decide the split2
c       call compb2d
       write(*,*)'2nd scattering'
       cmcount1=cmcount1+1
c       split2 = max0(min0(int(ew/ewcsv/1000),100000000),1)
c       split2 = max0(min0(int((511.d0/xnu)**2/1000000),1000000),1)
       ewcsv = ewcsv/split2
ccccccccccccccccccccccc
       
       do 220 ii=1,split2
       ew = ewcsv
       dcen = dcencsv
       xnu = xnucsv
       zpre = zprecsv
       rpre = rprecsv
       wmu = wmucsv
       phi = phicsv
       jph = jcsv
       kph = kcsv
       jgpsp = jgpspcsv
       jgplc = jgplccsv
       jgpmu = jgpmucsv
       
       ewold = ew
       call compb2d(i_gam)
c_______________________________________
       if(ew.gt.ewold*split2*60)then
          ewold = ewcsv/split3
          do 217, ii2=1,split3
215         continue
             ew = ewcsv/split3
             dcen = dcencsv
             xnu = xnucsv
             zpre = zprecsv
             rpre = rprecsv
             wmu = wmucsv
             phi = phicsv
             jph = jcsv
             kph = kcsv
             jgpsp = jgpspcsv
             jgplc = jgplccsv
             jgpmu = jgpmucsv
             call compb2d(i_gam)
            if(ew.le.ewold*split2*60)goto 215
            edep(jph,kph) = edep(jph,kph) + ew - ewold ! Xuhui 2/15/11

            E_IC(i_gam) = E_IC(i_gam) + ew - ewold
c            cmener = cmener + ew - ewold
            cmcount2 = cmcount2 +1
            if(phi.gt.2.d0*pi) phi = phi - 2.d0*pi
            if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
               Eta_switch = 1
            else
               Eta_switch = -1
            endif
            call imctrk2d
217       continue
c___________________________________________________________________
       else
         edep(jph,kph) = edep(jph,kph) + ew - ewold
         E_IC(i_gam) = E_IC(i_gam) + ew - ewold
c         cmener = cmener + ew - ewold
c       if (ew .lt. wtmin)goto 220 ! Xuhui 5/18/09
c      if the energy weight is below wtmin, kill the photon.
         cmcount2 = cmcount2 +1 ! Xuhui

         if(phi.gt.2.d0*pi) phi = phi - 2.d0*pi
         if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
            Eta_switch = 1
         else
            Eta_switch = -1
         endif

         call imctrk2d
       endif
220    continue
c
c       if (gone) then
c       gone=.false.
c       goto 900
c       endif
c
c      goto 100
      goto 900 ! Xuhui

 230  continue         
 900  continue
c
      return
      end    
c

       
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:26:40 EDT 2006
c version: 2
c Name: J. Finke
c Changed common block 'random'.     
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Fri Jun 16 12:31:21 EDT 2006
c version: 3
c Name: J. Finke
c seed is now stored in census files. imcfield 
c is now determinable.      
