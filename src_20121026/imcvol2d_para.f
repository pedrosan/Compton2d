       subroutine imcvol2d
c      This subroutine convert the radiative energy lose of the zone in this 
c      time step into photons, and write them into census file. 
c      ! comment by Xuhui Chen 7/8/08

       implicit none
       include 'mpif.h'
       include 'general.pa'
       include 'commonblock.f'
c
c     MPI variables
      integer status(MPI_STATUS_SIZE), 
     1        num_sent, end_signal, sender,
     1        num_zones       
      integer l, zone
c
c
c     MPI Initialization 
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      num_sent = 0
      end_signal = jmax*kmax+1
c
c
c     master part
      if(myid.eq.master) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      call vol_bcast
      num_zones = nr*nz
      call vol_create_job_type
c
c
c     This sends the first round of zones to the slaves for processing.
      do 500 l= 1, min(num_zones, (numprocs-1))
         zone = l
         call vol_send_job(l, zone)
         num_sent = num_sent + 1
 500  continue  
c
c     As slaves complete processing a zone, this recieves the results
c     and sends the next zone to the slaves.
      do 501 l = 1, num_zones
         call vol_recv_result(sender, zone)
         if(num_sent.lt.num_zones) then
            zone = num_sent+1
            call vol_send_job(sender, zone)
            num_sent = num_sent + 1
         else
            zone = num_sent
            call vol_send_end_signal(sender,zone)
         endif
 501  continue
c
c
c
      else if(myid.ne.master) then
c     beginning of slave part
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      call vol_bcast
      num_zones = nr*nz
      call vol_create_job_type
c
c     if there are more nodes than work skip this node.
      if(myid.gt.num_zones) goto 990
c
 991  continue
      call vol_recv_job(zone)
      if(zone.eq.end_signal) goto 990
      call vol_calc(zone)
      call vol_send_result(zone)
      goto 991
 990  continue
c
      endif
c     end of slave part
c
c
c
c     end of imcvol2d
      return
      end           
c
c
c
c
      subroutine vol_calc(zone)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer zone
c    
      integer i, jv, kv
      integer npsurf, isurfmu
      integer npvol, npick, ivol, more
      integer jmid, jbot, jtop
      integer m, n
c
      double precision fas(jmax, kmax)
      double precision esurf, ew_save
      double precision ewcd, f_thermal
      double precision rnum0, rnum
      double precision psi, x1, x2
      double precision hubot, hutop, evol, norwlk, xxmore, rnpick
      double precision f_outer, f_inn, f_lower, f_upper, delz
      double precision fibran

c
c
c
      call get_j_k(zone, jv, kv, nr)
      call initialize_rand(jv,kv)
            if(kv.eq.1) then
             endif
c
            ewcd = ewsv(jv,kv)
            ew_save = ewcd
            more = nsv(jv,kv)
                        
            f_thermal = Eloss_th(jv,kv)/Eloss_tot(jv,kv)
         
            if (jv.eq.1) then
               delz = z(1)
            else
               delz = z(jv) - z(jv-1)
            endif

            if (kv.eq.1) then
               f_inn  = (4.4d1/7.d0*rmin*delz)/zsurf(jv, kv)
               f_outer  = (4.4d1/7.d0*r(kv)*delz)/zsurf(jv, kv)
               f_upper = (2.2d1/7.d0 *(r(kv)**2
     1                        - rmin**2))/zsurf(jv, kv)
               f_lower = (2.2d1/7.d0 *(r(kv)**2
     1                        - rmin**2))/zsurf(jv, kv)
            else
               f_inn  = (4.4d1/7.d0*r(kv-1)*delz)/zsurf(jv, kv)
               f_outer  =  (4.4d1/7.d0*r(kv)*delz)/zsurf(jv, kv)
               f_upper = (2.2d1/7.d0 *(r(kv)**2
     1                        - r(kv-1)**2))/zsurf(jv, kv)
               f_lower = (2.2d1/7.d0 *(r(kv)**2
     1                        - r(kv-1)**2))/zsurf(jv, kv)
            endif
            
            f_outer = f_inn + f_outer
            f_upper = f_outer + f_upper
            f_lower = f_upper + f_lower
c
  50        format (' f_lower = ',e14.7)
            if (dabs(f_lower - 1.d0).gt.1.e-2) then
               write(4,50) f_lower
            endif
c                  
c      rseed = 8535           
            do 450 ivol = 1, more
               
               jph = jv
               kph = kv
               ew = ew_save
c               dcen = c_light*dt(1)*drand(0)
               dcen = c_light*dt(1)*fibran()
c               rnum = drand(0)
               rnum = fibran()
               if (rnum.lt.f_thermal) then
                   i = 0
c                   rnum = drand(0)    
                   rnum = fibran()     
 120               i = i+1
                   if (eps_th(i, jph, kph).lt.rnum .and. i.lt.n_vol) 
     1                 go to 120
c
                   if (i.lt.n_vol) then
c                      xnu = E_ph(i)
c     1                    + drand(0)*(E_ph(i+1) - E_ph(i))
                      xnu = E_ph(i)
     1                    + fibran()*(E_ph(i+1) - E_ph(i))
                   else
                      xnu = E_ph(i) 
                   endif
c
c                   rnum0 = drand(0)
                   rnum0 = fibran()
                   if (rnum0.lt.f_inn) then
c                       wmu = 2.d0*drand(0) - 1.d0
                       wmu = 2.d0*fibran() - 1.d0
                       if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                       if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
c
c                       x1 = drand(0)
c                       x2 = drand(0)
                       x1 = fibran()
                       x2 = fibran()
                       if (x1.lt.0.5) then
                          phi = 1.1d1/7.d0 + (1.1d1/7.d0)*x2
                          if (phi.lt.1.57079638d0) phi = 1.57079638d0
                       else 
                          phi = -1.1d1/7.d0 - (1.1d1/7.d0)*x2
                          if (phi.gt.-1.57079638d0) phi = -1.57079638d0
                       endif 

                       if (kph.eq.1) then
                          rpre = 1.00001*rmin
                       else
                          rpre = 1.00001*r(kph-1)
                       endif
                       if  (jph.eq.1) then          
                          zpre = z(1)*fibran()
                       else
                          zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                       endif

                   else if (rnum0.lt.f_outer) then
c
c                        wmu = 2.d0*drand(0) - 1.d0
c
                        wmu = 2.d0*fibran() - 1.d0
                        if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                        if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                        rpre = 0.999999*r(kph)
c
                        if (jph.eq.1) then          
                            zpre = z(1)*fibran()
                        else
                            zpre = z(jph-1)+fibran()*(z(jph)-z(jph-1))
                        endif
c                        phi = -1.1d1/7.d0 + 2.2d1/7.d0*drand(0)
                        phi = -1.1d1/7.d0 + 2.2d1/7.d0*fibran()
                        if (phi.lt.-1.5707963d0) phi = -1.57079063d0
                        if (phi.gt.1.5707963d0) phi = 1.5707963d0
c
c                        
                    else if (rnum0.lt.f_upper) then
c
                        wmu = fibran()
                        if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                        if (wmu.lt.0.d0) wmu = 0.d0
                        phi = 4.4d1/7.d0*fibran()
                        if (phi.gt.2.d0*pi) phi = 2.d0*pi
c
c                        wmu = drand(0)
c                        phi = 4.4d1/7.d0*drand(0)
c
c                        psi = drand(0)
                        psi = fibran()
                        Zpre= 0.999999*z(jph)    
                        if (kph.eq.1) then
                           rpre = dsqrt((rmin**2) + psi*(r(kph)**2
     1                                            - rmin**2))
                        else
                           rpre = dsqrt((r(kph-1)**2) + psi*(r(kph)**2 
     1                                               - r(kph-1)**2))
                        endif
c
                    else

c                        wmu = -drand(0)
c                        phi = 4.4d1/7.d0*drand(0)
c
                        wmu = -fibran()
                        phi = 4.4d1/7.d0*fibran()
                        if (wmu.gt.0.d0) wmu = 0.d0
                        if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                        if (phi.gt.2.d0*pi) phi = 2.d0*pi
c
                        if (jph.eq.1) then
                           Zpre = 1.000001*zmin
                           if (zpre.le.zmin) zpre = zmin + 1.d-6
                        else
                           Zpre = 1.000001*z(jph-1)
                           if (zpre.le.z(jph-1)) zpre = z(jph-1) + 1.d-6
                        endif
                          
c                        psi = drand(0)
                        psi = fibran()
                        if (kph.eq.1) then
                           rpre = dsqrt((rmin**2)+ psi*(r(kph)**2
     1                                            - rmin**2))
                        else
                           rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2
     1                                            - r(kph-1)**2))  
                        endif
                    endif             
              
                else
                    i = 0
c                    rnum = drand(0)
                    rnum = fibran()
  125               i  =  i + 1
                    if (eps_tot(i, jph, kph).lt.rnum .and. i.lt.n_vol) 
     1                 go to 125
c
                   if (i.lt.n_vol) then
c                      xnu = E_ph(i)
c     1                    + drand(0)*(E_ph(i+1) - E_ph(i))
                      xnu = E_ph(i)
     1                    + fibran()*(E_ph(i+1) - E_ph(i))
                   else
                      xnu = E_ph(i) 
                   endif
c
c                    wmu = 2.d0*drand(0) - 1.d0
c                    phi = 4.4d1/7.d0*drand(0)
c
                    wmu = 2.d0*fibran() - 1.d0
                    phi = 4.4d1/7.d0*fibran()
                    if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                    if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                    if (phi.gt.2.d0*pi) phi = 2.d0*pi
c
                    if (jph.eq.1) then          
                        zpre = z(1)*fibran()
                    else
                        zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                    endif
                    
c                    psi = drand(0)
                    psi = fibran()
                    if (kph.eq.1) then
                       rpre = dsqrt((rmin**2)+ psi*(r(kph)**2
     1                                        - rmin**2))
                    else
                       rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2
     1                                        - r(kph-1)**2))  
                    endif

                endif          
      
c
c     Determine photon group of new photon energy
c     in photon energy bin structure for spectrum
c                   
                jbot = 1
                jtop = nphtotal + 1
                hubot = 1.000001*hu(1)
                hutop =.999999*hu(jtop)
                if( xnu .ge. hutop ) then
                   jgpsp = nphtotal
                   go to 300
                endif
                if( xnu .le. hubot ) then
                   jgpsp = 0
                   go to 300
                endif
c
  140           continue
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
  160           continue
                jgpsp = jmid
c
  300           continue
c
c
c     Determine photon group of new photon energy
c   in photon energy bin structure for light curves
c
                jgplc = 0
                do 400 m = 1, nph_lc
                   if ((xnu.gt.Elcmin(m)).and.
     1                 (xnu.le.Elcmax(m))) then
                      jgplc = m
                      goto 600
                   endif
  400           continue
c
  600           continue
c
c
c    Determine angular bin of photon
c
                do 700 n = 1, nmu
                   if (wmu.le.mu(n)) then
                      jgpmu = n
                      goto 750
                   endif
 700            continue
                jgpmu = nmu
c
 750            continue
c
c
                call imctrk2d(-1)
c
c
 450        continue
            if(kv.eq.1) then
c             write(*,*) 'end jv=',jv,' rand=',fibran()
             endif
c
c
 800  format(' rpre = ',e14.7,', zpre = ',e14.7)
 801  format(' xnu = ',e14.7)
 802  format(' phi = ',e14.7,',  wmu = ',e14.7)
 803  format(' jgpsp = ',i3,',  jgplc = ',i2)
 804  format(' jgpmu = ',i2)
 805  format(' ew = ',e14.7)
 806  format('     F = ',e14.7,', xnu = ',e14.7)
 807  format('     F = ',e14.7,', Delta (F) = ',e14.7)
 808  format('  dcen = ',e14.7)
c
c
c
      return
      end
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun 13 13:30:44 EDT 2006
c version: 2
c Name: J. Finke
c Changed 'write' statments.      
