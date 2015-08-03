c
c       This subroutine calculates the emissivities
c       (cumulative probability distributions)
c       and absorption opacities due to thermal
c       bremsstrahlung, cyclotron, non-thermal
c       synchrotron, and pair annihilation
c       emission. (MB/03/Dec/99)
c
c
      subroutine volume_em(j, k, T_keV, ne_local, B, l_min)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      integer j, k
      double precision T_keV, ne_local, B, l_min
c
      integer i, m, n_harmonics, i_el, i_pos
c
      double precision mm, f_m
      double precision P(n_vol), P_sum, Theta, kappa_C
      double precision dE, T, nu, E, x, y, G_ff, sum, kappa_br
      double precision kappa_sy, kappa_cy, j_sy, j_br, j_cy
      double precision sum_k, a, g, dg1, sd, sd_k, p_1
      double precision N_e_nt, f_cy
      double precision sum_th, P_th(n_vol), tau_tot, j_th
      double precision E_m, nu_m, D_m, x_br, g_av, gamma_r, f_rz
      double precision K2, nu_c, v, nu_min, B_nu, f_rel_br
      double precision ge, gp, vdsigma, j_pa, sum1_pa, eps
      double precision nu_p, gamma_bar, McDonald
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i2, i3, i_ab, F_sync_fit_flag
      double precision F_sync_x(36), F_sync_t(36), F_sync,
     1       gamma0(num_nt), nu_b, face, Ub, facg, tt,
     1       eq43, eq13, expk43, expk13, ff, es, gamp(num_nt)
      real(8), parameter ::  sigmaT = 6.6524616d-25
      real(8), parameter ::  hplanck= 6.626075d-27
      real(8), parameter ::  emc2   = 8.187111059d-7
      real(8), parameter ::  ee     = 4.803d-10
      real(8), parameter ::  em     = 9.109d-28
c
c
c
      n_harmonics = 5
c
      nu_b = ee*B/(2*pi*em*c_light)
      Ub = B**2.d0/(8.d0*pi)
      face = 3.d0**1.5d0*sigmaT*c_light*Ub/(pi*nu_b)
      ! gamma*p (p==momentum), needed for the self-abs. frequency
      gamma0(:) = gnt(:)+1.d0
      gamp(:) = gamma0(:)*sqrt(gamma0(:)**2.d0-1.d0)

c      dE = dexp(dlog(5.d8)/dble(n_vol))
      dE = dexp(dlog(1.d20)/dble(n_vol)) ! Xuhui
      T = 1.16d7*T_keV
      Theta = T_keV/5.11d2
      kappa_C = 6.65d-25*ne_local
      if (Theta.lt.2.d-1) then
         K2 = 1.2533d0*dsqrt(Theta)*(1. + 1.875d0*Theta 
     1       + 8.2031d-1*(Theta**2) - 2.03d-1*(Theta**3))/dexp(Theta)
c      else
c         K2 = 2.*(Theta**2)
c      endif
      else
         K2 = McDonald(2.d0, (1./Theta))
      endif
      nu_c = 2.8d6*B
      nu_min = dble(n_harmonics)*nu_c
      nu_p = 9.d3*dsqrt(ne_local)
c
      if (Theta.gt.0.1) then
         f_rel_br = 1.41d0*dsqrt(Theta)*(dlog(2.d0*Theta) + 9.228d-1)
     1            - 1.d0
         f_rel_br = 1.d0 + (Theta**2)*f_rel_br/(1.d0 + (Theta**2))
         if (f_rel_br.lt.1.d0) f_rel_br = 1.d0
      else
         f_rel_br = 1.d0
      endif
c
      dg1 = 1.05d0
      p_1 = 1. - p_nth(j,k)
      if ((amxwl(j,k).gt.9.9999d-1).or.
     1    (gmax(j,k).lt.1.01*gmin(j,k))) then
         N_e_nt = 0.d0
      else if (dabs(p_1).gt.1.d-3) then
         N_e_nt = ne_local*(1. - amxwl(j,k))*p_1
     1           /(gmax(j,k)**p_1 - gmin(j,k)**p_1)
      else
         N_e_nt = ne_local*(1. - amxwl(j,k))/log(gmax(j,k)/gmin(j,k))
      endif
c
      P_sum = 0.
      sum_th = 0.
      Eloss_cy(j,k) = 0.
      Eloss_th(j,k) = 0.
c
      E = 1.d-10/dE ! Xuhui

      x_br = 7.353d1*T_keV
c
      g_av = gamma_bar(Theta)
c
      gamma_R = 2.1d-3*sqrt(ne_local)/(B*dsqrt(g_av))
      y = gamma_R/g_av
      if (y.lt.1.d2) then
         f_rz = dexp(-y)
      else
         f_rz = 0.
      endif
c
      i_ab = 0
      do 50 i = 1, n_vol
         E = E*dE
         E_ph(i) = E
         nu = 2.41487d17*E
c
c           Thermal bremsstrahlung
c
         y = E/T_keV
c
         if ((x_br.lt.1.d8).and.(x_br.gt.1.d-8).and.
     1       (y.gt.1.d-8).and.(y.lt.1.d8)) then
            if (y.ge.1.d0) then
               if (x_br.gt.1.d0) then
                  G_ff = 9.772d-1/dsqrt(y)
               else
                  if (x_br*y.lt.1.d0) then
                     G_ff = 1.d0
                  else
                     G_ff = 3.464d0/dsqrt(x_br*y)
                  endif
               endif
            else
               if (x_br.lt.1.) then
                  if (y.gt.dsqrt(x_br)) then
                     G_ff = 1.d0
                  else
                     G_ff = 2.7566d0*dlog(1.5802d1*x_br/y)
                  endif
               else
                  G_ff = 5.513d-1*dlog(6.9298d0/y)
               endif
            endif
         else
            G_ff = 0.
         endif
c
         if (y.lt.150.) then
            j_br = G_ff*dexp(-y)*1.6416d-20*(ne_local**2)*amxwl(j,k)
     1            *f_rel_br/dsqrt(T)
         else
            j_br = 0.
         endif
c
         kappa_br = 2.6d-44*(ne_local**2)*G_ff*amxwl(j,k)
     1             *(1. - dexp(-y))*f_rel_br/((T**1.5)*(E**3))
c
c           Non-thermal synchrotron
c
c         if ((nu.le.nu_p).or.(amxwl(j,k).gt.9.999d-1)) then
         if (nu.le.nu_p)then
            j_sy = 0.
            kappa_sy = 0.
            goto 30
         endif
         sum = 0.
         sum_k = 0.
         a = 5.75d10*1.5d0*E/B  ! a/g**2 is the 'x' in R&L book. Xuhui 4/2/09
c        1.5d0 factor is for the sine term in the eq (6.17c) of R&L.
         F_sync_fit_flag = 0
         if(F_sync_fit_flag.eq.1)then
ccccccccc         c             c               c             c
          g = gmin(j,k)
 20       y = a/(g**2) + g/gmax(j,k) + gamma_R/g
          if (y.lt.150.) then
             sd = dexp(-y)/(g**(.66666667 + p_nth(j,k)))
             sd_k = dexp(-y)/(g**(1.6666667 + p_nth(j,k)))
          else
             sd = 0.
             sd_k = 0.
          endif
          sum = sum + g*(dg1 - 1.)*sd
          sum_k = sum_k + g*(dg1 - 1.)*sd_k
          g = g*dg1
          if (g.lt.gmax(j,k)) goto 20
          j_sy = 6.34d-7*(nu**.3333333)*(B**.66666667)*sum*N_e_nt
          kappa_sy = 1.45d3*(B**.6666667)*(2. + p_nth(j,k))
     1              *N_e_nt*sum_k/(nu**1.6666667)

ccccccccc          c            c               c              c
         else
c          do i2=1,num_nt-1
c            gamma0(i2) = gnt(i2) + 1.d0
c            y = a/(gamma0(i2)*gamma0(i2))
c            if(y.le.0.001d0)then
c              F_sync = 2.7083d0*(y/2)**0.3333333
c            else if(y.ge.10.d0)then
c              F_sync = 1.2533d0*dexp(-y)*dsqrt(y)
c            else
c                do i3=1,36
c                  if(y.lt.F_sync_x(i3)) exit
c                enddo
c                F_sync =(y-F_sync_x(i3-1))/(F_sync_x(i3)-F_sync_x(i3-1))
c     1                 *(F_sync_t(i3)-F_sync_t(i3-1))+F_sync_t(i3-1)
c            endif
c
            do i2=1,num_nt-1
              facg = 3.d0*(gamma0(i2)**2.d0)*nu_b
              tt = nu/facg
              if( tt < 1.0d4 ) then
                 eq43 = expk43(tt)
                 eq13 = expk13(tt)
                 ff = tt*tt*(eq43*eq13 - 0.6*tt*(eq43-eq13)*(eq43+eq13))
                 es = face*ff*dexp(-2.d0*tt)
              else
                 es = 0.0
              end if
c            sd = f_nt(j,k,i2)*F_sync
c            do i3=1,num_nt-1
c               if(gnt(i3).ge.(gnt(i2)-E/5.11d2)) exit
c            enddo
c            sd_k=(f_nt(j,k,i3)*gamma0(i2)**2/gamma0(i3)**2-f_nt(j,k,i2))
c     1            *F_sync
c            if(i3.gt.2)sd_k = sd_k +
c     1      (1/gnt(i3)**2 - 1/(gnt(i2)-E/5.11d2)**2)/
c     1      (1/gnt(i3)**2 - 1/gnt(i3-1)**2)
c     1      *(f_nt(j,k,i3-1)*gamma0(i3)**2/gamma0(i3-1)**2-f_nt(j,k,i3))
c     1      *F_sync ! grid correction Xuhui

            sd = f_nt(j,k,i2)*es
            sd_k = gamp(i2)*es
            sum = sum + (gnt(i2+1)-gnt(i2))*sd
            sum_k = sum_k +(f_nt(j,k,i2)/gamp(i2)-f_nt(j,k,i2+1)/
     1              gamp(i2+1))*sd_k
          enddo
c          j_sy = 9.01d-6*sum*ne_local
c          kappa_sy =9.01d-6*sum_k*ne_local
c     1              /2.41487d17*5.400d45/nu/nu/nu ! Xuhui sync
          j_sy = sum*ne_local/(4.d0*pi)
          kappa_sy = sum_k*ne_local/(8.d0*pi*em*nu**2)


         endif
c
         if(kappa_sy .lt. 0.)then
            write(*,*)'kappa_sy < 0 at i =',i, kappa_sy
            kappa_sy = -1.d0*kappa_sy
         endif

c
c
c           Thermal cyclotron
c
 30      f_m = 1.d0
         j_cy = 0.
         kappa_cy = 0.
c
c        Sum first n harmonics
c
         if (nu.le.nu_p) then
            j_cy = 0.d0
            kappa_cy = 0.d0
            goto 26
         endif
c
         do 25 m = 1, n_harmonics
            mm = dble(m)
            f_m = f_m/(4.*mm)
            nu_m = mm*nu_c
            E_m = 4.14d-18*nu_m
            D_m = 7.07d-1*Theta*E_m
            x = ((E - E_m)/D_m)**2
            y = E_m/T_keV
c
            if (x.lt.50.) then
               f_cy = f_rz*dexp(-x)*ne_local*(B**2)*(Theta**(mm-1.5d0))
     1               *(mm + 1.d0)*f_m*(mm**(2.d0*mm + 1.d0))
               j_cy = j_cy + 8.46d-14*f_cy*(E**2)/(E_m**3)
               if (y.lt.150.) then
                  kappa_cy = kappa_cy + 5.705d33*(dexp(y) - 1.d0)
     1                                 *f_cy/(nu*(nu_m**3))
               else
                  x = y - x
                  if (x.gt.150.) then 
                     kappa_cy = 1.d70
                  else if (x.gt.-100.) then
                     kappa_cy = kappa_cy+f_rz*5.705d33*dexp(x)*ne_local
     1                         *(B**2)*(Theta**(mm - 1.5d0))*f_m
     2                         *(mm + 1.d0)*(mm**(2.d0*mm + 1.d0))
     3                         /(nu*(nu_m**3))
                  endif
               endif
            endif
  25     continue
c
c        Use Mahadevan, Narayan & Yi (1996) formula
c          for the higher harmonics (nu > n*nu_c)
c
         if (nu.gt.nu_min) then
            v = nu/(nu_c*(Theta**2))
            y = 4.5*v
            if (y.lt.1.d6) then
               j_cy = j_cy + 4.652d-12*ne_local*nu/(K2*(v**1.6666667d-1)
     1                      *dexp(y**3.33333d-1))
            endif
            y = E/T_keV
            if (y.gt.100.) then
               B_nu = 1.d-50
            else if (y.lt.1.d-6) then
               B_nu = 3.56d-30*(nu**3)/y
            else
               B_nu = 3.56d-30*(nu**3)/(dexp(E/T_keV) - 1.d0)
            endif
            if (y.lt.100.) kappa_cy = kappa_cy + j_cy/B_nu
         endif
  26     continue
c
c
c             Pair annihilation radiation
c
c
         j_pa = 0.d0
         if ((pair_switch.eq.0).or.(f_pair(j,k).lt.1.d-10)) goto 47
         eps = 1.957d-3*E
c
         do 45 i_el = 1, num_nt-1
            if (f_nt(j, k, i_el).lt.1.d-10) goto 45
            ge = gnt(i_el) + 1.d0
            sum1_pa = 0.d0
            do 40 i_pos = 1, num_nt-1
               if (n_pos(j, k, i_pos).lt.1.d-10) goto 40
               gp = gnt(i_pos) + 1.d0
c
               sum1_pa = sum1_pa + (gnt(i_pos+1) - gnt(i_pos))
     1                      *n_pos(j, k, i_pos)*vdsigma(eps, ge, gp)
  40        continue
            j_pa = j_pa + (gnt(i_el+1) - gnt(i_el)) 
     1                   *ne_local*f_nt(j, k, i_el)*sum1_pa
  45     continue
         j_pa = j_pa*eps*1.6d-9
c
c
c  47     kappa_tot(i, j, k) = kappa_br + kappa_sy + kappa_cy
  47     kappa_tot(i, j, k) = kappa_sy
c        Xuhui 6/1/09 deactivated br and cy absorption

c
c         if (kappa_tot(i, j, k).lt.(1.d1*kappa_C)) then
         if (kappa_tot(i, j, k).lt.dmax1(1.d0/l_min,
     1       1d1*kappa_C))then ! Xuhui 5/8/09
c            P_sum = P_sum + (j_br + j_sy + j_cy + j_pa)*E*(dE - 1.d0)
            P_sum = P_sum + j_sy*E*(dE - 1.d0)
c         Xuhui 6/1/09; deactivated any spectrum except synchrotron
            Eloss_cy(j, k) = Eloss_cy(j, k) + j_cy*E*(dE - 1.d0)
         else
            i_ab = i ! Xuhui 11/16/10
            x = E/T_keV
            tau_tot = kappa_tot(i, j, k)*l_min
            if (x.lt.1.d2) then
               j_th = 1.47d-47*(nu**3)/(dexp(x) - 1.d0)
            else
               j_th = 1.d-50
            endif
            if (tau_tot.lt.5.d1) 
     1         j_th = j_th*(1.d0 - dexp(-tau_tot))
            sum_th = sum_th + j_th*E*(dE - 1.d0)
            Eloss_th(j, k) = Eloss_th(j, k) + j_th*E*(dE - 1.d0)
         endif
         P(i) = P_sum
         P_th(i) = sum_th
c
  50  continue
      if(j.eq.1.and.k.eq.1.and.i_ab.gt.0)then
            write(*,*)'kappa_tot is big at E=',E_ph(i_ab),
     1               'kappa =',kappa_tot(i_ab,j,k)
      endif

c
      do 100 i = 1, n_vol
         if (P_sum.gt.1.d-50) then
            eps_tot(i, j, k) = P(i)/P_sum
         else
            eps_tot(i, j, k) = 0.
         endif
         if (sum_th.gt.1.d-50) then
            eps_th(i, j, k) = P_th(i)/sum_th
         else
            eps_th(i, j, k) = 0.
         endif
c
 100  continue
c
c
      return
      end
c
c
c        This subroutine calculates the gamma-gamma
c        pair production opacity at photon energy E
c
c
      double precision function kgg_calc(E, j, k)
      implicit none
      include 'general.pa'
      include 'commonblock.f'
      double precision E
      integer j, k
c

      integer i
      double precision sum_mu, mu_thr, mu_local, dmu, eps1, eps2, beta,f
c
      kgg_calc = 0.d0
      eps2 = 1.957d-3*E
c
      do 200 i = 1, n_gg-1
         eps1 = 1.957d-3*E_gg(i)
         if ((n_ph(i, j, k).gt.1.d-20).and.(eps1.gt.1.d0/eps2)) then
            sum_mu = 0.
            mu_thr = 1.d0 - 2.d0/(eps1*eps2)
            if (mu_thr.gt.1.d0) mu_thr = 1.d0
            dmu = .01*(1. + mu_thr)
            mu_local = 5.d-1*dmu - 1.d0
c
 100        continue
            beta = dsqrt(1.d0 - 2.d0/(eps1*eps2*(1.d0 - mu_local)))
            if ((beta.ge.1.d0).or.(beta.le.0.d0)) goto 150
            f = (1.d0 - beta**2)*((3.d0 - beta**4)
     1          *dlog((1.d0 + beta)/(1.d0 - beta))
     2          - 2.d0*beta*(2.d0 - beta**2))
            sum_mu = sum_mu + (1.d0 - mu_local)*f*dmu
 150        mu_local = mu_local + dmu
            if (mu_local.lt.mu_thr) goto 100
c
            kgg_calc = kgg_calc 
     1               + sum_mu*n_ph(i, j, k)*(E_gg(i+1) - E_gg(i))
         endif
 200  continue
c
      kgg_calc = kgg_calc*6.234d-26
      return
      end
c
c
c       This subroutine calculates the photon-energy-
c      differential cross section for pair annihilation
c
c
      double precision function vdsigma(eps, ge, gp)
      double precision eps, ge, gp
c
      double precision be, bp, gcm_u, gcm_l, f_vds
      double precision gcm2, gcmstar, gcmmax, eps_u, eps_l
c
      if ((ge.lt.1.000001d0).or.(gp.lt.1.0000001d0)) then
         vdsigma = 0.d0
         goto 900
      endif
      be = dsqrt(1.d0 - 1.d0/(ge**2)) + 1.d-10
      bp = dsqrt(1.d0 - 1.d0/(gp**2)) + 1.d-10
      eps_u = 5.d-1*(gp*(1.d0 + bp) + ge*(1.d0 + be))
      eps_l = 5.d-1*(gp*(1.d0 - bp) + ge*(1.d0 - be))
      if ((eps.le.eps_l).or.(eps.ge.eps_u)) then
         vdsigma = 0.d0
         goto 900
      endif
c
      gcm2 = 5.d-1*(1.d0 + ge*gp*(1.d0 - be*bp))
      if (gcm2.le.1.00001d0) then
         vdsigma = 0.d0
         goto 900
      endif
      gcm_l = dsqrt(gcm2)
c
      gcm2 = 5.d-1*(1.d0 + ge*gp*(1.d0 + be*bp))
      if (gcm2.le.1.d0) then
         vdsigma = 0.d0
         goto 900
      endif
      gcmmax = dsqrt(gcm2)
      gcm2 = eps*(ge + gp - eps)
      if (gcm2.le.1.00001d0) then
         vdsigma = 0.d0
         goto 900
      endif
      gcmstar = dsqrt(gcm2)
      gcm_u = dmin1(gcmstar, gcmmax)
c
      if (gcm_u.gt.(1.0001d0*gcm_l)) then
         vdsigma = 7.48d-15*(f_vds(gcm_u, ge, gp, eps) 
     1                     - f_vds(gcm_l, ge, gp, eps))
     2                     /(be*bp*((ge*gp)**2))
      else
         vdsigma = 0.d0
      endif
c
 900  return
      end
c
c
c
      double precision function f_vds(gcm, ge, gp, eps)
      double precision gcm, ge, gp, eps
c
      double precision H_pa, D
c
      D = (ge + gp)**2 - 4.*(gcm**2)
      if (D.le.1.d-20) then
         f_vds = 0.d0
         goto 900
      endif
c
      f_vds = dsqrt(D) + H_pa(gcm, ge, gp, eps) 
     1                 + H_pa(gcm, gp, ge, eps)
c
 900  return
      end
c
c
c
      double precision function H_pa(gcm, ge, gp, eps)
      double precision gcm, ge, gp, eps
c
      double precision c, d, u, I_pa, u2, gcms2, gstar
c
c
      c = (ge - eps)**2 - 1.d0
      d = ge*(gp + ge) + eps*(gp - ge)
c
      gcms2 = eps*(ge + gp - eps)
      if (gcms2.lt.1.00001d0) then
         H_pa = 0.d0
         goto 900
      endif
      gstar = dsqrt(gcms2)
      u2 = c*(gcm**2) + gcms2
      if (u2.lt.1.d-20) then
         H_pa = 0.d0
         goto 900
      endif
      u = dsqrt(u2)
c
      if (dabs(c).gt.1.d-8) then
         H_pa = (2.d0 + (1.d0 - gcms2)/c)*I_pa(c, gcm, gstar, u)
     1        + (1.d0/gcm - gcm/c + 5.d-1*gcm*(2.d0*c - d)/gcms2)/u
     2        + gcm*u/c
      else
         H_pa = (2.d0*(gcm**3)/3.d0 + 2.d0*gcm + 1.d0/gcm)/gstar
     1        + 5.d-1*(2.d0*(gcm**3)/3.d0 - d*gcm)/(gstar**3)
      endif
c
 900  return
      end
c
c
c
      double precision function I_pa(c, gcm, gcmstar, u)
      double precision c, gcm, gcmstar, u
c
      if (c.ge.1.d-8) then
         I_pa = dlog(gcm*dsqrt(c) + u)/dsqrt(c)
      else if (c.le.-1.d-8) then
         I_pa = asin(gcm*dsqrt(-c)/gcmstar)/dsqrt(-c)
      else
         I_pa = 0.d0
      endif
c
 900  return
      end
c
c
c
      double precision function gamma_bar(Theta)
      double precision Theta
c
      double precision sum, sd, t, ts, s, dt, K2, K3
      double precision McDonald
c
      if (Theta.lt.0.2) then
         gamma_bar = (1. + 4.375*Theta + 7.383*(Theta**2) 
     1              + 3.384*(Theta**3))/(1. + 1.875*Theta 
     2              + .8203*(Theta**2)) - Theta
         goto 900
      endif
c
      K2 = McDonald(2.d0, (1.d0/Theta))
      K3 = McDonald(3.d0, (1.d0/Theta))
c
      gamma_bar = K3/K2 - Theta
c
 900  continue
      if (gamma_bar.lt.1.d0) gamma_bar = 1.d0
c
      return
      end
c
c
c
      double precision function McDonald(nu, z)
      double precision nu, z
c
      double precision sum, sd, t, ts, s, dt, d, y, a,
     1                 GammaF
c
      sum = 0.d0
      t = 1.d0
      dt = 1.001d0
      d = dt - 1.d0
      s = 5.d-1*(1.d0 + dt)
      a = nu - 5.d-1
c
 100  ts = t*s
      y = z*ts
      if (y.lt.2.25d2) then
         sd = ((ts*ts - 1.d0)**a)/dexp(y)
      else
         sd = 0.d0
      endif
      sum = sum + d*t*sd
      t = t*dt
      if ((t.lt.2.d0).or.(sd.gt.1.d-8)) goto 100
c
      McDonald = dsqrt(3.14159265d0)*((5.d-1*z)**nu)*sum
     1          /GammaF(5.d-1 + nu)
c
      return
      end
c
c
c
c     This calculates the gamma function of xx.  It uses
c     the function gammln from "Numerical Recipies".  It is
c     called by McDonald.
      Function GammaF(xx)
      implicit none
      double precision xx, GammaF, gammln
      GammaF = exp(gammln(xx))
      return
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c  This function is from "Numerical Recipies in FORTRAN",
c  second edition, page207.  It calculates the ln of the
c  gamma function of xx.

      FUNCTION gammln(xx)
      implicit none
      DOUBLE PRECISION gammln, xx
      INTEGER j
      DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
      SAVE cof, stp
      DATA cof, stp/76.18009172947146d0,-86.50532032941677d0,
     1    24.01409824083091d0,-1.231739572450155d0, 
     2    .1208650973866179d-2,-.5395239384953d-5,
     3    2.5066282746310005d0/
       x=xx
       y=x
       tmp=x+5.5d0
       tmp=(x+0.5d0)*log(tmp)-tmp
       ser=1.000000000190015d0
       do 11 j=1,6
          y=y+1.d0
          ser=ser+cof(j)/y
 11    continue
       gammln=tmp+log(stp*ser/x)
       return
       END
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        double precision function expk13(t)

!-------------------------------------------------------------------------------
!       Returns exp(t)*(Modified Bessel function order 1/3)(t)
!-------------------------------------------------------------------------------

        implicit none
        include 'general.pa'

        !! save c1,c2
        double precision  z3, zs, z, t, poly
        double precision  c1, c2
        double precision  f1, f2

        data c1,c2 / 0.35502805 , 0.25881940 /

        if( t .le. 1.0 ) then 
           !--------------------------------------------------
           ! Small argument use Airy function expansion A&S p446
           !
           z3 = 1.5d0*t
           zs = z3**0.3333333333333333d0
           z  = zs*zs
           z3 = z3*z3

           f1 = 1.0d0 + z3/6.0d0*(1.0d0 + z3/30.0d0*(1.0d0 + z3/56.0d0))
           f2 = z*(1.0d0 + z3/12.0d0*(1.0d0 
     &          + z3/42.0d0*(1.0d0 + z3/90.0d0)))
           expk13 = exp(t)*pi*1.7320508d0/zs*(c1*f1 - c2*f2)

        else 
           !--------------------------------------------------
           ! Large arguments use asymptotic expansion
           !
           z   = 1.0d0/(72.0d0*t)
           poly= 1.0d0 - 5.0d0*z*(1.0d0 - 38.5d0*z)
           expk13 = dsqrt(0.5d0*pi/t)*poly
     &              /(1.0d0 + 1.0d0/(1.0d0 + 58.0d0*t*t))

        end if 

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        double precision function expk43(t)

!-------------------------------------------------------------------------------
!      Returns exp(t)**(modified Bessel function order 4/3)(t)
!-------------------------------------------------------------------------------
        
        implicit none
        include 'general.pa'

        real(8) :: z, poly, t

        if( t .le. 1.0 ) then 
           !--------------------------------------------------
           ! Small argument use fit.
           !
           poly   = 1.0d0 + t*(0.9757317d0 - 7.6790616d-2*t)
           expk43 = 0.44648975d0*(2.0d0/t)**1.333333333d0*poly

        else 
           !--------------------------------------------------
           ! Large arguments use asymptotic argument.
           !
           z    = 1.0/(72.0*t)
           poly = 1.0+55.0*z*(1.0-8.5*z)
           expk43 = dsqrt(0.5*pi/t)*poly*(1.0 + 1.0/(1.0 + 50.0*t*t))

        end if

        return
        end
c
cccccccccccccccccccccccccccccccccccccccccccccccccc 
c Tue Jun  6 23:56:12 EDT 2006
c version: 2
c Name: J. Finke
c Changed a few double constants with 'e' to 
c 'd'. Changed variable 'm' to an integer from 
c double.        
