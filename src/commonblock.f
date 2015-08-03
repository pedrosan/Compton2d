c     common blocks used in updated.f, imctrk2d.f, imcgen.f, xec2d.f, record.f
      integer ndxout
      integer cf_sentinel, nr, nz, ncycle, pair_switch
      integer ierr, myid, numprocs, master
      integer jph, kph ,jgpsp, jgplc, jgpmu
      integer nsurfu(kmax), nsurfl(kmax), nsurfo(jmax), nsurfi(jmax),
     1        nsv(jmax,kmax), npcen(jmax, kmax)
      integer ibufin(iucens), ibufout(iucens)
      integer rseed, seeds(jmax,kmax), zseeds(jmax), 
     1     rseeds(kmax), rand_switch
      integer T_const
      integer pfrl(30)
      integer incounter, outcounter,cmcount1,cmcount2,lkcount
      integer nunit_evt
      integer nmu
      integer nelectron(num_nt)
      integer nphreg, nphbins(nregmax)
      integer nphtotal, nph_lc
      integer cr_sent, upper_sent
      integer nst
      integer ntime
      integer ep_switch(jmax,kmax)
      integer dh_sentinel
      integer star_switch
      integer spec_switch
      integer inj_switch, inj_dis, g2var_switch, pick_sw
      integer FP_bcast_type, FP_Bcast_type2
      integer ti
      integer npsurfr(kmax)
      integer vol_job_type, integer_ans(int_ans_size)
      integer nfile
      integer npsurfz(jmax)
      integer randcounter
      integer idead
      integer split1, split2, split3, spl3_trg
      double precision r_flare, z_flare, t_flare, sigma_r, 
     1                 sigma_z, sigma_t, flare_amp
      double precision amxwl(jmax, kmax), gmin(jmax, kmax),
     1                 gmax(jmax, kmax), p_nth(jmax, kmax),
     2                 gbar_nth(jmax, kmax), N_nth(jmax, kmax)
      double precision dne_pa(jmax, kmax, num_nt), 
     1                 dnp_pa(jmax, kmax, num_nt)
      double precision dn_pp(jmax, kmax, num_nt)
      double precision q_turb(jmax, kmax), turb_lev(jmax,kmax)
      double precision ec_old(jmax,kmax)
      double precision n_pos(jmax,kmax,num_nt)
      double precision erini(jmax), erino(jmax), 
     1                 erinu(kmax), erinl(kmax)
      double precision erlki(jmax), erlko(jmax), 
     1                 erlku(kmax), erlkl(kmax)
      double precision ecens(jmax,kmax)
      double precision edep(jmax,kmax), prdep(jmax,kmax)
      double precision kappa_tot(n_vol, jmax, kmax),
     1                 eps_tot(n_vol, jmax, kmax), 
     2                 eps_th(n_vol, jmax, kmax), E_ph(n_vol)
      double precision Eloss_tot(jmax,kmax) ,Eloss_br(jmax,kmax),  
     1                 Eloss_cy(jmax,kmax), 
     2                 Eloss_sy(jmax, kmax), Eloss_th(jmax, kmax)
      double precision tea(jmax,kmax), tna(jmax,kmax),
     1                 B_field(jmax,kmax), n_e(jmax, kmax)
      double precision z(jmax), r(kmax), rmin, zmin, vol(jmax,kmax),
     1                 zsurf(jmax,kmax)
      double precision Asurfu(kmax), Asurfl(kmax), 
     1                 Asurfi(jmax), Asurfo(jmax)
      double precision time, dt(2), tstop, mcdt, t0(ntmax), t1(ntmax)
      double precision f_pair(jmax,kmax)
      double precision Pnt(jmax,kmax,num_nt), gnt(num_nt)
      double precision f_nt(jmax,kmax,num_nt), dg
      double precision E_IC(num_nt)
      double precision n_field(nphfield,jmax,kmax), E_field(nphfield)
      double precision xnu, wmu, phi, rpre, zpre, dcen, ew
      double precision dbufin(ducens), dbufout(ducens)
      double precision n_ph(n_gg, jmax, kmax), k_gg(n_gg, jmax, kmax), 
     1                 E_gg(n_gg)    
      double precision mu(nmumax) 
      double precision P_ref(n_ref, n_ref), W_abs(n_ref,n_ref),
     1                 E_ref(n_ref), F_ib, dE_in, ref_delay, 
     2                 f_ref, A_in, n_disk
      double precision Ed_abs(kmax), Ed_ref(kmax), Ed_in(kmax)
      double precision Ephmin(nregmax), Ephmax(nregmax)
      double precision hu(nphomax+1), Elcmin(nphlcmax), Elcmax(nphlcmax)
      double precision tbbu(kmax, ntmax), tbbl(kmax, ntmax), 
     1                 tbbi(jmax, ntmax), tbbo(jmax, ntmax)
      double precision Rstar, dist_star
      double precision inj_g1, inj_g2, inj_p, inj_t, inj_gg,
     1                 inj_sigma, inj_L, inj_v, g_bulk, pick_rate
      double precision r_esc, r_acc
      double precision R_blr, fr_blr, R_ir, fr_ir, R_disk, d_jet 
      double precision F_IC(num_nt, nphfield)
      double precision Te_new(jmax, kmax),dt_new, lnL, dT_max
      double precision E_tot_old, E_tot_new, hr_total,
     1                 hr_st_total, f_t_implicit, dr, dz
      double precision lib_dg_ce(num_nt,num_temp_max),
     1                 lib_dg_cp(num_nt, num_temp_max),
     2                 lib_disp_ce(num_nt, num_temp_max),
     3                 lib_disp_cp(num_nt, num_temp_max),
     4                 rate_gm1(num_nt), tea_min, tea_max,
     5                 tna_min, tna_max
      double precision ewsurfu(kmax), ewsurfl(kmax), ewsurfo(jmax), 
     1                 ewsurfi(jmax), ewsv(jmax, kmax)
      double precision double_ans(double_ans_size)
      double precision E_file(nfmax), I_file(nfmax), P_file(nfmax),
     1                 F_file(nfmax), alpha(nfmax), a1(nfmax), int_file
      double precision comp0(201), enxtab(13,66),
     1                 enx_nth(26,6,8,66)
      double precision randlist(randmax)
      double precision T_sum(jmax,kmax), time_sum
      double precision edout(nmumax, nphlcmax),
     1                 fout(nmumax, nphomax)
      double precision fac_old, Emiss_tot, Emiss_old
      double precision t_bound
      character *30 temp_file
      character *30 spname, phname, eventfile, lcname(nmumax)
      character *30 i_fname(jmax, ntmax), o_fname(jmax, ntmax),
     1              u_fname(kmax, ntmax), l_fname(kmax, ntmax)
      real etotal, etotal_old

      common / ndx / ndxout
      common / i_flare / cf_sentinel
      common / nc / ncycle
      common / izones / nz, nr
      common / ps / pair_switch
      ! MPI block not to be recorded
      common / MPI / ierr, myid, numprocs, master
      ! iphoton block not to be recorded
      common / iphoton / jph, kph, jgpsp, jgplc, jgpmu
      common / ph_numbers / nsurfu, nsurfl, nsurfo, nsurfi, nsv, 
     1                      npcen
      common / ibuffer / ibufin, ibufout
      common / random / rseed, seeds, zseeds, rseeds, rand_switch
      common / tc / T_const
      common / frl / pfrl
c     counter block not to be recorded
      common / counter / incounter,outcounter,cmcount1,cmcount2,lkcount
      common / event_pointer / nunit_evt
      common / imubins / nmu
      common / nel / nelectron
      common / ieb_setup / nphreg, nphbins
      common / ienergy / nphtotal, nph_lc
      common / icr / cr_sent, upper_sent
      common / phnumber / nst
      common / itimes / ntime
      common / izq / ep_switch
      common / idh / dh_sentinel
      common / star_sw / star_switch
      common / spec_sw / spec_switch
      common / injswi / inj_switch, inj_dis, g2var_switch, pick_sw
      common / Bcast / FP_bcast_type, FP_bcast_type2
      common / timeindex / ti
      common / imcsurfr / npsurfr
      common / vol_var_type  / vol_job_type, integer_ans
      common / inspi / nfile
      common / imcsurfz / npsurfz
      common / rcount / randcounter
      common / itrk / idead
      common / split / split1, split2, split3, spl3_trg

      common / d_flare / r_flare, z_flare, t_flare, sigma_r, 
     1                   sigma_z, sigma_t, flare_amp
      common / nontherm / amxwl, gmin, gmax, p_nth, gbar_nth, N_nth
      common / annihil / dne_pa, dnp_pa
      common / dnpp / dn_pp
      common / turbulence / q_turb, turb_lev
      common / ecold / ec_old
      common / n_pair / n_pos
      common / s_energies / erini, erino, erinu, erinl
      common / s_leakage / erlki, erlko, erlku, erlkl
      common / cens_energy / ecens
      common / deposition / edep, prdep
      common / vol_em / kappa_tot, eps_tot, eps_th, E_ph, Eloss_tot,
     1                  Eloss_br, Eloss_cy, Eloss_sy, Eloss_th
      common / zone_quantities / tea, tna, B_field, n_e
      common / zones / z, r, rmin, zmin, vol, Asurfu, Asurfl,
     1                 Asurfi, Asurfo, zsurf
      common / times / time, dt, tstop, mcdt, t0, t1
      common / fpair / f_pair
      common / Pnth / Pnt, gnt, f_nt, dg
      common / IC_track / E_IC 
      common / photonfield / n_field, E_field
c     7 photon variables do not need to be recorded
      common / photon / xnu, wmu, phi, rpre, zpre, dcen, ew 
      common / dbuffer / dbufin, dbufout
      common / gg_abs / n_ph, k_gg, E_gg
      common / mubins / mu
      common / reflection / P_ref, W_abs, E_ref, F_ib, dE_in,
     1                      ref_delay, f_ref, A_in, n_disk
      common / ddh / Ed_abs, Ed_ref, Ed_in
      common / eb_setup / Ephmin, Ephmax
      common / energy / hu, Elcmin, Elcmax
      common / bdtemp / tbbu, tbbl, tbbi, tbbo
      common / star / Rstar, dist_star
      common / inj / inj_g1, inj_g2, inj_p, inj_t, inj_gg,
     1               inj_sigma, inj_L, inj_v, g_bulk, pick_rate
      common / escacc / r_esc, r_acc
      common / extrad / R_blr, fr_blr, R_ir, fr_ir, R_disk, d_jet
      common / fic / F_IC
      common / d_update / Te_new, dt_new, lnL, dT_max
      common / to_fp_calc / E_tot_old, E_tot_new, hr_total, 
     1                      hr_st_total, f_t_implicit, dr, dz
c     these coulomb_sc not actually used
      common / coulomb_sc / lib_dg_ce, lib_disp_ce, lib_dg_cp, 
     1                      lib_disp_cp, rate_gm1, tea_min, tea_max,
     2                      tna_min, tna_max
      common / eweights / ewsurfu, ewsurfl, ewsurfi, ewsurfo, ewsv
      common / vol_var_type_dbl / double_ans
      common / insp / E_file, I_file, P_file, F_file, alpha, a1,int_file
      common / ctot / comp0, enxtab, enx_nth
      common / rlist / randlist
      common / T_sums / T_sum, time_sum
      common / outputs / edout, fout
      common / fac_o / fac_old, Emiss_tot, Emiss_old
      common / t_esc / t_bound

      common / tf / temp_file
      common / filenames / spname, phname, eventfile, lcname
      common / s_fnames / i_fname, o_fname, u_fname, l_fname

c     does not record elapse time. It is process dependent
      common / elapse_time / etotal, etotal_old

















