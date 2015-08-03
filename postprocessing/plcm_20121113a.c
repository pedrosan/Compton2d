#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define ch_max 20
#define st_max 1024
#define mu_max 10

#define qu(x) ((x)*(x))

/* ------------------------------------------------------------ */
double max(x, y)
double x, y;
{
   if (x > y) return (x);
         else return (y);
}

/* ------------------------------------------------------------ */
int digit(c)
char c;
{
   if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
       c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
        return (1);
   else return (0);
}

/* ------------------------------------------------------------ */
void add_dat(str)
char *str;
{
   int len;
   
   len = strlen(str);
   if (str[len-4] != '.') {
       str[len] = '.';
       str[len+1] = 'd';
       str[len+2] = 'a';
       str[len+3] = 't';
       str[len+4] = '\0';
   }
   return;
}

/* ------------------------------------------------------------ */
double dinput(x)
double x;
{
   char a[30];
   
   gets(a);
   if (strlen(a)) return (atof(a));
             else return (x);
}

/* ------------------------------------------------------------ */
int iinput(x)
int x;
{
   char a[30];
   
   gets(a);
   if (strlen(a)) return (atoi(a));
             else return (x);
}

/* ------------------------------------------------------------ */
void file_error(fname)
char *fname;
{
   printf ("\n\n File error: Can not open file '%s'! - Exit.\n\n", fname);
   exit(0);
}

/* ------------------------------------------------------------ */
main()
{
double E0[ch_max], E1[ch_max], t0[st_max], t1[st_max];
double time, E, ew, dt=7.e2, F[st_max][mu_max][ch_max];
double F1[st_max][mu_max][ch_max]; // GF
double F2[st_max][mu_max][ch_max]; // GF 
double sF[st_max][mu_max][ch_max]; // GF 
double E_lower=1.e-3, E_upper=3.e-3, dE;
double Pi = 3.1415926536;
double mu0[mu_max], mu1[mu_max];
double t_offset=0.0, t_max=7.e4, t_start=0.0, t_end=1e6;
double t_bound, r, rmax=1.e16, z, mu, phi, cdt;
double gam_bulk=33.,doppler, betta;

int k, l, m, n, factor=0, regions=7, reg, n_bins=0, n_r=1, n_mu=1;
int particle_sum[st_max][mu_max][ch_max], n_file, max_file=300;

char infile[30]="p001_evb.dat", outfile1[mu_max][50],  outfile2[mu_max][50], outfile3[50], a[30];

FILE *fpin, *fpout1, *fpout2, *fpout3, *fpext, *fpcheck;

printf ("\n\n Input from event file: [%s]  ", infile);
gets(a);
if( strlen(a) ) {
   if( digit(a[0]) ) {
      for (k=0; k<strlen(a); k++) infile[2+k] = a[k];
      infile[2+k] = '\0';
   } else {
      strcpy(infile, a);
   }
   add_dat(infile);
}

// if (!(fpin=fopen(infile, "r"))) file_error(InFile);

l = strlen(infile);

printf ("\n gam_bulk = [%e]  ", gam_bulk);
gam_bulk = dinput(gam_bulk);
printf ("\n rmax = [%e]  ", rmax);
rmax = dinput(rmax);

/*printf ("\n rmax = [%e]  ", rmax);
rmax = dinput(rmax);
printf (" zmax = [%e]  ", zmax);
zmax = dinput(zmax);*/

strcpy(outfile1[0], infile);
outfile1[0][0] = 'l';
outfile1[0][1] = 'c';
outfile1[0][3] = '7';

/*
outfile1[0][l-4] = '_';
outfile1[0][l-3] = 'e';
outfile1[0][l-2] = 'v';
outfile1[0][l-1] = 't';
outfile1[0][l] = '_';
outfile1[0][l+1] = '0';
outfile1[0][l+2] = '\0'; 
*/

outfile1[0][l-5] = '0';
add_dat(outfile1[0]);

printf ("\n Output file for mu-bin 0: [%s]  ", outfile1[0]);
gets(a);
if( strlen(a) ) strcpy(outfile1[0], a);
add_dat(outfile1[0]);
l = strlen(outfile1[0]);
strncpy(outfile2[0], outfile1[0],l-4);
strcat(outfile2[0],"_aux");
strcat(outfile2[0],".dat");
//add_dat(outfile2[0]);
strncpy(outfile3, outfile1[0],l-4);
strcat(outfile3,"_particles");
strcat(outfile3,".dat");

printf("\n");
printf(" output file name 1 : %s\n",outfile1[0]);
printf(" output file name 2 : %s\n",outfile2[0]);
printf(" output file name 3 : %s\n",outfile3);

printf ("\n Number of angular bins: [%1d]  ", n_mu);
n_mu = iinput(n_mu);
printf (" ==> %1d\n",n_mu);

if (n_mu > mu_max) { 
    printf ("\n Number of angular bins set to maximum allowed number = %d\n", mu_max);
    n_mu = mu_max;
}

mu0[0] = 0.99944;
mu1[0] = 0.99964;
for (k=0; k<n_mu; k++) {
   printf ("\n----------------------------------------\n Angular bin no. %1d:", k);
   if (k) {
      strcpy(outfile1[k], outfile1[0]);
      outfile1[k][l-5] = (char)(48 + k);
      strcpy(outfile2[k], outfile2[0]);
      outfile2[k][l-5] = (char)(48 + k);
   }
   printf("\n");
   printf(" PRE output file name 1 [%d]: %s\n",k,outfile1[k]);
   printf(" PRE output file name 2 [%d]: %s\n",k,outfile2[k]);
   
   printf ("\n Output file [%s]: ", outfile1[k]);
   gets(a);
   if( strlen(a) ) strcpy(outfile1[k], a);
   add_dat(outfile1[k]);
   //printf (" ==> %s\n",outfile1[k]);

   l = strlen(outfile1[k]);
   printf("\n");
   printf(" MID1 output file name 2 [%d]: %s (%d)\n",k,outfile2[k],l);
   strncpy(outfile2[k], outfile1[k],l-4);
   printf(" MID2 output file name 2 [%d]: %s (%d)\n",k,outfile2[k],l);
   //strcat(outfile2[k],"_aux");
   //strcat(outfile2[k],".dat");
   // add_dat(outfile2[k]);

   printf("\n");
   printf(" output file name 1 [%d]: %s\n",k,outfile1[k]);
   printf(" output file name 2 [%d]: %s\n",k,outfile2[k]);
   
   printf("\n");
   printf (" mu_min(%1d) = [%e]  ", k, mu0[k]);
   mu0[k] = dinput(mu0[k]);
   printf (" ==> %e\n",mu0[k]);

   printf (" mu_max(%1d) = [%e]  ", k, mu1[k]);
   mu1[k] = dinput(mu1[k]);
   printf (" ==> %e\n",mu1[k]);

   // if (k < (mu_max - 1)) mu0[k+1] = mu1[k+1] = mu1[k];
   // GF: do not understand this!
   if( k < (mu_max - 1) ) { 
      mu0[k+1] = mu1[k+1] = mu1[k];
      printf(" for k= %d : mu0[k+1]= %f , mu1[k+1]= %f , mu1[k]= %f\n",k,mu0[k+1],mu1[k+1],mu1[k]);
   }
}

printf ("\n Time step: dt = [%e] s  ", dt);
dt = dinput(dt);
printf (" Time offset: [%e] s  ", t_offset);
t_offset = dinput(t_offset);
printf (" Stop at t = [%e] s  ", t_max);
t_max = dinput(t_max);

t_start=t_offset; //GF
t_end  =t_max;    //GF

printf ("\n\n");
printf (" Time step: dt = %e s\n", dt);
printf (" Time offset   = %e s\n", t_start);
printf (" Stop at t     = %e s\n", t_end);

t_max-=t_offset;

t0[0] = 0.;
for (k=0; k<st_max-1; k++) {
   t1[k] = t0[k+1] = t0[k] + dt;
}
t1[st_max-1] = t0[st_max-1] + dt;

/* ------------------------------------------------------------ */
m1:
printf ("\n");
printf (" Preparing energy grid\n");
printf (" Number of energy regions: [%d]  ", regions);
regions = iinput(regions);
printf (" ==> %d\n",regions);
for (reg=0; reg<regions; reg++) {
    if(reg==0)
      {E_lower = 1e-3; E_upper = 3e-3;}
    else if(reg==1)
      {E_lower = 2; E_upper = 4;}
    else if(reg==2)
      {E_lower = 9; E_upper = 15;} 
    else if(reg==3)
      {E_lower = 15; E_upper = 20;} 
    else if(reg==4)
      {E_lower = 20; E_upper = 60;} 
    else if(reg==5)
      {E_lower = 5e5; E_upper = 5e7;}
    else if(reg==6)
      {E_lower = 1e9; E_upper = 1e10;}
    else
      {E_lower = E_upper;}

    printf ("\n\n   Region no. %d:", reg+1);
    printf ("\n   Lower boundary: [%e] keV  ", E_lower);
    E_lower = dinput(E_lower);
    printf ("   Upper boundary: [%e] keV  ", E_upper);
    E_upper = dinput(E_upper);
    printf (" ==> %e - %e keV\n",E_lower,E_upper);
    printf ("   Number of energy bins in this region: [%d]  ", n_r);
    n_r = iinput(n_r);
    printf (" ==> %d\n",n_r);

    if( (n_r + n_bins) > ch_max ) {
        printf ("\n Not more than 20 energy channels!");
        printf ("\n Re-specify energy grid! \n");
        goto m1;
    }

    printf ("   Logarithmic (0) or linear (1) bins? [0]  ");
    gets(a);
    if (a[0] == '1') { 
        dE = (E_upper - E_lower)/((double)(n_r));
        E0[n_bins] = E_lower;
        for (k=n_bins; k<(n_bins + n_r - 1); k++) {
           E1[k] = E0[k+1] = E0[k] + dE;
	}
        E1[n_bins + n_r - 1] = E0[n_bins + n_r - 1] + dE;
    } else {
        dE = exp(log(E_upper/E_lower)/((double)(n_r)));
        E0[n_bins] = E_lower;
        for (k=n_bins; k<(n_bins + n_r - 1); k++) {
           E1[k] = E0[k+1] = E0[k]*dE;
	}
        E1[n_bins + n_r - 1] = E0[n_bins + n_r - 1]*dE;
    }
    n_bins+=n_r;
}
// end of reading energy regions   

printf ("\n\n----------------------------------------\n");
printf ("\n Energy binning:\n");
for (k=0; k<n_bins; k++) {
   printf (" [%d] E0 = %e keV , E1 = %e keV\n", k, E0[k], E1[k]);
}
printf ("\n");

printf ("\n Angular binning:\n");
for (n=0; n<n_mu; n++) {
   printf (" [%d] mu_min = %e , mu_max = %e", n, mu0[n], mu1[n]);
}
printf ("\n");
printf ("----------------------------------------\n");
//printf ("\n number of input files: [%1d]", max_file);
//max_file = iinput(max_file);

// Initialize energy and particle-count arrays
for (k=0; k<st_max; k++) {
   for (l=0; l<n_bins; l++) {
      for (n=0; n<n_mu; n++) {
         F[k][n][l]  = 0.;
         F1[k][n][l] = 0.;
         F2[k][n][l] = 0.;
         particle_sum[k][n][l] =0;
      }
   }
}
m = 0;

// GF: if (!(fpext=fopen("TeV.dat", "w"))) file_error("TeV.dat");
// if (!(fpcheck=fopen("particles.dat", "w"))) file_error("particles.dat");
if( !(fpcheck=fopen(outfile3, "w")) ) file_error(outfile3);
for (l=0; l<n_bins; l++) {
   fprintf (fpcheck, "# energy range %i : %e %e (keV)\n", l+1, E0[l], E1[l]);
}
fprintf (fpcheck, "# time   : %e %e %e\n",t_start,t_end,dt);
fprintf (fpcheck, "# mu     : %f %f\n",mu0[0],mu1[0]);
// fprintf (fpcheck, "# factor : %i\n", factor);

fprintf(fpcheck,"#\n");
fprintf(fpcheck,"# time, E, ew, mu, r, z, k(time), n(mu), l(energy)\n");
fprintf(fpcheck,"#\n");

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// m2:
n_file = 1;
while(n_file <= max_file) {
  
   /*
   if(infile[3]=='9'){
   infile[2]++;
   infile[3]='0';}
   else
   infile[3]++;
   n_file++;
   */
  
   // if (n_file <= max_file)
   if (!(fpin=fopen(infile, "r"))) {
      do {
         if( infile[1]!='0' || infile[2]!='0' || infile[3]!='1' ) factor++;
         infile[7]++;
         infile[1]='0';
         infile[2]='0';
         infile[3]='1';
      }
      while( (!(fpin=fopen(infile,"r"))) && infile[7]<'l' );
         if(infile[7]>='l') break;
   } 
   //file_error(infile);
   printf ("\n Input from another event file? [%s] or n \n", infile);
  
   while ((a[0]=fgetc(fpin)) != EOF) {
  
      /* it reads one photon at a time */
      fseek(fpin, -1, 1);
      fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf\n", &t_bound, &E, &ew, &r, &z, &mu, &phi);
  
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /*    if (z > zmax*(1.-1.e-7)) continue;          */
      mu      = -mu;
      betta   = sqrt(1.-1./(gam_bulk*gam_bulk));
      doppler = gam_bulk*(1.+mu*betta);
      t_bound = (t_bound-betta*z*3.33333333e-11)/doppler;
      E       = E*doppler;
      ew      = ew*doppler;
      mu      = (mu+betta)/(1.+mu*betta);
      cdt     = z*mu/gam_bulk + sqrt(1.-mu*mu)*(rmax-r*cos(phi));
      time    = t_bound + 3.33333333e-11*cdt;
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
      if (m++ < 10) { 
         printf (" t = %e, E = %e, ew = %e\n", t_bound, E, ew);
         printf (" r = %e, z = %e, mu = %e, phi = %e\n", r, z, mu, phi);
         printf (" t_corr = %e\n", time);
      }
      
      /*  cdt = (rmax - r*cos(phi))*sqrt(1. - qu(mu)) + (zmax - z)*mu; */
  
      time-=t_offset;
      if (time < 0.) continue;
  
      /* GF: this is where it bins stuff up */
      for (k=0; k<st_max; k++) {
          // printf(" CHECK: doing k %3d\n",k);
          if (time >= t0[k] && time < t1[k]) break;
          // printf(" CHECK: time found %3d\n",k);
      }
        
      for (n=0; n<n_mu; n++) {
          // printf("    CHECK: doing mu %3d\n",n);
          if (mu >= mu0[n] && mu < mu1[n]) break;
          // printf(" CHECK: mu   found %3d\n",n);
      }
        
      // GF for (l=0; l<n_bins; l++)
      // GF   if (E >= E0[l] && E < E1[l]) break;
      //
      /* GF
      if(k == 50 && l == 6 && n == 0){
         fprintf(fpext,"%e %e %e %e %e %e %e\n", t_bound, E, ew, r, z, mu, phi);
      } */
      /* GF
      if (k < st_max && l < n_bins && n < mu_max) {
          F[k][n][l]+=ew;
          F2[k][n][l]+=(ew*ew);
          particle_sum[k][n][l]++;  
      }
      */ 
  
      // printf(" CHECK: time %d , mu %d\n",k,n);
      for (l=0; l<n_bins; l++) {
        if( E >= E0[l] && E < E1[l] ) {
          // if (k < st_max && l < n_bins && n < mu_max) 
          if (k < st_max && l < n_bins && n < n_mu) {
            F[k][n][l]+=ew;
            F1[k][n][l]+=ew;
            F2[k][n][l]+=(ew*ew);
            particle_sum[k][n][l]++;  
	    // GF: too ad-hoc to check the 7th array...
	    /* 
            if( l == 6 ) {
              //fprintf(fpcheck,"%e %e %e %e %e\n", t_bound, time, E, ew, mu);
              fprintf(fpcheck," %e %e %e %f  %4d %2d %2d  %f %f\n", time, E, ew, mu, k,n,l,mu0[n],mu1[n]);
            }
	    */ 
	    /*  TEMP FIX
            if (t1[k] <= (t_max+dt) ) {
               fprintf(fpcheck," %e %e %e %f %e %e  %4d %2d %2d\n", time, E, ew, mu, r, z, k, n, l);
	    }
	    */ 
          }
        }
      }
   }
   /* end of while on input photon file */
  
   fclose(fpin);
  
   if(infile[3]=='9') {
      if(infile[2]=='9') {
         infile[1]++;
         infile[2]='0';
      } else {
         infile[2]++;
      }
  
      infile[3]='0';
   } else {
      infile[3]++;
   }
  
   n_file++;
}
/* end of loop over input files */
// GF: fclose(fpext);
fclose(fpcheck);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Xuhui */

for (k=0; k<st_max; k++) {
  for (l=0; l<n_bins; l++) {
    for (n=0; n<n_mu; n++) {

        F[k][n][l]/=(dt*(mu1[n]-mu0[n])/2.);

        if( particle_sum[k][n][l] > 0 ) {
           F1[k][n][l]/= (double)particle_sum[k][n][l];                  // GF
           F2[k][n][l]/= (double)particle_sum[k][n][l];                  // GF
           sF[k][n][l] = sqrt( F2[k][n][l] - F1[k][n][l]*F1[k][n][l] );  // GF
        } else {
           F1[k][n][l] = 1e-20;                                          // GF
           F2[k][n][l] = 1e-20;                                          // GF
           sF[k][n][l] = 1e-20;                                          // GF
        }
        if( particle_sum[k][n][l] == 1 ) {
           sF[k][n][l] = 0.0; //GF
        } else {
        /*F[k][n][l]/=(dt*(E1[l] - E0[l])*(mu1[n]-mu0[n])/2.);*/
        /*F[k][n][l]/=(dt*(E1[l] - E0[l])); */
    }
  }
}

for (n=0; n<n_mu; n++) {
   if( !(fpout1=fopen(outfile1[n], "w")) ) file_error(outfile1[n]);
   if( !(fpout2=fopen(outfile2[n], "w")) ) file_error(outfile2[n]);

   // GF : reformat header later 
   for (l=0; l<n_bins; l++) {
      fprintf (fpout1, "#energy range%i: %e %e (keV)\n", l+1, E0[l], E1[l]);
      fprintf (fpout2, "#energy range%i: %e %e (keV)\n", l+1, E0[l], E1[l]);
   }
   fprintf (fpout1, "#time: %e %e %e\n",t_start,t_end,dt);
   fprintf (fpout1, "#angle: %f %f\n",mu0[0],mu1[0]);
   fprintf (fpout1, "#factor: %i\n", factor);
   fprintf (fpout1, "#   time      +----------------------------- luminosity ----------------------------+ +--------------------------- particle_sum ---------------------+ +---------------------------- var(ew)/<ew> -------------------------+\n");
   fprintf (fpout1, "#             |    1         2         3         4         5         6         7    | |       1        2        3        4        5        6        7| |    1         2         3         4         5         6         7  |\n");

   fprintf (fpout2, "#time: %e %e %e\n",t_start,t_end,dt);
   fprintf (fpout2, "#angle: %f %f\n",mu0[0],mu1[0]);
   fprintf (fpout2, "#factor: %i\n", factor);
   // fprintf (fpout2, "#   time     +----------------------------------------- <ew> -----------------------------------------+ +---------------------------------------- <ew^2> ----------------------------------------+ +---------------------------------------- var(ew) ---------------------------------------+\n");
   // fprintf (fpout2, "#            |     1            2            3            4            5            6            7    | |     1            2            3            4            5            6            7    | |     1            2            3            4            5            6            7    |\n");

   for (k=0; k<st_max; k++) {
      fprintf (fpout1, "%e", max(1.e-20, (t0[k]+t1[k])/2.));
      fprintf (fpout2, "%e", max(1.e-20, (t0[k]+t1[k])/2.));

      for (l=0; l<n_bins; l++) fprintf (fpout1, " %e", max(1.e-20, F[k][n][l]));
      for (l=0; l<n_bins; l++) fprintf (fpout1, " %d", particle_sum[k][n][l]);
      for (l=0; l<n_bins; l++) fprintf (fpout1, " %e", sF[k][n][l]/F1[k][n][l]);  // GF
      fprintf (fpout1, "\n");

      for (l=0; l<n_bins; l++) fprintf (fpout2, " %e", F1[k][n][l]);  // GF
      for (l=0; l<n_bins; l++) fprintf (fpout2, " %e", F2[k][n][l]);  // GF
      for (l=0; l<n_bins; l++) fprintf (fpout2, " %e", sF[k][n][l]);  // GF
      fprintf (fpout2, "\n");

      /* 
      fprintf (fpout1, "\n%e", t1[k]);
      for (l=0; l<n_bins; l++) fprintf (fpout1, " %e", max(1.e-20, F[k][n][l]));
      for (l=0; l<n_bins; l++) fprintf (fpout1, " %i", particle_sum[k][n][l]);
      fprintf (fpout1, "\n"); 
      */
      if (t1[k] > t_max) break;
   }
   fclose(fpout1);
   fclose(fpout2);
}
/* END of loop writing files */
}
/* END of main */

/* F is isotropic luminosity (ergs) per second, integreted over the entire energy band. */
