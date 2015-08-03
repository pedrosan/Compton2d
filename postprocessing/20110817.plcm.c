#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define ch_max 20
#define st_max 1024
#define mu_max 10

#define qu(x) ((x)*(x))



double max(x, y)
double x, y;
{
if (x > y) return (x);
else return (y);
}



int digit(c)
char c;
{
if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4'
 || c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
    return (1);
else return (0);
}



void add_dat(str)
char *str;
{
int len;

len = strlen(str);
if (str[len-4] != '.')
  { str[len] = '.';
    str[len+1] = 'd';
    str[len+2] = 'a';
    str[len+3] = 't';
    str[len+4] = '\0';
    }
return;
}


double dinput(x)
double x;
{
char a[30];

gets(a);
if (strlen(a)) return (atof(a));
else return (x);
}


int iinput(x)
int x;
{
char a[30];

gets(a);
if (strlen(a)) return (atoi(a));
else return (x);
}



void file_error(fname)
char *fname;
{
printf ("\n\n File error: Can not open file '%s'! - Exit.\n\n", fname);
exit(0);
}



main()
{
double E0[ch_max], E1[ch_max], t0[st_max], t1[st_max];
double time, E, ew, dt=7.e2, F[st_max][mu_max][ch_max];
double E_lower=1.e-3, E_upper=3.e-3, dE, t_offset=0.;
double Pi = 3.1415926536;
double mu0[mu_max], mu1[mu_max], tmax=7.e4;
double t_bound, r,rmax=1.e16, z, mu, phi, cdt;
double gam_bulk=33.,doppler,betta;
int k, l, m, n, factor=0, regions=7, reg, n_bins=0, n_r=1, n_mu=1;
int particle_sum[st_max][mu_max][ch_max], n_file, max_file=200;
char infile[30]="p001_evb.dat", outfile[mu_max][50], a[30];
FILE *fpin, *fpout, *fpext;

printf ("\n\n Input from event file: [%s]  ", infile);
gets(a);
if (strlen(a))
  { if (digit(a[0]))
      { for (k=0; k<strlen(a); ++k) infile[2+k] = a[k];
        infile[2+k] = '\0';
        }
    else strcpy(infile, a);
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

strcpy(outfile[0], infile);
outfile[0][0] = 'l';
outfile[0][1] = 'c';
outfile[0][3] = '7';
/*outfile[0][l-4] = '_';
outfile[0][l-3] = 'e';
outfile[0][l-2] = 'v';
outfile[0][l-1] = 't';
outfile[0][l] = '_';
outfile[0][l+1] = '0';
outfile[0][l+2] = '\0'; */
outfile[0][l-5] = '0';
add_dat(outfile[0]);

printf ("\n Output file for mu-bin 0: [%s]  ", outfile[0]);
gets(a);
if (strlen(a)) strcpy(outfile[0], a);
add_dat(outfile[0]);
l = strlen(outfile[0]);

printf ("\n Number of angular bins: [%1d]  ", n_mu);
n_mu = iinput(n_mu);

if (n_mu > mu_max) 
  { printf ("\n Number of angular bins set to maximum number = %d.\n", mu_max);
    n_mu = mu_max;
    }

mu0[0] = .99944;
mu1[0] = .99964;
for (k=0; k<n_mu; ++k)
  { if (k)
      { strcpy(outfile[k], outfile[0]);
        outfile[k][l-5] = (char)(48 + k);
        }
    printf ("\n Angular bin no. %1d:", k);
    printf ("\n\n Output file: [%s]  ", outfile[k]);
    gets(a);
    if (strlen(a)) strcpy(outfile[k], a);
    add_dat(outfile[k]);
    printf (" mu_min(%1d) = [%e]  ", k, mu0[k]);
    mu0[k] = dinput(mu0[k]);
    printf (" mu_max(%1d) = [%e]  ", k, mu1[k]);
    mu1[k] = dinput(mu1[k]);
    if (k < (mu_max - 1)) mu0[k+1] = mu1[k+1] = mu1[k];
    }

printf ("\n Time step: dt = [%e] s  ", dt);
dt = dinput(dt);
printf (" Time offset: [%e] s  ", t_offset);
t_offset = dinput(t_offset);
printf (" Stop at t = [%e] s  ", tmax);
tmax = dinput(tmax);

tmax-=t_offset;

t0[0] = 0.;
for (k=0; k<st_max-1; ++k) 
   t1[k] = t0[k+1] = t0[k] + dt;
t1[st_max-1] = t0[st_max-1] + dt;

m1:
printf ("\n Energy grid: \n");
printf ("\n Number of energy regions: [%d]  ", regions);
regions = iinput(regions);
for (reg=0; reg<regions; ++reg)
  { if(reg==0)
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
    printf ("   Region no. %d:", reg+1);
    printf ("\n   Lower boundary: [%e] keV  ", E_lower);
    E_lower = dinput(E_lower);
    printf ("   Upper boundary: [%e] keV  ", E_upper);
    E_upper = dinput(E_upper);
    printf ("   Number of energy bins in this region: [%d]  ", n_r);
    n_r = iinput(n_r);
    if ((n_r + n_bins) > ch_max)
      { printf ("\n Not more than 20 energy channels!");
        printf ("\n Re-specify energy grid! \n");
        goto m1;
        }
    printf ("   Logarithmic (0) or linear (1) bins? [0]  ");
    gets(a);
    if (a[0] == '1')
      { dE = (E_upper - E_lower)/((double)(n_r));
        E0[n_bins] = E_lower;
        for (k=n_bins; k<(n_bins + n_r - 1); ++k)
           E1[k] = E0[k+1] = E0[k] + dE;
        E1[n_bins + n_r - 1] = E0[n_bins + n_r - 1] + dE;
        }
    else
      { dE = exp(log(E_upper/E_lower)/((double)(n_r)));
        E0[n_bins] = E_lower;
        for (k=n_bins; k<(n_bins + n_r - 1); ++k)
           E1[k] = E0[k+1] = E0[k]*dE;
        E1[n_bins + n_r - 1] = E0[n_bins + n_r - 1]*dE;
        }
    n_bins+=n_r;
    }
   

printf ("\n Energy binning:\n");
for (k=0; k<n_bins; ++k)
   printf ("\n %d. E0 = %e keV, E1 = %e keV", k, E0[k], E1[k]);
printf ("\n");

printf ("\n Angular binning:\n");
for (n=0; n<n_mu; ++n)
   printf ("\n %d. mu_min = %e, mu_max = %e", n, mu0[n], mu1[n]);
printf ("\n");
//printf ("\n number of input files: [%1d]", max_file);
//max_file = iinput(max_file);

for (k=0; k<st_max; ++k)
  for (l=0; l<n_bins; ++l)
    for (n=0; n<n_mu; ++n)
      {F[k][n][l] = 0.;
       particle_sum[k][n][l] =0;}

m = 0;

      if (!(fpext=fopen("TeV.dat", "w"))) file_error("TeV.dat");

n_file = 1;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// m2:
while(n_file <= max_file)
{
  /*if(infile[3]=='9'){
  infile[2]++;
  infile[3]='0';}
  else
  infile[3]++;
  n_file++;*/
// if (n_file <= max_file)
  if (!(fpin=fopen(infile, "r"))){
    do{
     if(infile[1]!='0' || infile[2]!='0' || infile[3]!='1')factor++;
     infile[7]++;
     infile[1]='0';
     infile[2]='0';
     infile[3]='1';}
    while((!(fpin=fopen(infile,"r"))) && infile[7]<'l');
    if(infile[7]>='l') break;} //file_error(infile);
  printf ("\n Input from another event file? [%s] or n \n", infile);

  while ((a[0]=fgetc(fpin)) != EOF)
    { fseek(fpin, -1, 1);
      fscanf(fpin, "%lf %lf %lf %lf %lf %lf %lf\n", &t_bound, &E, &ew, &r, &z, &mu, &phi);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*    if (z > zmax*(1.-1.e-7)) continue;          */
    mu = -mu;
    betta = sqrt(1.-1./(gam_bulk*gam_bulk));
    doppler = gam_bulk*(1.+mu*betta);
    t_bound = (t_bound-betta*z*3.33333333e-11)/doppler;
    E = E*doppler;
    ew = ew*doppler;
    mu = (mu+betta)/(1.+mu*betta);
    cdt = z*mu/gam_bulk + sqrt(1.-mu*mu)*(rmax-r*cos(phi));
    time = t_bound + 3.33333333e-11*cdt;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

      if (m++ < 10) 
        { printf (" t = %e, E = %e, ew = %e\n", t_bound, E, ew);
          printf (" r = %e, z = %e, mu = %e, phi = %e\n", r, z, mu, phi);
          printf (" t_corr = %e\n", time);
          }
    
/*      cdt = (rmax - r*cos(phi))*sqrt(1. - qu(mu)) + (zmax - z)*mu; */

      time-=t_offset;
      if (time < 0.) continue;

      for (k=0; k<st_max; ++k)
        if (time >= t0[k] && time < t1[k]) break;
      for (l=0; l<n_bins; ++l)
        if (E >= E0[l] && E < E1[l]) break;
      for (n=0; n<n_mu; ++n)
        if (mu >= mu0[n] && mu < mu1[n]) break;
      if(k == 50 && l == 6 && n == 0){
         fprintf(fpext,"%e %e %e %e %e %e %e\n", t_bound, E, ew, r, z, mu, phi);}
      if (k < st_max && l < n_bins && n < mu_max) {
          F[k][n][l]+=ew;
          particle_sum[k][n][l]++;}
      }

  fclose(fpin);
  if(infile[3]=='9'){
   if(infile[2]=='9'){
     infile[1]++;
     infile[2]='0';}
   else
  infile[2]++;
  infile[3]='0';}
  else
  infile[3]++;
  n_file++;
}
fclose(fpext);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Xuhui */

for (k=0; k<st_max; ++k)
  for (l=0; l<n_bins; ++l)
    for (n=0; n<n_mu; ++n)
      F[k][n][l]/=(dt*(mu1[n]-mu0[n])/2.);
      /*F[k][n][l]/=(dt*(E1[l] - E0[l])*(mu1[n]-mu0[n])/2.);*/
      /*F[k][n][l]/=(dt*(E1[l] - E0[l])); */

for (n=0; n<n_mu; ++n)
  { if (!(fpout=fopen(outfile[n], "w"))) file_error(outfile[n]);
    for (l=0; l<n_bins; ++l) fprintf (fpout, "#energy range%i: %e %e (keV)\n", l, E0[l], E1[l]);
    fprintf (fpout, "#angle: %f %f\n",mu0[0],mu1[0]);
    fprintf (fpout, "#factor: %i\n", factor);
    for (k=0; k<st_max; ++k)
      { fprintf (fpout, "%e", max(1.e-20, (t0[k]+t1[k])/2.));
        for (l=0; l<n_bins; ++l) fprintf (fpout, " %e", max(1.e-20, F[k][n][l]));
        for (l=0; l<n_bins; ++l) fprintf (fpout, " %i", particle_sum[k][n][l]);
        fprintf (fpout, "\n");
/*        fprintf (fpout, "\n%e", t1[k]);
        for (l=0; l<n_bins; ++l) fprintf (fpout, " %e", max(1.e-20, F[k][n][l]));
        for (l=0; l<n_bins; ++l) fprintf (fpout, " %i", particle_sum[k][n][l]);
        fprintf (fpout, "\n"); */
        if (t1[k] > tmax) break;
        }
    fclose(fpout);
    }
}

/* F is
isotropic luminosity (ergs) per second, integreted over the entire energy band. */
