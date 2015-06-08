#ifndef DEFS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define DEFS_H

# include "structure_definitions.h"
const double m =1.0;
const double inv_mass =1.0/m;
const double dt=0.0001;
const int NrParticles=2;
const double r_cut  = pow(2.0,1.0/6.0);
const double r_cut2 = (r_cut)*(r_cut);
const double sigma =1, epsilon =1; 
const double r_min = pow(2.0,1.0/6.0)*sigma;
const double r_min2= r_min*r_min;
const double rs = 3.0*r_min/4.0; // saturation radius, below this potential is assumed linear and force remains constant, to prevent calculation of huge forces at extremely close contacts 
const double rs2=rs*rs;
const double eta=1;
const double mu = 1.0/(6.0*pi*eta*(sigma/2.0)); // mu - mobility, eta - viscosity, r-radius of particle suspensions
const double mu_sqrt=sqrt(mu);
const double kb=1.0 ;
const double T0=1.0 ; 
const double kbT_dt=sqrt(2.0*kb*T0*dt);

const double Lx=4.5, Ly=4.5, Lz=4.5;// , R_cut=2.5/*= 1.1225 */,R_shell = 0; // = 0.3775;
const double Volume =Lx*Ly*Lz;
const double Volume_inv = 1/ Volume;
const int cellx=(int) ceil(Lx/r_cut);
const int celly=(int) ceil(Ly/r_cut);
const int cellz=(int) ceil(Lz/r_cut);

const vctr3D box(Lx, Ly, Lz);
const vctr3D rbox(1.0/Lx, 1.0/Ly, 1.0/Lz);
const vctr3D havbox (Lx/2.0, Ly/2.0 , Lz/2.0);
const vctr3D shift(r_cut+r_min,r_cut+r_min,r_cut+r_min);

const double Fs=4*epsilon*(12.0*pow(sigma/rs,12.0)-6.0*pow(sigma/rs,6.0))/rs;
const double phis =4.0*epsilon*(pow(sigma/rs,12.0)-pow(sigma/rs,6.0));
const double phicut =4.0*epsilon*(pow(sigma/r_cut,12)-pow(sigma/r_cut,6.0));
const double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
const double sigma12 = sigma6*sigma6;

const mtrx3D Unit_diag(1,1,1);


#endif

