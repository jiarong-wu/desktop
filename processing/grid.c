/**
   Change the grid to cartesian so that we get smoother velocity profile output. */

// #include "grid/cartesian.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "curvature.h"
#include "vof.h"


void outfield (double time) {
  
  char fieldname[100], etaname[100];
  scalar pos[];
#if dimension == 2 
  coord G = {0.,1.}, Z = {0.,0.};
#else
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#endif
  position (f, pos, G, Z);    
  scalar FirOrder[];
  scalar SecOrder[];
  foreach()
    FirOrder[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
  foreach() 
    SecOrder[] = (FirOrder[0,1] - FirOrder[0,-1])/(2.*Delta);
  sprintf (fieldname, "field%g", time);
  FILE * ffield = fopen (fieldname, "w");
  output_field ({u.x, pos, FirOrder, SecOrder}, ffield, n = 2048); // linear = true
  fclose (ffield);
  sprintf (etaname, "eta%g", time);
  FILE * feta = fopen (etaname, "w");
  foreach(){
    if (interfacial (point, f)) 
      fprintf (feta, "%g,%g\n", x, pos[]);
  }
  fclose (feta);
  // return 0;
}

int NUMBER = 32; // Time gap in taking the snapshots
double TIME = 2; // Total number of run time
double snapshot_time; // Current snapshot time
double uemax = 0.0000001;
double femax = 0.0000001;
int LEVEL = 12;
int counting = 0;

/* All the constant that comes with the case, will be changed after input */
double ak = 0.05;
double BO = 200.;
double RE = 40000.;
double PRESSURE = 0;
double m = 5.;  // vary between 5 and 8
double B = 0.;
double Karman = 0.41;   // Karman universal turbulence constant
double UstarRATIO = 1; // Ratio between Ustar and c

double k_ = 2.*pi;
double h_ = 0.5;
double g_ = 1.;

double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4);

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);   
  if (argc > 5)
    m = atof(argv[5]);
  if (argc > 6)
    B = atof(argv[6]);
  if (argc > 7)
    UstarRATIO = atof(argv[7]);
  if (argc > 8)
    PRESSURE = atof(argv[8]);
  
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
  p[left] = dirichlet(PRESSURE);
  p[right] = dirichlet(0); 
  u.n[top] = dirichlet(0);
  u.t[top] = neumann(0);

  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

  for (int i = 0; i < (TIME*NUMBER); i++)
  {
    snapshot_time = i * 1 / (double)NUMBER;
    printf("Run starts!\n"); 
    run();
    printf("Run ends!\n"); 
  }

  return 0;
}

event init (i = 0)
{
  char targetname[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) 
  	printf("Not restored!\n");
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
  printf("Adapt for the %d time\n", i)
}

event end (i = 10) {
	printf("end\n");
}

event finalize (t = end) {
	outfield(snapshot_time);
}
