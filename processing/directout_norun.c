/**
   Run from dump with one or two level more refined grid. Direct output for more accurate interpolation.
   in python  */

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

/* Direct output on each grid point */
void outfield (double time) {
  
  char fieldname[100], etaname[100];
  scalar pos[];
#if dimension == 2 
  coord G = {0.,1.}, Z = {0.,0.};
#else
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
#endif
  position (f, pos, G, Z);    
  sprintf (fieldname, "field%g", time);
  FILE * ffield = fopen (fieldname, "w");
  fprintf(ffield, "x,y,u.x,u.y,p\n");
  // direct output of 5 fields
  foreach()
  {
    fprintf (ffield, "%g,%g,%g,%g,%g\n", x, y, u.x[], u.y[], p[]);
  }
  fclose (ffield);

  // direct output of the interface (x, eta)
  // also direct output of pressure, tau in cartesian coord, normal vector at the interface 
  sprintf (etaname, "eta%g", time);
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,eta,p,tau.x,tau.y,n.x,n.y\n");
  vector tau[];
  tensor SDeform = new tensor;

  foreach(){
    if (interfacial (point, f)){
      // Getting the local normal vector
      coord n = mycs (point, f);
      double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
      double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
      // double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
      double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
      // double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      // double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
      // double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
      // double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      SDeform.x.x[] = dudx;
      SDeform.x.y[] = 0.5*(dudy + dvdx);
      // double SDeformxz = 0.5*(dudz + dwdx);
      SDeform.y.x[] = SDeform.x.y[];
      SDeform.y.y[] = dvdy;
      double mu_eff = mu1/rho[]*f[] + mu2/rho[]*(1. - f[]); // compute effective viscosity
      tau.x[] = 2*mu_eff*(SDeform.x.x[]*n.x + SDeform.y.x[]*n.y);
      tau.y[] = 2*mu_eff*(SDeform.x.y[]*n.x + SDeform.y.y[]*n.y);
      // double SDeformyz = 0.5*(dvdz + dwdy);
      // double SDeformzx = SDeformxz;
      // double SDeformzy = SDeformyz;
      // double SDeformzz = dwdz; 
      fprintf (feta, "%g,%g,%g,%g,%g,%g,%g\n", x, pos[], p[], tau.x[], tau.y[], n.x, n.y);
    }
  }
  fclose (feta);
  // return 0;
}

/* Some info about the snapshot time and the refining criteria */
int NUMBER = 32; // Time gap in taking the snapshots
double TIME = 2; // Total number of run time
double snapshot_time; // Current snapshot time

/* Even no run, still need viscosity info for stress */
double RE = 40000.;
double RATIO = (1./850.);
double MURATIO = (17.4e-6/8.9e-4);


int main (int argc, char * argv[])
{
  if (argc > 1)
    NUMBER = atoi (argv[1]);
  if (argc > 2)
    TIME = atoi(argv[2]);
  if (argc < 3){
    printf("Missing parameter!\n");
    return 1;
  }
  
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  char targetname[100];

  for (int i = 0; i < (TIME*NUMBER); i++)
  {
    snapshot_time = i * 1 / (double)NUMBER;
    sprintf (targetname, "dump%g", snapshot_time);
    if (!restore (targetname)) 
      printf("Not restored!\n");
    printf("Processing %d!\n", (i+1)); 
    outfield(snapshot_time);
  }

  return 0;
}


