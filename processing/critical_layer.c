/**
   Main file for processing the dump files */


#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"
#include "curvature.h"
#include "vof.h"

int outfiled (int i);

int outfield (int i) {
  
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
  foreach() {
    FirOrder[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
  }
  foreach() {
    SecOrder[] = (FirOrder[0,1] - FirOrder[0,-1])/(2.*Delta);
  }
  sprintf (fieldname, "field%d", i);
  FILE * ffield = fopen (fieldname, "w");
  output_field ({u.x, pos, FirOrder, SecOrder}, ffield, n = 128); // linear = true
  fclose (ffield);
  sprintf (etaname, "eta%d", i);
  FILE * feta = fopen (etaname, "w");
  foreach(){
    if (interfacial (point, f)) 
      fprintf (feta, "%g,%g\n", x, pos[]);
  }
  fclose (feta);
  return 0;
}

int main ()
{
  for (int i = 0; i<16; i++) {
    char targetname[100];
    double t = 1/32.*i;
    sprintf (targetname, "dump%g", t);
    if (restore(targetname)){
      restore (targetname); //This brings back all the scalar fields f, u, etc
      outfield (i);
    }
  }    
  return 0;
}
