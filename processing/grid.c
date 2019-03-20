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


void outfield (int i) {
  
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
  // return 0;
}

double snapshot_time;
int INDEX;
double uemax = 0.00001;
double femax = 0.00001;
int LEVEL = 12;

int main (int argc, char * argv[])
{
	double PERIOD = 1/32.;
	if (argc > 1)
		snapshot_time = atof (argv[1]);
	INDEX = snapshot_time/PERIOD; 
	printf("INDEX=%d\n", INDEX);

		// restore (targetname); //This brings back all the scalar fields f, u, etc
	// outfield(INDEX);
	run();
    return 0;
}

int j = 0;

event init (t = 0) {
	char targetname[100];
	sprintf (targetname, "dump%g", snapshot_time);
	if (!restore(targetname)){
		printf("No dump file!\n");
	}
	fprintf(stdout, "%d\n", i);
}


event count (i++, t<20000) {
	j++;
	outfield(INDEX);
	if (j==5)
		return 1;	

}

// event something (t = end) {
// 	fprintf(stdout, "%g %d\n", t, i);

// 		outfield(INDEX);
// 		// fprintf(stdout, "%g %d\n", t, i);
// 		// return 0;
// }



// event end (t = 20000) {
// 	printf("toto\n");
// }

// event adapt (i++) {
//   adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
// }
