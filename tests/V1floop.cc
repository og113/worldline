/*
	V1floop
		test of speed of integrating V1 for floops, using gsl monte carlo integration,
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <ctime>
//#include <cmath>
#include <gsl/gsl_randist.h> 	// Distributions of random numbers
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "simple.h"

using namespace std;

/* 
	Computation of the integral,

      I = int dx.dy /( (x-y)^2 + a^2 )

   over floops, with x and y approximated by sums of sins and cosines.
*/

struct GSLMonteParams {
	size_t			N;
	double 			a;
	vector<double>  fcoeffs;
};

double g (double *k, size_t dim, void *params) {
	(void)(dim); /* avoid unused parameter warnings */
	struct GSLMonteParams *p = (struct GSLMonteParams *) params; 
	double denom = (p->a)*(p->a);
	double x0 = 0.0, x1 = 0.0, dx0 = 0.0, dx1 = 0.0;
	vector<double>& fcr = p->fcoeffs;
	for (size_t i=0; i<(p->N); i++) {
		x0 += fcr[2*i]*cos(2.0*M_PI*(i+1.0)*k[0]) + fcr[2*i+1]*sin(2.0*M_PI*(i+1.0)*k[0]);
		x1 += fcr[2*i]*cos(2.0*M_PI*(i+1.0)*k[1]) + fcr[2*i+1]*sin(2.0*M_PI*(i+1.0)*k[1]);
		dx0 += 2.0*pi*i*(-fcr[2*i]*sin(2.0*pi*(i+1.0)*k[0]) + fcr[2*i+1]*cos(2.0*pi*(i+1.0)*k[0]));
		dx1 += 2.0*pi*i*(-fcr[2*i]*sin(2.0*pi*(i+1.0)*k[1]) + fcr[2*i+1]*cos(2.0*pi*(i+1.0)*k[1]));
	}
	denom += pow(x0-x1,2);
	return dx0*dx1/denom;
}

void display_results (char *title, double result, double error) {
	printf ("%s ==================\n", title);
	printf ("result = % .6f\n", result);
	printf ("sigma  = % .6f\n", error);
}

int main () {
  double res, err;

  double tl[2] = { 0.0, 0.0 };
  double tu[2] = { 1.0, 1.0 };

  const gsl_rng_type *T;
  gsl_rng *r;
  
  size_t N = 32;
  double a = 0.1;
  vector<double> fcoeffs(2*N);
  
  uint Seed = time(NULL);
  gsl_rng* Generator = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(Generator,Seed);
  number sigma = SQRT2/pi;
  
  for (size_t i=0; i<N; i++) {
  	fcoeffs[2*i] = gsl_ran_gaussian (Generator, sigma); //gsl_ran_gaussian_ziggurat is another option
  	fcoeffs[2*i+1] = gsl_ran_gaussian (Generator, sigma); //gsl_ran_gaussian_ziggurat is another option
  }
  
  delete Generator;

  GSLMonteParams gslParams;
  gslParams.N = N;
  gslParams.a = a;
  gslParams.fcoeffs = fcoeffs;

  gsl_monte_function G = { &g, 2, &gslParams };

  size_t calls = 1e6;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  {
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&G, tl, tu, 2, calls, r, s, 
                               &res, &err);
    gsl_monte_plain_free (s);

    display_results ("plain", res, err);
  }

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
    gsl_monte_miser_integrate (&G, tl, tu, 2, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    display_results ("miser", res, err);
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    gsl_monte_vegas_integrate (&G, tl, tu, 2, 10000, r, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);

    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, tl, tu, 2, calls/5, r, s,
                                   &res, &err);
        printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }

  gsl_rng_free (r);

  return 0;
}
