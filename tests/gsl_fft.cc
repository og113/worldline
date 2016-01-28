/*
	quick program to test gsl_fft_real and gsl_fft_halfcomplex (the inverse).
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

using namespace std;

int main() {

int i, n=100, cutoff=10;
double data[n];

string xFile = "data/temp/gsl_fft_x.dat";
string xftFile = "data/temp/gsl_fft_xft.dat";

gsl_fft_real_wavetable * real;
gsl_fft_halfcomplex_wavetable * hc;
gsl_fft_real_workspace *work;

for (i=0; i<n; i++)
	data[i] = 0.0;
	
for (i=n/3; i<2*n/3; i++)
	data[i] = 1.0;
	
ofstream os(xFile.c_str());
for (i=0; i<n; i++)
	os << setw(24) << i << setw(24) << data[i] << endl;
os.close();

work = gsl_fft_real_workspace_alloc(n);
real = gsl_fft_real_wavetable_alloc(n);

gsl_fft_real_transform(data,1,n,real,work);

gsl_fft_real_wavetable_free(real);

for (i=cutoff+1; i<n; i++)
	data[i] = 0;
	
hc = gsl_fft_halfcomplex_wavetable_alloc(n);

gsl_fft_halfcomplex_inverse(data,1,n,hc,work);
gsl_fft_halfcomplex_wavetable_free(hc);

os.open(xftFile.c_str());
for (i=0; i<n; i++)
	os << setw(24) << i << setw(24) << data[i] << endl;
os.close();
	
gsl_fft_real_workspace_free(work);

return 0;
}
