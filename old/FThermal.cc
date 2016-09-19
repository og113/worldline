
// sinch
static number sinch(const number& x) {
	return (abs(x)<MIN_NUMBER? 1.0: sinh(x)/x);
}

// coshsinhc
/*static number coshsinhc(const number& x) {
	return (abs(x)<MIN_NUMBER? 0.0: cosh(x)/x - sinh(x)/pow(x,2));
}*/

// coshsinhcc
static number coshsinhcc(const number& x) {
	return (abs(x)<MIN_NUMBER? 1.0/3.0: cosh(x)/pow(x,2) - sinh(x)/pow(x,3));
}

// FThermal
static number FThermal(const number& r, const number& t, const number& beta, const number& a) {
	return sinch((2.0*r*PI)/beta)/(2.0*pow(beta,2)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*r*PI)/beta)));
}

// DFThermalDr - use Onr version instead as finite at r=0
/*static number DFThermalDr(const number& r, const number& t, const number& beta, const number& a) {
	return (PI*coshsinhc((2.0*PI*r)/beta))/\
 (pow(beta,3)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta))) \
 + (2.0*pow(PI,2)*r*pow(sinch((2.0*PI*r)/beta),2))/\
 (pow(beta,4)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
}*/

// DFThermalDrOnr
static number DFThermalDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (2.0*pow(PI,2)*coshsinhcc((2.0*PI*r)/beta))/\
 (pow(beta,4)*((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta))) \
 + (2.0*pow(PI,2)*pow(sinch((2.0*PI*r)/beta),2))/\
 (pow(beta,4)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
 }

// DFThermalDt
static number DFThermalDt(const number& r, const number& t, const number& beta, const number& a) {
	return (PI*sin((2.0*PI*t)/beta)*sinch((2.0*PI*r)/beta))/\
 (pow(beta,3)*pow((-2.0*pow(a,2)*pow(PI,2))/pow(beta,2) + cos((2.0*PI*t)/beta) - cosh((2.0*PI*r)/beta),2));
}

// DDFThermalDrDr
static number DDFThermalDrDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*(-((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 (3.0*cosh((2.0*PI*r)/beta) - 2.0*sinch((2.0*PI*r)/beta))*sinch((2.0*PI*r)/beta)) \
 + 8.0*pow(PI,2)*pow(r,2)*pow(sinch((2.0*PI*r)/beta),3) \
 + (pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),2)*
 (-2.0*coshsinhcc((2.0*PI*r)/beta) + sinch((2.0*PI*r)/beta)))/pow(beta,2)))/\
 pow(-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}

// DDFThermalDtDr - use Onr version instead as finite at r=0
/*static number DDFThermalDtDr(const number& r, const number& t, const number& beta, const number& a) {
	return  (2.0*pow(PI,2)*sin((2.0*PI*t)/beta)*((2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta))*
 coshsinhc((2.0*PI*r)/beta) - 4.0*PI*r*beta*pow(sinch((2.0*PI*r)/beta),2)))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}*/

 // DDFThermalDtDrOnr
static number DDFThermalDtDrOnr(const number& r, const number& t, const number& beta, const number& a) {
	return (4.0*pow(PI,3)*sin((2.0*PI*t)/beta)*((-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta))*\
 coshsinhcc((2.0*PI*r)/beta) + 2.0*pow(beta,2)*pow(sinch((2.0*PI*r)/beta),2)))/\
 (beta*pow(-2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cos((2.0*PI*t)/beta) - pow(beta,2)*cosh((2.0*PI*r)/beta),3));
}

// DDFThermalDtDt
static number DDFThermalDtDt(const number& r, const number& t, const number& beta, const number& a) {
	return  (pow(PI,2)*(-3.0*pow(beta,2) + pow(beta,2)*cos((4.0*PI*t)/beta) \
 + 2.0*cos((2.0*PI*t)/beta)*(2.0*pow(a,2)*pow(PI,2) + pow(beta,2)*cosh((2.0*PI*r)/beta)))*sinch((2.0*PI*r)/beta))/\
 pow(2.0*pow(a,2)*pow(PI,2) - pow(beta,2)*cos((2.0*PI*t)/beta) + pow(beta,2)*cosh((2.0*PI*r)/beta),3);
}
