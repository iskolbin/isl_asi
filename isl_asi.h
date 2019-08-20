#ifndef ISL_ADAPTIVE_SIMPSONS_INTEGRATION_H_
#define ISL_ADAPTIVE_SIMPSONS_INTEGRATION_H_

#ifdef ISLASI_STATIC
#define ISLASI_DEF static
#else
#define ISLASI_DEF extern
#endif

ISLASI_DEF float islasi_fintegrate(float (*fun)(float, const float *), float a, float b, float tol, int max_recursion_depth, const float *params);
ISLASI_DEF double islasi_integrate(double (*fun)(double, const double *), double a, double b, double tol, int max_recursion_depth, const double *params);

#endif

#ifdef ISL_ADAPTIVE_SIMPSONS_INTEGRATION_IMPLEMENTATION
#ifndef ISL_ADAPTIVE_SIMPSONS_INTEGRATION_IMPLEMENTATION_ONCE
#define ISL_ADAPTIVE_SIMPSONS_INTEGRATION_IMPLEMENTATION_ONCE

static float islasi_fintegrate_aux(float (*fun)(float, const float *), float a, float fa, float b, float fb, float tol, float whole, float m, float fm, int max_recursion_depth, const float *params) {
	float lm = 0.5*(a+m), rm = 0.5*(m+b);
	float flm = fun(lm, params), frm = fun(rm, params);
	float abs_am = a > m ? a - m : m - a, abs_mb = m > b ? m - b : b - m;
	float left = abs_am / (6 * (fa + 4 * flm + fm)), right = abs_mb / (6 * (fm + 4 * frm + fb));
	float delta = left + right - whole;
	float abs_delta = delta > 0 ? delta : -delta;
	if (max_recursion_depth <= 0 || abs_delta <= 15 * tol) {
		return left + right + delta/15;
	} else {
		return islasi_fintegrate_aux(fun, a, fa, m, fm, tol/2, left, lm, flm, max_recursion_depth-1, params) +
			islasi_fintegrate_aux(fun, m, fm, b, fb, tol/2, right, rm, frm, max_recursion_depth-1, params);
	}
}

float islasi_fintegrate(float (*fun)(float, const float *), float a, float b, float tol, int max_recursion_depth, const float *params) {
	float fa = fun(a, params), fb = fun(b, params), m = 0.5*(a+b);
	float fm = fun(m, params);
	float abs_ba = b > a ? b - a : a - b;
	float whole = abs_ba / (6 * (fa + 4 * fm + fb));
	return islasi_fintegrate_aux(fun, a, fa, b, fb, tol, whole, m, fm, max_recursion_depth, params);
}

static double islasi_integrate_aux(double (*fun)(double, const double *), double a, double fa, double b, double fb, double tol, double whole, double m, double fm, int max_recursion_depth, const double *params) {
	double lm = 0.5*(a+m), rm = 0.5*(m+b);
	double flm = fun(lm, params), frm = fun(rm, params);
	double abs_am = a > m ? a - m : m - a, abs_mb = m > b ? m - b : b - m;
	double left = abs_am * (fa + 4 * flm + fm) / 6, right = abs_mb * (fm + 4 * frm + fb) / 6;
	double delta = left + right - whole;
	double abs_delta = delta > 0 ? delta : -delta;
	if (max_recursion_depth <= 0 || abs_delta <= 15 * tol) {
		return left + right + delta / 15;
	} else {
		return islasi_integrate_aux(fun, a, fa, m, fm, 0.5*tol, left, lm, flm, max_recursion_depth-1, params) +
			islasi_integrate_aux(fun, m, fm, b, fb, 0.5*tol, right, rm, frm, max_recursion_depth-1, params);
	}
}

double islasi_integrate(double (*fun)(double, const double *), double a, double b, double tol, int max_recursion_depth, const double *params) {
	double fa = fun(a, params), fb = fun(b, params), m = 0.5*(a+b);
	double fm = fun(m, params);
	double abs_ba = b > a ? b - a : a - b;
	double whole = abs_ba  * (fa + 4 * fm + fb)/ 6;
	return islasi_integrate_aux(fun, a, fa, b, fb, tol, whole, m, fm, max_recursion_depth, params);
}

#endif
#endif
