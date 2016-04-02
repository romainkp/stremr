#include <Rcpp.h>
using namespace Rcpp;

// replacing .Call(stats:::"C_logit_linkinv", eta) with a C++ function
// to compile: sourceCpp("/Users/olegsofrygin/GoogleDrive/Network_TMLE/TMLENET_package/tmlenet/glm_binom_family.cpp")

#define DOUBLE_EPS     DBL_EPSILON

static const double THRESH = 30.;
static const double MTHRESH = -30.;
static const double INVEPS = 1/DOUBLE_EPS;

static double x_d_opx(double x) {return x/(1 + x);}

// [[Rcpp::export]]
NumericVector logit_linkinv(NumericVector eta) {
  int n = eta.size();
  NumericVector ans(n);
  for(int i = 0; i < n; ++i) {
    double etai = eta[i], tmp;
    tmp = (etai < MTHRESH) ? DOUBLE_EPS :
      ((etai > THRESH) ? INVEPS : exp(etai));
    ans[i] = x_d_opx(tmp);
  }
  return ans;
}

/*** R
library(microbenchmark)
x <- rnorm(n=1e5, mean=0, sd=1000)
microbenchmark(
  logitlinkinv(x),
  logit_linkinv_new(x)
)
*/
// Unit: microseconds
//                  expr     min       lq     mean  median       uq      max neval cld
//       logitlinkinv(x) 580.418 590.3630 770.1703 605.873 656.4225 1546.020   100   a
//  logit_linkinv_new(x) 550.359 555.3005 726.0105 563.787 630.5980 1559.106   100   a


// -----------------------------------------------------------------
// OLD C FUNCTIONS from r-source/src/library/stats/src/family.c
// -----------------------------------------------------------------
/**
 * Evaluate x/(1 - x). An inline function is used so that x is
 * evaluated once only.
 * @param x input in the range (0, 1)
 * @return x/(1 - x)
 */

/**
static R_INLINE double x_d_omx(double x) {
    if (x < 0 || x > 1)
  error(_("Value %g out of range (0, 1)"), x);
    return x/(1 - x);
}
SEXP logit_link(SEXP mu)
{
    int i, n = LENGTH(mu);
    SEXP ans = PROTECT(shallow_duplicate(mu));
    double *rans = REAL(ans), *rmu=REAL(mu);
    if (!n || !isReal(mu))
  error(_("Argument %s must be a nonempty numeric vector"), "mu");
    for (i = 0; i < n; i++)
  rans[i] = log(x_d_omx(rmu[i]));
    UNPROTECT(1);
    return ans;
}
*/