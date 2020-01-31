/**
 * \file statistics.cpp
 *
 * Contains implementations of statistical functions
 *
 */

#include "statistics.hpp"

namespace vg {


double median(std::vector<int> &v) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    int vn = v[n];
    if (v.size()%2 == 1) {
        return vn;
    } else {
        std::nth_element(v.begin(), v.begin()+n-1, v.end());
        return 0.5*(vn+v[n-1]);
    }
}

// from Python exmaple here:
// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
void wellford_update(size_t& count, double& mean, double& M2, double new_val) {
    ++count;
    double delta = new_val - mean;
    mean += delta / (double)count;
    double delta2 = new_val - mean;
    M2 += delta * delta2;
}

pair<double, double> wellford_mean_var(size_t count, double mean, double M2, bool sample_variance) {
    if (count == 0 || (sample_variance && count == 1)) {
        return make_pair(nan(""), nan(""));
    } else {
        return make_pair(mean, M2 / (double)(sample_variance ? count - 1 : count));
    }
}

double phi(double x1, double x2) {
    return (std::erf(x2/std::sqrt(2)) - std::erf(x1/std::sqrt(2)))/2;
}

// Modified from qnorm function in R source:
// https://svn.r-project.org/R/trunk/src/nmath/qnorm.c
double normal_inverse_cdf(double p) {
    assert(0.0 < p && p < 1.0);
    double q, r, val;
    
    q = p - 0.5;
    
    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     
     Produces the normal deviate Z corresponding to a given lower
     tail area of P; Z is accurate to about 1 part in 10**16.
     
     (original fortran code used PARAMETER(..) for the coefficients
     and provided hash codes for checking them...)
     */
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
        q * (((((((r * 2509.0809287301226727 +
                   33430.575583588128105) * r + 67265.770927008700853) * r +
                 45921.953931549871457) * r + 13731.693765509461125) * r +
               1971.5909503065514427) * r + 133.14166789178437745) * r +
             3.387132872796366608)
        / (((((((r * 5226.495278852854561 +
                 28729.085735721942674) * r + 39307.89580009271061) * r +
               21213.794301586595867) * r + 5394.1960214247511077) * r +
             687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */
        
        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = 1.0 - p;
        else
            r = p;
        
        r = sqrt(- log(r));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
        
        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177) *
                       r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734)
            / (((((((r *
                     1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                   .14810397642748007459) * r + .68976733498510000455) *
                 r + 1.6763848301838038494) * r +
                2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                      r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
            / (((((((r *
                     2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                    r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + .0148753612908506148525)
                 * r + .13692988092273580531) * r +
                .59983220655588793769) * r + 1.);
        }
        
        if(q < 0.0)
            val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return val;
}

// https://stackoverflow.com/a/19039500/238609
double slope(const std::vector<double>& x, const std::vector<double>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}

//https://stats.stackexchange.com/a/7459/14524
// returns alpha parameter of zipf distribution
double fit_zipf(const vector<double>& y) {
    // assume input is log-scaled
    // fit a log-log model
    assert(y.size());
    vector<double> ly(y.size());
    for (int i = 0; i < ly.size(); ++i) {
        //cerr << y[i] << " ";
        ly[i] = log(y[i]);
    }
    //cerr << endl;
    vector<double> lx(y.size());
    for (int i = 1; i <= lx.size(); ++i) {
        lx[i-1] = log(i);
    }
    return -slope(lx, ly);
}

double fit_max_exponential(const vector<double>& x, double N, double tolerance) {
    
    // Fit S for a fixed N with the density of the maximum of N exponential variables
    //
    //   NS exp(-Sx) (1 - exp(-Sx))^(N - 1)
    //
    // where S is the scale
    
    double x_sum = 0;
    double x_max = numeric_limits<double>::lowest();
    for (const double& val : x) {
        x_sum += val;
    }
    
    // compute the log of the 1st and 2nd derivatives for the log likelihood (split up by positive and negative summands)
    // we have to do it this wonky way because the exponentiated numbers get very large and cause overflow otherwise
    
    double log_deriv_neg_part = log(x_sum);
    
    function<double(double)> log_deriv_pos_part = [&](double scale) {
        double accumulator = numeric_limits<double>::lowest();
        for (const double& val : x) {
            accumulator = add_log(accumulator, log(val) - scale * val - log(1.0 - exp(-scale * val)));
        }
        accumulator += log(N - 1.0);
        return add_log(accumulator, log(num_simulations / scale));
    };
    
    function<double(double)> log_deriv2_neg_part = [&](double scale) {
        double accumulator = numeric_limits<double>::lowest();
        for (const double& val : x) {
            accumulator = add_log(accumulator, 2.0 * log(val) - scale * val - 2.0 * log(1.0 - exp(-scale * length)));
        }
        accumulator += log(N - 1.0);
        return add_log(accumulator, log(num_simulations / (scale * scale)));
    };
    
    // use Newton's method to find the MLE
    double scale = 1.0 / x_max;
    double prev_scale = scale * (1.0 + 10.0 * tolerance);
    while (abs(prev_scale / scale - 1.0) > tolerance) {
        prev_scale = scale;
        double log_d2 = log_deriv2_neg_part(scale);
        double log_d_pos = log_deriv_pos_part(scale);
        double log_d_neg = log_deriv_neg_part;
        // determine if the value of the 1st deriv is positive or negative, and compute the
        // whole ratio to the 2nd deriv from the positive and negative parts accordingly
        if (log_d_pos > log_d_neg) {
            scale += exp(subtract_log(log_d_pos, log_d_neg) - log_d2);
        }
        else {
            scale -= exp(subtract_log(log_d_neg, log_d_pos) - log_d2);
        }
    }
    return scale;
}

}
