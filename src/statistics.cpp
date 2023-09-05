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

double Phi(double x) {
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

// Modified from qnorm function in R source:
// https://svn.r-project.org/R/trunk/src/nmath/qnorm.c
double Phi_inv(double p) {
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

double lognormal_pdf(double x, double mu, double sigma) {
    const static double root_2pi = sqrt(2.0 * 3.14159265358979323846);
    double density;
    if (x > 0.0) {
        double z = (log(x) - mu) / sigma;
        density = exp(-z * z / 2.0) / (sigma * x * root_2pi);
    }
    else {
        density = 0.0;
    }
    return density;
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

double fit_fixed_shape_max_exponential(const vector<double>& x, double shape, double tolerance) {
    
    // Fit S for a fixed N with the density of the maximum of N exponential variables
    //
    //   NS exp(-Sx) (1 - exp(-Sx))^(N - 1)
    //
    // where S is the rate
    // where N is the shape
    
    double x_sum = 0;
    double x_max = numeric_limits<double>::lowest();
    for (const double& val : x) {
        x_sum += val;
        x_max = max(x_max, val);
    }
    
    // compute the log of the 1st and 2nd derivatives for the log likelihood (split up by positive and negative summands)
    // we have to do it this wonky way because the exponentiated numbers get very large and cause overflow otherwise
    
    double log_deriv_neg_part = log(x_sum);
    
    function<double(double)> log_deriv_pos_part = [&](double rate) {
        double accumulator = numeric_limits<double>::lowest();
        for (const double& val : x) {
            if (val > 0.0) {
                // should always be > 0, but just so we don't blow up on some very small graphs
                accumulator = add_log(accumulator, log(val) - rate * val - log(1.0 - exp(-rate * val)));
            }
        }
        accumulator += log(shape - 1.0);
        return add_log(accumulator, log(x.size() / rate));
    };
    
    function<double(double)> log_deriv2_neg_part = [&](double rate) {
        double accumulator = numeric_limits<double>::lowest();
        for (const double& val : x) {
            if (val > 0.0) {
                // should always be > 0, but just so we don't blow up on some very small graphs
                accumulator = add_log(accumulator, 2.0 * log(val) - rate * val - 2.0 * log(1.0 - exp(-rate * val)));
            }
        }
        accumulator += log(shape - 1.0);
        return add_log(accumulator, log(x.size() / (rate * rate)));
    };
    
    // set a maximum so this doesn't get in an infinite loop even when numerical issues
    // prevent convergence
    size_t max_iters = 1000;
    size_t iter = 0;
    
    // use Newton's method to find the MLE
    double rate = 1.0 / x_max;
    double prev_rate = rate * (1.0 + 10.0 * tolerance);
    while (abs(prev_rate / rate - 1.0) > tolerance && iter < max_iters) {
        prev_rate = rate;
        double log_d2 = log_deriv2_neg_part(rate);
        double log_d_pos = log_deriv_pos_part(rate);
        double log_d_neg = log_deriv_neg_part;
        // determine if the value of the 1st deriv is positive or negative, and compute the
        // whole ratio to the 2nd deriv from the positive and negative parts accordingly
        if (log_d_pos > log_d_neg) {
            rate += exp(subtract_log(log_d_pos, log_d_neg) - log_d2);
        }
        else {
            rate -= exp(subtract_log(log_d_neg, log_d_pos) - log_d2);
        }
        ++iter;
    }
    return rate;
}


double fit_fixed_rate_max_exponential(const vector<double>& x, double rate, double tolerance) {
        
    // Fit N for a fixed S with the density of the maximum of N exponential variables
    //
    //   NS exp(-Sx) (1 - exp(-Sx))^(N - 1)
    //
    // where S is the rate
    // where N is the shape
    
    function<double(double)> log_likelihood = [&](double shape) {
        return max_exponential_log_likelihood(x, rate, shape);
    };
    // expand interval until we find a region where the likelihood is decreasing with
    // shape increasing
    double max_shape = 1.0;
    double max_shape_likelihood = log_likelihood(max_shape);
    double prev_max_shape_likelihood = max_shape_likelihood - 1.0;
    while (prev_max_shape_likelihood <= max_shape_likelihood) {
        prev_max_shape_likelihood = max_shape_likelihood;
        max_shape *= 2.0;
        max_shape_likelihood = log_likelihood(max_shape);
    }
    
    // use golden section search to find the maximum
    return golden_section_search(log_likelihood, 0.0, max_shape, tolerance);
}

pair<double, double> fit_max_exponential(const vector<double>& x,
                                         double tolerance) {

    // set a maximum so this doesn't get in an infinite loop even when numerical issues
    // prevent convergence
    size_t max_iters = 1000;
    size_t iter = 0;
    
    // alternate maximizing shape and rate until convergence
    double shape = 1.0;
    double rate = fit_fixed_shape_max_exponential(x, shape, tolerance / 2.0);
    double prev_shape = shape + 10.0 * tolerance;
    double prev_rate = rate + 10.0 * tolerance;
    while ((abs(prev_rate / rate - 1.0) > tolerance / 2.0
            || abs(prev_shape / shape - 1.0) > tolerance / 2.0)
           && iter < max_iters) {
        prev_shape = shape;
        prev_rate = rate;
        
        shape = fit_fixed_rate_max_exponential(x, rate, tolerance / 2.0);
        rate = fit_fixed_shape_max_exponential(x, shape, tolerance / 2.0);
        
        ++iter;
    }
    
    return pair<double, double>(rate, shape);
}

//tuple<double, double, double> fit_offset_max_exponential(const vector<double>& x,
//                                                         const function<double(double)>& shape_prior,
//                                                         double tolerance) {
//
//    // the max log likelihood of the data for a fixed location parameter
//    function<double(double)> fit_log_likelihood = [&](double loc) {
//        vector<double> x_offset(x.size());
//        for (size_t i = 0; i < x.size(); ++i) {
//            x_offset[i] = x[i] - loc;
//        }
//        pair<double, double> params = fit_max_exponential(x_offset);
//        return max_exponential_log_likelihood(x, params.first, params.second, loc) + log(shape_prior(shape));
//    };
//
//    // the maximum value of location so that all data points are in the support
//    double max_loc = *min_element(x.begin(), x.end());
//    // search with exponentially expanding windows backward to find the window
//    // that contains the highest likelihood MLE for the location
//    double min_loc = max_loc - 1.0;
//    double log_likelihood = numeric_limits<double>::lowest();
//    double probe_log_likelihood = fit_log_likelihood(min_loc);
//    while (probe_log_likelihood > log_likelihood) {
//        log_likelihood = probe_log_likelihood;
//        double probe_loc = max_loc - 2.0 * (max_loc - min_loc);
//        probe_log_likelihood = fit_log_likelihood(probe_loc);
//        min_loc = probe_loc;
//    }
//
//    // find the MLE location
//    double location = golden_section_search(fit_log_likelihood, min_loc, max_loc, tolerance);
//
//    // fit the scale and shape given the locatino
//    vector<double> x_offset(x.size());
//    for (size_t i = 0; i < x.size(); ++i) {
//        x_offset[i] = x[i] - location;
//    }
//    auto params = fit_max_exponential(x_offset);
//
//    return make_tuple(params.first, params.second, location);
//}

double max_exponential_log_likelihood(const vector<double>& x, double rate, double shape,
                                      double location) {
    double accumulator_1 = 0.0;
    double accumulator_2 = 0.0;
    for (const double& val : x) {
        if (val <= location) {
            // this should be -inf, but doing this avoids some numerical problems
            continue;
        }
        accumulator_1 += log(1.0 - exp(-rate * (val - location)));
        accumulator_2 += (val - location);
    }
    return x.size() * log(rate * shape) - rate * accumulator_2 + (shape - 1.0) * accumulator_1;
}

pair<double, double> fit_weibull(const vector<double>& x) {
    // Method adapted from Datsiou & Overend (2018) Weibull parameter estimation and
    // goodness-of-fit for glass strength data
    
    assert(x.size() >= 3);
    
    vector<double> x_local = x;
    sort(x_local.begin(), x_local.end());
    
    // regress the transformed ordered data points against the inverse CDF
    vector<vector<double>> X(x_local.size() - 1, vector<double>(2, 1.0));
    vector<double> y(X.size());
    for (size_t i = 1; i < x_local.size(); ++i) {
        X[i - 1][1] = log(x_local[i]);
        y[i - 1] = log(-log(1.0 - double(i) / double(x.size())));
    }
    vector<double> coefs = regress(X, y);
    
    // convert the coefficients into the parameters
    return make_pair(exp(-coefs[0] / coefs[1]), coefs[1]);
}

tuple<double, double, double> fit_offset_weibull(const vector<double>& x,
                                                 double tolerance) {
    
    // the max log likelihood of the data for a fixed location parameter
    function<double(double)> fit_log_likelihood = [&](double loc) {
        vector<double> x_offset(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            x_offset[i] = x[i] - loc;
        }
        pair<double, double> params = fit_weibull(x_offset);
        return weibull_log_likelihood(x, params.first, params.second, loc);
    };
    
    // the maximum value of location so that all data points are in the support
    double max_loc = *min_element(x.begin(), x.end());
    
    // search with exponentially expanding windows backward to find the window
    // that contains the highest likelihood MLE for the location
    double min_loc = max_loc - 1.0;
    double log_likelihood = numeric_limits<double>::lowest();
    double probe_log_likelihood = fit_log_likelihood(min_loc);
    while (probe_log_likelihood > log_likelihood) {
        log_likelihood = probe_log_likelihood;
        double probe_loc = max_loc - 2.0 * (max_loc - min_loc);
        probe_log_likelihood = fit_log_likelihood(probe_loc);
        min_loc = probe_loc;
    }
    
    // find the MLE location
    double location = golden_section_search(fit_log_likelihood, min_loc, max_loc, tolerance);
    
    // fit the scale and shape given the locatino
    vector<double> x_offset(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        x_offset[i] = x[i] - location;
    }
    auto params = fit_weibull(x_offset);
    
    return make_tuple(params.first, params.second, location);
}

double weibull_log_likelihood(const vector<double>& x, double scale, double shape,
                              double location) {
    double sum_1 = 0.0, sum_2 = 0.0;
    for (const double& val : x) {
        sum_1 += log(val - location);
        sum_2 += pow((val - location) / scale, shape);
    }
    return x.size() * (log(shape) - shape * log(scale)) + (shape - 1.0) * sum_1 - sum_2;
}

double golden_section_search(const function<double(double)>& f, double x_min, double x_max,
                             double tolerance) {
     
    const static double inv_phi = (sqrt(5.0) - 1.0) / 2.0;
    
    // the number of steps needed to achieve the required precision (precalculating avoids
    // fiddly floating point issues on the breakout condition)
    size_t steps = size_t(ceil(log(tolerance / (x_max - x_min)) / log(inv_phi)));
    
    // the two interior points we will evaluate the function at
    double x_lo = x_min + inv_phi * inv_phi * (x_max - x_min);
    double x_hi = x_min + inv_phi * (x_max - x_min);
    
    // the function value at the two interior points
    double f_lo = f(x_lo);
    double f_hi = f(x_hi);
    
    for (size_t step = 0; step < steps; ++step) {
        if (f_lo < f_hi) {
            // there is a max in one of the right two sections
            x_min = x_lo;
            x_lo = x_hi;
            x_hi = x_min + inv_phi * (x_max - x_min);
            f_lo = f_hi;
            f_hi = f(x_hi);
        }
        else {
            // there is a max in one of the left two sections
            x_max = x_hi;
            x_hi = x_lo;
            x_lo = x_min + inv_phi * inv_phi * (x_max - x_min);
            f_hi = f_lo;
            f_lo = f(x_lo);
        }
    }
    
    // return the midpoint of the interval we narrowed down to
    if (f_lo > f_hi) {
        return (x_min + x_hi) / 2.0;
    }
    else {
        return (x_lo + x_max) / 2.0;
    }
}

double phred_to_prob(uint8_t phred) {
    // Use a statically initialized lookup table
    static std::vector<double> prob_by_phred([](void) -> std::vector<double> {
        std::vector<double> to_return;
        to_return.reserve((int)numeric_limits<uint8_t>::max() + 1);
        for (int i = 0; i <= numeric_limits<uint8_t>::max(); i++) {
            to_return.push_back(phred_to_prob((double) i));
        }
        return to_return;
    }());
    
    // Look up in it
    return prob_by_phred[phred];
}

double phred_for_at_least_one(size_t p, size_t n) {

    /**
     * Assume that we have n <= MAX_AT_LEAST_ONE_EVENTS independent events with probability p each.
     * Let x be the AT_LEAST_ONE_PRECISION most significant bits of p. Then
     *
     *   phred_at_least_one[(n << AT_LEAST_ONE_PRECISION) + x]
     *
     * is an approximate phred score of at least one event occurring.
     *
     * We exploit the magical thread-safety of static local initialization to
     * fill this in exactly once when needed.
     */
    static std::vector<double> phred_at_least_one([](void) -> std::vector<double> {
        // Initialize phred_at_least_one by copying from the result of this function.
        std::vector<double> to_return;
        size_t values = static_cast<size_t>(1) << AT_LEAST_ONE_PRECISION;
        to_return.resize((MAX_AT_LEAST_ONE_EVENTS + 1) * values, 0.0);
        for (size_t n = 1; n <= MAX_AT_LEAST_ONE_EVENTS; n++) {
            for (size_t p = 0; p < values; p++) {
                // Because each p represents a range of probabilities, we choose a value
                // in the middle for the approximation.
                double probability = (2 * p + 1) / (2.0 * values);
                // Phred for at least one out of n.
                to_return[(n << AT_LEAST_ONE_PRECISION) + p] = prob_to_phred(1.0 - std::pow(1.0 - probability, n));
            }
        }
        return to_return;
    }());
    
    // Make sure we don't go out of bounds.
    assert(n <= MAX_AT_LEAST_ONE_EVENTS);
    
    p >>= 8 * sizeof(size_t) - AT_LEAST_ONE_PRECISION;
    return phred_at_least_one[(n << AT_LEAST_ONE_PRECISION) + p];
}

// This is just like phred_for_at_least_one but we don't prob_to_phred
// TODO: combine the code somehow?
double prob_for_at_least_one(size_t p, size_t n) {

    /**
     * Assume that we have n <= MAX_AT_LEAST_ONE_EVENTS independent events with probability p each.
     * Let x be the AT_LEAST_ONE_PRECISION most significant bits of p. Then
     *
     *   prob_at_least_one[(n << AT_LEAST_ONE_PRECISION) + x]
     *
     * is an approximate probability of at least one event occurring.
     *
     * We exploit the magical thread-safety of static local initialization to
     * fill this in exactly once when needed.
     */
    static std::vector<double> prob_at_least_one([](void) -> std::vector<double> {
        // Initialize prob_at_least_one by copying from the result of this function.
        std::vector<double> to_return;
        size_t values = static_cast<size_t>(1) << AT_LEAST_ONE_PRECISION;
        to_return.resize((MAX_AT_LEAST_ONE_EVENTS + 1) * values, 0.0);
        for (size_t n = 1; n <= MAX_AT_LEAST_ONE_EVENTS; n++) {
            for (size_t p = 0; p < values; p++) {
                // Because each p represents a range of probabilities, we choose a value
                // in the middle for the approximation.
                double probability = (2 * p + 1) / (2.0 * values);
                // Prob for at least one out of n.
                to_return[(n << AT_LEAST_ONE_PRECISION) + p] = 1.0 - std::pow(1.0 - probability, n);
            }
        }
        return to_return;
    }());
    
    // Make sure we don't go out of bounds.
    assert(n <= MAX_AT_LEAST_ONE_EVENTS);
    
    p >>= 8 * sizeof(size_t) - AT_LEAST_ONE_PRECISION;
    return prob_at_least_one[(n << AT_LEAST_ONE_PRECISION) + p];
}

vector<vector<double>> transpose(const vector<vector<double>>& A) {
    vector<vector<double>> AT(A.front().size());
    for (size_t i = 0; i < AT.size(); ++i) {
        AT[i].resize(A.size());
        for (size_t j = 0; j < A.size(); ++j) {
            AT[i][j] = A[j][i];
        }
    }
    return AT;
}

vector<vector<double>> matrix_multiply(const vector<vector<double>>& A,
                                       const vector<vector<double>>& B) {
    assert(A.front().size() == B.size());
    
    vector<vector<double>> AB(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        AB[i].resize(B.front().size(), 0.0);
        for (size_t j = 0; j < B.front().size(); ++j) {
            for (size_t k = 0; k < B.size(); ++k) {
                AB[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return AB;
}

vector<double> matrix_multiply(const vector<vector<double>>& A,
                               const vector<double>& b) {
    assert(A.front().size() == b.size());
    
    vector<double> Ab(A.size(), 0.0);
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A.front().size(); ++j) {
            Ab[i] += A[i][j] * b[j];
        }
    }
    return Ab;
}

vector<vector<double>> matrix_invert(const vector<vector<double>>& A) {
    
    // invert by Gaussian elimination
    
    assert(A.front().size() == A.size());
    
    vector<vector<double>> A_inv(A.size());
    
    for (size_t i = 0; i < A.size(); ++i) {
        A_inv[i].resize(A.size(), 0.0);
        A_inv[i][i] = 1.0;
    }
    
    // a non-const local copy
    auto A_loc = A;
    
    // forward loop, make upper triangular
    
    for (int64_t i = 0; i < A_loc.size(); ++i) {
        int64_t ii = i;
        while (A_loc[ii][i] == 0.0 && ii < A_loc.size()) {
            ++ii;
        }
        if (ii == A_loc.size()) {
            std::runtime_error("error: matrix is not invertible!");
            
        }
        swap(A_loc[i],A_loc[ii]);
        swap(A_inv[i], A_inv[ii]);
                
        // make the diagonal entry 1
        double factor = A_loc[i][i];
        for (int64_t j = 0; j < A_loc.size(); ++j) {
            A_loc[i][j] /= factor;
            A_inv[i][j] /= factor;
        }
        
        // make the off diagonals in one column 0's
        for (ii = i + 1; ii < A_loc.size(); ++ii) {
            factor = A_loc[ii][i];
            for (size_t j = 0; j < A_loc.size(); ++j) {
                A_loc[ii][j] -= factor * A_loc[i][j];
                A_inv[ii][j] -= factor * A_inv[i][j];
            }
            
        }
    }
    
    // backward loop, make identity
    
    for (int64_t i = A_loc.size() - 1; i >= 0; --i) {
        // make the off diagonals in one column 0's
        for (int64_t ii = i - 1; ii >= 0; --ii) {
            double factor = A_loc[ii][i];
            for (size_t j = 0; j < A_loc.size(); ++j) {
                A_loc[ii][j] -= factor * A_loc[i][j];
                A_inv[ii][j] -= factor * A_inv[i][j];
            }
        }
    }
    
    return A_inv;
}

vector<double> regress(const vector<vector<double>>& X, vector<double>& y) {
    auto X_t = transpose(X);
    return matrix_multiply(matrix_multiply(matrix_invert(matrix_multiply(X_t, X)), X_t), y);
}

}
