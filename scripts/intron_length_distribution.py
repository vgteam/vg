#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import math
import numpy as np
import scipy.optimize as opt
import scipy.stats as stats
import scipy.special as special
from collections import Counter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
    
def parse_intron_length_distr(filepath, label, unique = False):
    
    observed_introns = set()
    lengths = Counter()
    with open(filepath) as f:
        curr_tx = None
        prev_end = None
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split("\t")
            if tokens[2] != "exon":
                continue
            chrom = tokens[0]
            start = None
            end = None
            strand = (tokens[6] == "+")
            if strand:
                start = int(tokens[3])
                end = int(tokens[4])
            else:
                start = int(tokens[4])
                end = int(tokens[3])
            annotations = tokens[8].strip().split(";")
            tx_id = None
            for annotation in annotations:
                annotation = annotation.replace("\"", "").strip()
                if annotation.startswith(label):
                    tx_id = annotation.split()[1]
                    break
            assert(tx_id is not None)
            intron = (chrom, strand, prev_end, start)
            if tx_id == curr_tx and (not unique or intron not in observed_introns):
                # inclusive indexing on interval
                lengths[abs(start - prev_end) - 1] += 1
                observed_introns.add(intron)
            else:
                curr_tx = tx_id
            prev_end = end
    return lengths
    
def frechet_log_likelihood(x, a, s, m):
    if x < m:
        return -sys.float_info.max * np.ones_like(x)
    z = (x - m) / s
    return np.log(a / s) - (a + 1.0) * np.log(z) - np.power(z, -a)

def frechet_pdf_raw(x, a, s, m):
    f1 = np.log(a / s)
    f2 = (-a - 1) * np.log((x - m) / s)
    f3 = -np.power((x - m) / s, -a)
    return np.exp(f1 + f2 + f3)

def frechet_pdf(x, a, s, m):
    if np.isscalar(x):
        if x > m:
            return frechet_pdf_raw(x, a, s, m)
        else:
            return 0.0
    else:
        p = np.zeros_like(x)
        mask = x > m
        x_pos = x[mask]
        p[mask] = frechet_pdf_raw(x_pos, a, s, m)
        return p

def frechet_dshape(x, a, s, m):
    z = (x - m) / s
    return 1.0 / a + np.log(z) * (np.power(z, -a) - 1.0)

def frechet_dscale(x, a, s, m):
    return (a / s) * (1.0 - np.power((x - m) / s, -a))

def frechet_dlocation(x, a, s, m):
    return (1.0 + a - np.power((x - m) / s, -a)) / (x - m)

def frechet_mean(a, s, m):
    return m + s * special.gamma(1.0 - 1.0 / a)

def frechet_variance(a, s):
    k1 = special.gamma(1.0 - 1.0 / a)
    return s * s * (special.gamma(1.0 - 2.0 / a) - k1 * k1)
 
def frechet_skewness(a):
    assert(a > 3.0)
    k1 = special.gamma(1.0 - 1.0 / a)
    k2 = special.gamma(1.0 - 2.0 / a)
    k3 = special.gamma(1.0 - 3.0 / a)
    return (k3 - 3.0 * k1 * k2 + 2.0 * k1 * k1 * k1) / np.power(k2 - k1 * k1, 1.5)

# exponential search for method of moments
def find_shape(skewness, tol = 1e-8):
    assert(skewness > 0.0)
    lo = 3.0 # skewness is not defined for shape <= 3
    step = 1.0
    while frechet_skewness(lo + step) > skewness:
        lo += step
        step *= 2
    hi = lo + step
    while hi - lo > tol:
        mid = 0.5 * (hi + lo)
        if frechet_skewness(mid) > skewness:
            lo = mid
        else:
            hi = mid
    return 0.5 * (hi + lo)

def normal_pdf(x, s, m):
    z = (x - m) / s
    return np.exp(-z * z / 2.0) / (s * 2.5066282746310002) # root(2pi)


# fit frechet distribution with method of moments
def frechet_mom(mean, variance, skewness):
    # skewness depends only on shape
    shape = find_shape(skewness)
    # variance depends only on shape and scale
    k1 = math.gamma(1.0 - 1.0 / shape)
    k2 = math.gamma(1.0 - 2.0 / shape)
    scale = (variance / (k2 - k1 * k1))**0.5
    # location depends on all three parameters
    location = mean - scale * k1
    return shape, scale, location


# returns tuple of (frechet params, mixture weights)
def frechet_M_step(vals, counts, Z, params):
    
    assert(len(vals) == len(counts))
    assert(len(vals) == Z.shape[0])
    assert(Z.shape[1] * 3 == len(params))
    
    
    # mixture weights have an easy closed form solution
    m0 = np.sum(Z * counts[:,np.newaxis], 0)
    mix_weights = m0 / np.sum(m0)     
    
    # # get the cumulants
    # k1 = np.sum(Z * (vals * counts)[:,np.newaxis], 0) / m0
    # k2 = np.sum(Z * (np.power(vals, 2) * counts)[:,np.newaxis], 0) / m0
    # k3 = np.sum(Z * (np.power(vals, 3) * counts)[:,np.newaxis], 0) / m0
    # # compute central moments from the cumulants
    # mean = k1
    # var = k2 - np.power(k1, 2)
    # skew = (k3 - 3 * mean * var - np.power(mean, 3)) / np.power(var, 1.5)
    # # choose initial values for maximization
    # init_params = np.zeros_like(params)
    # for j in range(Z.shape[1]):
    #     if skew[j] > 0.0:
    #         # use method of moments to choose good initial values
    #         init_params[3*j:3*j+3] = np.array(frechet_mom(mean[j], var[j], skew[j]))
    #     else:
    #         # it can happen that sample skew is negative by random
    #         # sampling. if so, the method of moments is undefined, so 
    #         # we just use previous value as the starting position
    #         init_params[3*j:3*j+3] = params[3*j:3*j+3]
    init_params = params
    
    # function to compute negative log likelihood
    def neg_expected_log_likelihood(par):
        ll = 0.0
        for j in range(Z.shape[1]):
            a, s, m = par[3*j:3*j+3]
            for i in range(Z.shape[0]):
                v = vals[i] 
                if v <= m:
                    continue
                lx = frechet_log_likelihood(v, a, s, m)
                ll += counts[i] * Z[i,j] * lx
        return -ll   
    
    # function to compute negative gradient
    def neg_gradient(par):
        grad = np.zeros_like(par)
        for j in range(Z.shape[1]):
            a, s, m = par[3*j:3*j+3]
            for i in range(Z.shape[0]):
                v = vals[i] 
                if v <= m:
                    continue
                C = Z[i,j] * counts[i]
                grad[3*j] += C * frechet_dshape(v, a, s, m)
                grad[3*j+1] += C * frechet_dscale(v, a, s, m)
                grad[3*j+2] += C * frechet_dlocation(v, a, s, m)
        return -grad
    
    
    bounds = []
    for k in range(len(params)):
        if k % 3 != 2:
            # make sure the scale and shape are positive
            bounds.append((1e-12, None))
        else:
            # keep m from leaving the support
            j = k // 3
            max_m = min(vals[i] for i in range(len(vals)) if Z[i,j] != 0.0) - 1e-12
            bounds.append((None, max_m))
    
    # TODO: the documentation in scipy says that both nelder-mead and full
    # memory BGFS should be able to take bounds, but i can only do it with
    # L-BGFS-B
    
    # # improve the initial guess with nelder-mead, which supposedly is less
    # # likely to be fooled by shallow regions
    # res = opt.minimize(fun = neg_expected_log_likelihood,
    #                    x0 = init_params,
    #                    bounds = bounds,
    #                    method = "Nelder-Mead",
    #                    options = {"maxiter":25})
    
    # use gradient-based quasi-newton method to refine optimum
    res = opt.minimize(fun = neg_expected_log_likelihood,
                       x0 = init_params,
                       jac = neg_gradient,
                       bounds = bounds,
                       method = "L-BFGS-B") 
    return res.x, mix_weights

def frechet_E_step(vals, Z, params, mix_weights, trunc):
    for i in range(Z.shape[0]):
        row = np.array([frechet_pdf(vals[i], params[3*j], params[3*j+1], params[3*j+2]) for j in range(Z.shape[1])])
        Z[i,:] = row / np.sum(row)

def init_Z_matrix(num_vals, num_comps, asymm):
    
    # randomly initialize the assignment matrix
    Z = np.ones((num_vals, num_comps))
    if num_comps > 1:
        for i in range(Z.shape[0]):
            # make a dirichlet shape that biases smaller lengths to be assigned to
            # earlier components
            # TODO: this is pretty hacky, but it should break symmetry...
            dir_shape = np.ones(num_comps)
            if Z.shape[0] > 1:
                for j in range(num_comps):
                    i_frac = i / (Z.shape[0] - 1.0)
                    j_frac = j / (num_comps - 1.0)
                    if j_frac == 0.0:
                        dir_shape[j] += asymm * (1.0 - i_frac)
                    elif j_frac == 1.0:
                        dir_shape[j] += asymm * i_frac
                    elif i_frac <= j_frac:
                        dir_shape[j] += asymm * (i_frac / j_frac)
                    else:
                        dir_shape[j] += asymm * ((1.0 - i_frac) / (1.0 - j_frac))
            Z[i,:] = stats.dirichlet.rvs(dir_shape)[0]
    return Z

def truncate_Z_matrix(Z, trunc_ratio):
    # clip very unlikely assignments to 0 to release the location parameter
    # to give them 0 probability
    for i in range(Z.shape[0]):
        M = np.max(Z[i,:])
        Z[i, Z[i,:] < M * trunc_ratio] = 0.0
        Z[i,:] /= np.sum(Z[i,:])

def fit_frechet_mixture(counter, num_comps, max_iters = 100, asymm = 2.0, tol = 1e-3, verbose = False):
    
    vals = sorted(counter)
    counts = np.array([float(counter[v]) for v in vals])
    vals = np.array([float(v) for v in vals])
    
    # randomly initialize the assignment matrix
    
    
    mix_weights = np.ones(num_comps) / num_comps
    params = np.zeros(3 * num_comps)
    for i in range(num_comps):
        # shape and scale should be positive
        params[3 * i] = 1.0
        params[3 * i + 1] = 1.0
        
    # fit a log-normal model initialize the Z matrix (it is much more stable
    # and easy to fit)
    if verbose:
        print("initial lognormal EM", file = sys.stderr)
    Z = fit_log_normal_mixture(counter, num_comps, max_iters, tol, asymm, return_Z = True, verbose = verbose)
    
    if verbose:
        print("main frechet EM", file = sys.stderr)
    for it in range(max_iters):
        if verbose and (it + 1) % 100 == 0:
            print("EM iter {}".format(it + 1), file = sys.stderr)
        truncate_Z_matrix(Z, 1e-5)
        prev_params = params
        prev_mix_weights = mix_weights
        
        # estimated the parameters based on the conditional likelihood matrix
        params, mix_weights = frechet_M_step(vals, counts, Z, params)
        

        
        # check for convergence
        converged = True
        err = 0.0
        for j in range(num_comps):
            err1 = abs(prev_mix_weights[j] - mix_weights[j])
            err2 = abs(params[3*i] - prev_params[3*i]) / prev_params[3*i]
            err3 = abs(params[3*i+1] - prev_params[3*i+1]) / prev_params[3*i+1]
            err4 = abs(params[3*i+2] - prev_params[3*i+2])
            err = max(err, err1, err2, err3, err4)
            if (err1 > tol or err2 > tol or err3 > tol or err4 > tol):
                converged = False
        if verbose and (it + 1) % 100 == 0:
            print("\tmax param velocity {}".format(err), file=sys.stderr)
        if converged:
            break
        
        # recompute the conditional likelihood matrix with new params
        frechet_E_step(vals, Z, params, mix_weights, 1e-10)
        
    return params, mix_weights
        

def normal_E_step(vals, Z, params, mix_weights):
    assert(len(mix_weights) == Z.shape[1])
    assert(len(mix_weights) * 2 == len(params))
    assert(len(vals) == Z.shape[0])
    
    for i in range(Z.shape[0]):
        row = np.array([mix_weights[j] * normal_pdf(vals[i], params[2*j], params[2*j+1]) for j in range(Z.shape[1])])
        Z[i,:] = row / np.sum(row)
        
def normal_M_step(vals, counts, Z):
    
    m0 = np.sum(Z * counts[:,np.newaxis], 0)
    m1 = np.sum(Z * (vals * counts)[:,np.newaxis], 0) / m0
    m2 = np.sqrt(np.sum(((Z * np.power(np.ones_like(m0) * vals[:,np.newaxis] - m1, 2)) * counts[:,np.newaxis]), 0) / m0)
    
    mix_weights = m0 / np.sum(m0)
    params = np.array([m2[i//2] if i % 2 == 0 else m1[i//2] for i in range(2*Z.shape[1])])
    
    return params, mix_weights



def fit_log_normal_mixture(counter, num_comps, max_iters = 100, tol = 1e-3, asymm = 2.0, return_Z = False, verbose = False):
    assert(0 not in counter)
    
    vals = sorted(counter)
    counts = np.array([float(counter[v]) for v in vals])
    vals = np.log(np.array(vals))
    
    Z = init_Z_matrix(len(vals), num_comps, asymm)
    
    params = np.array([1.0 if i % 2 == 0 else 0.0 for i in range(2 * num_comps)])
    mix_weights = np.ones(num_comps) / num_comps
    
    for it in range(max_iters):
        if verbose and (it + 1) % 100 == 0:
            print("EM iter {}".format(it + 1), file = sys.stderr)
        # print(Z)
        prev_params = params
        prev_mix_weights = mix_weights
        
        # estimated the parameters based on the conditional likelihood matrix
        params, mix_weights = normal_M_step(vals, counts, Z)
        
        # check for convergence
        converged = True
        err = 0.0
        for j in range(num_comps):
            err1 = abs(prev_mix_weights[j] - mix_weights[j])
            err2 = abs(params[2*j] - prev_params[2*j]) / prev_params[2*j]
            err3 = abs(params[2*j+1] - prev_params[2*j+1])
            err = max(err, err1, err2, err3)
            if (err1 > tol or err2 > tol or err3 > tol):
                converged = False
        if verbose and (it + 1) % 100 == 0:
            print("\tmax param velocity {}".format(err), file = sys.stderr)
        if converged:
            break
        
        # recompute the conditional likelihood matrix with new params
        normal_E_step(vals, Z, params, mix_weights)
    
    if return_Z:
        return Z
    else:
        return params, mix_weights
 

def frechet_mixture_pdf(x, params, mix_weights):
    d = 0.0
    for j in range(len(mix_weights)):
        d += mix_weights[j] * frechet_pdf(x, params[3*j], params[3*j+1], params[3*j+2])
    return d

def normal_mixture_pdf(x, params, mix_weights):
    d = 0.0
    for j in range(len(mix_weights)):
        d += mix_weights[j] * normal_pdf(x, params[2*j], params[2*j+1])
    return d

def plot_frechet_log_likelihood(counter, x_min, x_max, y_min, y_max, grid, a, s, m):
    vals = sorted(counter)
    counts = np.array([float(counter[v]) for v in vals])
    vals = np.array([float(v) for v in vals])
    
    X = np.array([np.linspace(x_min, x_max, grid) for i in range(grid)])
    Y = np.array([[y for j in range(grid)] for y in np.linspace(y_min, y_max, grid)])
    Z = []
    for i in range(grid):
        Z.append([])
        for j in range(grid):
            sh, sc, l = None, None, None
            if a == "x":
                sh = X[i,j]
            elif a == "y":
                sh = Y[i,j]
            else:
                sh = a
            if s == "x":
                sc = X[i,j]
            elif s == "y":
                sc = Y[i,j]
            else:
                sc = s
            if m == "x":
                l = X[i,j]
            elif m == "y":
                l = Y[i,j]
            else:
                l = m
            Z[-1].append(sum(c * frechet_log_likelihood(v, sh, sc, l) for v,c in zip(vals, counts)))
    
    Z = np.array(Z)
    
    elev = 75
    azim = 225
    
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(elev=elev, azim=azim)
    ax.plot_wireframe(X, Y, Z)
    plt.show()
    
    fig = plt.contour(X, Y, Z)
    plt.show()

def log_flatten(counter):
    x = np.zeros(sum(counter.values()) - counter[0])
    i = 0
    for l in counter:
        if l == 0:
            continue
        c = counter[l]
        x[i:i+c] = np.log10(l)
        i += c
    return x

def plot_length_distribution(counter):

    x = log_flatten(counter)
    
    plt.hist(x, bins = 100)

def plot_log_normal_mixture(params, mix_weights, log_density = False, num_sigmas = 4.0, log_x = True, xlim = None):
    
    log_xmin, log_xmax = None, None
    for i in range(len(mix_weights)):
        s = params[2 * i]
        m = params[2 * i + 1]
        
        lo = m - num_sigmas * s
        hi = m + num_sigmas * s
        if log_xmin is None or lo < log_xmin:
            log_xmin = lo
        if log_xmax is None or hi > log_xmax:
            log_xmax = hi
    
    
    logx = np.linspace(log_xmin, log_xmax, 500)
    x = np.exp(logx)
    log10x = logx / np.log(10.)
    
    y = np.zeros_like(logx)
    for i in range(len(logx)):
        d = normal_mixture_pdf(logx[i], params, mix_weights) / x[i]
        if log_density:
            y[i] = np.log(d)
        else:
            y[i] = d
    
    if log_density:
        y -= np.max(y)
    
    if log_x:
        plt.plot(log10x, y)
        ticks, labels = plt.xticks()
        plt.xticks(ticks, [str(int(round(t))) for t in 10**np.array(ticks)])
    else:
        plt.plot(x, y)
        
    if xlim is not None:
        plt.xlim(xlim)

def plot_frechet_mixture(params, mix_weights):
    
    log_xmin = 0
    log_xmax = 1
    for i in range(len(mix_weights)):
        a = params[3*i]
        s = params[3*i+1]
        m = params[3*i+2]
        if a > 2:
            ev = m + s * np.gamma(1.0 - 1.0 / a)
            sd = s * np.sqrt(np.gamma(1.0 - 2.0 / a) - np.power(np.gamma(1.0 - 1.0 / a), 2))
            hi = ev + 8 * sd
        else:
            mode = m + s * pow(a / (1.0 + a), 1.0 / a)
            hi = m + 10 * (mode - m)
        
        log_xmax = max(log_xmax, hi)
      
    logx = np.linspace(log_xmin, log_xmax, 500)
    log10x = logx / np.log(10.)
    y = np.zeros_like(logx)
    for i in range(len(logx)):
        y[i] = frechet_mixture_pdf(np.exp(logx[i]), params, mix_weights)
            
    plt.plot(log10x, y)


def log_normal_BIC(counter, params, mix_weights):
    vals = sorted(counter)
    counts = np.array([float(counter[v]) for v in vals])
    vals = np.array([float(v) for v in vals])
    log_vals = np.log(vals)
    # one parameter isn't free in the weights
    p = len(params) + len(mix_weights) - 1
    N = np.sum(counts)
    
    log_likelihood = 0.0
    for i in range(len(vals)):
        likelihood = 0.0
        for j in range(len(mix_weights)):
            likelihood += mix_weights[j] * normal_pdf(log_vals[i], params[2*j], params[2*j+1]) / vals[i]
        log_likelihood += counts[i] * np.log(likelihood)
    
    return float(p) * np.log(N) - 2.0 * log_likelihood
    

def main():
    parser = argparse.ArgumentParser(description='Estimate the intron length distribution from transcript annotations')
    
    parser.add_argument('-g', '--gtf', required = True, type = str, metavar = "FILE",
                        help = 'GTF/GFF file containing transcript annotations')
    parser.add_argument('-o', '--out', required = True, type = str, metavar = "FILE",
                        help = 'where to save distribution')
    parser.add_argument('-l', '--label', required = False, type = str, metavar = "STR",
                        default = 'transcript_id',
                        help = 'label of transcripts in GTF (default: transcript_id)')
    
    args = parser.parse_args()
    print("parsing transcript annotations from {}...".format(args.gtf), file = sys.stderr)
    lengths = parse_intron_length_distr(args.gtf, args.label)
    
    training_results = {}
    for num_comps in range(1, 6):
        print("training model with {} components...".format(num_comps), file = sys.stderr)
        params, weights = fit_log_normal_mixture(lengths, 
                                                 num_comps, 
                                                 max_iters = 500 * num_comps,
                                                 asymm = 5 * num_comps,
                                                 tol = 1e-4)
        training_results[num_comps] = params, weights
        
        
    best_BIC = None
    best_num_comps = None
    for num_comps in sorted(training_results):
        params, weights = training_results[num_comps]
        BIC = log_normal_BIC(lengths, params, weights)
        print("{} component model achieves BIC {}".format(num_comps, BIC), file = sys.stderr)
        if best_BIC is None or BIC < best_BIC:
            best_BIC = BIC
            best_num_comps = num_comps
                    
    print("outputting {} component model to {}".format(best_num_comps, args.out), file = sys.stderr)
    params, weights = training_results[best_num_comps]
    with open(args.out, 'w') as f:
        print(str(best_num_comps), file = f)
        for param in list(weights) + list(params):
            print(str(param), file = f)
  

if __name__ == "__main__":
    main()