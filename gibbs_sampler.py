#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import gammaln
from scipy import sparse as sp
from scipy.stats import beta
import itertools

def log_dirichlet_pdf(x, alpha):
    return float(np.sum((alpha - 1) * np.log(x)) - (np.sum(gammaln(alpha)) - gammaln(np.sum(alpha))))

def _proj_floor_and_renorm(p, eps):
    p = np.asarray(p, dtype=float).copy()
    m = p < eps
    if m.any():
        p[m] = eps
        s_not = p[~m].sum()
        if s_not <= 0:
            p = np.ones_like(p) / p.size
        else:
            p[~m] *= (1.0 - eps * m.sum()) / s_not
    else:
        s = p.sum()
        if s <= 0:
            p = np.ones_like(p) / p.size
        else:
            p /= s
    return p

def effective_sample_size(x, max_lag=20000):
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = x.size
    if n < 20:
        return float(n)
    fftlen = 1 << (int(np.ceil(np.log2(2 * n - 1))))
    fx = np.fft.rfft(x, fftlen)
    ac = np.fft.irfft(fx * np.conjugate(fx))[:n]
    ac = np.real(ac)
    ac /= max(ac[0], 1e-20)
    rho = ac[1:]
    neg = np.where(rho < 0)[0]
    cutoff = neg[0] if neg.size > 0 else min(n - 1, max_lag)
    tau = 1 + 2 * np.sum(rho[:cutoff])
    ess = n / max(tau, 1e-6)
    return float(max(1.0, ess))

def eps_from_H(H, cap=1e-4, safety=100):
    e = min(cap, 1.0 / (safety * H))
    return min(e, 0.49 / H)

def _colsum_Z(C):
    Z = np.asarray(C.sum(axis=0)).ravel().astype(float)
    Z[Z == 0.0] = 1.0   
    return Z

def _preprocess_patterns(patterns, alleles_per_site):
    compiled = []
    for pat in patterns:
        req = []
        for j, code in enumerate(pat):
            if code == '-':
                continue
            allele_char = alleles_per_site[j][0 if code == '0' else 1]
            req.append((j, allele_char))
        compiled.append(req)
    return compiled

def build_compat_idx(haplotypes, compiled_patterns):
    idx = []
    for req in compiled_patterns:
        ok = []
        for h, hap in enumerate(haplotypes):
            good = True
            for (j, allele_char) in req:
                if hap[j] != allele_char:
                    good = False
                    break
            if good:
                ok.append(h)
        idx.append(np.array(ok, dtype=int))
    return idx

def build_compat_matrix(haplotypes, compiled_patterns, counts=None):
    rows, cols, data = [], [], []
    for r, req in enumerate(compiled_patterns):
        for h, hap in enumerate(haplotypes):
            good = True
            for (j, allele_char) in req:
                if hap[j] != allele_char:
                    good = False
                    break
            if good:
                rows.append(r); cols.append(h); data.append(1.0)
    R = len(compiled_patterns); H = len(haplotypes)
    C = sp.csr_matrix((data, (rows, cols)), shape=(R, H), dtype=np.float64)
    if counts is not None:
        bad = (counts > 0) & (np.array(C.sum(axis=1)).ravel() == 0)
        if bad.any():
            where = np.where(bad)[0]
            print(f"[WARN] {where.size} patterns with positive counts have no compatible hap: e.g., r={where[:5]}")
    return C

def compute_log_likelihood_from_C(p, C, counts, Z=None, eps=1e-300):
    if Z is None:
        Z = _colsum_Z(C)
    q = C.dot(p / Z) + eps
    return float(np.dot(counts, np.log(q)))

def prepare_model(var_table, prior_table, enc_table, eps_floor=1e-4, rng_seed=42):
    np.random.seed(rng_seed)
    var_df = pd.read_csv(var_table, sep="\t", dtype={
        "Site": str, "REF": str, "ALT": str,
        "REF_coverage": int, "ALT_coverage": int,
        "Depth": int, "ALT_freq": float
    })
    var_df["pos"] = var_df["Site"].str.split(":").str[1].astype(int)
    loci = {
        row.pos: {
            row.REF: row.REF_coverage / row.Depth if row.Depth > 0 else 0.5,
            row.ALT: row.ALT_coverage / row.Depth if row.Depth > 0 else 0.5
        }
        for row in var_df.itertuples()
    }
    positions = sorted(loci.keys())
    alleles_per_site = [list(loci[pos].keys()) for pos in positions]
    alt_base_per_site = {row.pos: row.ALT for row in var_df.itertuples()}
    cov_map = {row.pos: (int(row.ALT_coverage), int(row.Depth)) for row in var_df.itertuples()}
    site_alt_counts = np.array([cov_map[p][0] for p in positions], dtype=int)
    site_depths     = np.array([cov_map[p][1] for p in positions], dtype=int)
    haplotypes = [''.join(h) for h in itertools.product(*alleles_per_site)]
    H = len(haplotypes)
    good_guess = np.array([
        np.prod([loci[pos][base] for pos, base in zip(positions, hap)])
        for hap in haplotypes
    ], dtype=float)
    good_guess = (good_guess / good_guess.sum()) if good_guess.sum() > 0 else np.ones(H) / H
    uniform_guess = np.ones(H) / H
    enc_df = pd.read_csv(enc_table, sep="\t")
    read_patterns = enc_df["pattern"].astype(str).str.replace("/", "", regex=False).tolist()
    pattern_counts = enc_df["read_count"].astype(int).to_numpy()
    prior_uniform = np.ones(H, dtype=float) / H
    prior_df = pd.read_csv(prior_table, sep="\t")
    prior_map = dict(zip(prior_df["Haplotype"].astype(str), prior_df["Frequency"].astype(float)))
    raw_emp = np.array([prior_map.get(h, 0.0) for h in haplotypes], dtype=float)
    prior_empirical = _proj_floor_and_renorm(raw_emp, eps_floor)
    prior_full = prior_empirical
    
    compiled_patterns = _preprocess_patterns(read_patterns, alleles_per_site)
    compat_idx = build_compat_idx(haplotypes, compiled_patterns)
    C = build_compat_matrix(haplotypes, compiled_patterns, counts=pattern_counts)
    Z = _colsum_Z(C)
    return {
        "haplotypes": haplotypes,
        "H": H,
        "prior_full": prior_full,
        "prior_uniform": prior_uniform,
        "prior_empirical": prior_empirical,
        "good_guess": good_guess,
        "uniform_guess": uniform_guess,
        "read_patterns": read_patterns,
        "pattern_counts": pattern_counts,
        "compiled_patterns": compiled_patterns,
        "compat_idx": compat_idx,
        "C": C,
        "Z": Z,
        "positions": positions,
        "alt_base_per_site": alt_base_per_site,
        "site_alt_counts": site_alt_counts,
        "site_depths": site_depths,
    }

def mcmc_collapsed_gibbs(
    init_p, prior_full, counts, compat_idx,
    C=None, alpha0=1e4, n_iter=100000, burnin=20000, thin=1, eps=1e-4,
    record_traces=True, rng_seed=None
):
    if rng_seed is not None:
        np.random.seed(rng_seed)
    H = prior_full.size
    p = _proj_floor_and_renorm(init_p, eps)
    alpha_prior = alpha0 * prior_full

    Z_local = None
    if C is not None:
        Z_local = _colsum_Z(C)
        def _ll(pvec):
            return compute_log_likelihood_from_C(pvec, C, counts, Z=Z_local)
    else:
        def _ll(pvec):
            ll = 0.0
            for r, k in enumerate(counts):
                if k <= 0:
                    continue
                ids = compat_idx[r]
                if ids.size == 0:
                    return -np.inf
                q = pvec[ids].sum()
                if q <= 0:
                    return -np.inf
                ll += k * np.log(q)
            return float(ll)

    trace_ll, trace_pr, trace_post = [], [], []
    samples = []
    suff = np.zeros(H, dtype=np.int64)

    for it in range(n_iter):
        suff[:] = 0
        for r, k in enumerate(counts):
            if k <= 0:
                continue
            ids = compat_idx[r]
            if ids.size == 0:
                continue
            probs = p[ids] if Z_local is None else (p[ids] / Z_local[ids])
            s = probs.sum()
            if (not np.isfinite(s)) or (s <= 0):
                probs = np.ones_like(probs) / probs.size
            else:
                probs = probs / s
            y = np.random.multinomial(int(k), probs)
            np.add.at(suff, ids, y)

        alpha_post = alpha_prior + suff + eps
        p = np.random.dirichlet(alpha_post)
        p = _proj_floor_and_renorm(p, eps)

        if record_traces:
            ll = _ll(p)
            pr = log_dirichlet_pdf(p, alpha_prior)
            trace_ll.append(ll)
            trace_pr.append(pr)
            trace_post.append(ll + pr)

        if it >= burnin and ((it - burnin) % thin == 0):
            samples.append(p.copy())

    out = {
        "samples": np.vstack(samples) if len(samples) > 0 else np.empty((0, H)),
        "trace_ll": np.array(trace_ll),
        "trace_pr": np.array(trace_pr),
        "trace_post": np.array(trace_post)
    }
    return out

def run_post_traces(init_vec, n_chains, n_iter, burnin, thin, seed_base,
                    prior_full, counts, compat_idx, C, alpha0, eps):
    traces = []
    for i in range(n_chains):
        out = mcmc_collapsed_gibbs(
            init_p=init_vec, prior_full=prior_full, counts=counts,
            compat_idx=compat_idx, C=C, alpha0=alpha0,
            n_iter=n_iter, burnin=burnin, thin=thin, eps=eps,
            record_traces=True, rng_seed=seed_base + i
        )
        traces.append(out["trace_post"][burnin:])
    return traces

def run_ll_traces(init_vec, n_chains, n_iter, burnin, thin, seed_base,
                  prior_full, counts, compat_idx, C, alpha0, eps):
    traces = []
    for i in range(n_chains):
        out = mcmc_collapsed_gibbs(
            init_p=init_vec, prior_full=prior_full, counts=counts,
            compat_idx=compat_idx, C=C, alpha0=alpha0,
            n_iter=n_iter, burnin=burnin, thin=thin, eps=eps,
            record_traces=True, rng_seed=seed_base + i
        )
        traces.append(out["trace_ll"][burnin:])
    return traces

def plot_trace_list(traces, title, ylabel, burnin, legend_prefix="Chain", absolute_index=False):
    plt.figure(figsize=(10, 5))
    for i, tr in enumerate(traces, 1):
        if absolute_index:
            x = np.arange(burnin, burnin + len(tr))
            xlabel = "Iteration"
        else:
            x = np.arange(len(tr))
            xlabel = f"Iteration (after {burnin:,} burn-in)"
        plt.plot(x, tr, alpha=0.6, label=f"{legend_prefix} {i}")
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel)
    plt.legend(loc="upper right"); plt.tight_layout(); plt.show()

def plot_topk_posteriors(samples, haplotypes, top_k=10):
    post_means = samples.mean(axis=0)
    post_vars = samples.var(axis=0)
    K = min(top_k, len(haplotypes))
    idx_top = np.argsort(post_means)[::-1][:K]
    idx_top5 = idx_top[:min(5, K)]
    idx_next5 = idx_top[min(5, K):K]
    plt.figure(figsize=(8, 5))
    for i in idx_top5:
        plt.hist(samples[:, i], bins=50, density=True, alpha=0.6,
                 label=f"{haplotypes[i]}  (μ={post_means[i]:.3f})")
    plt.xlabel(r"Posterior sample of haplotype frequency, $p_h$")
    plt.ylabel("Estimated posterior density")
    plt.title("Posterior density (Top haplotypes)")
    plt.legend(title="Haplotype (mean)", frameon=False)
    plt.tight_layout(); plt.show()
    if len(idx_next5) > 0:
        plt.figure(figsize=(8, 5))
        for i in idx_next5:
            plt.hist(samples[:, i], bins=50, density=True, alpha=0.6,
                     label=f"{haplotypes[i]} (μ={post_means[i]:.3f})")
        plt.xlabel("Frequency $p_h$"); plt.ylabel("Density")
        plt.title("Posterior density (Next haplotypes)")
        plt.legend(); plt.tight_layout(); plt.show()
    print("\n[TOP {} posterior means]".format(K))
    for rank, i in enumerate(idx_top, 1):
        print(f"{rank:2d}. {haplotypes[i]}  mean={post_means[i]:.4f}  var={post_vars[i]:.4e}")

def fit_two_inits(model, alpha0, n_iter=100000, burnin=20000, thin=1, eps=1e-4,
                  seed_data_main=24680, seed_uni_main=13579):
    H = model["H"]
    prior_full = model["prior_full"]
    good_guess = model["good_guess"]
    counts = model["pattern_counts"]
    compat_idx = model["compat_idx"]
    C = model["C"]
    out_data = mcmc_collapsed_gibbs(
        init_p=good_guess, prior_full=prior_full, counts=counts,
        compat_idx=compat_idx, C=C, alpha0=alpha0,
        n_iter=n_iter, burnin=burnin, thin=thin, eps=eps,
        record_traces=True, rng_seed=seed_data_main
    )
    out_uni = mcmc_collapsed_gibbs(
        init_p=np.ones(H) / H, prior_full=prior_full, counts=counts,
        compat_idx=compat_idx, C=C, alpha0=alpha0,
        n_iter=n_iter, burnin=burnin, thin=thin, eps=eps,
        record_traces=True, rng_seed=seed_uni_main
    )
    return out_data, out_uni

def fit_two_inits_and_plots(model, alpha0, n_iter=100000, burnin=20000, thin=1, eps=1e-4,
                            n_chains_plot=5, seed_data=24680, seed_uni=13579):
    haplotypes = model["haplotypes"]
    H = model["H"]
    prior_full = model["prior_full"]
    good_guess = model["good_guess"]
    counts = model["pattern_counts"]
    compat_idx = model["compat_idx"]
    C = model["C"]
    out_data, out_uni = fit_two_inits(
        model, alpha0, n_iter=n_iter, burnin=burnin, thin=thin, eps=eps,
        seed_data_main=seed_data, seed_uni_main=seed_uni
    )
    iters = np.arange(burnin, len(out_data["trace_post"]))
    plt.figure(figsize=(10, 5))
    plt.plot(iters, out_data["trace_post"][burnin:], alpha=0.6, label="Data-based Init")
    plt.plot(iters, out_uni["trace_post"][burnin:], alpha=0.6, label="Uniform Init")
    plt.title("Log Posterior Trace (Collapsed Gibbs, after burn-in)")
    plt.xlabel("Iteration"); plt.ylabel("Log Posterior"); plt.legend(loc="upper right")
    plt.tight_layout(); plt.show()
    iters = np.arange(burnin, len(out_data["trace_ll"]))
    plt.figure(figsize=(10, 5))
    plt.plot(iters, out_data["trace_ll"][burnin:], alpha=0.6, label="Data-based Init")
    plt.plot(iters, out_uni["trace_ll"][burnin:], alpha=0.6, label="Uniform Init")
    plt.title("Log Likelihood Trace (Collapsed Gibbs, after burn-in)")
    plt.xlabel("Iteration"); plt.ylabel("Log Likelihood"); plt.legend(loc="upper right")
    plt.tight_layout(); plt.show()
    traces_post_data = run_post_traces(
        init_vec=good_guess, n_chains=n_chains_plot, n_iter=n_iter, burnin=burnin, thin=thin,
        seed_base=seed_data, prior_full=prior_full, counts=counts, compat_idx=compat_idx, C=C,
        alpha0=alpha0, eps=eps
    )
    traces_post_uni = run_post_traces(
        init_vec=np.ones(H)/H, n_chains=n_chains_plot, n_iter=n_iter, burnin=burnin, thin=thin,
        seed_base=seed_uni, prior_full=prior_full, counts=counts, compat_idx=compat_idx, C=C,
        alpha0=alpha0, eps=eps
    )
    traces_ll_data = run_ll_traces(
        init_vec=good_guess, n_chains=n_chains_plot, n_iter=n_iter, burnin=burnin, thin=thin,
        seed_base=seed_data, prior_full=prior_full, counts=counts, compat_idx=compat_idx, C=C,
        alpha0=alpha0, eps=eps
    )
    traces_ll_uni = run_ll_traces(
        init_vec=np.ones(H)/H, n_chains=n_chains_plot, n_iter=n_iter, burnin=burnin, thin=thin,
        seed_base=seed_uni, prior_full=prior_full, counts=counts, compat_idx=compat_idx, C=C,
        alpha0=alpha0, eps=eps
    )
    plot_trace_list(traces_post_data, "Log Posterior Trace (Data-based init)", "Log Posterior", burnin)
    plot_trace_list(traces_post_uni,  "Log Posterior Trace (Uniform init)",   "Log Posterior", burnin)
    plot_trace_list(traces_ll_data,   "Log Likelihood Trace (Data-based init)", "Log Likelihood", burnin)
    plot_trace_list(traces_ll_uni,    "Log Likelihood Trace (Uniform init)",    "Log Likelihood", burnin)
    plot_topk_posteriors(out_data["samples"], haplotypes, top_k=10)
    return out_data, out_uni

def plot_topk_errorbars_two_inits(out_data, out_uni, haplotypes, prior_full,
                                  top_k=20, jitter=0.15):
    samp_data = out_data["samples"]
    samp_uni  = out_uni ["samples"]
    mean_data = samp_data.mean(axis=0)
    mean_uni  = samp_uni .mean(axis=0)
    lo_data   = np.quantile(samp_data, 0.025, axis=0)
    hi_data   = np.quantile(samp_data, 0.975, axis=0)
    lo_uni    = np.quantile(samp_uni , 0.025, axis=0)
    hi_uni    = np.quantile(samp_uni , 0.975, axis=0)
    K = min(top_k, mean_data.size)
    idx_top = np.argsort(mean_data)[::-1][:K]
    y_pos   = np.arange(K)
    plt.figure(figsize=(8, 6))
    plt.errorbar(
        mean_data[idx_top], y_pos,
        xerr=[mean_data[idx_top] - lo_data[idx_top],
              hi_data[idx_top]  - mean_data[idx_top]],
        fmt='o', markersize=6, capsize=3,
        color='midnightblue', label='Data-based init'
    )
    plt.errorbar(
        mean_uni[idx_top], y_pos + jitter,
        xerr=[mean_uni[idx_top] - lo_uni[idx_top],
              hi_uni[idx_top]  - mean_uni[idx_top]],
        fmt='o', markersize=6, capsize=3,
        mfc='white', mec='darkorange', ecolor='darkorange',
        label='Uniform init'
    )
    plt.scatter(
        prior_full[idx_top], y_pos - jitter,
        marker='x', s=50, color='crimson',
        label='Prior (floored & renorm)'
    )
    plt.yticks(y_pos, [haplotypes[i] for i in idx_top])
    plt.gca().invert_yaxis()
    plt.xlabel("Posterior mean frequency ± 95% credible interval")
    plt.title(f"Top-{K} haplotypes: posterior estimates under two initializations")
    plt.legend(frameon=False, loc='lower right')
    plt.tight_layout(); plt.show()

def plot_top5_posterior_side_by_side(out_data, out_uni, haplotypes,
                                     good_guess, uniform_guess,
                                     top_k=5, bins=80):
    samples_data = out_data["samples"]
    samples_uni  = out_uni ["samples"]
    post_means_data = samples_data.mean(axis=0)
    K = min(top_k, post_means_data.size)
    idx_top = np.argsort(post_means_data)[::-1][:K]
    colors = plt.cm.tab10(np.arange(K) % 10)
    from matplotlib import gridspec
    fig = plt.figure(figsize=(14, 5))
    gs  = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace=0.25)
    ax0 = fig.add_subplot(gs[0])
    for j, i in enumerate(idx_top):
        ax0.hist(samples_data[:, i], bins=bins, density=True,
                 alpha=0.65, color=colors[j],
                 label=f"{haplotypes[i]} (μ={post_means_data[i]:.3f})")
        ax0.axvline(good_guess[i], color=colors[j], ls=':', lw=2.0, alpha=0.9)
    ax0.set_xlabel(r"Frequency $p_h$"); ax0.set_ylabel("Density")
    ax0.set_title("(a) Posterior distributions – Data-based init", fontsize=13)
    ax0.legend(frameon=True, fontsize=8)
    ax1 = fig.add_subplot(gs[1], sharey=ax0)
    post_means_uni = samples_uni.mean(axis=0)
    for j, i in enumerate(idx_top):
        ax1.hist(samples_uni[:, i], bins=bins, density=True,
                 alpha=0.65, color=colors[j],
                 label=f"{haplotypes[i]} (μ={post_means_uni[i]:.3f})")
        ax1.axvline(uniform_guess[i], color=colors[j], ls=':', lw=2.0, alpha=0.9)
    ax1.set_xlabel(r"Frequency $p_h$")
    ax1.set_title("(b) Posterior distributions – Uniform init", fontsize=13)
    ax1.legend(frameon=True, fontsize=8)
    x_max = max(samples_data[:, idx_top].max(), samples_uni[:, idx_top].max()) * 1.05
    for ax in (ax0, ax1):
        ax.set_xlim(0, x_max)
    fig.suptitle("Top-5 haplotypes: posterior density under two initializations",
                 fontsize=15, y=1.03)
    plt.tight_layout(); plt.show()

def r_hat(chains_1d):
    m, n = chains_1d.shape
    chain_means = chains_1d.mean(axis=1)
    grand_mean  = chain_means.mean()
    B = n * ((chain_means - grand_mean) ** 2).sum() / (m - 1)
    W = (chains_1d.var(axis=1, ddof=1)).mean()
    var_hat = ((n - 1) / n) * W + B / n
    return float(np.sqrt(var_hat / W))

def collect_sample_chains(init_vec,
                          n_chains,
                          n_iter,
                          burnin,
                          thin,
                          seed_base,
                          *,
                          prior_full, counts, compat_idx, C, alpha0, eps):
    chains = []
    for i in range(n_chains):
        out = mcmc_collapsed_gibbs(
            init_p        = init_vec,
            prior_full    = prior_full,
            counts        = counts,
            compat_idx    = compat_idx,
            C             = C,
            alpha0        = alpha0,
            n_iter        = n_iter,
            burnin        = burnin,
            thin          = thin,
            eps           = eps,
            record_traces = False,
            rng_seed      = seed_base + i
        )
        chains.append(out["samples"])
    return chains

def compute_rhat_table_topk(out_data, out_uni, model, alpha0,
                            top_k = 20,
                            n_chains = 5,
                            n_iter = 100_000,
                            burnin = 20_000,
                            thin = 1,
                            eps = 1e-4,
                            seed_data = 51510,
                            seed_uni = 62620):
    haplotypes   = model["haplotypes"]
    H            = model["H"]
    prior_full   = model["prior_full"]
    good_guess   = model["good_guess"]
    counts       = model["pattern_counts"]
    compat_idx   = model["compat_idx"]
    C            = model["C"]
    mean_data = out_data["samples"].mean(axis=0)
    K = min(top_k, H)
    idx_topK = np.argsort(mean_data)[::-1][:K]
    chains_data = collect_sample_chains(
        init_vec   = good_guess,
        n_chains   = n_chains,
        n_iter     = n_iter,
        burnin     = burnin,
        thin       = thin,
        seed_base  = seed_data,
        prior_full = prior_full,
        counts     = counts,
        compat_idx = compat_idx,
        C          = C,
        alpha0     = alpha0,
        eps        = eps
    )
    chains_uni = collect_sample_chains(
        init_vec   = np.ones(H)/H,
        n_chains   = n_chains,
        n_iter     = n_iter,
        burnin     = burnin,
        thin       = thin,
        seed_base  = seed_uni,
        prior_full = prior_full,
        counts     = counts,
        compat_idx = compat_idx,
        C          = C,
        alpha0     = alpha0,
        eps        = eps
    )
    rows = []
    for idx in idx_topK:
        mat_data = np.stack([c[:, idx] for c in chains_data])
        mat_uni  = np.stack([c[:, idx] for c in chains_uni ])
        rows.append({
            "Haplotype":             haplotypes[idx],
            "R̂ (Data-based init)":   r_hat(mat_data),
            "R̂ (Uniform init)":     r_hat(mat_uni),
        })
    df = pd.DataFrame(rows)
    df = df.set_index("Haplotype").loc[[haplotypes[i] for i in idx_topK]].reset_index()
    return df

def _summary_mean_var_for_chains(chains, idx, haplotypes, ddof = 1):
    rows = []
    for j in idx:
        all_samples = np.concatenate([c[:, j] for c in chains], axis=0)
        rows.append({
            "Haplotype":  haplotypes[j],
            "Post. mean": float(all_samples.mean()),
            "Post. var":  float(all_samples.var(ddof=ddof)),
        })
    return pd.DataFrame(rows)

def compute_mean_var_table_topk(out_data, model, alpha0,
                                top_k = 20,
                                n_chains = 5,
                                n_iter = 100_000,
                                burnin = 20_000,
                                thin = 1,
                                eps = 1e-4,
                                seed_data = 97531,
                                seed_uni = 86420,
                                chains_data = None,
                                chains_uni  = None,
                                idx_topK    = None):
    haplotypes   = model["haplotypes"]
    H            = model["H"]
    prior_full   = model["prior_full"]
    good_guess   = model["good_guess"]
    counts       = model["pattern_counts"]
    compat_idx   = model["compat_idx"]
    C            = model["C"]
    if idx_topK is None:
        mean_data = out_data["samples"].mean(axis=0)
        K = min(top_k, H)
        idx_topK = np.argsort(mean_data)[::-1][:K]
    else:
        idx_topK = np.asarray(idx_topK, dtype=int)
    if chains_data is None:
        chains_data = collect_sample_chains(
            init_vec   = good_guess,
            n_chains   = n_chains,
            n_iter     = n_iter,
            burnin     = burnin,
            thin       = thin,
            seed_base  = seed_data,
            prior_full = prior_full,
            counts     = counts,
            compat_idx = compat_idx,
            C          = C,
            alpha0     = alpha0,
            eps        = eps
        )
    if chains_uni is None:
        chains_uni = collect_sample_chains(
            init_vec   = np.ones(H)/H,
            n_chains   = n_chains,
            n_iter     = n_iter,
            burnin     = burnin,
            thin       = thin,
            seed_base  = seed_uni,
            prior_full = prior_full,
            counts     = counts,
            compat_idx = compat_idx,
            C          = C,
            alpha0     = alpha0,
            eps        = eps
        )
    summ_emp = _summary_mean_var_for_chains(chains_data, idx_topK, haplotypes, ddof=1) \
                 .rename(columns={"Post. mean": "Mean (Data-based)", "Post. var": "Var (Data-based)"})
    summ_uni = _summary_mean_var_for_chains(chains_uni,  idx_topK, haplotypes, ddof=1) \
                 .rename(columns={"Post. mean": "Mean (Uniform)", "Post. var": "Var (Uniform)"})
    summary_df = summ_emp.merge(summ_uni, on="Haplotype", how="inner") \
                         .set_index("Haplotype") \
                         .loc[[haplotypes[i] for i in idx_topK]] \
                         .reset_index()
    return summary_df

def _select_site_indices(model, sites=None, max_sites=12, order="pos_asc"):
    positions = np.array(model["positions"])

    if sites is None:
        k_vec = model["site_alt_counts"]
        idx = np.argsort(k_vec)[::-1][:min(max_sites, len(k_vec))]
    else:
        pos2idx = {pos: i for i, pos in enumerate(positions)}
        idx = [pos2idx[p] for p in sites if p in pos2idx]
        idx = np.array(idx[:min(max_sites, len(idx))], dtype=int)

    if idx.size > 0:
        if order in ("pos_asc", "pos"):
            idx = idx[np.argsort(positions[idx])]
        elif order == "pos_desc":
            idx = idx[np.argsort(positions[idx])[::-1]]

    return np.asarray(idx, dtype=int)

def _alt_indicator_rows(haplotypes, positions, alt_base_map, idx_sel):
    H = len(haplotypes)
    A = np.zeros((len(idx_sel), H), dtype=float)
    for r, j in enumerate(idx_sel):
        pos = positions[j]
        alt = alt_base_map[pos]
        A[r, :] = np.fromiter((1.0 if h[j] == alt else 0.0 for h in haplotypes), float, count=H)
    return A

def _binom_ci_clopper_pearson(k, n, alpha=0.05):
    if n <= 0:
        return np.nan, np.nan
    lo = 0.0 if k == 0 else beta.ppf(alpha/2, k, n-k+1)
    hi = 1.0 if k == n else beta.ppf(1-alpha/2, k+1, n-k)
    return float(lo), float(hi)

def run_chain_alt_traces_full_iters(init_p, model, alpha0, idx_sel, n_iter=100_000, eps=1e-4, rng_seed=None):
    if rng_seed is not None:
        np.random.seed(rng_seed)
    H             = model["H"]
    haplotypes    = model["haplotypes"]
    compat_idx    = model["compat_idx"]
    counts        = model["pattern_counts"].astype(int)
    prior_full    = model["prior_full"]
    positions     = model["positions"]
    alt_base_map  = model["alt_base_per_site"]

    Z = model["Z"]               

    A_sel = _alt_indicator_rows(haplotypes, positions, alt_base_map, idx_sel)
    p = _proj_floor_and_renorm(np.asarray(init_p, dtype=float), eps)
    alpha_prior = alpha0 * prior_full
    J_sel = len(idx_sel)
    f_tr  = np.zeros((n_iter, J_sel), dtype=float)
    suff = np.zeros(H, dtype=np.int64)

    for it in range(n_iter):
        suff[:] = 0
        for r, k in enumerate(counts):
            if k <= 0:
                continue
            ids = compat_idx[r]
            if ids.size == 0:
                continue
            probs = p[ids] / Z[ids]
            s = probs.sum()
            probs = np.ones_like(probs)/probs.size if (s <= 0 or not np.isfinite(s)) else (probs/s)
            y = np.random.multinomial(int(k), probs)
            np.add.at(suff, ids, y)

        alpha_post = alpha_prior + suff + eps
        p = np.random.dirichlet(alpha_post)
        p = _proj_floor_and_renorm(p, eps)
        f_tr[it, :] = A_sel.dot(p)

    return f_tr

def _moving_avg(y, w):
    if not w or w <= 1 or w >= y.shape[0]:
        return None
    ker = np.ones(int(w))/int(w)
    return np.convolve(y, ker, mode="valid")

def plot_alt_traces_full_iterations(
    model, alpha0,
    n_iter=100_000, burnin=20_000, eps=1e-4,
    sites=None, max_sites=12,
    seed_emp=24680, seed_uni=13579,
    ci_color="#BEBEBE",
    color_emp="#7FA6C9",
    color_uni="#9B8EC1",
    line_lw=1.0,
    smooth_window=0,
    zoom_to_ci=True,
    ylim_pad=0.03,
    min_height=0.08,
    row_height=2.2,
    fig_width=10
):
    idx_sel = _select_site_indices(model, sites=sites, max_sites=max_sites, order="pos_asc")
    if len(idx_sel) == 0:
        print("[INFO] No sites selected."); return
    H          = model["H"]
    positions  = np.array(model["positions"])
    k_vec      = model["site_alt_counts"].astype(int)
    n_vec      = model["site_depths"].astype(int)
    col_map    = {j: i for i, j in enumerate(idx_sel)}
    f_emp = run_chain_alt_traces_full_iters(
        init_p=model["good_guess"], model=model, alpha0=alpha0,
        idx_sel=idx_sel, n_iter=n_iter, eps=eps, rng_seed=seed_emp
    )
    f_uni = run_chain_alt_traces_full_iters(
        init_p=np.ones(H)/H, model=model, alpha0=alpha0,
        idx_sel=idx_sel, n_iter=n_iter, eps=eps, rng_seed=seed_uni
    )
    def _panel(f_mat, title_prefix, trace_color):
        n_rows = len(idx_sel)
        fig, axes = plt.subplots(
            n_rows, 1, figsize=(fig_width, max(2.0, row_height*n_rows)), squeeze=False, sharex=True
        )
        axes = axes.ravel()
        x = np.arange(n_iter)
        for ax, j in zip(axes, idx_sel):
            pos = positions[j]
            k, n = int(k_vec[j]), int(n_vec[j])
            obs  = (k/n) if n > 0 else np.nan
            lo, hi = _binom_ci_clopper_pearson(k, n, alpha=0.05)
            c     = col_map[j]
            ax.fill_between([0, n_iter-1], lo, hi, color=ci_color, alpha=0.50, linewidth=0)
            if np.isfinite(obs):
                ax.axhline(obs, ls="--", lw=1.0, color="black", alpha=0.9)
            ax.axvline(burnin, ls="--", lw=1.0, color="gray", alpha=0.9)
            ax.plot(x, f_mat[:, c], lw=line_lw, color=trace_color, alpha=0.95)
            sm = _moving_avg(f_mat[:, c], smooth_window)
            if sm is not None:
                xs = np.arange(sm.size) + (smooth_window - 1) / 2.0
                ax.plot(xs, sm, lw=1.2, color=trace_color, alpha=0.9)
            if zoom_to_ci and np.isfinite(lo) and np.isfinite(hi):
                qlo = np.quantile(f_mat[:, c], 0.01)
                qhi = np.quantile(f_mat[:, c], 0.99)
                y0 = max(0.0, min(lo, qlo) - ylim_pad)
                y1 = min(1.0, max(hi, qhi) + ylim_pad)
                if (y1 - y0) < min_height:
                    mid = 0.5*(y0+y1)
                    y0, y1 = max(0.0, mid-min_height/2), min(1.0, mid+min_height/2)
                ax.set_ylim(y0, y1)
            else:
                ax.set_ylim(0, 1)
            ax.set_xlim(0, n_iter-1)
            ax.set_title(f"{title_prefix}: pos {pos}", fontsize=11)
            ax.set_xlabel("Iteration"); ax.set_ylabel("ALT frequency")
            ax.grid(alpha=0.15)
        plt.tight_layout(rect=[0, 0.05, 1, 1]); plt.show()
    _panel(f_emp, "Data-based init", color_emp)
    _panel(f_uni, "Uniform init",  color_uni)

def plot_top20_haplotype_traces(
    out_data, out_uni, model,
    good_guess=None, uniform_guess=None,
    burnin=20_000,
    nrow=4, ncol=5,
    thin=1,
    zoom=True, pad=0.05, min_span=0.08, q_lo=0.01, q_hi=0.99,
    cmap="tab20",
    lw=0.9, alpha=0.85,
    init_ls=":", init_lw=1.2,
    title_emp="Trace of Top 20 haplotype frequencies – Data-based init",
    title_uni="Trace of Top 20 haplotype frequencies – Uniform init"
):
    haplotypes = model["haplotypes"]
    H = len(haplotypes)
    if good_guess is None:
        good_guess = model["good_guess"]
    if uniform_guess is None:
        uniform_guess = model["uniform_guess"]
    post_means_emp = out_data["samples"].mean(axis=0)
    idx_top20 = np.argsort(post_means_emp)[::-1][:min(20, H)]
    colors20 = plt.cm.get_cmap(cmap, 20)(np.arange(20))
    def _auto_ylim(y):
        lo = np.quantile(y, q_lo); hi = np.quantile(y, q_hi)
        span = hi - lo
        if span < min_span:
            mid = 0.5*(lo+hi); lo, hi = mid - min_span/2, mid + min_span/2
        lo = max(0.0, lo - pad); hi = min(1.0, hi + pad)
        return float(lo), float(hi)
    def _trace_grid(samples, init_vec, title):
        iters = np.arange(samples.shape[0])[::thin]
        fig, axes = plt.subplots(nrow, ncol, figsize=(18, 11), sharex=True, sharey=False)
        axes = axes.ravel()
        for k, hidx in enumerate(idx_top20):
            ax = axes[k]
            y = samples[::thin, hidx]
            ax.plot(iters, y, color=colors20[k], lw=lw, alpha=alpha)
            ax.axhline(init_vec[hidx], ls=init_ls, lw=init_lw, color=colors20[k], alpha=0.95)
            ax.axvline(min(burnin, iters[-1]), ls="--", lw=1.2, color="grey", alpha=0.9)
            ax.set_title(haplotypes[hidx], fontsize=12); ax.tick_params(labelsize=10)
            if zoom: ax.set_ylim(*_auto_ylim(y))
            else:    ax.set_ylim(0, 1)
        for ax in axes[len(idx_top20):]: ax.axis("off")
        for ax in axes[:len(idx_top20)]: ax.set_xlim(iters[0], iters[-1])
        fig.suptitle(title, fontsize=20, y=0.97)
        fig.text(0.5, 0.04, "Iteration", ha="center", fontsize=13)
        fig.text(0.06, 0.5, r"Frequency $p_h$", va="center", rotation="vertical", fontsize=13)
        plt.tight_layout(rect=[0.08, 0.06, 1, 0.95]); plt.show()
    _trace_grid(out_data["samples"], good_guess, title_emp)
    _trace_grid(out_uni["samples"],  uniform_guess, title_uni)


#%% [CFG] Paths & seeds
np.random.seed(42)

VAR_TABLE   = "variants_final_table.tsv"
PRIOR_TABLE = "consensus_all_haplotypes.tsv"
ENC_TABLE   = "encoding_counts.tsv"

#%% [STEP 1] Prepare model
model = prepare_model(VAR_TABLE, PRIOR_TABLE, ENC_TABLE, eps_floor=1e-4, rng_seed=42)
print("H =", model["H"], "| patterns =", len(model["read_patterns"]), "| reads =", int(model["pattern_counts"].sum()))

#%% [STEP 2] Choose alpha0 (prior carries ω of total evidence)

model["prior_full"] = model["prior_empirical"]      # informative prior
#model["prior_full"] = model["prior_uniform"]      # uniform prior

H        = model["H"]
w_full   = model["prior_full"]
C        = model["C"]
counts   = model["pattern_counts"]
omega = 0.10
row_nonempty = (np.array(C.sum(axis=1)).ravel() > 0)
N_eff = int(np.asarray(counts)[row_nonempty].sum())

alpha0 = (omega / (1.0 - omega)) * N_eff
print(f"[alpha0-omega] ω={omega:.3f}, N_eff={N_eff:,} -> alpha0={alpha0:.1f}")
print(f"[prior type] {'uniform' if np.allclose(w_full, np.ones(H)/H) else 'empirical'}")

alpha0 = H
EPS = eps_from_H(H)
print(f"[INFO] adaptive eps = {EPS:.3e}  (H={H})")

#%% [STEP 3] Fit two inits + base plots
out_data, out_uni = fit_two_inits_and_plots(
    model, alpha0, n_iter=100000, burnin=20000, thin=1, eps=EPS,
    n_chains_plot=5, seed_data=97531, seed_uni=86420
)

#%% [STEP 4] Extra LL trace overlays
ll_traces_data = run_ll_traces(model["good_guess"], 5, 100_000, 20_000, 1, 97531,
                               model["prior_full"], model["pattern_counts"], model["compat_idx"], model["C"], alpha0, EPS)
ll_traces_uni  = run_ll_traces(np.ones(model["H"])/model["H"], 5, 100_000, 20_000, 1, 86420,
                               model["prior_full"], model["pattern_counts"], model["compat_idx"], model["C"], alpha0, EPS)
plot_trace_list(ll_traces_data, "Log Likelihood Trace (Data-based init)", "Log Likelihood", burnin=20_000, absolute_index=True)
plot_trace_list(ll_traces_uni,  "Log Likelihood Trace (Uniform init)",    "Log Likelihood", burnin=20_000, absolute_index=True)

#%% [STEP 5] Top-k errorbars
plot_topk_errorbars_two_inits(out_data, out_uni, model["haplotypes"], model["prior_full"], top_k=20)

#%% [STEP 6] Side-by-side posterior histograms (Top-5)
plot_top5_posterior_side_by_side(
    out_data, out_uni,
    haplotypes=model["haplotypes"],
    good_guess=model["good_guess"],
    uniform_guess=model["uniform_guess"],
    top_k=5,
    bins=80
)

#%% [STEP 7] R-hat table (Top-20 by data-based posterior mean)
rhat_df = compute_rhat_table_topk(
    out_data=out_data,
    out_uni=out_uni,
    model=model,
    alpha0=alpha0,
    top_k=20,
    n_chains=5,
    n_iter=100_000,
    burnin=20_000,
    thin=1,
    eps=EPS,
    seed_data=97531,
    seed_uni=86420
)
pd.set_option("display.float_format", "{:.3f}".format)
print("\nGelman–Rubin diagnostics (top-20 haplotypes)\n")
print(rhat_df.to_string(index=False))

#%% [STEP 8] Posterior mean/variance table (Top-20; reuse same order if needed)
summary_df = compute_mean_var_table_topk(
    out_data=out_data,
    model=model,
    alpha0=alpha0,
    top_k=20,
    n_chains=5,
    n_iter=100_000,
    burnin=20_000,
    thin=1,
    eps=EPS,
    seed_data=97531,
    seed_uni=86420
)
pd.set_option("display.float_format", "{:.4e}".format)
print("\nPosterior mean & variance (top-20 haplotypes)\n")
print(summary_df.to_string(index=False))

#%% [STEP 9] ALT frequency trace panels per site
plot_alt_traces_full_iterations(
    model, alpha0,
    n_iter=100_000, burnin=20_000,
    max_sites=12,
    ci_color="#6E6E6E",
    color_emp="#7FA6C9",
    color_uni="#9B8EC1",
    zoom_to_ci=True,
    ylim_pad=0.03,
    min_height=0.08,
    row_height=1.5, fig_width=14,
    smooth_window=0,
    eps=EPS
)
#%% [STEP 10] Top-20 haplotype trace panels
plot_top20_haplotype_traces(
    out_data, out_uni, model,
    good_guess=model["good_guess"],
    uniform_guess=model["uniform_guess"],
    burnin=20_000,
    nrow=4, ncol=5,
    thin=1,
    zoom=True, pad=0.05, min_span=0.08, q_lo=0.01, q_hi=0.99,
    cmap="tab20",
    lw=0.9, alpha=0.85
)
