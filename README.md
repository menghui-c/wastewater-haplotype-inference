# 🧬 Wastewater Haplotype Inference
**Bayesian haplotype inference using short-read viral sequencing data from wastewater samples**

![illustrative diagram](images/simple_illu.png)

---

## 📖 Overview

This project implements a Bayesian framework to infer viral haplotype frequencies from deep-sequenced SARS-CoV-2 wastewater data. We first run two workflows to precompute the fixed inputs for inference: (1) read based SNV and encoding pipeline and (2) consensus panel prior pipeline. Then a collapsed Gibbs sampler is used to estimate the haplotype frequencies.

---

## 🧠 Model Description

### Data and notation

* A Deep sequenced wastewater sample with $J$ biallelic sites is indexed by $j=1,\dots,J$.
* Each site has REF/ALT alleles; a haplotype $h \in \lbrace 1,\ldots,H \rbrace$ is a length $J$ allele vector; with biallelic sites, $H=2^J$.
* Short reads are grouped into **read patterns** $r = 1, 2, \ldots, R$; each pattern is a length $J$ string over `{0,1,-}`(`0` = REF, `1` = ALT,  = not covered).
* Let $k_r$ be the number of reads with pattern $r$; $\sum_{r=1}^R k_r = N$. Write $\mathbf{k}=(k_1,\dots,k_R)^\top$.

---

### Compatibility and pattern probabilities

Let $C \in \lbrace 0,1 \rbrace^{R\times H}$ be the compatibility matrix: 

$$
C_{rh} =
\begin{cases}
1, & \text{if pattern } r \text{ is compatible with haplotype } h, \\
0, & \text{otherwise.}
\end{cases}
$$

Define

$$
A_{rh}=\frac{C_{rh}}{Z_h},\qquad
Z_h=\sum_{r=1}^R C_{rh},\qquad
\sum_{r=1}^R A_{rh}=1
$$

for each $h$.

Let $\mathbf p=(p_1,\dots,p_H)^\top$ denote haplotype frequencies on the $H$-simplex.

The pattern distribution is

$$
\mathbf q = A \mathbf p
$$

$$
q_r = \sum_{h=1}^H A_{rh} p_h
$$

---

### Read-pattern likelihood

Given $\mathbf p$, the pattern counts follow a multinomial model:

$$
\mathbf K|\mathbf p \sim \mathrm{Multinomial}(N;\mathbf q),
$$

and the log-likelihood (without constants) is:

$$
\ell(\mathbf p;\mathbf k) = \sum_{r=1}^R k_r \log q_r = \mathbf k^\top \log(A\mathbf p).
$$

---

### Dirichlet prior on haplotype frequencies

We place a Dirichlet prior on $\mathbf p$:

$$
\mathbf p \sim \mathrm{Dirichlet}(\boldsymbol\alpha),
\quad
\boldsymbol\alpha = \alpha_0\\mathbf w.
$$

- **Informative prior:** $\mathbf w$ is derived from normalized GISAID haplotype frequencies.  
  To avoid zeros, a small floor $\varepsilon$ is applied before renormalization. 

- **Non-informative prior:** $\mathbf w = \tfrac{1}{H}\mathbf 1$ and $\alpha_0 = H$.

---

### Collapsed Gibbs sampling via data augmentation

We introduce allocation counts $y_{rh}$: number of reads with pattern $r$ that are attributed to haplotype $h$.  
They satisfy $\sum_h y_{rh}=k_r$ for each r.
Our model alternates between these two steps:

1. **Allocation step:**
   
$$
(y_{r1},\dots,y_{rH})|\mathbf p
\sim \mathrm{Multinomial}\\big(k_r;\ w_{r1},\dots,w_{rH}\big),
\quad
w_{rh} = \frac{A_{rh}p_h}{q_r}.
$$

3. **Update step:**
   
$$
\mathbf p|\mathbf y \sim 
\mathrm{Dirichlet}(\alpha_1+n_1,\dots,\alpha_H+n_H),
\quad n_h = \sum_r y_{rh}.
$$

This results in a **collapsed Gibbs sampler** with Dirichlet–multinomial conjugacy.

---

### Initialization, constraints, reproducibility

There are two parallel chains:
- **Data-based initialization:** product over sites: $\hat f_j$ for ALT allele frequency, $1-\hat f_j$ for REF; 
- **Uniform initialization:** $\mathbf p^{(0)} = \tfrac{1}{H}\mathbf 1$.  

At each iteration, we set:

$$
p_h \leftarrow \max(p_h,\varepsilon), \quad 
\text{renormalize s.t. } \sum_h p_h=1.
$$

The typical configurations are:
$n_{\text{iter}}=10^5$, burn-in $B=2\times10^4$, no thinning, fixed seeds.

---

### Posterior diagnostics

We compute:
- Log-posterior traces
- Gelman--Rubin statistic $\hat R$ convergence diagnostic
- Posterior mean $\pm$ 95% credible intervals for top-\(K\) haplotypes (ranked by posterior mean frequency)
- ALT allele posterior frequency traces

---

## ⚙️ Code structure

```text
pipelines/
├─ read_based_pipeline.sh            # End-to-end: Read-based SNV and Encoding pipeline
├─ consensus_based_pipeline.sh       # End-to-end: Consensus Panel Prior pipeline

read_based:                         # used by read_based_pipeline.sh
├─ get_snvs.py                       # call SNVs (LoFreq)
├─ filter_lofreq_vcf.py              # filter by VAF/RAF/DP
├─ extract_snv_table.py              # build SNVs summary table
├─ read_encoding.py                  # encode read patterns {0,1,-}
├── compute_haplotype_freqs_allele.py       # compute data-based initial haplotype frequencies

consensus_panel:                     # used by consensus_based_pipeline.sh
├─ consensus_summarize_variant_alleles.py  # summarize per-site allele counts from panel
├─ consensus_haplotype_freq.py              # compute panel-based haplotype frequencies as priors

inference:                         
├─ gibbs_sampler.py                  # whole collapsed Gibbs sampler code for haplotype inference
```

---

## 🧩 Workflow

![Pipeline diagram](images/workflow.jpg)

**Step 1.** Reads alignment → call SNVs → summarize read patterns  
**Step 2.** Consensus genomes alignment → compute panel imformative prior  
**Step 3.** Run Gibbs sampler to infer posterior haplotype frequencies  
**Step 4.** Summarize diagnostics and credible intervals  

---

## 🚀 Run the pipeline

```bash
# Step 1. Run read-based pipeline

Paired-end (recommended)
./read_based_pipeline.sh \
  -r reference.fasta \
  -1 short_read_sequencing_R1.fastq \
  -2 short_read_sequencing_R2.fastq \
  --primer-bed primers.bed \
  -t 8 \
  --min-cov 30 \
  --min-alt-freq 0.10 \
  --min-ref-freq 0.10 \
  --outdir results/read_based_output

Single-file interleaved
./read_based_pipeline.sh \
  -r reference.fasta \
  -q short_read_sequencing.fastq \
  --interleaved \
  --primer-bed primers.bed \
  -t 8 \
  --outdir results/read_based_output

Single-end
./read_based_pipeline.sh \
  -r reference.fasta \
  -q short_read_sequencing.fastq \
  -t 8 \
  --outdir results/read_based_output

# Step 2. Run consensus-panel pipeline
chmod +x consensus_based_pipeline.sh
./consensus_based_pipeline.sh \
  -r reference.fasta \
  -s samples.fasta \
  --readdir results/read_based_output \
  [-t threads] \
  [-o outdir]

# Step 3. Perform Gibbs sampling
```

---

## 📊 Dependencies
	•	Python ≥ 3.11
	•	NumPy, SciPy, Pandas, Matplotlib
	•	BWA, MAFFT, LoFreq
	•	Bash ≥ 4.0
	
---

## 🧑‍💻 Implementation details

All the calculations are vectorized using sparse matrix operations in NumPy and SciPy. The compatibility matrix (C) is stored in CSR format, and the log-likelihood terms are efficiently calculated through sparse-dense multiplication.

---

## 📬 Contact

For questions or collaboration:
Prof. Katia Koelle
📧 katia.koelle@emory.edu


