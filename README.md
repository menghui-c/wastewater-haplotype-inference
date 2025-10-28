# 🧬 Wastewater Haplotype Inference
**Bayesian haplotype inference using short-read viral sequencing data from wastewater samples**

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
- **Data-based initialization:** uses empirical sitewise ALT frequencies $\hat f_j$.  
- **Uniform initialization:** $\mathbf p^{(0)} = \tfrac{1}{H}\mathbf 1$.  

Each iteration enforces:

$$
p_h \leftarrow \max(p_h,\varepsilon), \quad 
\text{renormalize s.t. } \sum_h p_h=1.
$$

Typical configuration:  
$n_{\text{iter}}=10^5$, burn-in $B=2\times10^4$, no thinning, fixed seeds.

---

### Posterior diagnostics

We compute:
- Log-likelihood / log-posterior traces
- $\hat R$ convergence diagnostic
- Posterior mean $\pm$ 95% credible intervals for top-$K$ haplotypes
- ALT allele posterior projections
  $$
  f_j^{(t)} = \sum_h a_{hj} p_h^{(t)}
  $$

---

## ⚙️ Pipeline structure
├── read_based_pipeline.sh         # End-to-end SNV extraction and pattern encoding
├── consensus_based_pipeline.sh    # Consensus panel prior generation
├── get_snvs.py                    # SNV extraction via LoFreq
├── filter_lofreq_vcf.py           # Variant filtering (VAF, RAF, DP)
├── extract_snv_table.py           # Generate variant summary table
├── read_encoding.py               # Encode read patterns {0,1,-}
├── consensus_summarize_variant_alleles.py  # Summarize consensus allele frequencies
├── consensus_haplotype_freq.py    # Compute panel-derived haplotype priors
├── compute_haplotype_freqs_allele.py       # Aggregate allele-level frequency estimates
├── gibbs_sampler.py               # Collapsed Gibbs sampler for haplotype inference

---

## 🧩 Workflow

![Pipeline diagram](images/pipeline_workflow.png)

**Step 1.** Map reads with BWA → call SNVs with LoFreq → summarize read patterns  
**Step 2.** Align consensus genomes with MAFFT → compute panel prior  
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
```bash
chmod +x consensus_based_pipeline.sh
./consensus_based_pipeline.sh \
  -r reference.fasta \
  -s samples.fasta \
  --readdir path/to/read_based_output \
  [-t threads] \
  [-o outdir]

# Step 3. Perform Gibbs sampling

## 📊 Dependencies
	•	Python ≥ 3.11
	•	NumPy, SciPy, Pandas, Matplotlib
	•	BWA, MAFFT, LoFreq
	•	Bash ≥ 4.0

## 🧑‍💻 Implementation details

All computations are vectorized with NumPy and sparse matrix operations in SciPy.
The compatibility matrix (C) is stored in CSR format, and the log-likelihood term is efficiently computed via sparse–dense multiplication.

## 📚 References
	•	Bernardo, J. M., & Smith, A. F. M. (1994). Bayesian Theory.
	•	Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior distributions by data augmentation. JASA.
	•	Liu, J. S. (1994). Collapsed Gibbs sampling and other variance-reduction techniques. JASA.
	•	Morita et al. (2008). Determining the effective sample size of a Dirichlet prior.
	•	Elbe & Buckland-Merrett (2017). GISAID: Global initiative on sharing all influenza data.
	•	Katoh & Standley (2013). MAFFT multiple sequence alignment software.
	•	Li & Durbin (2009). BWA: Burrows–Wheeler Aligner.
	•	Wilm et al. (2012). LoFreq: sensitive variant calling for low-frequency variants.

## 📬 Contact

For questions or collaboration:
Menghui (Mona) Chen
📧 menghui.chen@emory.edu


