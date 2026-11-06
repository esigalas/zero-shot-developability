# 🧬 Universal Zero-Shot Protein & Antibody Developability Pipeline

An end-to-end, two-phase computational protein engineering platform for designing, scoring, and triaging aggregation-resistant, highly developable antibodies and proteins **without requiring prior wet-lab training data on the target scaffold**.

The pipeline strictly decouples **variant generation (Phase 1)** from **biophysical and evolutionary evaluation (Phase 2)**, combining 3D inverse folding, structure-based ΔΔG scanning, and protein/antibody language models with non-compensatory biophysical hard gates and Jaccard diversity pruning.

---

## ⚙️ Installation & Environment Setup

Follow these steps to set up the dual conda environments, clone external physics engines, and configure your local workspace paths.

### 1. Clone the Main Repository
First, clone this repository to your local system and change into the project directory:
```bash
git clone git@github.com:esigalas/zero-shot-developability.git
cd zero-shot-developability
```

### 2. Set Up Dual Conda Environments
To optimize memory footprints, library compatibility, and isolated execution, the pipeline splits model dependencies into two distinct conda environments:

*   **Environment 1: The Main Pipeline Environment (`environment_main.yml`)**
    This acts as the principal environment where you launch the Streamlit app and execute all structural scoring scripts (except for ESM-C DMS and scoring).
    ```bash
    # Create the environment from the config file
    conda env create -n zero-shot-env -f environment_main.yml
    
    # Activate the principal environment
    conda activate zero-shot-env

    # Install antifold python module separately with no dependencies
    pip install antifold==0.2.3 --no-deps

    # 1. Alias 'filter_backbone' in biotite for ESM-IF1 and AntiFold compatibility
    python -c "import biotite.structure as bs, pathlib; p = pathlib.Path(bs.__file__); txt = p.read_text(); p.write_text(txt + '\n# Backwards compatibility for AntiFold and ESM-IF1\nfilter_backbone = filter_peptide_backbone\n') if 'filter_backbone = filter_peptide_backbone' not in txt else None; print('✅ Biotite patched!')"

    # 2. Pre-download AntiFold weights to bypass the auto-download bug in antifold >= 1.0.1
    python -c "import os, urllib.request, antifold; d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(antifold.__file__))), 'models'); os.makedirs(d, exist_ok=True); urllib.request.urlretrieve('[https://opig.stats.ox.ac.uk/data/downloads/AntiFold/models/model.pt](https://opig.stats.ox.ac.uk/data/downloads/AntiFold/models/model.pt)', os.path.join(d, 'model.pt')); print('✅ AntiFold weights ready!')"
    ```

*   **Environment 2: Isolated ESM-C Environment (`environment_esmc.yml`)**
    This environment is dedicated to running ESM-C predictions (which operate in isolated subprocess calls under `run_esmc_dms.py` and `run_esmc_scoring.py` to prevent library and memory conflicts with PyRosetta).
    ```bash
    # Create the isolated ESM-C environment
    conda env create -n esmc_env -f environment_esmc.yml
    ```

---

### 3. Clone External Dependencies
To run Phase 1 generative models and Phase 2 scoring metrics, you must clone three external third-party repositories into your workspace directory. Run these clone commands in your workspace:

1.  **ThermoMPNN** (Phase 1 generative ΔΔG scanning model):
    ```bash
    git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git
    ```
2.  **ProteinMPNN** (Phase 1 GNN inverse-folding model):
    ```bash
    git clone https://github.com/dauparas/ProteinMPNN.git
    ```
3.  **JanusDDG** (Phase 2 sequence-level thermodynamic scoring model):
    ```bash
    git clone git@github.com:compbiomed-unito/JanusDDG.git
    ```

---

### 4. Configure Local Workspace Paths in `pipeline.py`
Open `pipeline.py` and configure the static variables at the top of the file to point to your local installation directories and conda environment names:

```python
# =====================================================================
# GLOBAL PATHS & ENVIRONMENT CONFIGURATION (Edit These For Your Setup)
# =====================================================================
PROTEIN_MPNN_SCRIPT = "/home/esigalas/workspace/ProteinMPNN/protein_mpnn_run.py"
THERMO_MPNN_SCRIPT  = "/home/esigalas/workspace/ThermoMPNN/analysis/custom_inference.py"
JANUS_SRC_PATH      = "/home/esigalas/workspace/JanusDDG/src"
JANUS_MODELS_DIR    = "/home/esigalas/workspace/JanusDDG/models"
ESMC_CONDA_ENV      = "esmc_env"
```

*Ensure these absolute paths match where you cloned the repositories in Step 3.*

---

### 5. Launch the Web App
Once your environments are built and your paths are configured, activate your main pipeline environment and start the Streamlit interface:

```bash
# Make sure you are in your main environment
conda activate zero-shot-env

# Start the web app
streamlit run app.py
```

---

## ✨ Key Features & Supported Formats

* **Universal Multi-Chain Format Support:**
  * `scFv (Single-chain Fv)`
  * `Standard Monoclonal Antibody (IgG / Fab)`
  * `scFv-Fc Dimer`
  * `VHH-Fc / Single-Domain Dimer` (Nanobodies)
  * `Bispecific Fab-scFv` & `Bispecific IgG (Asymmetric Heavy/Light Chains)` (up to 4–6 chains)
  * `Monomeric Enzyme / Non-Antibody` & `Custom Multi-Chain Complex`
* **6-Model Orthogonal Generative Ensemble (Phase 1):** Pools single-point and multi-mutant combinations across `AntiFold`, `ProteinMPNN (Diverse)`, `ProteinMPNN (Soluble)`, `ThermoMPNN`, `ESM-C (600M)`, and `AbLang2-Paired`.
* **Unified IMGT & Disulfide Masking:** Automatically segments `VH`, `VL`, and `VHH` chains via `abnumber` (IMGT scheme) to lock user-selected CDR/FR regions, freezes native cysteines to protect disulfide bonds, and enforces custom 1-indexed residue masks across all generators.
* **16-Metric Biophysical & Evolutionary Triage (Phase 2):** Evaluates 3D packing (`PyRosetta`), 1D folding stability (`JanusDDG`), `VH–VL` interface affinity (`BINDPRED`), electrostatics (`pI`, `Net Charge`, `VH–VL Charge Asymmetry`), intrinsic solubility (`CamSol`, `GRAVY`), and 5 evolutionary/antibody language models (`ESM-2 650M`, `ESM-C 600M`, `AntiBERTy`, `AbLang2-Paired`, `ProtBERT`).
* **Agnostic Gated + Invariant Scoring:** Prevents 1D language model "consensus inflation" in multi-mutant designs using one-sided downside folding hinges and three non-compensatory biophysical hard gates.
* **Mutational & Positional Jaccard Diversity Triage:** Eliminates redundant "combo clones" in the final wet-lab synthesis shortlist while diagnosing cross-generator structural hotspot convergence.

---

## 🏗️ Pipeline Architecture

```text
Input FASTA (N Chains)
   │
   ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ PHASE 1: STRUCTURE PREPARATION & ORTHOGONAL GENERATION (pipeline.py)    │
│  1. Unified Masking: IMGT CDR/FR parsing (abnumber) + Cys lock          │
│  2. 3D Folding: ESMFold v1 Multimer (+512 chain-break offset)           │
│  3. 3D Refinement: Constrained PyRosetta FastRelax (ref2015)            │
│  4. 6-Model Ensemble Generation:                                        │
│     • AntiFold (3D Antibody/Nanobody Inverse Folding)                   │
│     • ProteinMPNN - Diverse & Soluble (3D Backbone GNNs)                │
│     • ThermoMPNN (3D ΔΔG Single + Combo Scanner)                        │
│     • ESM-C 600M (1D Evolutionary Masked Marginals, Isolated Env)       │
│     • AbLang2-Paired (1D Paired VH-VL Germline Scanner)                 │
│  5. Post-Filtering: Mask enforcement, mutation cap (1..K), &            │
│     cross-model deduplication (`+` Consensus tagging)                   │
└─────────────────────────────────────────────────────────────────────────┘
   │
   ├──► Export `{run_name}_camsol_upload.fasta` ──► Cambridge CamSol Server (Optional)
   │
   ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ PHASE 2: ZERO-SHOT BIOPHYSICAL & EVOLUTIONARY SCORING (pipeline.py)     │
│  • Fail-Fast CamSol ID & Sequence Integrity Verification                │
│  • 3D & 1D Thermodynamics: PyRosetta ΔΔG, JanusDDG ΔΔG                  │
│  • Domain Interface & Colloidal Physics: BINDPRED Kd, Charge Asymmetry, │
│    Net Charge (pH 7.4), Isoelectric Point (pI & |pI - 7.4|), GRAVY      │
│  • 5 Protein/Antibody LMs: ESM-2 (650M), ESM-C (600M), AntiBERTy,       │
│    AbLang2-Paired, ProtBERT (+ optional PoET-2 API)                     │
│  • Ranking: Agnostic Gated + Invariant Score & Format-Weighted Score    │
└─────────────────────────────────────────────────────────────────────────┘
   │
   ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ INTERACTIVE TRIAGE & WET-LAB SHORTLISTING (app.py & Triage Scripts)     │
│  • Interactive Streamlit 3D Viewer (py3Dmol), Multi-Gate Filter,        │
│    4D Developability Landscape, & Normalized Radar Footprints           │
│  • Greedy Jaccard Diversity Pruning (Exact Substitution & Position)     │
│  • In-Silico nDCG@20 Regime Diagnostics & 16-Metric Scorecard Heatmaps  │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 🔬 Phase 1: Structure Preparation & Orthogonal Generative Models

### 1. Multi-Chain Structure Prediction & Constrained Relaxation

1. **ESMFold Multimer Folding:** Concatenates all input chains into a single pass through `facebook/esmfold_v1` (`torch.bfloat16`), injecting a cumulative `+512` positional index offset at each chain boundary so the folding trunk models distinct polypeptide chains accurately. Output coordinates are automatically re-indexed into standard per-chain PDB numbering (`A`, `B`, `C`...) with `TER` cards.
2. **Constrained PyRosetta `FastRelax`:** Relaxes the predicted structure under the `ref2015` full-atom energy function with backbone coordinate constraints (`-relax:constrain_relax_to_start_coords -ex1 -ex2 -use_input_sc`) to relieve local side-chain clashes without distorting the global fold.

### 2. The Six Generative Models

Relying on a single generative architecture introduces systematic bias: 3D inverse-folding GNNs frequently mutate buried/vernier residues at the expense of germline grammar, whereas 1D language models gravitate toward germline consensus residues that can introduce 3D steric clashes. Phase 1 resolves this by sampling across six orthogonal generators:

| Generator Model | Modality & Architecture | Sampling & Combination Strategy | What It Captures in Theory |
| :--- | :--- | :--- | :--- |
| **`AntiFold`** | 3D Antibody/Nanobody Inverse Folding (`antifold.main`) | Samples full sequences conditioned on the relaxed PDB across mutable IMGT regions (`FWH/L`, `CDRH/L`); supports both paired `VH/VL` and `--nanobody_chain` modes. | Optimizes 3D rotamer packing while respecting immunoglobulin framework and loop geometry. |
| **`ProteinMPNN (Diverse)`** | 3D Backbone Graph Neural Network (`protein_mpnn_run.py`) | Samples conditioned on Cα/N/C/O coordinates using `fixed_positions.jsonl`, `--sampling_temp`, and Gaussian `--backbone_noise` (0.02 Å). | Identifies sequences that maximize 3D conformational stability on the target backbone. |
| **`ProteinMPNN (Soluble)`** | 3D Surface-Reweighted Backbone GNN (`--use_soluble_model`) | Identical coordinate conditioning to `ProteinMPNN`, using weights trained strictly on soluble proteins. | Replaces solvent-exposed hydrophobic patches with polar or charged side chains to prevent surface sticking. |
| **`ThermoMPNN`** | 3D ΔΔG Point-Mutation Scanner (`custom_inference.py`) | Filters all single mutations with ΔΔG_pred < cutoff (default -0.3 kcal/mol); builds **50% top single mutants** and **50% multi-mutant combinations** (2..K mutations across top 30 positions). | Directly targets local thermodynamic unfolding bottlenecks via transfer-learned Megascale ΔΔG predictions. |
| **`ESM-C (600M)`** | 1D General Protein Language Model (`biohub/ESMC-600M-hf`) | Runs in an isolated Conda environment (`esmc_env` via `run_esmc_dms.py`) to rank masked-marginal point mutations; builds **50% single mutants + 50% combos**. | Captures deep evolutionary conservation and secondary-structure propensities across natural proteins. |
| **`AbLang2-Paired`** | 1D Paired Antibody Language Model (`ablang2-paired`) | Scans each unmasked position with `*` in `mode='probability'`, retaining mutations where P(mut) > P(WT); builds **50% single mutants + 50% combos**. | Captures co-evolutionary `VH–VL` pairing grammar and human germline residue frequencies. |

### 3. Post-Filtering, Mutation Cap & Cross-Model Consensus Tracking

* **Strict Mask & Cap Enforcement (`enforce_mutation_cap_and_filters`):** Every raw sequence is checked against the global 1-indexed frozen position set and banned amino acid list (`omit_aas`, default `["C"]`). If a multi-mutant sequence from `AntiFold`, `ProteinMPNN`, or `SolubleMPNN` exceeds `max_muts_per_seq` (*K*), excess mutations are uniformly subsampled down to `[1..K]` so all models span a comparable mutational load.
* **Cross-Model Consensus Detection:** Identical sequences proposed independently by multiple generators are collapsed into a single record and tagged with `+` in their sequence ID (e.g., `ThermoMPNN_W165L_rank3+ESMC`). In Phase 2, these variants are grouped as **`Consensus (Multi-Model)`** and consistently exhibit the highest gate-pass rates.

---

## 📊 Phase 2: 16-Metric Biophysical & Evolutionary Scoring Engine

Phase 2 (`score_universal_variants`) evaluates the pooled FASTA library alongside the unmutated Wild-Type (`{run_name}_WT`) across 16 structural, colloidal, and evolutionary metrics. Smart metric-level caching automatically merges existing raw calculations if a run is resumed or updated with a new CamSol file.

### Complete Metric Reference & Theoretical Interpretation

| Category | Metric Column (`CSV`) | UI Display Name | Favorable Direction | Theoretical Property Captured |
| :--- | :--- | :--- | :---: | :--- |
| **3D Structure & Interface** | `pyrosetta_ddg` | Stability ΔΔG (PyRosetta) | **Lower (↓)**<br>*(Negative = Stabilizing)* | Atomic van der Waals steric clashes, H-bonding, and rotamer strain after side-chain repacking (`mutate_residue` under `ref2015`) relative to WT (REU). |
| | `bindpred_kd` | Affinity log₁₀ Kd (BINDPRED) | **Lower (↓)**<br>*(Tighter Pairing)* | `VH–VL` domain interface binding affinity predicted via CatBoost over mean-pooled ESM-2 layer-33 representations. Prevents domain "breathing" that exposes hydrophobic cores. |
| **1D Folding Stability** | `janus_ddg` | Stability ΔΔG (Janus) | **Higher (↑)**<br>*(Positive = Stabilizing)* | Bidirectional cross-attention transformer prediction of sequence-level folding ΔΔG (kcal/mol), providing a coordinate-free check against unfolding. |
| **Electrostatics & Colloidal Repulsion** | `pI` & `pi_delta` | Isoelectric Point (pI) & \|pI − 7.4\| | **Higher `pi_delta` (↑)** | Distance between the protein's isoelectric point and physiological buffer (pH 7.4). Larger \|pI − 7.4\| increases net monopole repulsion (DLVO stability). |
| | `charge_asymmetry` | VH-VL Charge Delta | **Lower (↓)**<br>*(Weaker Dipole)* | Mean absolute charge difference \|Q_VH(7.4) − Q_VL(7.4)\| across paired Fv arms. High asymmetry creates an electrostatic dipole that drives head-to-tail self-association and viscosity. |
| | `net_charge` | Net Charge (pH 7.4) | Context-Dependent | Total polypeptide charge at pH 7.4 computed via Henderson-Hasselbalch pKa summation (`Bio.SeqUtils.ProtParam`). |
| **Intrinsic Solubility & Hydropathy** | `camsol_score` | Intrinsic Solubility (CamSol) | **Higher (↑)** | Sliding-window physicochemical propensity for intrinsic solubility and resistance to β-sheet aggregation (mapped from Cambridge CamSol server `.txt` output). |
| | `gravy` | Hydrophobicity Index (GRAVY) | **Lower (↓)**<br>*(More Negative)* | Grand Average of Hydropathicity across the full sequence using the Kyte-Doolittle scale; lower values indicate reduced bulk hydrophobic sticking. |
| **1D Evolutionary & Antibody PLMs** | `esm_fitness` | ΔLog-Likelihood (ESM-2 General) | **Higher (↑)** | Sum of masked-marginal log-odds Σ [log P(mut_i) − log P(WT_i)] from `esm2_t33_650M_UR50D`. Penalizes disruption of conserved structural motifs. |
| | `esmc_fitness` | ΔLog-Likelihood (ESM-C) | **Higher (↑)** | Next-generation evolutionary log-likelihood delta computed via `biohub/ESMC-600M-hf` in `esmc_env` (`run_esmc_scoring.py`). |
| | `antiberty_fitness` | ΔLog-Likelihood (AntiBERTy) | **Higher (↑)** | Per-domain antibody masked-marginal log-odds from `jeffruffolo/AntiBERTy` evaluated independently on each `VH`, `VL`, or `VHH` chain. |
| | `ablang_fitness` | ΔLog-Likelihood (AbLang2 Paired) | **Higher (↑)** | Paired heavy-light germline pseudo-log-likelihood delta (PLL_mut − PLL_WT) from `ablang2-paired`. Because `pseudo_log_likelihood` returns negative cross-entropy loss (≤ 0), values closer to 0 (positive ΔPLL) indicate superior germline pairing (`invert=False`). |
| | `protbert_fitness` | ΔLog-Likelihood (ProtBERT) | **Higher (↑)** | Masked-marginal log-likelihood delta from `Rostlab/prot_bert`, providing an orthogonal general-protein evolutionary check. |
| | `poet2_fitness` | ΔLog-Likelihood (PoET-2) | **Higher (↑)** | *(Optional)* Homology/MSA-prompted generative log-likelihood delta via the OpenProtein.AI API. |
| **Parsimony & Composite** | `num_mutations` | Number of Mutations | **Lower (↓)** | Total amino acid substitutions relative to WT. Acts as a parsimony tie-breaker to minimize epistatic risk. |
| | `agnostic_gated_score` | Agnostic Gated Score | **Higher (↑)** | Primary zero-shot Z-score ranking metric combining the 5-term Invariant Core, 3-term Downside Hinge, and non-compensatory Hard Gate penalties. |
| | `global_composite_score` | Global Composite Score (Legacy) | **Higher (↑)** | Format-specific weighted average (`FORMAT_WEIGHTS`) of min-max normalized `[0, 1]` metrics retained for historical comparison. |

---

## 🛡️ The Agnostic Gated + Invariant Score

Standard linear composite scores fail in multi-mutant antibody engineering because generative models can stack 4–8 germline consensus mutations to inflate 1D language model scores (`esm_fitness`) while introducing 3D steric clashes that degrade wet-lab aggregation (`%HMW`).

The **Agnostic Gated + Invariant Score** resolves this using a **2 : 1 : 1 architecture** calculated on within-campaign Z-scores ($Z_x$) of deltas relative to the Wild-Type baseline ($\Delta x = x_{\text{mut}} - x_{\text{WT}}$):

### 1. Component A: The Invariant Physical Core (Weight Budget = `1.00`)

Uses Dawes-style rank-proportional staircase weights (`0.30, 0.25, 0.20, 0.15, 0.10`) across the five physical metrics whose anti-aggregation directionality is invariant across antibody scaffolds:

$$
\text{Core}_i = -0.30 \cdot Z_{\text{PyRosetta}}(i) - 0.25 \cdot Z_{\Delta\text{BINDPRED}}(i) - 0.20 \cdot Z_{\Delta\text{Asym}}(i) + 0.15 \cdot Z_{\Delta\text{pI-delta}}(i) - 0.10 \cdot Z_{N_{\text{mut}}}(i)
$$

### 2. Component B: One-Sided Downside Folding Hinge (Weight Budget = `0.50`)

Applies $h_x = \min(0.0, \Delta x)$ to `esm_fitness`, `janus_ddg`, and `camsol_score` before Z-scoring. Variants scoring above WT receive a capped hinge of `0.0` (upside-agnostic to consensus inflation), whereas variants dropping below WT are proportionally penalized:

$$
\text{Hinge}_i = +0.25 \cdot Z_{h(\text{ESM})}(i) + 0.15 \cdot Z_{h(\text{Janus})}(i) + 0.10 \cdot Z_{h(\text{CamSol})}(i)
$$

### 3. Component C: Non-Compensatory Biophysical Hard Gates (`-0.50` Penalty Each)

Three boolean failure gates flag severe structural or evolutionary defects:

* **Gate 1 (`3D_Clash`):** `(Z_PyRosetta > 0.67 AND ΔΔG_PyRosetta > 0.0 REU) OR (ΔΔG_PyRosetta > +5.0 REU)`
* **Gate 2 (`Interface_Kd`):** `has_paired_fv AND (Z_ΔBINDPRED > 0.67 AND Δlog10(Kd) > 0.0)`
* **Gate 3 (`1D_Unfolding`):** `(ΔESM-2 < -10.0) OR (ΔΔG_Janus < -0.5 kcal/mol)`

$$
\text{Score}_{\text{Agnostic}}(i) = \text{Core}_i + \text{Hinge}_i - 0.50 \cdot \left( g_{\text{clash}}(i) + g_{\text{interface}}(i) + g_{\text{unfold}}(i) \right)
$$

All candidates passing all three gates (`passed_all_gates == True`, `gate_status == "PASS"`) are ranked ahead of gate-failing variants.

---

## 🖥️ Interactive Triage & Jaccard Diversity Diagnostics

### 1. Streamlit Web Application (`app.py`)

* **Tab 1 (Generation):** Interactive format selector, IMGT region freezing presets, sampling hyperparameters, live multi-chain FASTA validator, CamSol export button, and an embedded `py3Dmol` 3D structure viewer with RCSB-style residue hover tooltips and chain/surface toggles.
* **Tab 2 (Evaluation & Multi-Gate Triage):**
  * **Fail-Fast CamSol Validation:** Verifies ID overlap and exact amino acid sequence integrity before running PyRosetta; automatically removes mismatched uploads and disables WT-relative checkboxes if the WT baseline is `NaN`.
  * **Interactive Anti-Aggregation Multi-Gate Filter:** Live sliders (`max_rosetta = +10.0 REU`, `min_janus = -0.5 kcal/mol`, `min_ablang = -1.0 ΔLL`, `min_esm = -10.0 ΔLL`) and WT-relative checkboxes (`GRAVY <= WT`, `CamSol >= WT`, `BINDPRED <= WT`, `Charge Dipole <= WT`).
  * **4D Developability Scatter Landscape & Normalized Radar Footprint:** Interactive Plotly charts with WT baseline crosshairs and `[0, 1]` min-max normalized multi-spoke footprints.

