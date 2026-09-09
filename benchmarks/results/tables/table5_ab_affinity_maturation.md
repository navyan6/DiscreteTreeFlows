# Paper Table 5 — Antibody affinity maturation forecasting

**Updated:** 2026-08-26. Numbers from `t5_summaries_betty.json` / `table6_ab_track_c.json` (do not invent).  
Screenshot draft mixed ε levels and some OAS vs pathogen cells — **use the paste table below**.

## Experiment summary

**Question.** Can TreeSBM generate plausible human IGH maturation trees from early clonal sequences, or does it only capture generic protein-family drift? The assay checks CDR-targeted mutation recovery, SHM load, terminal diversity, and neighborhood coverage — not RF (topology is forced for CTMC/TreeSBM rows).

**Data.** Rodriguez AIRR-seq race production IGH parental–child pairs (PCP), DASM-processed:  
`…/rodriguez-airr-seq-race-prod-NoWinCheck_igh_pcp_2024-11-12_MASKED_NI_noN_no-naive.csv.gz`  
Trees: `antibody_benchmark/data/processed/benchmark_trees.jsonl`.  
**Holdout:** **82** clonal families (`ab_family_ids_82.txt`; sample-igg-SC-* / W-*). CDR site masks from PCP codon coords (IMGT CDR).

**Split / protocol (Track C).** These 82 families are a **frozen eval set**, not the OAS train split. For each family: condition on the **observed early/root sequence**, paint **N=20** rollouts on the **observed Newick + branch lengths** (forced topo for Neutral / Thrifty / TreeSBM / pLM; CoSiNE = unguided Gillespie on observed edges; AR = free topology from an N16 ARTreeFormer pool + JTT). Score up to **K=100** generated leaf AAs vs true leaves.

### Fairness / training provenance (re-checked)

| Method | Ab-trained? | What it actually is |
|---|---|---|
| **Neutral SHM** | No params | JC69 NT CTMC — domain-shaped prior, not learned on Rod.82 |
| **Thrifty** (additive) | **Yes** (SHM rates) | `ThriftyHumV0.2-59` — Ab context SHM; **not** the Neutral row |
| **DASM+Thrifty** (additive) | **Yes** | Ab DASM + Thrifty |
| **CoSiNE** | **Yes** | `cosine_dasm.ckpt` (HF `thematrixmaster/cosine`) — Ab DASM CTMC |
| **AR tree-edit** | **No** | Viral ARTreeFormer N16 pool + shared **JTT**; not Ab-trained |
| **pLM prior** (additive) | **No** | ESM-2 8M general protein LM |
| **TreeSBM (main T5 row)** | **No** | `checkpoints/best.pt` ep237 val=16.79 — **pathogen** multi-virus ckpt (`full.yaml`); OOD on Ab |
| **TreeSBM OAS v1** (additive) | **Yes** | `ab_oas_1m_v1` donor-disjoint OAS clones → Rod.82 |
| **TreeSBM ab_dasm** (additive) | **Yes*** | `ab_dasm_v1` on non-Rod DASM trees (*loss collapsed; weak*) |

**Verdict:** the screenshot/main table is **not a matched-domain bake-off**. CoSiNE (and Thrifty/DASM) are Ab SHM models; the printed TreeSBM row is a **viral** ckpt. Even the fairer Ab-trained OAS row is **worse** on CDR (0.150 vs pathogen 0.278) — so domain mismatch is real, but Ab fine-tune with ESM Q0 did not close the SHM gap.

**Metrics (paper columns).**
| Metric | Definition | Better |
|---|---|---|
| CDR mut. recall | Fraction of true root→leaf CDR substitutions recovered in gen leaves | ↑ |
| SHM load error | Error between gen vs true SHM-count distributions on terminals | ↓ |
| Terminal diversity error | \|mean pairwise Hamming(gen leaves) − mean pairwise Hamming(true leaves)\| | ↓ |
| Coverage@100 | Absolute Hamming **ε=5** coverage of true leaves by the K=100 gen pool | ↑ |

Viral tables use ε=2; Ab paper Coverage@100 is **ε=5** (neighborhood is wider under SHM).

**Methods in the main table.**
| Row | Process |
|---|---|
| Neutral SHM | JC69 independent-site NT CTMC on observed topo+BL (**not** Thrifty/S5F) |
| Autoregressive tree-edit | Adapted ARTreeFormer + JTT; free topo |
| CTMC (CoSiNE) | CoSiNE DASM ckpt; Gillespie on observed edges |
| TreeSBM | Pathogen `best.pt` on Rod.82 (forced topo+BL) |

---

## Paste-ready (corrects screenshot)

| Method | CDR mut. recall ↑ | SHM load error ↓ | Terminal diversity error ↓ | Coverage@100 ↑ |
|---|---:|---:|---:|---:|
| Neutral SHM model | **0.688** | 0.054 | 17.435 | 0.059 |
| Autoregressive tree-edit model | 0.000 | 0.117 | 19.142 | **0.081** |
| CTMC model (CoSiNE) | 0.658 | **0.037** | **3.019** | 0.076 |
| TreeSBM (pathogen ckpt) | 0.278 | 0.055 | 5.104 | **0.084** |

**Screenshot fixes:** CoSiNE CDR is **0.658** (not 0.678); AR Cov@100 is **0.081** (ε=5), not 0.0417 (that was ε=3); TreeSBM terminal-div is **5.10** (pathogen), not 9.09 (that tracks OAS `treesbm_ab_oas` ≈9.89).

### Additive context (not main T5 rows)

| Method | CDR | SHM err | Term-div err | Cov@ε5 |
|---|---:|---:|---:|---:|
| Thrifty (≠ Neutral) | 0.746 | 0.037 | 12.97 | 0.058 |
| DASM+Thrifty | 0.702 | 0.028 | 2.71 | 0.085 |
| pLM prior (ESM-2 8M) | 0.297 | 0.032 | 2.63 | 0.076 |
| TreeSBM OAS v1 (`treesbm_ab_oas`) | 0.150 | 0.077 | 9.89 | 0.081 |

**Takeaway.** SHM-aware CTMCs (Neutral JC69 / CoSiNE / Thrifty) dominate **CDR mut. recall**. Pathogen TreeSBM wins Ab **Coverage@ε5** slightly but under-recovers CDR mutations — expected with an ESM conservation Q0. OAS-trained TreeSBM does not close the CDR gap without an SHM Q0 (v2–v4 still ≪ Neutral).

Full ε grid + jobs: [`table6_ab_track_c.md`](table6_ab_track_c.md). Tracks A/B/C: [`antibody_tracks.md`](antibody_tracks.md).
