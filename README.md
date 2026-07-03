# TreeSBM: Tree-Valued Schrödinger Bridge Matching

Generative modeling of protein evolutionary trees from root sequences.

## Setup

### Option A: Conda (recommended)
```bash
conda env create -f environment.yml
conda activate treesbm
```

### Option B: Pip (with system MAFFT + FastTree)
```bash
pip install -r requirements.txt
brew install mafft fasttree  # macOS
# or
apt install mafft fasttree   # Ubuntu/Debian
```

## Preprocessing

Prepare your data with ancestral sequence reconstruction:
```bash
python scripts/preprocess_trees.py data/processed
```

## Training

```bash
python train.py --config configs/protein_family.yaml
```

## Documentation

- `PLAN.md` — full implementation plan covering all phases
- `scripts/` — data preprocessing and training scripts
- `src/` — core model implementation
- `dataloaders/` — data loading and batching
- `configs/` — YAML configuration files
- `tests/` — unit tests
