# Betty submit (VOC threat panel)

## Blocker
Laptop SSH ControlMaster expired; Duo required again. Scripts+cases are ready locally but **not yet synced**.

## On laptop (after `kinit` + Duo SSH)
```bash
cd ~/Documents/GitHub/DiscreteTreeFlows
bash scripts/betty_voc_threat_go.sh
```

This will:
1. Re-open SSH (approve Duo)
2. Ensure `data/covid/{test,val,train}` → lab DiscreteTreeFlows_data (already fixed once)
3. Rsync VOC scripts + primary cases
4. Submit **9 TreeSBM** jobs (`slurm_voc_threat_panel.sh`)
5. Submit **9 NeutralBD baseline** jobs (`slurm_voc_threat_baselines.sh`)

Dry-run only:
```bash
DRY=1 bash scripts/betty_voc_threat_go.sh
```

## Primary cases (9)
| VOC | split | group |
|-----|-------|------:|
| Gamma | test | 3 |
| Delta | test | 7 |
| Delta | test | 11 |
| Beta | train | 302 |
| Delta | train | 191 |
| Omicron_BA1 | train | 309 |
| Alpha | train | 323 |
| Lambda | train | 289 |
| Mu | train | 34 |

## After jobs finish
```bash
scp nnori@login.betty.parcc.upenn.edu:~/DiscreteTreeFlows/results/voc_threat_panel/submit_ids.json results/voc_threat_panel/
scp nnori@login.betty.parcc.upenn.edu:~/DiscreteTreeFlows/results/voc_threat_panel/baseline_submit_ids.json results/voc_threat_panel/
rsync -avP nnori@login.betty.parcc.upenn.edu:~/DiscreteTreeFlows/results/voc_threat_panel/cases/ results/voc_threat_panel/cases/
```
