#!/bin/bash
# Watch CPU tree pipelines; on COMPLETED, verify groups and submit mig-max precompute.
# nohup bash scripts/_ops/watch_tree_pipes_then_precompute.sh &
set -uo pipefail

export SLURM_CONF=/cm/shared/apps/slurm/etc/slurm/slurm.conf
SQUEUE=/cm/local/apps/slurm/current/bin/squeue
SACCT=/cm/local/apps/slurm/current/bin/sacct
SBATCH=/cm/local/apps/slurm/current/bin/sbatch
export PATH=/cm/local/apps/slurm/current/bin:/vast/home/n/nnori/.conda/envs/treesbm/bin:/usr/bin:/bin

cd ~/DiscreteTreeFlows
mkdir -p logs scripts/_ops
LOG=logs/watch_tree_pipes_precompute.log
STATE=scripts/_ops/watch_state.env
touch "$STATE"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

PIPES=(
  "7350452|data/h1n1_temporal|566|h1n1_temporal"
  "7350453|data/covid_temporal|1280|covid_temporal"
  "7350454|data/covid_cladeholdout|1280|covid_cladeholdout"
)

job_state() {
  local jid="$1"
  local st
  st=$("$SACCT" -j "$jid" -n -X -o State 2>/dev/null | head -1 | tr -d ' ')
  if [[ -z "$st" ]]; then
    st=$("$SQUEUE" -j "$jid" -h -o %T 2>/dev/null | head -1 | tr -d ' ')
  fi
  echo "${st:-UNKNOWN}"
}

verify_groups() {
  local root="$1"
  local ok=1
  for split in train val test; do
    local d="$root/$split"
    local n_anc n_nwk n_both
    n_anc=$(ls "$d"/group_*_anc_aa.fasta 2>/dev/null | wc -l | tr -d ' ')
    n_nwk=$(ls "$d"/group_*_rooted.nwk 2>/dev/null | wc -l | tr -d ' ')
    n_both=$(python3 - <<PY
from pathlib import Path
p=Path("$d")
anc={int(x.stem.split("_")[1]) for x in p.glob("group_*_anc_aa.fasta")}
nwk={int(x.stem.split("_")[1]) for x in p.glob("group_*_rooted.nwk")}
print(len(anc & nwk))
PY
)
    log "  verify $d: anc=$n_anc rooted=$n_nwk both=$n_both"
    if [[ "$n_both" -lt 1 ]]; then
      log "  FAIL: no complete groups in $d"
      ok=0
    fi
  done
  return $(( 1 - ok ))
}

already_submitted() {
  local label="$1"
  grep -q "^PRECOMPUTE_${label}=" "$STATE" 2>/dev/null
}

mark_submitted() {
  local label="$1" jid="$2"
  grep -v "^PRECOMPUTE_${label}=" "$STATE" > "${STATE}.tmp" 2>/dev/null || true
  echo "PRECOMPUTE_${label}=$jid" >> "${STATE}.tmp"
  mv "${STATE}.tmp" "$STATE"
}

submit_precompute() {
  local root="$1" msl="$2" label="$3"
  local jname
  jname=$(echo "pre_${label}" | cut -c1-18)
  local jid
  jid=$("$SBATCH" --qos=mig-max --parsable \
    --job-name="$jname" \
    --export=ALL,DATA_ROOT="$root",MAX_SEQ_LEN="$msl",BATCH_SIZE=8 \
    scripts/_ops/slurm_precompute_dataset.sh)
  log "SUBMITTED precompute $label -> job $jid (DATA_ROOT=$root MAX_SEQ_LEN=$msl)"
  mark_submitted "$label" "$jid"
}

poll_once() {
  local all_done=1
  for entry in "${PIPES[@]}"; do
    IFS='|' read -r jid root msl label <<< "$entry"
    local st
    st=$(job_state "$jid")
    log "pipe $jid ($label): $st"
    case "$st" in
      COMPLETED)
        if already_submitted "$label"; then
          log "  precompute already recorded: $(grep "^PRECOMPUTE_${label}=" "$STATE")"
        else
          if verify_groups "$root"; then
            submit_precompute "$root" "$msl" "$label" || log "  ERROR submitting precompute for $label"
          else
            log "  SKIP submit: verification failed for $label"
            # treat as not fully handled so watcher can retry later if files appear
            all_done=0
          fi
        fi
        ;;
      FAILED|TIMEOUT|CANCELLED|NODE_FAIL|OUT_OF_MEMORY|BOOT_FAIL|DEADLINE)
        log "  TERMINAL FAILURE for $label ($st) — not submitting precompute"
        ;;
      RUNNING|PENDING|COMPLETING|CONFIGURING)
        all_done=0
        local lf n_anc
        lf=$(ls -t logs/*_${jid}.log 2>/dev/null | head -1 || true)
        if [[ -n "${lf:-}" ]]; then
          n_anc=$(grep -c 'stopped after translate' "$lf" 2>/dev/null || echo 0)
          log "  progress: translate_done~$n_anc log=$(basename "$lf")"
        fi
        ;;
      *)
        all_done=0
        log "  unexpected state $st (squeue/sacct may be flaky)"
        ;;
    esac
  done
  return $(( 1 - all_done ))
}

log "=== watcher start pid=$$ ==="
log "squeue=$SQUEUE sacct=$SACCT"
MAX_HOURS="${MAX_HOURS:-20}"
INTERVAL_SEC="${INTERVAL_SEC:-180}"
END_TS=$(( $(date +%s) + MAX_HOURS * 3600 ))

while true; do
  if poll_once; then
    log "All pipelines handled. Exiting."
    break
  fi
  now=$(date +%s)
  if [[ $now -ge $END_TS ]]; then
    log "Reached MAX_HOURS=$MAX_HOURS without all done."
    cat "$STATE" | tee -a "$LOG" || true
    log "Re-run: bash scripts/_ops/watch_tree_pipes_then_precompute.sh"
    break
  fi
  sleep "$INTERVAL_SEC"
done

log "=== watcher end ==="
