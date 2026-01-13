#!/bin/bash
#==============================================================================
# DynamicBind - Submit All Target Simulations
#
# This script submits GaMD simulations for all target proteins to Polaris.
# Manages job dependencies and resource allocation across 5,000 node-hours.
#
# Usage:
#   ./submit_all_targets.sh [--test] [--dry-run]
#
# Options:
#   --test     Submit short test jobs (1 hour each)
#   --dry-run  Print commands without submitting
#==============================================================================

set -e

# Parse arguments
TEST_MODE=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            TEST_MODE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Project configuration
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
SCRIPT_DIR="${PROJECT_DIR}/scripts/polaris"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p ${LOG_DIR}

echo "=================================================="
echo "DynamicBind - Batch Job Submission"
echo "=================================================="
echo "Project Directory: ${PROJECT_DIR}"
echo "Test Mode: ${TEST_MODE}"
echo "Dry Run: ${DRY_RUN}"
echo "=================================================="

# Define targets and their PDB IDs
# Format: TARGET:PDB_ID:PRIORITY:NODE_HOURS
declare -a JOBS=(
    # c-Myc - HIGHEST PRIORITY (70% of cancers)
    "cmyc:6G6K:1:100"
    "cmyc:1NKP:2:80"
    "cmyc:2A93:2:80"

    # p53 - HIGH PRIORITY (50% of cancers)
    "p53:2AHI:1:100"
    "p53:2FEJ:2:80"

    # Alpha-Synuclein - HIGH PRIORITY (Parkinson's)
    "alpha_synuclein:8A9L:1:80"
    "alpha_synuclein:2KKW:1:80"
    "alpha_synuclein:8FPT:2:60"

    # Amyloid-Beta - MEDIUM PRIORITY (Alzheimer's)
    "abeta42:2NAO:1:80"
    "abeta42:1Z0Q:2:60"

    # Tau K18 - MEDIUM PRIORITY (Alzheimer's)
    "tau_k18:2MZ7:2:60"
)

# Track total hours
TOTAL_HOURS=0

# Submit function
submit_job() {
    local TARGET=$1
    local PDB_ID=$2
    local PRIORITY=$3
    local HOURS=$4

    if [ "$TEST_MODE" = true ]; then
        HOURS=1
        WALLTIME="01:00:00"
        SCRIPT="submit_gamd_test.pbs"
    else
        WALLTIME="24:00:00"
        SCRIPT="submit_gamd.pbs"
    fi

    TOTAL_HOURS=$((TOTAL_HOURS + HOURS))

    echo ""
    echo "Submitting: ${TARGET}/${PDB_ID}"
    echo "  Priority: ${PRIORITY}"
    echo "  Hours: ${HOURS}"
    echo "  Total Allocated: ${TOTAL_HOURS}"

    if [ "${TOTAL_HOURS}" -gt 5000 ]; then
        echo "  WARNING: Exceeding 5000 node-hour allocation!"
    fi

    CMD="qsub -v TARGET=${TARGET},PDB_ID=${PDB_ID} ${SCRIPT_DIR}/${SCRIPT}"

    if [ "$DRY_RUN" = true ]; then
        echo "  [DRY RUN] ${CMD}"
    else
        JOB_ID=$(${CMD})
        echo "  Submitted: ${JOB_ID}"
        echo "${TARGET},${PDB_ID},${JOB_ID},${HOURS}" >> ${LOG_DIR}/submitted_jobs.csv
    fi
}

# Initialize job log
if [ "$DRY_RUN" = false ]; then
    echo "target,pdb_id,job_id,hours" > ${LOG_DIR}/submitted_jobs.csv
fi

# Sort by priority and submit
echo ""
echo "Submitting Priority 1 jobs first..."
for job in "${JOBS[@]}"; do
    IFS=':' read -r target pdb_id priority hours <<< "$job"
    if [ "$priority" = "1" ]; then
        submit_job "$target" "$pdb_id" "$priority" "$hours"
    fi
done

echo ""
echo "Submitting Priority 2 jobs..."
for job in "${JOBS[@]}"; do
    IFS=':' read -r target pdb_id priority hours <<< "$job"
    if [ "$priority" = "2" ]; then
        submit_job "$target" "$pdb_id" "$priority" "$hours"
    fi
done

echo ""
echo "=================================================="
echo "Submission Complete"
echo "Total Hours Requested: ${TOTAL_HOURS}"
echo "Remaining Allocation: $((5000 - TOTAL_HOURS))"
echo "=================================================="

if [ "$DRY_RUN" = false ]; then
    echo "Job log: ${LOG_DIR}/submitted_jobs.csv"
fi
