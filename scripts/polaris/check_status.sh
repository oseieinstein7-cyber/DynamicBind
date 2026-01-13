#!/bin/bash
#==============================================================================
# DynamicBind - Check Job Status
#
# Monitor running and completed jobs on Polaris.
#
# Usage:
#   ./check_status.sh
#==============================================================================

echo "=================================================="
echo "DynamicBind - Job Status"
echo "=================================================="
echo ""

# Check queue status
echo "=== Current Jobs ==="
qstat -u $USER

echo ""
echo "=== Recent Job History ==="
qstat -x -u $USER | head -20

echo ""
echo "=== Allocation Usage ==="
# Show project allocation status (command varies by site)
echo "Check ALCF website for current allocation balance"

echo ""
echo "=== Disk Usage ==="
du -sh /eagle/projects/dynamicbind/* 2>/dev/null || echo "Check /eagle/projects/dynamicbind/ manually"

echo ""
echo "=== Recent Log Files ==="
PROJECT_DIR="/eagle/projects/dynamicbind/DynamicBind"
ls -lt ${PROJECT_DIR}/logs/*.log 2>/dev/null | head -10 || echo "No log files found"

echo ""
echo "=== Completed Trajectories ==="
find ${PROJECT_DIR}/data/trajectories -name "production_traj.dcd" -ls 2>/dev/null | head -10 || echo "No completed trajectories yet"

echo ""
echo "=================================================="
