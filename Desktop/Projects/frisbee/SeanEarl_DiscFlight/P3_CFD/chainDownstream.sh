#!/bin/bash
# Wait for the CFD follow-up to finish, then run the whole downstream pipeline.
cd "$(dirname "$0")"
until grep -q "FOLLOW-UP COMPLETE" cases/followup.log 2>/dev/null; do sleep 20; done
echo "### CFD chain done; starting downstream ###"
bash ../runDownstream.sh 20 150
