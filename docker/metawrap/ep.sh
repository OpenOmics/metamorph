#!/bin/bash --login
. /opt/conda/etc/profile.d/conda.sh && conda activate metawrap-env
export PATH="/home/metaWRAP/bin:$PATH"
exec "$@"