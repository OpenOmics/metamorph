#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh && conda activate metawrap-env
source ~/.bashrc
exec "$@"