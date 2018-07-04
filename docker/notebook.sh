#!/bin/sh
cd /build/ATOM/benchmark
exec /usr/local/bin/jupyter notebook --no-browser --allow-root --ip=0.0.0.0 --NotebookApp.token='' >>/var/log/notebook.log 2>&1 &
