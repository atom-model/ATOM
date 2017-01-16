#!/bin/sh
cd /workspace
exec /usr/local/bin/jupyter notebook --no-browser --ip=0.0.0.0 --NotebookApp.token='' >>/var/log/notebook.log 2>&1
