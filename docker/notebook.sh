#!/bin/sh
cd /home/ATOM/
exec /usr/local/bin/jupyter notebook --no-browser --allow-root --ip=0.0.0.0 --NotebookApp.token=''
