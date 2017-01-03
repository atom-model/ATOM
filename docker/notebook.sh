#!/bin/sh
exec /usr/bin/jupyter notebook >>/var/log/notebook.log 2>&1

