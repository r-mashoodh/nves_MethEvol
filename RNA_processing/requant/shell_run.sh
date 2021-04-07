#!/bin/bash

for pid in $(pidof -x test.sh); do
    if [ $pid != $$ ]; then
        echo "[$(date)] : test.sh : Process is already running with PID $pid"
        exit 1
    fi
done
