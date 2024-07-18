#!/bin/bash
make "$1"

if [ "$1" != "debug" ]; then
    ./main.out  # use a debugger to launch it instead
fi
