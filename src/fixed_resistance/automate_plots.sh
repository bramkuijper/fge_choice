#!/usr/bin/env bash
find . -iname "output_*"
find . -iname "output_*" -print0 | xargs -0 -P8 -I % ./plot_output.r %
