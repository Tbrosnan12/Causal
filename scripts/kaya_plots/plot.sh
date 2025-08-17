#!/bin/bash

csv=$1

pdf=$(math -script plot.m $csv)
xdg-open "pdfs/$pdf" >/dev/null 2>&1 &
