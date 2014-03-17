#!/bin/bash

ffmpeg -qscale 5 -r 30 -b 9600 -i $1_%d.png $2
