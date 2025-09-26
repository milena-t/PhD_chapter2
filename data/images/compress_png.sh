#!/bin/sh

exiftool -all= -overwrite_original $1 && optipng -o7 -strip all $1