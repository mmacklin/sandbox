#!/bin/sh
ffmpeg -r 60 -f image2 -i ./dump/frame%d.tga -threads 0 -vcodec libx264 -vpre slow -crf 20 $1
