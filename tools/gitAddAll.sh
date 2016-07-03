#!/bin/bash

# program to add all the usual suspects to git

# adding pics
git add pics/*/*.png
git add pics/*/*.svg
git add pics/*/*/*.png
git add pics/*/*/*.svg

# adding .cc and .h files
git add mpisrc/*.cc
git add csrc/*.cc
git add src/*.cc
git add tests/*.cc
git add include/*.h

# adding .sh files
git add tools/*.sh

# adding notes
git add notes/*.txt

# adding reports
git add reports/*

# adding gp files
git add gp/*.gp

# adding mout files
git add mout/*.txt

# adding database outputs
git add db/nr/plots/*.csv
git add db/nr/plots/*/*.csv

# adding results
git add results/*/*.csv
