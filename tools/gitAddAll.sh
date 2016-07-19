#!/bin/bash

# program to add all the usual suspects to git

# adding pics
git add --all pics/*/*.png
git add --all pics/*/*.svg
git add --all pics/*/*/*.png
git add --all pics/*/*/*.svg

# adding .cc and .h files
git add --all mpisrc/*.cc
git add --all csrc/*.cc
git add --all src/*.cc
git add --all tests/*.cc
git add --all include/*.h

# adding .sh files
git add --all tools/*.sh

# adding notes
git add --all notes/*.txt

# adding reports
git add --all reports/*

# adding gp files
git add --all gp/*.gp

# adding mout files
git add --all mout/*.txt

# adding database stuff
git add --all db/nr/*.sqlite
git add --all db/nr/nr*.csv
git add --all db/nr/plots/*.csv
git add --all db/nr/plots/*/*.csv

# adding results
git add --all results/*/*.csv

