#!/bin/bash

if [ -z "$1" ]
  then
    echo "must supply input file"
else
	echo "changing FThermal to GThermal"
	sed -i 's/FThermal/GThermal/g' $1
fi
