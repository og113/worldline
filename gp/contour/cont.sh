#!/bin/bash
gawk -v d=$2 -v w=$3 -v os=$4 'function abs(x) { return (x>=0?x:-x) }
    {
            if($0~/# Contour/) nr=0
            if(nr==int(os+w/2) && d==0) {i++; a[i]=$1; b[i]=$2; c[i]=$3;}
            if(abs(nr-os-w/2)>w/2 && d==1) print $0
            nr++
    }
    END {   if(d==0) {
                    for(j=1;j<=i;j++)
                    printf "set label %d \"%g\" at %g, %g centre front\n", j, c[j], a[j], b[j]
            }
    }' $1
