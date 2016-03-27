#!/bin/bash
# program to simplify results for semiclassical calculation by reducing the number of parameters by one

i1='results/nr/nrmain_cosmos.dat'
i2='results/nr/nrmain_cosmos_2.dat'
i3='results/nr/nrmain_cosmos_3.dat'
i4='results/nr/nrmain_cosmos_4.dat'
i5='results/nr/nrmain_cosmos_5.dat'
i6='results/nr/nrmain_cosmos_6.dat'
i7='results/nr/nrmain_cosmos_7.dat'
i7l='results/nr/nrmain_laptop_7.dat'
i8='results/nr/nrmain_cosmos_8.dat'

o1='db/nr/nr1.csv'
o2='db/nr/nr2.csv'
o3='db/nr/nr3.csv'
o4='db/nr/nr4.csv'
o5='db/nr/nr5.csv'
o6='db/nr/nr6.csv'
o7='db/nr/nr7.csv'
o7l='db/nr/nr7l.csv'
o8='db/nr/nr8.csv'

o='db/nr/nr.csv'

bash toCsv.sh $i1
bash toCsv.sh $i2
bash toCsv.sh $i3
bash toCsv.sh $i4
bash toCsv.sh $i5
bash toCsv.sh $i6
bash toCsv.sh $i7
bash toCsv.sh $i7l
bash toCsv.sh $i8

m1=$(echo "$i1" | sed 's/\.\([a-z]*\)$/.csv/');
m2=$(echo "$i2" | sed 's/\.\([a-z]*\)$/.csv/');
m3=$(echo "$i3" | sed 's/\.\([a-z]*\)$/.csv/');
m4=$(echo "$i4" | sed 's/\.\([a-z]*\)$/.csv/');
m5=$(echo "$i5" | sed 's/\.\([a-z]*\)$/.csv/');
m6=$(echo "$i6" | sed 's/\.\([a-z]*\)$/.csv/');
m7=$(echo "$i7" | sed 's/\.\([a-z]*\)$/.csv/');
m7l=$(echo "$i7l" | sed 's/\.\([a-z]*\)$/.csv/');
m8=$(echo "$i8" | sed 's/\.\([a-z]*\)$/.csv/');

i1=${i1##*/}
i2=${i2##*/}
i3=${i3##*/}
i4=${i4##*/}
i5=${i5##*/}
i6=${i6##*/}
i7=${i7##*/}
i7l=${i7l##*/}
i8=${i8##*/}

sp='\([[:space:]]\+\)' # a space
spn='\([[:space:]]*\)' # a space, or nor
nm='\([0-9.]\+\)' # a number

echo "reduced results:"

# for nrmain_cosmos_8.csv
awk -F "," '{print $1","$2","$3","$4","$5**3*$6","$7","$8","$9","$10","$11","$5*$6*$12","$13",,,,,,"$14","$15","$16","$17","$18","$19","$20","$21","$i8}' $m8  > $o8
echo $o8

# for nrmain_laptop_7.csv
awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$4*$5*$10","$11",,,,,,"$12","$13","$14","$15","$16","$17","$18","$19","$i7l}' $m7l  > $o7l
echo $o7l

# same for nrmain_cosmos_7.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$4*$5*$10","$11",,,,,,"$12","$13","$14","$15","$16","$17","$18","$19","$i7}' $m7  > $o7
echo $o7

# for nrmain_cosmos_6.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$4*$5*$10","$11",,,,,,"$12","$13","$14","$15","$16","$17","$18",,"$i6}' $m6  > $o6
echo $o6

# for nrmain_cosmos_5.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5","$6","$7","$8",,"$9","$4*$5*$10","$11",,,,,,"$12","$13","$14","$15","$16","$17",,,"$i5}' $m5  > $o5
echo $o5

# for nrmain_cosmos_4.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5",,"$6","$7",,"$8","$4*$5*$9","$10",,,,,,"$11","$12","$13","$14","$15","$16",,,"$i4}' $m4  > $o4
echo $o4

# for nrmain_cosmos_3.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5",,"$6",1.0,,"$7","$4*$5*$8",,,,,,,"$9","$10",,,"$11","$12",,,"$i3}' $m3  > $o3
echo $o3

# for nrmain_cosmos_2.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5",,"$6",1.0,,"$7","$4*$5*$8",,,,,,,"$9","$10",,,"$11",,,,"$i2}' $m2  > $o2
echo $o2

# for nrmain_cosmos.dat
awk -F "," '{print $1","$2",0,"$3","$4**3*$5",0,"$6",1,,"$7","$10",,"$8","$9",,,,"$11",,,,,,,,"$i1}' $m1  > $o1
echo $o1

cat $o1 $o2 $o3 $o4 $o5 $o6 $o7 $o7l $o8 > $o
echo "combined into:"
echo $o
