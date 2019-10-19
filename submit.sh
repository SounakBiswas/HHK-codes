#! /bin/sh
SEED=11
RANDOM=$SEED
for num in `seq 1 10`
do
for j2 in  0.0000 0.2000 1.0000
do
for size in  8
do
for temp in  0.2000 0.4000 0.6000 0.8000 1.0000
do
sed -e "s/#LX#/${size}/g" gen_global.h > foo1.h 
sed -e "s/#J2#/${j2}/g" foo1.h > foo2.h
sed -e "s/#T#/${temp}/g" foo2.h > foo3.h
sed -e "s/#NUM#/${num}/g" foo3.h > foo4.h
sed -e "s/#SEED#/${RANDOM}/g" foo4.h > global.h
rm foo1.h foo2.h foo3.h foo4.h foo4.h

sed -e "s/#LX#/${size}/g" gen_compile.h > foo1.h 
sed -e "s/#J2#/${j2}/g" foo1.h > foo2.h
sed -e "s/#T#/${temp}/g" foo2.h > foo3.h
sed -e "s/#NUM#/${num}/g" foo3.h > compile.sh


rm foo1.h foo2.h foo3.h

chmod +x compile.sh
./compile.sh 

addqueue -c "5-6 hours,hhkjobs,testing" -n 1 -m 2 ./lx${size}t${temp}j2${j2}.out


done
done
done
done
