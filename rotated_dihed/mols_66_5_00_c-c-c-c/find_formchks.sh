home=$(pwd)
module load ccs/gaussian/g16-A.03/g16-sandybridge
mkdir $home/formchks

for m in mol*/
do
  cd ${m}
  formchk *.chk
  cp *.fchk  $home/formchks
  cd ..
done

