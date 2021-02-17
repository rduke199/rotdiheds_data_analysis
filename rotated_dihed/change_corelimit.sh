home=$(pwd)

for m in mol*/
do
  cd ${m}
  for d in mol*deg*/;
  do
    cd ${d} || break 
    echo ${d}
    sed -i "s|%nproc=1|%nproc=8|g" *.gjf
    cd ..
  done
  cd ..
done
