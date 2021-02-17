home=$(pwd)

for mol in mol*/
do
  cd ${mol} || continue
  for d in mol*/;
  do
    cd ${d} || continue 
#    for f in *.chk.chk; do echo ${f} & done
    for f in *.chk.chk; do mv -- "$f" "${f%.chk.chk}.chk" & done
    cd ..
  done
  cd ..
done

