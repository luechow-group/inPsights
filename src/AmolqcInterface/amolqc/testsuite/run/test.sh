#! /bin/bash
echo ""
echo ""
echo "                 AMOLQC tests"
echo ""
FAILED=0
PAR=0
GEN=0
export REF=ref
export AMOLQC=${PWD%testsuite/run}
export AMOLQCRUN=$AMOLQC/src/amolqc
if [ "$1" == "PARALLEL" ] ;then
export PAR=1
echo " Parallel test(2 CPUs)"
fi
if [ "$1" == "PARALLELGEN" ] ;then
export PAR=1
export GEN=1
echo " Parallel test(2 CPUs)"
fi
if [ "$1" == "GEN" ] ;then
export GEN=1
fi
if [ $PAR == 1 ] ;then
  export AMOLQCRUN="mpirun -n 2 $AMOLQC/src/amolqc"
  export REF=ref.par
fi


echo " Current version informations: "
echo "    instalation folder = "  $AMOLQC
a=`cat ../../version`
echo "    version            = " $a
if [ $GEN == 1 ] ; then
  echo " updating the references ..."
  echo "version: " $a >reference
  echo "date:" `date` >>reference
fi
echo ""
echo " Reference informations: "
a=`grep "version" reference`
echo "    version            = " ${a#version:}
a=`grep "date" reference`
echo "    date               = " ${a#date:}
echo ""
fcount=0
n=1
for dir in */; do fcount=$[$fcount+1] ;done

for dir in */
do
   cd $dir
   echo " starting test" $n "out of" $fcount
   source RUN
   n=$[$n+1]
   echo ""
   cd ..
done

if [ $FAILED == 0 ] ; then
  >&2 echo -e " All tests \e[32msuccessfully\e[39m completed."
else
  >&2 echo  -e " At least one test \e[31mfailed \e[39m please check the output carefully."
fi
