#/usr/bin/env bash

gfortran parametersData.F90 -o xpara
./xpara
rm xpara

cd ../1_picks/
gfortran dataSelectSourceGather.F90 -o xgather
./xgather
rm xgather

gfortran receiverGather.F90 -o xrec
./xrec
rm xrec

mv evtstaminmax nevtstaray traveltimeReceiverGathers sources receivers ../../
rm earthquakesOnlyInResearchRegion selectedData src-rec

cd ../0_commandCenter
rm blockSelection earthquakeRegion phasewanted stationRegion timePeriod
