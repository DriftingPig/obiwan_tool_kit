#!/bin/bash -l
DR5ELGs="/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/splitdata/weighed_uniform/data_step4_subset"
ObiwanELGs="/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/splitdata/weighed_uniform/uniform_step5_subset"
dat=".dat"
OutputELGs="/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/BinHist/weighed_uniform/BinHist"
D="D"
R="R"
File_Num=20
Nproc=32
function PushQue {    
Que="$Que $1"
Nrun=$(($Nrun+1))
}
function GenQue {     
OldQue=$Que
Que=""; Nrun=0
for PID in $OldQue; do
if [[ -d /proc/$PID ]]; then
PushQue $PID
fi
done
}
function ChkQue {     
OldQue=$Que
for PID in $OldQue; do
if [[ ! -d /proc/$PID ]] ; then
GenQue; break
fi
done
}

i=0
./ObiwanCorr $DR5ELGs${i}$dat $DR5ELGs${i}$dat 1 $OutputELGs$D${i}$D${i}$dat "w" 0 ${i} ${i}
for((i=1; i<$File_Num; i++)); do #DD auto
	echo $DR5ELGs${i}$dat
	echo $OutputELGs$D${i}$D${i}$dat
	./ObiwanCorr $DR5ELGs${i}$dat $DR5ELGs${i}$dat 1 $OutputELGs$D${i}$D${i}$dat "a" 0 ${i} ${i}&
	PID=$!
	PushQue $PID
	while [[ $Nrun -ge $Nproc ]]; do
		ChkQue
		sleep 1
	done
done

for((i=0; i<$File_Num; i++)); do #DR cross
	for((j=0; j<$File_Num; j++)); do
		./ObiwanCorr $DR5ELGs${i}$dat $ObiwanELGs${j}$dat 0 $OutputELGs$D${i}$R${j}$dat "a" 1 ${i} ${j} &
		PID=$!
		PushQue $PID
	        while [[ $Nrun -ge $Nproc ]]; do
			ChkQue
			sleep 1
		done
	done	
done

for((i=0; i<$File_Num; i++)); do #RR auto
	./ObiwanCorr $ObiwanELGs${i}$dat $ObiwanELGs${i}$dat 1 $OutputELGs$R${i}$R${i}$dat "a" 2 ${i} ${i} &
	PID=$!
	PushQue $PID
	while [[ $Nrun -ge $Nproc ]]; do
		ChkQue
		sleep 1
	done
done

for((i=0; i<$File_Num; i++)); do #DD cross
	for((j=i+1; j<$File_Num; j++)); do
		./ObiwanCorr $DR5ELGs${i}$dat $DR5ELGs${j}$dat 0 $OutputELGs$D${i}$D${j}$dat "a" 0 ${i} ${j} &
		PID=$!
		PushQue $PID
		while [[ $Nrun -ge $Nproc ]]; do
			ChkQue
			sleep 1
		done
	done
done

for((i=0; i<$File_Num; i++)); do #RR cross
	for((j=i+1; j<$File_Num; j++)); do
		./ObiwanCorr $ObiwanELGs${i}$dat $ObiwanELGs${j}$dat 0 $OutputELGs$R${i}$R${j}$dat "a" 2 ${i} ${j} &
		PID=$!
		PushQue $PID
		while [[ $Nrun -ge $Nproc ]]; do
			ChkQue
			sleep 1
		done
	done
done

wait
