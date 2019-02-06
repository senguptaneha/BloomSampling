#./BSTtesting M nSets inputFile k T desAcc

M=100000
MS=100K
for a in 0.5 0.6 0.7 0.8 0.9 0.99
do
	./BSTtesting $M 100 ../SimulatedDatasets/${MS}/Uniform10K_100 3 0 $a
done

M=1000000
MS=1M
for a in 0.5 0.6 0.7 0.8 0.9 0.99
do
	./BSTtesting $M 100 ../SimulatedDatasets/${MS}/Uniform10K_100 3 0 $a
done

M=10000000
MS=10M
for a in 0.5 0.6 0.7 0.8 0.9 0.99
do
	./BSTtesting $M 100 ../SimulatedDatasets/${MS}/Uniform10K_100 3 0 $a
done
