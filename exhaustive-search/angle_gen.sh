rm run_file
for power in 0.001 0.0001; do
#for power in 1 0.1 0.01 0.001; do
for j in $(seq 1 1 100);
do
angle=`awk 'BEGIN { printf("%.5f\n", (rand()-0.5) * 3.1415926535897) }'`
#echo $angle, $power
	echo "./process $angle $power >> data/out_a${angle}_e${power}" >> run_file
done
done
chmod +x run_file
