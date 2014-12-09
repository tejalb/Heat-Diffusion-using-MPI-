#!/bin/bash
mpirun -np $1 ./heat_mpi $2
rm -rf map_mpi_$1_$2.txt
for ((i = 0; i < $1; i++));do cat map_mpi$i.txt >> map_mpi_$1_$2.txt;done

gnuplot<<-EOF
set term png
set output 'plot.png'
plot 'map_mpi_$1_$2.txt' using 1:2:3 with image
EOF


rm -rf map_mpi*
