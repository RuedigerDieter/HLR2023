
#!/bin/sh
make clean
make

mpirun -n 1 ./partdiff 1 1 64 2 2 100 > 1.txt
mpirun -n 2 ./partdiff 1 1 64 2 2 100 > 2.txt
mpirun -n 3 ./partdiff 1 1 64 2 2 100 > 3.txt
mpirun -n 4 ./partdiff 1 1 64 2 2 100 > 4.txt
mpirun -n 5 ./partdiff 1 1 64 2 2 100 > 5.txt
mpirun -n 6 ./partdiff 1 1 64 2 2 100 > 6.txt

