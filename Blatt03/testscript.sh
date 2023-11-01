#!/bin/bash
./partdiff-seq 1 2 64 1 2 10240 | grep Berechnungszeit >> 10240.txt
./partdiff-seq 1 2 64 2 2 5120 | grep Berechnungszeit >> 5120.txt
./partdiff-seq 1 2 64 1 2 10240 | grep Berechnungszeit >> 10240.txt
./partdiff-seq 1 2 64 2 2 5120 | grep Berechnungszeit >> 5120.txt
./partdiff-seq 1 2 64 1 2 10240 | grep Berechnungszeit >> 10240.txt
./partdiff-seq 1 2 64 2 2 5120 | grep Berechnungszeit >> 5120.txt
./partdiff-seq 1 2 64 1 2 10240 | grep Berechnungszeit >> 10240.txt
./partdiff-seq 1 2 64 2 2 5120 | grep Berechnungszeit >> 5120.txt
./partdiff-seq 1 2 64 1 2 10240 | grep Berechnungszeit >> 10240.txt
./partdiff-seq 1 2 64 2 2 5120 | grep Berechnungszeit >> 5120.txt
./partdiff-seq 1 2 64 1 2 10240 | grep Berechnungszeit >> 10240.txt
./partdiff-seq 1 2 64 2 2 5120 | grep Berechnungszeit >> 5120.txt
echo done >> 10240.txt
echo done >> 5120.txt