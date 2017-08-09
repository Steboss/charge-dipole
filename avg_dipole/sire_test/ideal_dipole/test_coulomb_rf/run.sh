#!/bin/bash


cd 30angstrom 

~/sire_LRC_coul/sire.app/bin/python rotation.py system.prmtop system.rst7

wait

cd ../50angstrom

~/sire_LRC_coul/sire.app/bin/python rotation.py system.prmtop system.rst7
wait

cd ../80angstrom

~/sire_LRC_coul/sire.app/bin/python rotation.py system.prmtop system.rst7
wait

cd ../
