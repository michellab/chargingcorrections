cd output
cd lambda-0.000
~/software/sire.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 0.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-0.000.dat 2> /dev/null

cd ..
cd lambda-1.000
~/software/sire.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 1.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-1.000.dat 2> /dev/null
cd ..
cd ..

wait

# utility script to get final LJ correction term
python output/parselj.py output/freenrg-LJCOR-lam-0.000.dat output/freenrg-LJCOR-lam-1.000.dat > freenrg-LJCOR.dat
