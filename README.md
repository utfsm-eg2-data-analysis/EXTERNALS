EXTERNALS
=========

Package that allows the calculation of inclusive electron radiative corrections in deep inelastic regime.

**Author:** Dave Gaskell.

### Structure

1. `INP/clasd2.inp` - main input file that points to the other required inputs.

2. `RUNPLAN/clas_kin.inp` - this file contains the kinematics at which to calculate the cross section and radiative corrections. Note that it's pretty sensitive to formatting.

3. `TARG/targ.D2tuna` - this file contains information about the target (Z, A, geometry, model to use, etc.). Right now, the only variables that matter are Z and A since the program only calculates the internal corrections. The model choice is also hardwired for now to use F1F209 from Peter Bosted.

4. `OUT/clasd2_details.out` - this gives a more detailed output than the summary table above.

5. To run, you just use the script `run_extern` in the top level directory.
The usage is: `./run_extern <input-file>`, leave the ".inp" off the master input file. 

### Compilation

To compile on the UTFSM cluster:

```
source set_env.sh
make
```

To compile on your local machine:

```
sudo apt install cernlib
make
```

### Execution

The repository is prepared to run the code for each different target:

```
./run_extern clasd2
./run_extern clasC12
./run_extern clasFe56
./run_extern clasPb208
```

Output files will be located in `OUT/`.

### References

1. **L. W. Mo, Y. S. Tsai.** Radiative corrections to elastic and inelastic $ep$ and $\mu p$ scattering. **Rev. of Mod. Phys. 41, 1 (1969)**

2. **Y. S. Tsai.** Radiative corrections to electron scatterings. **SLAC-PUB-848 (1971)**