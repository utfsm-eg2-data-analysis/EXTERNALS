EXTERNALS
===========

## Structure

1. `INP/clasd2.inp` - this file is the "master" input file that basically just points to the other required inputs.

2. `RUNPLAN/clas_kin.inp` - this file contains the kinematics at which to calculate the cross section and RC. Note that it's pretty sensitive to formatting.

3. `TARG/targ.D2tuna` - this file contains a lot of info about the target (Z, A, geometry, model to use, etc.). 
   Right now, the only variables that matter are Z and A since the program only calculates the internal corrections. 
   The model choice is also hardwired for now to use **F1F209** from **Peter Bosted**.

4. `OUT/clasd2_details.out` - this gives more detailed output than the summary table above.

The `Coulomb/` directory is supplementary, contains example of modified `externals_all.f` - the same code `rc_Coulomb.C` is in `RC/` and in `Coulomb/RC_CC`.

## Compilation

To compile at **UTFSM Cluster**, run:

1. `source set_env.sh`

2. `make`

## Execution

To run on Deuterium target:
```
./run_extern clasd2
```

To run on Lead target:
```
./run_extern clasPb208
```

## References

1. *L. W. Mo, Y. S. Tsai.* Radiative corrections to elastic and inelastic ep and µp scattering. Rev. of Mod. Phys. **41**, 1 (1969)

2. *Y. S. Tsai.* Radiative corrections to electron scatterings. **SLAC-PUB-848** (1971)

Also, there is an excellent lecture prepared by [**Taisiya Mineeva**](mailto:taya.mineeva@gmail.com) at `QED-RC.pdf`.
