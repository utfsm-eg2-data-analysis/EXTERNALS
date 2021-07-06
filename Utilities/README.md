EXTERNALS/Utilities
===================

1. Modify `include/Binning.hxx` with your binning of preference. It must be 2D binning.

2. Execute `bin/GetBinning`. A file named `binning.csv` will be created.

3. Execute `bin/GetCentroids` for each target (with the option `-t<TARGET>`). Files named `centroids_<TARGET>.txt` will be created. Here, `<TARGET>` can be `C`, `Fe`, `Pb` or `D`.

4. Copy the centroids text files into the directory `RUNPLAN/`. There, modify line 24 of `make_plan.py` to read one of the centroids files. Then execute,

    ```
    python make_plan.py > clas_kin.inp
    ```

5. After compiling `EXTERNALS`, execute it for each target. Here, `<TARGET*>` can be `C12`, `Fe56`, `Pb208` or `d2`.

    ```
    ./run_extern clas<TARGET*>
    ```

6. Make sure to repeat steps **4.** and **5.** for each target.
