# ezFCF 
## v.1.2 (formerly known as ezSpectrum) 2022
Authors:  Pawel Wojcik, Samer Gozem, Vadim Mozhayskiy, and Anna Krylov
If you are using mac type `make -f make.mac`.

The program is available for download at iopenshell.usc.edu/downloads/

### Compile
To compile go to the `ezspectrum_code` directory and type `make`.

ezFCF uses Armadillo which, in turn, has the optionally dependencies: OpenBLAS, 
LAPACK, ARPACK, and SuperLU.

It's recommended to install these dependecines using the system package 
manager, for example
```bash
apt install libopenblas-openmp-dev libarmadillo-dev
```

If you decide to install armadillo manually and in a non-default location (i.e.
not in `/usr/`), make sure to adjust the include, -I, and library, -L,
locations in the makefile.

If you use ezFCF on a compute cluster you should be able to load most of the
optional dependecines as modules. After you load this modules, make sure that
the corresponding include, and library directories are visible to the compiler
and linker, as they might not be added by default to the `CPATH`,
`LIBRARY_PATH`, and `LD_LIBRARY_PATH`.

### Linux (GNU) installation

On linux (GNU make)
```bash
make
make install
``` 
will install the program and its data under `/usr/local`. The global install
requires root privileges. Local installation is possible using the `prefix`
variable, e.g.
```bash
make prefix=~/.local
```
Using the `prefix` should be connected with setting the environmental variable 
`EZFCF_ROOT`, e.g., by adding in your `.bashrc` or submit scripts line like
```bash
export EZFCF_ROOT=~/.local
```

Uninstall also needs to use the prefix
```
make uninstall prefix=~/.local
```

