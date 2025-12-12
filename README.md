# Library of C functions

## Author

- Zdenek Futera

- University of South Bohemia  
      Faculty of Science  
      Branisovska 1760  
      370 07 Ceske Budejovice  
      Czech Republic

- <zfutera@prf.jcu.cz>


## Content

The library is divided to 4 parts:

1. #### General functions used in most of the C codes (*CMN*)

    * Command line argumens (*arg.h*)
    * General 2D and 3D dynamically allocated arrays (*array.h*)
    * Configuration files (*config.h*)
    * Directories (*dir.h*)
    * Files (*file.h*)
    * Linear linked lists (*list.h*)
    * Math functions (*math.h*)
    * Matrices (*matrix.h*)
    * Interactive menus (*menu.h*)
    * Error and warning messages (*message.h*)
    * Pipes (*pipe.h*)
    * Printing (*print.h*)
    * Progress bars (*progress.h*)
    * Quaternions (*quaternion.h*)
    * FIFO queues based on the linear linked lists (*queue.h*)
    * Random numbers (*random.h*)
    * LIFO stacks based on the linear linked lists (*stack.h*)
    * Dynamically allocated character strings (*string.h*)
    * Time and date (*time.h*)
    * Trees (*tree.h*)
    * Data types (*types.h*)
    * Physical units (*units.h*)
    * Vectors (*vector.h*)

2. #### Quantum chemistry

    * Basis sets (*basis.h*)
    * Gaussian type orbitals and their integrals (*gto.h*)
    * Pseudopotentials (*pseudo.h*)

3. #### Molecular structures and formats

    * Amino acids (*aacid.h*)
    * Atomic properties (*atom.h*)
    * Crystallographic cells (*cell.h*)
    * Inter-particle distances (*distance.h*)
    * Internal coordinates (*icoord.h*)
    * Molecular formats (*molec.h*)
    * Periodic boundary conditions (*pbc.h*)

4. #### External program formats

    * Amber molecular-mechanical simulation package (*amber.h*)
    * Bibliographic TeX formats (*bibtex.h*)
    * CP2K program package (*cp2k.h*)
    * CPMD density-functional-theory code (*cpmd.h*)
    * DL_POLY molecular-mechanical simulation package (*dlpoly.h*)
    * Gaussian quantum-chemistry code (*gauss.h*)
    * Gromacs molecular-mechanical simulation package (*gromacs.h*)
    * LAMMPS simulation package (*lammps.h*)
    * Molden visualizaton program (*molden.h*)


## Compilation

The library is compiled by GNU compiler (`gcc`) using the provided Makefile.

To compile the library type

```
    make
```

The compiled static libraries (`libcmm.a`, `libqmc.a`, `libmol.a`, and `libprg.a`) can be found in the automatically-created `lib` directory.


## Usage

To compile a C code *program.c* with the library type

```
    gcc -Wall program.c -o program.exe
        -Ilibcutils/include -Llibcutils/lib
        -lprg -lmol -lqmc -lcmn -lm
```


## License

The library is disributed under the [MIT license](LICENSE.md).
