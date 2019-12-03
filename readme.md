# Code Revamp

`2017-02-10 10:19` I have been grapling with the problem of not being able to reproduce the correct PMF for the ICAM system. I thought the problem may be due to insufficient sampling and in order to check this I simulated the system with two different biasing potentials and analyzed the results. I found that the results do not change. 

However, the code had the wrong expression for the `max_AB_ANT_distance=2*antigen_length`. This is not correct since this calculation is only performed between the the tip of the antigen and the tip of the antibody. 


`2017-02-02 10:27` A couple of bugs has been fixed in `nanocarrier.cpp`

```
- `creat_antibody()`: The comparison variable has been changed from `num_ab` to `num_abtype`
- `parse_abparameters(int index)`: The arguments for `parsval1, parsval2, parsval3` has been changed from 2 to `this->num_abtype`.
```

`2016-06-23 10:28` The code is now setup in a modular fashion with each component having their respective classes. The system members have three classes `nanocarrier.cpp` for NCs, `receptor.cpp` for antigens and `membrane.cpp` for the membrane. They interact with each other through specific call.

## `nanocarrier.cpp`

- A private class `_NANOCARRIER` which will be instantiated once for each nanocarrier.
- A public class `_nanocarrier` which will be used for performing all sorts of operations.

### setting up the system

1) specify the total number of nanocarriers in the system as

- `num_nanocarriers xx` in `CC_INPUT_FILES/system_parameters.inp`

# Using multiple particles in the same code

`2016-02-28`

- An additional routine `__module_datastruct:compute_nc_shadow_vertices_local(ncnum,vert)` has been introduced to locally determine if a moved vertex in subroutine `__module_mcsmoves::Vertex_move_biased()` is still within the shadow of an NC or not. This routine is organized as follows:

```
start routine
 -> compute dx,dy,dz for the given nanocarrier and membrane vertex, and compute the distance
 -> if  (outside cutoff) and (already  in shadow) : rearrange the shadow vertices and remove vert
 -> else if (inside cutoff) and (already not in shadow) : add vert to the shadow vertices
 -> else (inside cutoff) and (already in shadow) : update only pbcx, pbcy
end routine
```

## ------------------------------------------------------------------------------

`2016-02-26`

- The `datawriter.cpp` routine has been modified now to write one single system file that contains all the required information about the NC, MEMB, and Antigens. **However, writing all the objects as a single file `datawriter.cpp::VTK_system()` allows one only to visualize the system but does not allows to visualize the various scalars and vectors stored inside the xml tag `<Piece>`.**
- How to fix the multiple shadow problem? **Each vertex can be associated independently with an NC, and can contribute to the order parameter for different NCs. This has been incorporated by changing**
- ver(i)%shadownc -> ver(i)%shadownc(1:nnanocarrier)
- Does each particle get its shadow correctly? **Yes, each NC has its unique shadow, but the shadow regions may overlap. The overlap doe not affects the shadow calculation in any manner.**

  - In the present implementation a vertex can only remember one shadow it is associated with. This has to be changed since the ICAM molecules at a given region can have equal probability to bind to the ligands on any of the NCs.

- Does the particles get initialized correctly? **Yes, the code is equipped to do so**

## ------------------------------------------------------------------------------

# General Instruction for compilation and running

The makefile can generate two differnt executables `free_ener`, `free_ener_fixed`, and `MEMB_MC`, which respectively perform NC binding in the presence of undulating membranes, NC binding in the presence of fixed membranes, and equilibration of the membrane without the nanocarrier. Makefiles have been provided

```bash
make fe/fixed/memb
```

The executables are setup to run from inside the `RUNDIR`. The provided bash script can automatically do this for you. The runscript should be run as

```
bash runscript.sh free_ener 8 debug
```

Here `free_ener` is the executable to be run on `8` processors in the `debug` mode. In this mode, the code writes a certain number of statistical quantities to the screen which can be used to monitor the status the nanocarrier binding. This is not recommended for production runs which should be performed in the `silent` mode (If an option is not specified silent mode is enforced automatically).

# 24 Sep 15: Notes on the cutoff lengths used in the MC moves

The binding of the nanocarrier is very sensitive to the values of the various cutoff distances employed in the Monte Carlo simulations.

In our case, a number of parameters have been defined in the header of `movement.cpp`, which contains `void _movement :: bond_formation(int timestep)`. Given below are the parameters, their values and the rational for choosing them.

- Distance between the center of mass of the nanocarrier and any point on the membrane:

  `max_nc_memb_distance = pow(sqrt(2.0)*soft_radius,2);`

  This length is used to identify nanocarriers that have at least one vertex on the membrane within this cutoff

- The next distance is to identify all antigens on the membrane that can possibly form bonds with the antibodies on the nanocarrier:

  `max_nc_anttip_distance = pow(soft_radius+2.0*data.getmax_bond_stretch(),2);`

- The maximum distance between an antibody and the tip of an antigen is given by:

  `max_AB_ANT_distance = 2*pow(ant_size,2)*onemSqrt2;`

- The maximum stretch ($\Delta r$) of the Bell bond is set by the relation $lsnk_BT= \frac{1}{2}k_{\rm spring} \Delta r^2$ and this is defined in the code as:

  `Dreaction = pow(data.getmax_bond_stretch(),2);`

All the results I have were generated using these parameters. I have tested alternative values (smaller cutoffs) for these variables and found that the resulting PMF and energies from Thermodynamic integration were significantly different (and also unphysical) from those reported in our articles.
