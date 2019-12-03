# Possible extensions

## Computational side

- parallelize (MPI) the routines that handle the nanocarrier. This will allow for linear scaling of the code with increase in the number of particles. The parallelization may also be required to build cellular scale models for mechanotransduction and cell spreading.
- include an option to make the membrane into a closed geometry
- anisotropic particles
- flexible receptors (how does change in receptor entropy upon binding promotes cell motility.

  #### Scientific side

- appropriate models for the cortical and bulk actin

- models for the microtubules
- models for feedback from signaling networks
- models for cell-ECM interactions
- molecular study of the diffusion of the cytosolic side of the protein
- molecular study of how curvature remodeling protein organize in response to an applied stress (tether pulling)
- Talk to Irfan about how you can get access to the large datasets that would allow you to learn and contribute to aspects of data mining.

## Knowledge side

- biophysics and soft matter
- non-equilibrium hydrodynamics
- cell biology
- computations
- active matter -

## parallelization (2017-02-10 13:52)

- `OPENMP` and `CUDA` based parallelization for NC translation and rotation moves
- `MPI` parallelization for handling individual nanoparticles, bonding
- spatial decomposition of the membrane substrate
-
