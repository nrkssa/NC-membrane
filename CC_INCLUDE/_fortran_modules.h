#ifndef __FORTRAN_MODULES_H
#define __FORTRAN_MODULES_H
extern "C" void __module_mcsmoves_MOD_membrane_montecarlo();                         // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_mcsmoves_MOD_membrane_montecarlo_biased();                  // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_mcsmoves_MOD_link_flip();                                   // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_mcsmoves_MOD_link_flip_biased();                            // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_mcsmoves_MOD_vertex_move();                                 // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_mcsmoves_MOD_vertex_move_biased();                          // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_mcsmoves_MOD_antigen_diffusion_on_triangle();               // declare the fortran subroutine performing the Monte Carlo moves for the membrane
extern "C" void __module_makesurface_MOD_makesurface();
extern "C" void __module_curvcalc_MOD_compute_householdermatrix(int *, char *);      // Retrieve the HouseHolder Matrix for bond formation and breakage
extern "C" void __module_curvcalc_MOD_compute_antigen_orientation(int *);            // compute the polar and azimuthal angle for a chosen antigen
extern "C" void __module_writedata_MOD_write_membrane_xyz(int *,int *,int *);        // print membrane coordinates to a file
extern "C" void __module_writedata_MOD_write_membrane_area(int *, int *, int *);     // print membrane area to a file
extern "C" void __module_writedata_MOD_membrane_dump1(int *);                        // Dump the membrane files 
extern "C" void __module_writedata_MOD_membrane_dump2(int *, int *);                 // Dump the membrane files 
extern "C" void __module_mcsmoves_MOD_print_triangle_vert(int *);                  
extern "C" void __module_mcsmoves_MOD_print_image_vert(int *);                  
extern "C" int __module_mcsmoves_MOD_selfavoidance_allantigens();
extern "C" void __module_datastruct_MOD_interactive_shell();
extern "C" void __module_mcsmoves_MOD_antigen_hopping_on_vertex();
extern "C" void __module_mcsmoves_MOD_compute_flexed_base_tip_coord(int *);          // computes the base and tip coordinates of an antigen
extern "C" double __module_datastruct_MOD_compute_mean_membrane_z();                 // Compute mean height (z-value) of the membrane      
extern "C" void __module_datastruct_MOD_compute_mean_membrane_curvature(int *);      // Compute mean membrane curvature under the nanocarrier
extern "C" void __module_datastruct_MOD_store_restore_shadowvertices(int *,int *);   // Stores and Restores shadow vertices during NC move 
extern "C" void __module_datastruct_MOD_compute_nc_shadow_vertices(int *); 
extern "C" void __module_datastruct_MOD_ncshadow_r(int *);                           // compute mean r in the ncshadow.       
extern "C" double __module_datastruct_MOD_total_membrane_energy();
extern "C" double __module_datastruct_MOD_compute_closest_vertex(int*, double*, double*);
extern "C" double __module_datastruct_MOD_get_membrane_area();                       // function to get the membrane area
#endif
