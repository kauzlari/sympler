# Define source files to be compiled
set(basic
basic/cell.cpp
basic/colour_pair.cpp
basic/controller.cpp
basic/data_format.cpp
basic/energy_pump.cpp
basic/function_arbitrary.cpp
basic/function_fixed.cpp
basic/function_pair.cpp
basic/function_pair_fixed.cpp
basic/function_particle.cpp
basic/function.cpp
basic/general.cpp
basic/help.cpp
basic/help_system.cpp
basic/inlet_cell.cpp
basic/linked_list_creator.cpp
basic/manager_cell.cpp
basic/MArray2D.cpp
basic/node.cpp
basic/node_many_children.cpp
basic/node_one_child.cpp
basic/pair_creator.cpp
basic/pairdist.cpp
basic/pair_list.cpp
basic/particle.cpp
basic/particle_list.cpp
basic/phase.cpp
basic/property_list.cpp
basic/random.cpp
basic/simulation.cpp
basic/threads.cpp
basic/verlet_creator.cpp
basic/vtk_support.cpp
basic/write_restart_file.cpp
)

set(boundary
boundary/boundary.cpp
boundary/boundary_arbitrary.cpp
boundary/boundary_bubble_jet.cpp
boundary/boundary_couette.cpp
boundary/boundary_diffusor.cpp
boundary/boundary_cuboid.cpp
boundary/boundary_obstacle.cpp
boundary/boundary_step.cpp
boundary/boundary_stl.cpp
boundary/boundary_stl_periodic.cpp
boundary/boundary_with_inlet.cpp
)

set(calculator
calculator/calc.cpp
)

set(force
force/connector_basic.cpp
force/gen_connector.cpp
force/gen_f.cpp
force/gen_f_solo.cpp
force/gen_quintet.cpp
force/gen_triplet.cpp
force/forces_support.cpp
force/f_pair_arbitrary_BC.cpp
force/f_particle.cpp
force/f_wall_dpd.cpp
force/f_angular.cpp
force/f_pair_arbitrary.cpp
force/f_particle_scalar.cpp
force/f_wall_repulsion.cpp
force/f_centrifugal.cpp
force/f_pair_arbitrary_wf.cpp
force/f_particle_tensor.cpp
force/f_with_rng.cpp
force/f_coriolis.cpp
force/f_pair.cpp
force/f_particle_vector.cpp
force/f_curvature.cpp
force/f_pair_scalar.cpp
force/f_particle_vector_matrix.cpp
force/f_dpd.cpp
force/f_pair_tensor.cpp
force/f_particle_vector_rand_matrix.cpp
force/f_dpde.cpp
force/f_pair_vector.cpp
force/f_particle_vels.cpp
force/f_generic.cpp
force/f_pair_vels_BC.cpp
force/f_particle_vels_matrix.cpp
force/f_kinetic.cpp
force/f_pair_vels.cpp
force/f_rand.cpp
force/f_mdpd.cpp
force/f_pair_vels_wf.cpp
force/f_rand_var_noise.cpp
force/f_pair_wf.cpp
force/f_specific.cpp
force/ie_heat_conduction.cpp
force/ie_wall_heat_conduction.cpp
force/lennard_jones.cpp
force/lennard_jones_vc.cpp
)


set(function_parser
function_parser/binary_operators.cpp
function_parser/fp_scalar.cpp
function_parser/fp_tensor.cpp
function_parser/fp_vector.cpp
function_parser/function_compiler.cpp
function_parser/function_node.cpp
function_parser/function_parser.cpp
function_parser/typed_value.cpp
function_parser/unary_functions.cpp
function_parser/unary_operators.cpp
function_parser/variant.cpp
)

set(geometry
geometry/geometric_algorithms.cpp
geometry/geometric_primitives.cpp
geometry/vertex_list.cpp
geometry/wall.cpp
geometry/wall_container.cpp
geometry/wall_triangle.cpp
geometry/stl/stl_ascii.cpp
geometry/stl/stl_binary.cpp
)

set(integrator
integrator/integrator.cpp
integrator/integrator_energy.cpp
integrator/integrator_static.cpp
integrator/integrator_velocity_verlet_disp.cpp
integrator/integrator_lse.cpp
integrator/integrator_static_lse.cpp
integrator/integrator_velocity_verlet_disp_x.cpp
integrator/integrator_omelyan.cpp
integrator/integrator_tensor.cpp
integrator/integrator_velocity_verlet_pressure.cpp
integrator/integrator_omelyan_NR.cpp
integrator/integrator_tensor_lambda.cpp
integrator/integrator_velocity_verlet_symbols.cpp
integrator/integrator_position.cpp
integrator/integrator_pos_vel_step2.cpp
integrator/integrator_vector.cpp
integrator/integrator_velocity_verlet_x.cpp
integrator/integrator_scalar.cpp
integrator/integrator_scalar_step2.cpp
integrator/integrator_shift.cpp
integrator/integrator_vector_lambda.cpp
integrator/integrator_velocity_verlet_x_full.cpp
integrator/integrator_scalar_lambda.cpp
integrator/integrator_velocity_verlet.cpp
integrator/integrator_ISPH_const_rho_FDRHS.cpp
integrator/integrator_ISPH_const_rho.cpp
)

set(meter
meter/grid_averager_structured.cpp
meter/grid_averager.cpp
meter/grid_averager_circular.cpp
meter/grid_meter.cpp
meter/meter.cpp
meter/meter_autocorrelation_vector.cpp
meter/meter_autocorrelation_vector_t.cpp
meter/meter_average.cpp
meter/meter_bonded_crosscorrelation_vector.cpp
meter/meter_distribution.cpp
meter/meter_grouped.cpp
meter/meter_live.cpp
meter/meter_live_coloured.cpp
meter/meter_live_geomview.cpp
meter/meter_pair_distribution.cpp
meter/meter_pair_distribution_with_walls.cpp
meter/meter_pos_vel.cpp
meter/meter_relative_velocity.cpp
meter/grid_meter/gm_array_attr.cpp
meter/grid_meter/gm_density.cpp
meter/grid_meter/gm_ekin.cpp
meter/grid_meter/gm_local_density.cpp
meter/grid_meter/gm_pressure.cpp
meter/grid_meter/gm_scalar.cpp
meter/grid_meter/gm_scalar_times_scalar_sum.cpp
meter/grid_meter/gm_scalar_fourier_cos.cpp
meter/grid_meter/gm_scalar_fourier_sin.cpp
meter/grid_meter/gm_scalar_func.cpp
meter/grid_meter/gm_scalar_sum.cpp
meter/grid_meter/gm_self_diffusion.cpp
meter/grid_meter/gm_surface_tension.cpp
meter/grid_meter/gm_temperature.cpp
meter/grid_meter/gm_temperature_ie.cpp
meter/grid_meter/gm_vector.cpp
meter/grid_meter/gm_velocity.cpp
meter/grid_meter/gm_velocity_moment.cpp
meter/grid_meter/gm_volume_fraction.cpp
)

set(particle_creator
particle_creator/particle_creator.cpp
particle_creator/pc_connector.cpp
particle_creator/pc_file.cpp
particle_creator/pc_free.cpp
particle_creator/pc_free_f.cpp
particle_creator/pc_free_pcalc.cpp
particle_creator/pc_inlet.cpp
particle_creator/pc_lattice.cpp
particle_creator/pc_lattice_bcc.cpp
particle_creator/pc_lattice_frozen.cpp
particle_creator/pc_lse.cpp
particle_creator/pc_random.cpp
particle_creator/pc_static.cpp
particle_creator/pc_tube.cpp
particle_creator/pc_wall.cpp
particle_creator/pc_wall_textured.cpp
particle_creator/pc_with_rng.cpp
)

set(postprocessor
postprocessor/first_harmonic.cpp
postprocessor/linear_regression.cpp
postprocessor/output.cpp
postprocessor/output_dcd.cpp
postprocessor/output_file.cpp
postprocessor/output_mathematica.cpp
postprocessor/output_pdb.cpp
postprocessor/output_vtk.cpp
postprocessor/postprocessor.cpp
)

set(reflector
reflector/reflector.cpp
reflector/reflector_bounce_back.cpp
reflector/reflector_stochastic.cpp
reflector/reflector_thermalize_ie.cpp
reflector/reflector_mirror.cpp
reflector/reflector_thermalize.cpp
reflector/reflector_with_rng.cpp
)

set(callable
callable/apply_scalar_field.cpp
callable/apply_vector_field.cpp
callable/apply_vector_field_file.cpp
callable/apply_vel_field.cpp
callable/callable.cpp
callable/cbl_pair_particle_arbitrary.cpp
callable/cbl_pair_particle_tensor.cpp
callable/cbl_pair_particle_vector.cpp
callable/shift_particle.cpp
callable/triangle_CGMD_interpolation2points.cpp
callable/triangle_CGMD_interpolation.cpp
callable/thermostat.cpp
callable/thermostat_energy_rescaling.cpp
callable/thermostat_la_ec.cpp
callable/thermostat_wall_vel.cpp
callable/thermostat_la_with_weight_ec.cpp
callable/thermostat_wall_vel_ec.cpp
callable/thermostat_peters.cpp
callable/thermostat_wall_vel_iso.cpp
callable/thermostat_peters_ec.cpp
callable/thermostat_with_rng.cpp
callable/thermostat_peters_iso.cpp
callable/thermostat_la.cpp
callable/thermostat_vels.cpp
callable/vel_constraints.cpp
)

set(symbol
symbol/symbol.cpp
symbol/val_calculator.cpp
symbol/freesteam/pca_iapws-if97_1var.cpp
symbol/freesteam/pca_iapws-if97_2var.cpp
symbol/freesteam/pca_iapws-if97_cp.cpp
symbol/freesteam/pca_iapws-if97_eta.cpp
symbol/freesteam/pca_iapws-if97_kappa.cpp
symbol/freesteam/pca_iapws-if97_p.cpp
symbol/freesteam/pca_iapws-if97_rho.cpp
symbol/freesteam/pca_iapws-if97.cpp
symbol/particle_cache/particle_cache.cpp
symbol/particle_cache/particle_cache_arbitrary.cpp
symbol/particle_cache/particle_rand_norm_scalar.cpp
symbol/particle_cache/particle_rand_norm_vector.cpp
symbol/particle_cache/particle_scalar.cpp
symbol/particle_cache/particle_tensor.cpp
symbol/particle_cache/particle_vector.cpp
symbol/particle_cache/particle_vels.cpp
symbol/particle_cache/particle_cache_arb_rng.cpp
symbol/particle_cache/pca_2nd_sph_deriv_corr.cpp
symbol/particle_cache/pca_density_0oc.cpp
symbol/particle_cache/pca_density_self_contribution.cpp
symbol/particle_cache/pca_eigensystem.cpp
symbol/particle_cache/pca_energy_entropy.cpp
symbol/particle_cache/pca_matrix_inverse.cpp
symbol/particle_cache/pca_volume_self_contribution.cpp
symbol/particle_cache/symbol_f_particle_arbitrary.cpp
symbol/particle_cache/symbol_f_particle_scalar.cpp
symbol/particle_cache/symbol_f_particle_vels.cpp
symbol/particle_cache/transfer_particle_vector.cpp
symbol/val_calculator_pair/bonded_pair_scalar.cpp
symbol/val_calculator_pair/bonded_pair_vector.cpp
symbol/val_calculator_pair/bonded_pair_arbitrary.cpp
symbol/val_calculator_pair/pair_tensor.cpp
symbol/val_calculator_pair/pair_arbitrary.cpp
symbol/val_calculator_pair/pair_scalar.cpp
symbol/val_calculator_pair/pair_vector.cpp
symbol/val_calculator_pair/pair_rand_arbitrary.cpp
symbol/val_calculator_pair/pair_rand_scalar.cpp
symbol/val_calculator_pair/pair_rand_tensor.cpp
symbol/val_calculator_pair/pair_rand_vector.cpp
symbol/val_calculator_pair/val_calculator_pair.cpp
symbol/val_calculator_pair/val_calculator_BC.cpp
symbol/val_calculator_pair/val_calculator_dirichlet_BC_arbitrary.cpp
symbol/val_calculator_pair/val_calculator_dirichlet_BC_scalar.cpp
symbol/val_calculator_pair/val_calculator_dirichlet_BC_vector.cpp
symbol/val_calculator_pair/val_calculator_dirichlet_BC_vels.cpp
symbol/val_calculator_pair/val_calculator_kernel.cpp
symbol/val_calculator_pair/val_calculator_neg_dkernel_divr.cpp
symbol/val_calculator_part/bonded_pair_particle_arbitrary.cpp
symbol/val_calculator_part/bonded_pair_particle_calc.cpp
symbol/val_calculator_part/bonded_pair_particle_one_noise_vector.cpp
symbol/val_calculator_part/bonded_pair_particle_scalar.cpp
symbol/val_calculator_part/bonded_pair_particle_tensor_noise_vector.cpp
symbol/val_calculator_part/bonded_pair_particle_vector.cpp
symbol/val_calculator_part/pair_particle_scalar.cpp
symbol/val_calculator_part/pair_particle_tensor.cpp
symbol/val_calculator_part/pair_particle_vector.cpp
symbol/val_calculator_part/val_calculator_arbitrary.cpp
symbol/val_calculator_part/val_calculator_rho.cpp
symbol/val_calculator_part/val_calculator_volume.cpp
symbol/val_calculator_part/val_calculator_part.cpp
symbol/val_calculator_part/non_bonded_pair_particle_calculator.cpp
symbol/other/triplet_calculator.cpp
symbol/other/triplet_calc_central_part_scalar.cpp
symbol/other/triplet_calc_angular_f.cpp
symbol/other/triplet_calc_angular_falpha.cpp
symbol/other/triplet_calc_angular_dt2f.cpp
symbol/other/quintet_calculator.cpp
symbol/other/quintet_calc_curvature_f.cpp
symbol/other/quintet_calc_curvature.cpp
)

set(tools
tools/vtk2profile.cpp)

set(weighting_function
weighting_function/weighting_function.cpp
weighting_function/wf_input.cpp
weighting_function/wf_lucy_0th_order.cpp
weighting_function/wf_lucy.cpp
weighting_function/wf_linear.cpp
weighting_function/wf_lucy_1st_order.cpp
weighting_function/wf_square.cpp
)


# Create library of source files
add_library(basic ${basic})
add_library(boundary ${boundary})
add_library(calculator ${calculator})
add_library(force ${force})
add_library(function_parser ${function_parser})
add_library(geometry ${geometry})
add_library(integrator ${integrator})
add_library(meter ${meter})
add_library(particle_creator ${particle_creator})
add_library(postprocessor ${postprocessor})
add_library(reflector ${reflector})
add_library(callable ${callable})
add_library(symbol ${symbol})
add_library(weighting_function ${weighting_function})
target_link_libraries(geometry basic boundary calculator force function_parser integrator meter particle_creator postprocessor reflector callable symbol weighting_function)
set(src_lib basic boundary calculator force function_parser geometry integrator meter particle_creator postprocessor reflector callable symbol weighting_function)


# Create executable sympler
add_executable(sympler
main.cpp
)

# Link libraries to executable sympler
target_link_libraries(sympler
${src_lib}
${link_lib}
)


# Create executable vtk2profile
add_executable(vtk2profile
tools/vtk2profile.cpp
)

# Link libraries to executable vtk2profile
target_link_libraries(vtk2profile
${src_lib}
)


