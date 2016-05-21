
noinst_LTLIBRARIES = libbasic.la

libbasic_la_SOURCES = \
	bonded_pair_arbitrary.cpp \
	bonded_pair_particle_arbitrary.cpp \
	bonded_pair_particle_calc.cpp \
	bonded_pair_particle_one_noise_vector.cpp \
	bonded_pair_particle_scalar.cpp \
	bonded_pair_particle_tensor_noise_vector.cpp \
	bonded_pair_particle_vector.cpp \
	bonded_pair_scalar.cpp \
	bonded_pair_vector.cpp \
	boundary.cpp \
	callable.cpp \
	cell.cpp \
	colour_pair.cpp \
	controller.cpp \
	data_format.cpp \
	energy_pump.cpp \
	function.cpp \
	function_arbitrary.cpp \
	function_fixed.cpp \
	function_pair.cpp \
	function_pair_fixed.cpp \
	function_particle.cpp \
	gen_f.cpp \
	general.cpp \
	help.cpp \
	help_system.cpp \
	inlet_cell.cpp \
	integrator.cpp \
	linked_list_creator.cpp \
	lin_solve_tensor_vector.cpp \
	manager_cell.cpp \
	MArray2D.cpp \
	meter.cpp \
	node.cpp \
	node_many_children.cpp \
	node_one_child.cpp \
	pair_creator.cpp \
	pairdist.cpp \
	pair_list.cpp \
	pair_particle_scalar.cpp \
	pair_particle_tensor.cpp \
	pair_particle_vector.cpp \
	pair_rand_arbitrary.cpp \
	pair_rand_scalar.cpp \
	pair_rand_tensor.cpp \
	pair_rand_vector.cpp \
	pair_arbitrary.cpp\
	pair_scalar.cpp\
	pair_tensor.cpp\
	pair_vector.cpp\
	particle.cpp \
	particle_cache.cpp \
	particle_cache_arbitrary.cpp \
	particle_cache_arb_rng.cpp \
	particle_creator.cpp \
	particle_list.cpp \
	particle_rand_norm_scalar.cpp \
	particle_rand_norm_vector.cpp \
	particle_scalar.cpp \
	particle_tensor.cpp \
	particle_vector.cpp \
	pca_eigensystem.cpp \
	pca_density_0oc.cpp \
	pca_density_self_contribution.cpp \
	pca_energy_entropy.cpp \
	pca_matrix_inverse.cpp \
	pca_random.cpp \
	pca_volume_self_contribution.cpp \
	phase.cpp \
	postprocessor.cpp \
	property_list.cpp \
	quintet_calculator.cpp \
	quintet_calc_curvature.cpp \
	quintet_calc_curvature_f.cpp \
	random.cpp \
	reflector.cpp \
	simulation.cpp \
	symbol.cpp \
	threads.cpp \
	transfer_particle_vector.cpp \
	triangle_CGMD_interpolation.cpp \
	triangle_CGMD_interpolation2points.cpp \
	triplet_calculator.cpp \
	triplet_calc_angular_dt2f.cpp \
	triplet_calc_angular_f.cpp \
	triplet_calc_angular_falpha.cpp \
	vel_constraints.cpp\
	triplet_calc_central_part_scalar.cpp \
	val_calculator.cpp \
	val_calculator_arbitrary.cpp \
	val_calculator_BC.cpp \
	val_calculator_dirichletBC.cpp \
	val_calculator_dirichlet_BC_vels.cpp \
	val_calculator_kernel.cpp \
	val_calculator_neg_dkernel_divr.cpp \
	val_calculator_rho.cpp \
	val_calculator_volume.cpp \
	verlet_creator.cpp \
	vtk_support.cpp \
	weighting_function.cpp \
	write_restart_file.cpp 

libbasic_la_CXXFLAGS = @GSL_CFLAGS@ @SDL_CFLAGS@ @XML_CFLAGS@

INCLUDES = -I$(top_srcdir)/src/basic/head -I$(top_srcdir)/src/calculator/head -I$(top_srcdir)/src/function_parser/head -I$(top_srcdir)/src/geometry/head -I$(top_srcdir)/src/integrator/head -I$(top_srcdir)/src/meter/head -I$(top_srcdir)/src/particle_creator/head -I$(top_srcdir)/src/reflector/head
