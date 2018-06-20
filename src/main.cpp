/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
 * and others authors stated in the AUTHORS file in the top-level 
 * source directory.
 *
 * SYMPLER is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SYMPLER is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SYMPLER.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Please cite the research papers on SYMPLER in your own publications. 
 * Check out the PUBLICATIONS file in the top-level source directory.
 *
 * You are very welcome to contribute extensions to the code. Please do 
 * so by making a pull request on https://github.com/kauzlari/sympler
 * 
 */


//#include "svnversion.h"
//#ifdef SVN_REV
//#define FULL_VERSION "(revision "SVN_REV")"
//#else
//#define FULL_VERSION "(revision $Revision$)"
//#endif

#include <ctime>
#include <math.h>

#include "simulation.h"
#include "data_format.h"
#include "help_system.h"
#include "reflector.h"

#ifdef HAVE_SDL
  // We need to do this, otherwise SDL doesn't work...
  #include <SDL/SDL.h>
#endif

#ifdef _MPI
  #include <mpi.h>
#endif

using namespace std;

void link_all_properly();

void sayHello();

int main(int argc, char *argv[])
{

  sayHello();

//	srand(17003);
	srand((unsigned)time(NULL));
//  link_all_properly();

	int exitValue = 0;

//#ifdef _MPI
//  MPI_Init(&argc,&argv);
//  MPI_Comm_rank(MPI_COMM_WORLD,&thread_id);
//  MPI_Comm_size(MPI_COMM_WORLD,&n_procs);
//#endif

  DataFormat::alignDataFor(DATA_ALIGNMENT);

  // still the solution from good old subversion times
//   cout << "SYMPLER version " << SVNREV << endl;

  if (argc >= 2) {
    if (argc == 2 && !strcmp(argv[1], "--help"))
      main_help();
    else if (argc == 3)
      help(argv[2]);
    else if (argc == 2) {
      try {
        Simulation aSimulation;
        time_t t0, t1;
        aSimulation.readWithArg(argc, argv);
	clock_t tc = clock();
        time(&t0);
        aSimulation.run();
	MSG_DEBUG("main", "Simulation::run FINISHED");

        time(&t1);

//        MSG_INFO("main", "time: " << static_cast<double>(clock())/CLOCKS_PER_SEC << " seconds");
        MSG_INFO("main", "Simulation took " << difftime(t1, t0) << " seconds (time-function)");
        MSG_INFO("main", "Simulation took " << (clock() - tc)/(double) CLOCKS_PER_SEC << " seconds (clock-function)");

	MSG_INFO("main", "Total number of reflector hits: " << Reflector::s_n_hits);
      } catch (gError &err) {
        cout << "The following ERROR occured: " << endl 
	     << err.message() << endl 
	     << "---------------------" << endl
	     << "Aborted with ERROR-message." << endl
	     << "The ERROR-message might have been lengthy, so please search for the first keyword 'ERROR' occurring above." << endl;
        exitValue = err.exitValue();
      }
    } else
      syntax();
  } else
    syntax();


//#ifdef _MPI
//  MPI::Finalize();
//#endif


  return exitValue;
}

void sayHello() {
  cout << endl << "SYMPLER: SYMbolic ParticLE simulatoR" << endl
       << "Copyright 2002-2018, David Kauzlaric and " << endl
       << "other authors listed in the AUTHORS file." << endl
       << "This program comes with ABSOLUTELY NO WARRANTY;" << endl
       << "for details see the LICENSE file." << endl
       << "This is free software, and you are welcome to redistribute it" << endl
       << "under certain conditions; for details see the LICENSE file." << endl << endl;
}


/* --- Ugly, ugly, ugly --- */
// here we have to include those classes for which the SmartEnum database is used

#include "boundary_bubble_jet.h"
#include "boundary_couette.h"
#include "boundary_cuboid.h"
#include "boundary_diffusor.h"
#include "boundary_obstacle.h"
#include "boundary_stl.h"
#include "boundary_stl_periodic.h"
#include "boundary_step.h"

#include "pair_creator.h"
#include "linked_list_creator.h"
#include "verlet_creator.h"
// #include "vl_yao_creator.h"

#include "connector_basic.h"

#include "ie_heat_conduction.h"
#include "ie_wall_heat_conduction.h"
#include "f_angular.h"
#include "f_curvature.h"
#include "f_centrifugal.h"
#include "f_coriolis.h"
#include "f_dpd.h"
#include "f_dpde.h"
#include "f_generic.h"
// #include "f_grav.h"
#include "f_kinetic.h"
#include "f_mdpd.h"
#include "f_pair_scalar.h"
#include "f_pair_tensor.h"
#include "f_pair_vector.h"
#include "f_pair_vels.h"
#include "f_pair_vels_BC.h"
#include "f_pair_vels_wf.h"
#include "f_particle_scalar.h"
#include "f_particle_vector.h"
#include "f_particle_vector_matrix.h"
#include "f_particle_vector_rand_matrix.h"
#include "f_particle_vels.h"
#include "f_particle_vels_matrix.h"
#include "f_particle_tensor.h"
#include "f_rand.h"
#include "f_specific.h"
#include "f_wall_dpd.h"
#include "f_wall_repulsion.h"
#include "lennard_jones.h"
#include "lennard_jones_vc.h"
// #include "shear_migration.h"
// #include "viscosity_migration.h"
#include "gen_connector.h"
#include "gen_triplet.h"
#include "gen_quintet.h"
#include "gm_array_attr.h"
#include "gm_density.h"
#include "gm_ekin.h"
#include "gm_local_density.h"
#include "gm_pressure.h"
#include "gm_scalar.h"
#include "gm_scalar_fourier_cos.h"
#include "gm_scalar_fourier_sin.h"
#include "gm_scalar_func.h"
#include "gm_scalar_sum.h"
#include "gm_scalar_times_scalar_sum.h"
#include "gm_self_diffusion.h"
// #include "gm_shear_rate.h"
#include "gm_surface_tension.h"
#include "gm_temperature.h"
#include "gm_temperature_ie.h"
#include "gm_vector.h"
#include "gm_velocity.h"
#include "gm_velocity_moment.h"
#include "gm_volume_fraction.h"
#include "gm_vector.h"

// Integrators
#include "integrator_energy.h"
#include "integrator_ISPH_const_rho_FDRHS.h"
#include "integrator_lse.h"
#include "integrator_omelyan.h"
#include "integrator_omelyan_NR.h"
#include "integrator_position.h"
#include "integrator_scalar.h"
#include "integrator_scalar_lambda.h"
#include "integrator_static.h"
#include "integrator_tensor.h"
#include "integrator_velocity_verlet_disp.h"
#include "integrator_velocity_verlet_disp_x.h"
#include "integrator_tensor_lambda.h"
#include "integrator_vector.h"
#include "integrator_vector_lambda.h"
#include "integrator_velocity_verlet.h"
#include "integrator_velocity_verlet_pressure.h"
#include "integrator_velocity_verlet_symbols.h"
#include "integrator_velocity_verlet_x.h"
#include "integrator_velocity_verlet_x_full.h"
#include "integrator_tensor.h"
#include "integrator_static_lse.h"

#include "meter_autocorrelation_vector.h"
#include "meter_autocorrelation_vector_t.h"
#include "meter_average.h"
#include "meter_bonded_crosscorrelation_vector.h"
#include "meter_distribution.h"
// trying to get rid of those Meters
// #include "meter_e_kin.h"
// #include "meter_e_pot.h"
// #include "meter_e_total.h"
#include "meter_live.h"
#include "meter_live_coloured.h"
#include "meter_live_geomview.h"
#include "meter_pair_distribution.h"
#include "meter_pair_distribution_with_walls.h"
#include "meter_pos_vel.h"
#include "meter_relative_velocity.h"
// #include "meter_visco.h"
#include "grid_averager_circular.h"
#include "grid_averager_structured.h"

//#include "first_harmonic.h"
//#include "linear_regression.h"
#include "output_dcd.h"
#include "output_pdb.h"
#include "output_vtk.h"
#include "output_file.h"
#include "output_mathematica.h"

#include "pc_connector.h"
#include "pc_file.h"
#include "pc_inlet.h"
#include "pc_lattice.h"
#include "pc_lattice_bcc.h"
#include "pc_lattice_frozen.h"
#include "pc_lse.h"
#include "pc_random.h"
#include "pc_tube.h"
#include "pc_static.h"
#include "pc_wall.h"
#include "pc_wall_textured.h"

// ParticleCaches
#include "pca_eigensystem.h"
#include "pca_density_0oc.h"
#include "pca_iapws-if97_cp.h"
#include "pca_iapws-if97_eta.h"
#include "pca_iapws-if97_kappa.h"
#include "pca_iapws-if97_p.h"
#include "pca_iapws-if97_rho.h"
#include "pca_matrix_inverse.h"
#include "particle_rand_norm_scalar.h"
#include "particle_rand_norm_vector.h"
#include "particle_scalar.h"
#include "particle_tensor.h"
#include "particle_vector.h"
#include "particle_vels.h"
#include "symbol_f_particle_scalar.h"
#include "symbol_f_particle_vels.h"
// commented out because of complaint about redefinition
// #include "pca_density_self_contribution.h"
// #include "pca_energy_entropy.h"
// #include "pca_shear_self_contribution.h"
// #include "pca_volume_self_contribution.h"

#include "reflector_bounce_back.h"
#include "reflector_mirror.h"
#include "reflector_stochastic.h"
#include "reflector_thermalize.h"
#include "reflector_thermalize_ie.h"


// Callables
#include "apply_vector_field.h"
#include "apply_vector_field_file.h"
#include "apply_vel_field.h"
#include "cbl_pair_particle_tensor.h"
#include "cbl_pair_particle_vector.h"
#include "energy_pump.h"
#include "shift_particle.h"
#include "triangle_CGMD_interpolation.h"
#include "triangle_CGMD_interpolation2points.h"
#include "thermostat_energy_rescaling.h"
#include "thermostat_la.h"
#include "thermostat_la_ec.h"
#include "thermostat_la_with_weight_ec.h"
#include "thermostat_peters_ec.h"
#include "thermostat_peters_iso.h"
#include "thermostat_vels.h"
#include "thermostat_wall_vel.h"
#include "thermostat_wall_vel_ec.h"
#include "thermostat_wall_vel_iso.h"
#include "write_restart_file.h"

#include "transfer_particle_vector.h"

#ifdef HAVE_JAMA_JAMA_LU_H
  #include "vel_constraints.h"
#endif
// ValCalculators
#include "bonded_pair_particle_one_noise_vector.h"
#include "bonded_pair_particle_scalar.h"
#include "bonded_pair_particle_tensor_noise_vector.h"
#include "bonded_pair_particle_vector.h"
#include "bonded_pair_scalar.h"
#include "bonded_pair_vector.h"
#include "pair_particle_scalar.h"
#include "pair_particle_tensor.h"
#include "pair_particle_vector.h"
#include "pair_rand_scalar.h"
#include "pair_rand_tensor.h"
#include "pair_rand_vector.h"
#include "pair_scalar.h"
#include "pair_tensor.h"
#include "pair_vector.h"
#include "triplet_calc_angular_dt2f.h"
#include "triplet_calc_angular_f.h"
#include "triplet_calc_angular_falpha.h"
#include "triplet_calc_central_part_scalar.h"
#include "quintet_calc_curvature.h"
#include "quintet_calc_curvature_f.h"
// #include "val_calculator_dirichlet_BC_scalar.h"
#include "val_calculator_dirichlet_BC_scalar.h"
#include "val_calculator_dirichlet_BC_vels.h"
#include "val_calculator_kernel.h"
#include "val_calculator_neg_dkernel_divr.h"
#include "val_calculator_r_i.h"
#include "val_calculator_r_i6.h"
#include "val_calculator_rho.h"
// commented out because buggy
// #include "val_calculator_symmetry_BC_scalar.h"
#include "val_calculator_volume.h"

#include "wf_input.h"
#include "wf_linear.h"
#include "wf_lucy.h"
#include "wf_lucy_0th_order.h"
#include "wf_lucy_1st_order.h"
#include "wf_square.h"
//#include "wf_tabulated.h"


void link_all_properly()
{
  new BoundaryBubbleJet(NULL);
  new BoundaryCouette(NULL);
  new BoundaryCuboid(NULL);
  new BoundaryDiffusor(NULL);
  new BoundaryObstacle(NULL);
  new BoundarySTL(NULL);
  new BoundaryStlPeriodic(NULL);
  new BoundaryStep(NULL);

  new PairCreator(NULL);
  new LinkedListCreator(NULL);
  new VerletCreator(NULL);
//   new VLYaoCreator(NULL);

  new ConnectBasic(NULL);

  new IEHeatConduction(NULL);
  new IEWallHeatConduction(NULL);
  new FAngular(NULL);
  new FCurvature(NULL);
  new FCentrifugal(NULL);
  new FCoriolis(NULL);
  new FDPD(NULL);
  new FDPDE(NULL);
  new FGeneric(NULL);
//   new Fgrav(NULL);
  new FKinetic(NULL);
  new FMDPD(NULL);
  new FPairScalar(NULL);
  new FPairTensor(NULL);
  new FPairVector(NULL);
  new FPairVels(NULL);
  new FPairVelsBC(NULL);
  new FPairVelsWF(NULL);
  new FParticleScalar(NULL);
  new FParticleVector(NULL);
  new FParticleVectorMatrix(NULL);
  new FParticleVectorRandMatrix(NULL);
  new FParticleVels(NULL);
  new FParticleVelsMatrix(NULL);
  new FParticleTensor(NULL);
  new Frand(NULL);
  new Fspecific(NULL);
  new FWallDPD(NULL);
  new FWallRepulsion(NULL);
  new LJ(NULL);
  new LJVC(NULL);
//   new ShearMigration(NULL);
//   new ViscosityMigration(NULL);

#ifdef WITH_ARRAY_TYPES
  new GridMeterArrayAttr(NULL);
#endif
  new GridMeterDensity(NULL);
  new GridMeterEKin(NULL);
  new GridMeterLocalDensity(NULL);
  new GridMeterPressure(NULL);
  new GridMeterScalar(NULL);
  new GridMeterScalarFourierCos(NULL);
  new GridMeterScalarFourierSin(NULL);
  new GridMeterScalarFunc(NULL);
  new GridMeterScalarSum(NULL);
  new GridMeterScalarTimesScalarSum(NULL);
  new GridMeterSelfDiffusion(NULL);
//   new GridMeterShearRate(NULL);
  new GridMeterSurfaceTension(NULL);
  new GridMeterTemperature(NULL);
  new GridMeterTemperatureIE(NULL);
  new GridMeterVector(NULL);
  new GridMeterVelocity(NULL);
  new GridMeterVelocityMoment(NULL);
  new GridMeterVolumeFraction(NULL);
  new GridMeterVector(NULL);

  // Integrators
  new IntegratorEnergy(NULL);
  new IntegratorISPHconstRhoFDRHS(NULL);
  new IntegratorOmelyan(NULL);
  new IntegratorOmelyanNR(NULL);
#ifdef WITH_ARRAY_TYPES
#ifdef HAVE_JAMA_JAMA_LU_H
  new IntegratorLSE(NULL);
#endif
#endif
  new IntegratorPosition(NULL);
  new IntegratorScalar(NULL);
  new IntegratorScalarLambda(NULL);
  new IntegratorStatic(NULL);
#ifdef WITH_ARRAY_TYPES
#ifdef HAVE_JAMA_JAMA_LU_H
  new IntegratorStaticLSE(NULL);
#endif
#endif
  new IntegratorVector(NULL);
  new IntegratorVectorLambda(NULL);
  new IntegratorTensor(NULL);
  new IntegratorVelocityVerletDisp(NULL);
  new IntegratorVelocityVerletDispX(NULL);
  new IntegratorTensorLambda(NULL);
  new IntegratorVelocityVerlet(NULL);
#ifdef HAVE_JAMA_JAMA_LU_H
  new IntegratorVelocityVerletPressure(NULL);
#endif
  new IntegratorVelocityVerletSymbols(NULL);
  new IntegratorVelocityVerletX(NULL);
  new IntegratorVelocityVerletXFull(NULL);


  new MeterAutocorrelationVector(NULL);
  new MeterAutocorrelationVectorT(NULL);
  new MeterAverage(NULL);
  new MeterBondedCrosscorrelationVector(NULL);
  new MeterDistribution(NULL);
#ifdef HAVE_SDL
  new MeterLive(NULL);
  new MeterLiveColoured(NULL);
#endif
  new MeterLiveGeomview(NULL);
  new MeterPairDistribution(NULL);
  new MeterPairDistributionWithWalls(NULL);
  new MeterPosVel(NULL);
  new MeterRelativeVelocity(NULL);
  new GridAveragerCircular(NULL);
  new GridAveragerStructured(NULL);

  new OutputDCD(NULL, NULL);
  new OutputPDB(NULL, NULL);
  new OutputVTK(NULL, NULL);
  new OutputFile(NULL, NULL);
  new OutputMathematica(NULL, NULL);

  // ParticleChaches
  new ParticleCacheDensity0Oc(NULL);
  new ParticleCacheDensitySelfContribution(NULL);
//   new ParticleCacheShearSelfContribution(NULL);
  new ParticleCacheVolumeSelfContribution(NULL);
  new ParticleRandNormScalar(NULL);
  new ParticleRandNormVector(NULL);
  new ParticleScalar(NULL);
  new ParticleTensor(NULL);
  new ParticleVector(NULL);
  new PCacheIAPWSIF97Cp(NULL);
  new PCacheIAPWSIF97eta(NULL);
  new PCacheIAPWSIF97kappa(NULL);
  new PCacheIAPWSIF97p(NULL);
  new PCacheIAPWSIF97rho(NULL);
  new ParticleVels(NULL);
  new PCaEigensystem(NULL);
  new PCaMatrixInverse(NULL);
  new SymbolFParticleScalar(NULL);
  new SymbolFParticleVels(NULL);


  new PairRandScalar(NULL);
  new PairRandTensor(NULL);
  new PairRandVector(NULL);
  new PairScalar(NULL);
  new PairTensor(NULL);
  new PairVector(NULL);
  
  
  new ParticleConnectorFile(NULL);
  new ParticleCreatorFile(NULL);
  new ParticleCreatorInlet(NULL);
  new ParticleCreatorLattice(NULL);
  new ParticleCreatorLatticeBCC(NULL);
  new ParticleCreatorLatticeFrozen(NULL);
#ifdef WITH_ARRAY_TYPES
  new ParticleCreatorLSE(NULL);
#endif
  new ParticleCreatorRandom(NULL);
  new ParticleCreatorTube(NULL);
  new ParticleCreatorStatic(NULL);
  new ParticleCreatorWall(NULL);
  new ParticleCreatorWallTextured(NULL);

  new ReflectorBounceBack(NULL);
  new ReflectorMirror(NULL);
  new ReflectorStochastic(NULL);
  new ReflectorThermalize(NULL);
  new ReflectorThermalizeInternalEnergy(NULL);

  // Callables
  new ApplyVectorField(NULL);
  new ApplyVectorFieldFile(NULL);
  new ApplyVelField(NULL);
  new CblPairParticleTensor(NULL);
  new CblPairParticleVector(NULL);
  new EnergyPump(NULL);
  new ShiftParticle(NULL);
  new TriangleCGMDInterpolation(NULL);
  new TriangleCGMDInterpolation2Points(NULL);
  new ThermostatEnergyRescaling(NULL);
  new ThermostatLA(NULL);
  new ThermostatLAEnergyConserving(NULL);
  new ThermostatLAWithWeightEC(NULL);
  new ThermostatVels(NULL);
  new ThermostatPetersEnergyConserving(NULL);
  new ThermostatWallVelEC(NULL);
  new ThermostatPetersIso(NULL);
  new ThermostatWallVelIso(NULL);
  new WriteRestartFile(NULL);

  new TransferParticleVector(NULL);

#ifdef HAVE_JAMA_JAMA_LU_H
  new VelConstraints(NULL);
#endif
  //  ValCalculators
  new BondedPairParticleOneNoiseVector(NULL);
  new BondedPairParticleScalar(NULL);
  new BondedPairParticleTensorNoiseVector(NULL);
  new BondedPairParticleVector(NULL);
  new BondedPairScalar(NULL);
  new BondedPairVector(NULL);
  new PairParticleScalar(NULL);
  new PairParticleTensor(NULL);
  new PairParticleVector(NULL);
  new TripletCalcAngularDt2F(NULL);
  new TripletCalcAngularF(NULL);
  new TripletCalcAngularFalpha(NULL);
  new TripletCalcCentralPartScalar(NULL);
  new QuintetCalcCurvatureF(NULL);
  new QuintetCalcCurvature(NULL);
//   new ValCalculatorDirichletBCScalar(NULL);
  new ValCalculatorDirichletBCScalar(NULL);
  new ValCalculatorDirichletBCVels(NULL);
  new ValCalculatorKernel(NULL);
  new ValCalculatorNegDKernelDivr(NULL);
  new ValCalculatorRi(NULL);
  new ValCalculatorRi6(NULL);
  new ValCalculatorRho(NULL);
  // commented out because trying to get rid of it
//   new ValCalculatorShear(NULL);
//   new ValCalculatorShearX(NULL);
// commented out because buggy
//   new ValCalculatorSymmetryBCScalar(NULL);
  new ValCalculatorVolume(NULL);

  //  new TabulatedWeightingFunctionWithWall(NULL);
  new InputWF(NULL);
  new Linear(NULL);
  new Lucy(NULL);
  new LucyWithWall0thOrder(NULL);
  new LucyWithWall1stOrder(NULL);
  new Square(NULL);
}




/*#include "geometric_primitives.h"

int main(int argc, char *argv[])
{
    cuboid_t c;
    point_t a, b;

    c.corner1.x = 0;
    c.corner1.y = 0;
    c.corner1.z = 0;

    c.corner2.x = 1;
    c.corner2.y = 1;
    c.corner2.z = 1;

    a.x = -1;
    a.y = -1;
    a.z = -1;

    b.x = 0.5;
    b.y = 0.5;
    b.z = 0.5;

    cout << c.intersects(a, b) << endl;
    cout << c.intersects(b, a) << endl;

    return 0;
    }*/
