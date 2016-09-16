/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
 * David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
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


#include <fstream>

#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"
#include "triangle_CGMD_interpolation.h"

using namespace std;

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

const Callable_Register<TriangleCGMDInterpolation> triangle_CGMD_interpolation("TriangleCGMDInterpolation");

TriangleCGMDInterpolation::TriangleCGMDInterpolation(Simulation* sim)
  : Callable(sim), m_working_mat(NULL), m_inverse(NULL), m_permutation(NULL), m_zeroStored(false) 
{
  init();
}

TriangleCGMDInterpolation::~TriangleCGMDInterpolation()
{
  if(m_working_mat) gsl_matrix_free(m_working_mat);
  if(m_inverse) gsl_matrix_free(m_inverse);
  if(m_permutation) gsl_permutation_free(m_permutation);
}

void TriangleCGMDInterpolation::init()
{
  m_properties.setClassName("TriangleCGMDInterpolation");
  
  m_properties.setDescription("A TriangleCGMDInterpolation interpolates atomic displacements and velocities on the nodes of a user-defined triangular mesh. Linear CGMD interpolation is used as described in Rudd&Broughton, Phys. Rev. B 72 144104 \n"
 "* ASSUMPTION1: all original node-positions are INSIDE the box \n"
 "* ASSUMPTION2: the element sides are always smaller than box/2 in any direction \n"
"* ASSUMPTION3: all elements lie in (x,y)-plane (is partially used, e.g., for shape functions and element area)");
  
  STRINGPCINF
    (name, m_filename,
     "Name of the file containing the nodal indices of the elements. It is assumed that these indices indicate slots of \a Particle s of species specified by the attribute 'nodeSpecies'");

  m_filename = "undefined";

  INTPC
    (interval, m_interval, 0,
     "Call interval.");

  m_interval = 1;

  INTPC
    (activateAt, m_activate_at, -1,
     "Timestep at which the module should be activated.");

  m_activate_at = 0;

  INTPC
    (deactivateAt, m_deactivate_at, -2,
     "Timestep at which the module should be deactivated. -1 means never deactivate.");

  m_deactivate_at = -1;

  STRINGPC
    (nodeSpecies, m_nodeSpecies,"Species of the particles acting as nodes of the elements.");

  m_nodeSpecies = "undefined";

  STRINGPC
    (pSpecies, m_pSpecies, "Species of the particles delivering the data used for interpolation by the nodes.");

  m_pSpecies = "undefined";

  BOOLPC(periodic, m_periodic, "Should the module apply periodic boundary conditions?");

  m_periodic = true;

  STRINGPC 
      (targetSlot, m_targetSlotName,
       "Name for the symbol of the \"source particle\" storing the slot of the \"target particle\".");

  m_targetSlotName = "undefined";

  STRINGPC 
      (displacement, m_displacementName,
       "Name for the symbol of the node storing the non-locally interpolated displacement.");

  m_displacementName = "undefined";

  STRINGPC 
      (momentum1, m_momentum1Name,
       "Name for the symbol of the node storing the non-locally interpolated momentum.");

  m_momentum1Name = "undefined";

//   STRINGPC 
//       (momentum2, m_momentum2Name,
//        "Name for the symbol of the node storing the locally interpolated momentum.");

//   m_momentum2Name = "undefined";

  STRINGPC 
      (force1, m_force1Name,
       "Name for the symbol of the node storing the non-locally interpolated force.");

  m_force1Name = "undefined";

//   STRINGPC 
//       (force2, m_force2Name,
//        "Name for the symbol of the node storing the locally interpolated force.");

//   m_force2Name = "undefined";

  BOOLPC 
      (transient, m_transient,
       "If set to \"yes\", the interpolation matrix will not be computed as a function of the initial (assumed as equilibrium) atomic positions but as a function of the actual positions at every time step. This is of course more time-consuming because it requires the corresponding matrix inversion at every time step and not just once at the beginning.");

  m_transient = false;

  BOOLPC 
      (writeNN, m_writeNN,
       "Should the matrix of shape functions be written to a file?");

  m_writeNN = false;

  INTPC
    (writeNNevery, m_writeNNevery, 0,
     "Writing interval for matrix of shape functions. Notice: Only the calls of this module are counted, NOT every time step! If transient = \"no\", there is only one call. Inactive if 'writeNN = \"no\"'");

  m_writeNNevery = 1;

  STRINGPC(NNname, m_NNname, "File name for storing the matrix of shape functions in binary format. A counter and the ending '.bin' are added automatically. Inactive if 'writeNN = \"no\"'");

  m_NNname = "NNmatrix";

  BOOLPC 
      (writeF, m_write_f,
       "Should the CGMD interpolation matrix be written to a file?");

  m_write_f = false;

  INTPC
    (writeFEvery, m_write_fEvery, 0,
     "Writing interval for CGMD interpolation matrix. Notice: Only the calls of this module are counted, NOT every time step! If transient = \"no\", there is only one call. Inactive if 'writeF = \"no\"'");

  m_write_fEvery = 1;

  STRINGPC(fName, m_fName, "File name for storing the cgmd interpolation matrix in binary format. A counter and the ending '.bin' are added automatically. Inactive if 'writeF = \"no\"'");

  m_fName = "fMatrix";


} 

void TriangleCGMDInterpolation::setup()
{
  m_NNstepCounter = 0;
  m_NNfileCounter = 0;

  m_fStepCounter = 0;
  m_fFileCounter = 0;

  m_active = false;

  if(m_filename == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "Attribute 'nameInputFile' has value \"undefined\".");

  if(m_nodeSpecies == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "Attribute 'nodeSpecies' has value \"undefined\".");

  m_nodeColour = M_MANAGER->getColour(m_nodeSpecies);

  if(m_pSpecies == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "Attribute 'pSpecies' has value \"undefined\".");

  m_pColour = M_MANAGER->getColour(m_pSpecies);

  if(m_targetSlotName == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "'targetSlot' has value \"undefined\"");

  if(m_displacementName == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "'displacement' has value \"undefined\"");

  if(m_momentum1Name == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "'momentum1' has value \"undefined\"");

//   if(m_momentum2Name == "undefined")
//     throw gError("TriangleCGMDInterpolation::setup", "'momentum2' has value \"undefined\"");

  if(m_force1Name == "undefined")
    throw gError("TriangleCGMDInterpolation::setup", "'force1' has value \"undefined\"");

//   if(m_force2Name == "undefined")
//     throw gError("TriangleCGMDInterpolation::setup", "'force2' has value \"undefined\"");

  if(Particle::s_tag_format[m_pColour].attrExists(m_targetSlotName)) {
    m_targetSlotOffset = Particle::s_tag_format[m_pColour].indexOf(m_targetSlotName, DataFormat::INT);
    m_targetSlotOffset = Particle::s_tag_format[m_pColour].offsetByIndex(m_targetSlotOffset);
  }
  // if the target-slot symbol doesn't yet exist, we add it
  else {
    // This MUST be persistent !!!
    // Currently (2010-08-28) we create the symbol here, and then, no module writes it except ParticleCreatorFile once at the beginning
    // but if we use accidentally a m_targetSlotName which does not appear in the file containing the particles, this goes wrong! Hence there is a consistency check in TriangleCGMDInterpolation::setupAfterParticleCreation 
    m_targetSlotOffset = Particle::s_tag_format[m_pColour].addAttribute(m_targetSlotName, DataFormat::INT, true, m_targetSlotName).offset;
  MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for " << m_targetSlotName << " = " << m_targetSlotOffset);
  }

  if(Particle::s_tag_format[m_nodeColour].attrExists(m_displacementName)) {
    throw gError("TriangleCGMDInterpolation::setup", "Symbol " + m_displacementName  + " for species " + m_nodeSpecies + " already exists!");
  }
  // if the symbol doesn't yet exist, we add it
  else {
    // We make it persistent and initialise to zero right before computation
    m_displacementOffset = Particle::s_tag_format[m_nodeColour].addAttribute(m_displacementName, DataFormat::POINT, true, m_displacementName).offset;
  MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for " << m_displacementName << " = " << m_displacementOffset);
  }

  if(Particle::s_tag_format[m_nodeColour].attrExists(m_momentum1Name)) {
    throw gError("TriangleCGMDInterpolation::setup", "Symbol " + m_momentum1Name  + " for species " + m_nodeSpecies + " already exists!");
  }
  // if the symbol doesn't yet exist, we add it
  else {
    // We make it persistent and initialise to zero right before computation
    m_momentum1Offset = Particle::s_tag_format[m_nodeColour].addAttribute(m_momentum1Name, DataFormat::POINT, true, m_momentum1Name).offset;
  MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for " << m_momentum1Name << " = " << m_momentum1Offset);
  }

//   if(Particle::s_tag_format[m_nodeColour].attrExists(m_momentum2Name)) {
//     throw gError("TriangleCGMDInterpolation::setup", "Symbol " + m_momentum2Name  + " for species " + m_nodeSpecies + " already exists!");
//   }
//   // if the symbol doesn't yet exist, we add it
//   else {
//     // We make it persistent and initialise to zero right before computation
//     m_momentum2Offset = Particle::s_tag_format[m_nodeColour].addAttribute(m_momentum2Name, DataFormat::POINT, true, m_momentum2Name).offset;
//   MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for " << m_momentum2Name << " = " << m_momentum2Offset);
//   }

  if(Particle::s_tag_format[m_nodeColour].attrExists(m_force1Name)) {
    throw gError("TriangleCGMDInterpolation::setup", "Symbol " + m_force1Name  + " for species " + m_nodeSpecies + " already exists!");
  }
  // if the symbol doesn't yet exist, we add it
  else {
    // We make it persistent and initialise to zero right before computation
    m_force1Offset = Particle::s_tag_format[m_nodeColour].addAttribute(m_force1Name, DataFormat::POINT, true, m_force1Name).offset;
  MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for " << m_force1Name << " = " << m_force1Offset);
  }

//   if(Particle::s_tag_format[m_nodeColour].attrExists(m_force2Name)) {
//     throw gError("TriangleCGMDInterpolation::setup", "Symbol " + m_force2Name  + " for species " + m_nodeSpecies + " already exists!");
//   }
//   // if the symbol doesn't yet exist, we add it
//   else {
//     // We make it persistent and initialise to zero right before computation
//     m_force2Offset = Particle::s_tag_format[m_nodeColour].addAttribute(m_force2Name, DataFormat::POINT, true, m_force2Name).offset;
//   MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for " << m_force2Name << " = " << m_force2Offset);
//   }

  size_t __attrCounter = 0;
 
  while(Particle::s_tag_format[m_pColour].attrExists("__myElement" + ObjToString(__attrCounter))) {
    //     throw gError("TriangleCGMDInterpolation::setup", "Attribute with name \"__myElement\" already created by other module. Cannot create twice. Avoid creating a symbol with this name in the input file!");
    ++__attrCounter;
  }
  // This MUST be persistent !!! If we need a change we replace the old value by a new one
  m_elementSlotOffset = Particle::s_tag_format[m_pColour].addAttribute("__myElement" + ObjToString(__attrCounter), DataFormat::INT, true, "__myElement").offset;
  MSG_DEBUG("TriangleCGMDInterpolation::setup", "Offset for \"__myElement\"" << __attrCounter << " = " << m_targetSlotOffset);

  M_CONTROLLER->registerForSetupAfterParticleCreation(this);
}

void TriangleCGMDInterpolation::setupAfterParticleCreation() 
{
//   size_t pColour = M_MANAGER->getColour(m_pSpecies);
  size_t nOfParticles = M_PHASE->particles(m_pColour).size();

//   size_t nodeColour = M_MANAGER->getColour(m_nodeSpecies);
  size_t nOfNodes = M_PHASE->particles(m_nodeColour).size();

  m_boxSize = M_PHASE->boundary()->boundingBox().size();

  // consistency check for m_targetSlotName
  bool different = false;
  bool first = true;
  size_t newParticleNode, oldParticleNode;
//   oldParticleNode = newParticleNode;
  FOR_EACH_FREE_PARTICLE_C
    (M_PHASE, m_pColour,     
     newParticleNode = __iSLFE->tag.intByOffset(m_targetSlotOffset);
     if(!first) // check does not make sense in first round of loop
       if(newParticleNode != oldParticleNode) different = true;
     oldParticleNode = newParticleNode;
     if(first) first = false;
     );
  if(!different) 
    throw gError("TriangleCGMDInterpolation::setupAfterParticleCreation", "All particles have the same targetSlot. Currently I consider this as an application error and abort. Check for example the attribute 'targetSlot'!"); 

  // START:-------------allocation and initialisation for later----------------

  m_elementsOfNode.resize(nOfNodes);

  // resize number of lines  
  m_shapeMatrix.resize(nOfNodes);
  // resize number of columns  
  for(vector<vector<double> >::iterator m = m_shapeMatrix.begin(); m != m_shapeMatrix.end(); ++m)
    m->resize(nOfNodes);

  // resize number of lines  
  m_shapeF.resize(nOfParticles);
  // resize number of columns  
  for(vector<vector<double> >::iterator m = m_shapeF.begin(); m != m_shapeF.end(); ++m)
    m->resize(nOfNodes);

  // currently not needed
//   // resize number of lines  
//   m_shapeF_0.resize(nOfParticles);
//   // resize number of columns  
//   for(vector<vector<double> >::iterator m = m_shapeF_0.begin(); m != m_shapeF_0.end(); ++m)
//     m->resize(nOfNodes);

  // resize number of lines  
  m_f.resize(nOfParticles);
  // resize number of columns  
  for(vector<vector<double> >::iterator m = m_f.begin(); m != m_f.end(); ++m)
    m->resize(nOfNodes);

  // currently not needed
//   // resize number of lines  
//   m_f_0.resize(nOfParticles);
//   // resize number of columns  
//   for(vector<vector<double> >::iterator m = m_f_0.begin(); m != m_f_0.end(); ++m)
//     m->resize(nOfNodes);

  //   // resize number of lines  
  //   m_f_diff.resize(nOfNodes);
  //   // resize number of columns  
  //   for(vector<vector<double> >::iterator m = m_f_diff.begin(); m != m_f_diff.end(); ++m)
  //     m->resize(nOfParticles);
  
  // allocations and initialisations for matrix inversion with GSL
  m_working_mat = gsl_matrix_alloc(nOfNodes, nOfNodes);
  m_inverse = gsl_matrix_alloc(nOfNodes, nOfNodes);
  // this by hand setting is only necessary if you do not allocate memory for m_inverse but assigne another memory adress to its double*; Then, m_inverse would also have to be be of non-pointer type gsl_matrix !!!
//   m_inverse.size1 = m_inverse.size2 = m_inverse.tda = nOfNodes;
//   m_inverse.owner = 0;
//   m_inverse.data = NULL;
  m_permutation = gsl_permutation_alloc(nOfNodes);

  // resize array storing the initial positions
  m_rInit.resize(nOfParticles);
  // store the initial positions
  FOR_EACH_FREE_PARTICLE_C
    (M_PHASE, m_pColour,     
     m_rInit[__iSLFE->mySlot] = __iSLFE->r;
     );

  // END:-------------allocation and initialisation for later----------------

  // START:-------- read the element indices from the file -------------------------
  vector<size_t> element;
  element.resize(3);


  ifstream elementFile(m_filename.c_str());
  if (!elementFile)
    throw gError("TriangleCGMDInterpolation::setupAfterParticleCreation"+FILE_INFO,
		 "Error opening file '" + m_filename + "'.");
  
  while 
    (/*(elementFile >> skipws >> element[0] >> skipws >> element[1] >> skipws >> element[2]) && */
     !elementFile.eof()) {
    try {
      elementFile >> skipws >> element[0] >> skipws >> element[1] >> skipws >> element[2];

    }
    catch(ifstream::failure e) {

      cout << "Problem: " << strerror(errno) << endl;

      throw gError("ParticleConnectorFile::setupAfterParticleCreation::"+FILE_INFO,
		   "Error reading file '" + m_filename + "': " + e.what());
      
    }
    m_cornerIndices.push_back(element);
    size_t lastElement = m_cornerIndices.size()-1;
    // tell each node to which elements it belongs to
    m_elementsOfNode[element[0]].push_back(lastElement);
    m_elementsOfNode[element[1]].push_back(lastElement);
    m_elementsOfNode[element[2]].push_back(lastElement);

  } // end of while

  //-------- START: DEBUGGING --------------------

//   size_t debug_nOfNodes = m_elementsOfNode.size();
//   MSG_DEBUG("TriangleCGMDInterpolation::setupAfterParticleCreation", "Finished reading element nodes: nOfNodes = " << debug_nOfNodes); 
//   MSG_DEBUG("TriangleCGMDInterpolation::setupAfterParticleCreation", "This is in my array of elements:");
//   for(size_t element = 0; element < m_cornerIndices.size(); ++element) {
//     for(size_t node = 0; node < 3; ++node) {
//       cout << m_cornerIndices.at(element).at(node) << " ";
//       if(m_cornerIndices.at(element).at(node) > debug_nOfNodes)
// 	throw gError("TriangleCGMDInterpolation::setupAfterParticleCreation", "That's bad! nOfNodes = " + ObjToString(debug_nOfNodes));
//     }
//     cout << endl;
//   }
//   MSG_DEBUG("TriangleCGMDInterpolation::setupAfterParticleCreation", "again especially m_cornerIndices.at(element=864).at(node) = " << m_cornerIndices.at(864).at(0) << ", " << m_cornerIndices.at(864).at(1) << ", " << m_cornerIndices.at(864).at(2) );

//   MSG_DEBUG("TriangleCGMDInterpolation::setupAfterParticleCreation", "NofNodes according to array m_elementsOfNode = " << debug_nOfNodes);
//   for(size_t node = 0; node < debug_nOfNodes; ++node) {
//     size_t debug_nOfElements = m_elementsOfNode[node].size();
//     MSG_DEBUG("TriangleCGMDInterpolation::setupAfterParticleCreation", "nOfElements for node " << node << " = " << debug_nOfElements);
//     MSG_DEBUG("TriangleCGMDInterpolation::setupAfterParticleCreation", "elements of node " << node << " : ");
//     for(size_t element = 0; element < debug_nOfElements; ++element) {
//       cout << m_elementsOfNode.at(node).at(element) << " ";
//     }
//     cout << endl;
//   }

  //-------- END: DEBUGGING --------------------

  // END:-------- read the element indices from the file -------------------------


  //-------- START: initialisation of element corners and sides-------------------
  // this is very similar to what WallTriangle is doing. But for now, I decided to do everything locally without additional objects from other (old or new) classes

  size_t nOfElements = m_cornerIndices.size();

  // m_corners,... contain three point_t per element
  m_corners.resize(nOfElements);
  m_periodicCorners.resize(nOfElements);
  m_sides.resize(nOfElements);
  m_side_normals.resize(nOfElements);
  m_areasTwice.resize(nOfElements);

  for(size_t el = 0; el < nOfElements; ++el) {
    // resizing to three point_t per element
    m_corners[el].resize(3);
    m_periodicCorners[el].resize(3);
    m_sides[el].resize(3);
    m_side_normals[el].resize(3);
    // set corners; correction for periodicity comes below
    for (size_t corner = 0; corner < 3; ++corner) {
      m_corners[el][corner] = m_periodicCorners[el][corner] = ((M_PHASE->particles(m_nodeColour))[m_cornerIndices[el][corner]]).r;
    }
    // compute sides; correction for periodicity comes below
    for (size_t side = 0; side < 3; ++side) {
      m_sides[el][side] = m_periodicCorners[el][(side+1)%3] - m_periodicCorners[el][side];
    }
  }

  if(m_periodic) {
    for(size_t el = 0; el < nOfElements; ++el) {
      
      
      // START: PERIODICITY -----------------------------------------
      
      // ASSUMPTION1: all original node-positions are INSIDE the box
      // ASSUMPTION2: the element sides are always smaller than box/2. in any direction
      
      // for a triangle the following is true:
      // per direction, always two nodes are on the same side and one is periodically shifted
      // Hence, we can decide for each direction separately which particle's dir-coordinate to shift to where
      
      // bool toBeShifted[3]; // currently not needed
    
      for (size_t dir = 0; dir < SPACE_DIMS; ++dir) {

// 	MSG_DEBUG
// 	  ("TriangleCGMDInterpolation::setupAfterParticleCreation()", "START: periodicity check for element nodes for dir " << dir << endl 
// 	   << "corners[dir] = " << m_corners[el][0][dir] << ", " << m_corners[el][1][dir] << ", " << m_corners[el][2][dir] << endl
// 	   << "sides[dir] = " << m_sides[el][0][dir] << ", " << m_sides[el][1][dir] << ", " << m_sides[el][2][dir] << endl
// 	   << "box[dir] = " << m_boxSize[dir] << endl
// 	   << "element = " << el << endl
// 	   );

	// starting check with node 0
	if(m_sides[el][0][dir] > m_boxSize[dir]/2.) {
	  // then for sides 1+2 the SAME periodicity cannot be true
	  // take care! the m_sides[el] are defined in a cyclic way !!!
	  assert(m_sides[el][1][dir] < m_boxSize[dir]/2.);
	  assert(m_sides[el][2][dir] < m_boxSize[dir]/2.);
	  // cyclic !!!
	  if(m_sides[el][1][dir] < -m_boxSize[dir]/2.) {
	    // then node 1 must be shifted
	    m_corners[el][1][dir] -= m_boxSize[dir];
	    m_sides[el][0][dir] = m_corners[el][1][dir] - m_corners[el][0][dir];
	    m_sides[el][1][dir] = m_corners[el][2][dir] - m_corners[el][1][dir];
	    // and side 2 should be OK
	    assert(m_sides[el][2][dir] > -m_boxSize[dir]/2.);
	  }
	  // cyclic !!!
	  else if(m_sides[el][2][dir] < -m_boxSize[dir]/2.) {
	    // then node 0 must be shifted
	    m_corners[el][0][dir] += m_boxSize[dir];
	    m_sides[el][0][dir] = m_corners[el][1][dir] - m_corners[el][0][dir];
	    m_sides[el][2][dir] = m_corners[el][0][dir] - m_corners[el][2][dir];
	    // and side 1 should be OK
	    assert(m_sides[el][1][dir] > -m_boxSize[dir]/2.);
	  }
	}
	// "first else"
	else { 
	  if(m_sides[el][0][dir] < -m_boxSize[dir]/2.) {
	    // then for sides 1+2 the SAME periodicity cannot be true
	    // take care! the m_sides[el] are defined in a cyclic way !!!
	    assert(m_sides[el][1][dir] > -m_boxSize[dir]/2.);
	    assert(m_sides[el][2][dir] > -m_boxSize[dir]/2.);
	    // cyclic !!!
	    if(m_sides[el][1][dir] > m_boxSize[dir]/2.) {
	      // then node 1 must be shifted
	      m_corners[el][1][dir] += m_boxSize[dir];
	      m_sides[el][0][dir] = m_corners[el][1][dir] - m_corners[el][0][dir];
	      m_sides[el][1][dir] = m_corners[el][2][dir] - m_corners[el][1][dir];
	      // and side 2 should be OK
	      assert(m_sides[el][2][dir] < m_boxSize[dir]/2.);
	    }
	    // cyclic !!!
	    else if(m_sides[el][2][dir] > m_boxSize[dir]/2.) {
	      // then node 0 must be shifted
	      m_corners[el][0][dir] -= m_boxSize[dir];
	      m_sides[el][0][dir] = m_corners[el][1][dir] - m_corners[el][0][dir];
	      m_sides[el][2][dir] = m_corners[el][0][dir] - m_corners[el][2][dir];
	      // and side 1 should be OK
	      assert(m_sides[el][1][dir] < m_boxSize[dir]/2.);
	    }
	  }
	  // if we now come to "else if", we just have to check in which direction node 2 must be shifted because side 0 is OK
	  // side 1 points from node 1 to 2
	  // "second else"
	  else { 
	    if(m_sides[el][1][dir] > m_boxSize[dir]/2.) {
	      // m_sides[el][2][dir] should be large as well because side 0 is OK
	      // cyclic !!!
	      assert(m_sides[el][2][dir] < -m_boxSize[dir]/2.);
	      m_corners[el][2][dir] -= m_boxSize[dir];
	      // correct the two affected sides
	      m_sides[el][1][dir] = m_corners[el][2][dir] - m_corners[el][1][dir];
	      m_sides[el][2][dir] = m_corners[el][0][dir] - m_corners[el][2][dir];
	    }
	    else if(m_sides[el][1][dir] < -m_boxSize[dir]/2.) {
	      // m_sides[el][2][dir] should be large as well because side 0 is OK
	      // cyclic !!!
	      assert(m_sides[el][2][dir] > m_boxSize[dir]/2.);
	      m_corners[el][2][dir] += m_boxSize[dir];
	      // correct the two affected sides
	      m_sides[el][1][dir] = m_corners[el][2][dir] - m_corners[el][1][dir];
	      m_sides[el][2][dir] = m_corners[el][0][dir] - m_corners[el][2][dir];
	    }

	  } // end of "second else"
	      
	} // end of "first else"
		  
      } // end of dir-loop

      // END: PERIODICITY -----------------------------------------

    } // end loop over elements
  } // of if(m_periodic)

MSG_DEBUG
  ("TriangleCGMDInterpolation::setupAfterParticleCreation()", "END: periodicity of elements");

  for(size_t el = 0; el < nOfElements; ++el) {
    /* if done as follows, the surface normal ensures that the SIDE normals always point INTO the triangle 
       SUCCESFULLY CHECKED: NO MATTER WHETHER TRIANGLE-POINTS ARE ORDERED CLOCCK-WISE !!!! 
    */
    
    point_t surface_normal = m_sides[el][0].cross(m_sides[el][1]);
    surface_normal /= surface_normal.abs();
    
    /* not needed
       m_ndotr = m_corners[0]*m_surface_normal;
    */	
    
    /* m_side_normals point INTO the triangle. */
    for (size_t side = 0; side < 3; ++side) {
      m_side_normals[el][side] = surface_normal.cross(m_sides[el][side]);
      m_side_normals[el][side] /= m_side_normals[el][side].abs();
    }
    
    // compute two times the area (2xArea) of each element, 
    // which is in Mathematica notation: 
    // Det[{{1, x0, y0}, {1, x1, y1}, {1, x2, y2}}]
    // =  +x0 y1 -x1 y0 +x1 y2 -x2 y1 +x2 y0 -x0 y2   
    
    // Here we need ASSUMPTION3: all elements lie in (x,y)-plane
    m_areasTwice[el] = 
       m_corners[el][0][0] * m_corners[el][1][1] // +x0 y1
      -m_corners[el][0][0] * m_corners[el][2][1] // -x0 y2 
      -m_corners[el][1][0] * m_corners[el][0][1] // -x1 y0 
      +m_corners[el][1][0] * m_corners[el][2][1] // +x1 y2 
      -m_corners[el][2][0] * m_corners[el][1][1] // -x2 y1 
      +m_corners[el][2][0] * m_corners[el][0][1];// +x2 y0
  }

  //-------- END: initialisation of element corners and sides-------------------


  // currently not needed
//   // initialise m_f_0 and m_shapeF_0 to zero before first computation; rest is initialised in call()
//   for(vector<vector<double> >::iterator m = m_f_0.begin(); m != m_f_0.end(); ++m)
//     for(vector<double>::iterator n = m->begin(); n != m->end(); ++n)
//       (*n) = 0;

//   for(vector<vector<double> >::iterator m = m_shapeF_0.begin(); m != m_shapeF_0.end(); ++m)
//     for(vector<double>::iterator n = m->begin(); n != m->end(); ++n)
//       (*n) = 0;

  // computations before first timestep (which is also "0", so don't get confused)
  call(0);

  // the following is now (2011-03-21) in call()

//   for(size_t i = 0; i < nOfParticles; ++i)
//     for(size_t m = 0; m < nOfNodes; ++m) {
//       m_f_0[i][m] = m_f[i][m];
//       m_shapeF_0[i][m] = m_shapeF[i][m];
//     }

//   interpolate();
}

void TriangleCGMDInterpolation::call(size_t timestep)
{

  // !!! we assume that the initial configuration is the 
  // equilibrium configuration, no matter whether some 
  // equilibration comes first !!!
  if(!m_zeroStored) {
    // 		  MSG_DEBUG("TriangleCGMDInterpolation::call", "!m_zeroStored case");
    compute();
    m_zeroStored = true;
  }
  
  if (m_active) {
// 	MSG_DEBUG("TriangleCGMDInterpolation::call", "m_active == true case");
    if ((int) timestep == m_deactivate_at) {
      m_active = false;
// 	  MSG_DEBUG("TriangleCGMDInterpolation::call", "NOW DEACTIVATED");
    } 
    else {
      m_step--;
      if (!m_step) {
	if(m_transient) {
	  compute();
	}
	// always interpolate
	interpolate();
	m_step = m_interval;
      } // end if (!m_step)
    } // end else of if ((int) timestep == m_deactivate_at)
  } // end if(m_active) 
  else if ((int) timestep == m_activate_at) {
    m_active = true;
// 	MSG_DEBUG("TriangleCGMDInterpolation::call", "now ACTIVATING, m_step before = " << m_step);
    m_step = m_interval;    
// 	MSG_DEBUG("TriangleCGMDInterpolation::call", "m_step after = " << m_step);
    if(m_transient) {
//       MSG_DEBUG("TriangleCGMDInterpolation::call", "m_transient case");
      compute();
    }
    // always interpolate
    interpolate();
  } // end of if((int) timestep == m_activate_at)
}


void TriangleCGMDInterpolation::compute()
{
  Phase* phase = M_PHASE;

  size_t nOfParticles = phase->particles(m_pColour).size();
  size_t nOfNodes = phase->particles(m_nodeColour).size();

//   ParticleList& particles = phase->particles(m_pColour);
//   ParticleList& nodes = phase->particles(m_nodeColour);

  // initialisation to zero

  for(vector<vector<double> >::iterator m = m_f.begin(); m != m_f.end(); ++m)
    for(vector<double>::iterator n = m->begin(); n != m->end(); ++n)
      (*n) = 0;

  for(vector<vector<double> >::iterator m = m_shapeMatrix.begin(); m != m_shapeMatrix.end(); ++m)
    for(vector<double>::iterator n = m->begin(); n != m->end(); ++n)
      (*n) = 0;

  for(vector<vector<double> >::iterator m = m_shapeF.begin(); m != m_shapeF.end(); ++m)
    for(vector<double>::iterator n = m->begin(); n != m->end(); ++n)
      (*n) = 0;


//--------------START: loop over particles -------------------------- 

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_pColour,
     // the node this particle __iSLFE belongs to
     size_t particleNode = __iSLFE->tag.intByOffset(m_targetSlotOffset);

     // find (FIRST!!! avoid to find more due to COMPUT. GEOMETRY!!!) element el_i where __iSLFE is in
     // loop over all elements containing this node
     vector<int>* elements = &(m_elementsOfNode[particleNode]);
     vector<int>::iterator endEl = elements->end();
     for(vector<int>::iterator elit = elements->begin(); elit != endEl; ++elit) {
       // index of current element
       size_t el = *elit;

//         MSG_DEBUG("TriangleCGMDInterpolation::compute", "Starting element-search at element" << el);

       //---START check for periodic BCs, i.e. whether the particle-pos 
       // needs to be shifted periodically inside the element
       point_t testPos = __iSLFE->r;
       
       // find the corner in this element corresponding to this particleNode
       for(size_t corner = 0; corner < 3; ++corner) {
	 if(particleNode == m_cornerIndices[el][corner]) {
	   // check separately for each direction
	   for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
	     // Is the distance to the particle's (possibly shifted) corner
	     // which corresponds to particleNode periodic? 
	     // If yes, shift particle position.
	     if(testPos[dir] - m_corners[el][corner][dir] > m_boxSize[dir]/2)
	       testPos[dir] -= m_boxSize[dir];
	     if(testPos[dir] - m_corners[el][corner][dir] < -m_boxSize[dir]/2)
	       testPos[dir] += m_boxSize[dir];
	   }
	 } // end if
       } // end loop over corners
       
       //---END: check for periodic BCs

       // now comes the real check whether the particle is in this element
       bool isInside = true;

       /* If the normal product of the side normal with the intersection point
       becomes smaller than zero the point is outside. */
       for (int side = 0; side < 3; ++side) {
         if (m_side_normals[el][side] * (testPos - m_corners[el][side]) < c_wt_dist_eps)
         {
           isInside = false;
	   break; // should break the loop over sides
         }
       }
//          MSG_DEBUG("TriangleCGMDInterpolation::compute", "isInside-check finished with isInside = " << isInside);

       if(isInside) {
	 // store index of element for later
	 __iSLFE->tag.intByOffset(m_elementSlotOffset) = el;

	 // compute contribution to matrix N_mn += N_mi*Nni
         // only the product for nodes m, n of the same element (where i is in)
         // are different from zero 
	 for(size_t node_m = 0; node_m < 3; ++node_m) {
	   size_t global_m = m_cornerIndices[el][node_m];

	   // compute the matrix of linear interpolation values

	   m_shapeF[__iSLFE->mySlot][global_m] = shapeF(el, node_m, testPos); 

// 	   //----START:might vanish soon (2011-03-16)-------------
// 	   // compute here the locally interpolated momentum
// 	   nodes[global_m].tag.pointByOffset(m_momentum2Offset) +=
// 	     // currently (2011-03-16) all masses=1!!! 
// 	     shapeF(el, node_m, testPos) * __iSLFE->v;
// 	   // compute here the locally interpolated force
// 	   nodes[global_m].tag.pointByOffset(m_force2Offset) +=
// 	     // currently (2011-03-16) all masses=1!!! 
// 	     shapeF(el, node_m, testPos) * __iSLFE->force[forceIndex];
// 	   //----END:might vanish soon (2011-03-16)-------------

	   for(size_t node_n = 0; node_n < 3; ++node_n) {
	     size_t global_n = m_cornerIndices[el][node_n];
	     m_shapeMatrix[global_m][global_n] += 
// 	       shapeF(el, node_m, testPos) 
	       m_shapeF[__iSLFE->mySlot][global_m]
	       // the second factor might not yet have been stored 
	       // in m_shapeF, so we should compute it by calling 
	       // the function
	       * shapeF(el, node_n, testPos);
	   }
	 } // end of for(size_t node_m = 0;...

	 // should break the loop over elements because we found the right one
	 break;
       } // end if(isInside)

       vector<int>::iterator debug_elit = elit;
       ++debug_elit;
       if(debug_elit == endEl)
	 throw gError("TriangleCGMDInterpolation::compute", "last element not matching! None found!");

     } // end of loop over elements

     );
  //--------------END: loop over particles -------------------------- 

  // write m_shapeMatrix to file
  if(m_writeNN) {
    if(m_NNstepCounter == 0) {
      double* shapeMatrixCopy = new double[nOfNodes*nOfNodes];
      for(size_t i = 0; i < nOfNodes; ++i) {
	for(size_t j = 0; j < nOfNodes; ++j) {
	  shapeMatrixCopy[j+i*nOfNodes] = m_shapeMatrix[i][j];
	}}
      ofstream o;
      string fileName = string(m_NNname + ObjToString(m_NNfileCounter) + ".bin");
      o.open(fileName.c_str(), ios::out | ios::binary);
      o.write((char*) shapeMatrixCopy, nOfNodes*nOfNodes*sizeof(double));
      o.close();
    }
    ++m_NNstepCounter;
    if(m_NNstepCounter == m_writeNNevery) 
      m_NNstepCounter = 0;
    ++m_NNfileCounter;
  }


//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "Matrix Nmunu:START:");
//   for(size_t i = 0; i < m_shapeMatrix.size(); ++i) {
//     for(size_t j = 0; j < m_shapeMatrix[i].size(); ++j)
//       cout << m_shapeMatrix[i][j] << " ";
//     cout << endl;
//   }
//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "Matrix Nmunu:END");

//   time_t t0, t1;
//   std::time(&t0);

  // -------- START: matrix inverse ------------------

  for(size_t i = 0; i < nOfNodes; ++i) {
    for(size_t j = 0; j < nOfNodes; ++j) {
      GSL_MATRIX_SET(m_working_mat, i, j, m_shapeMatrix[i][j]);
      // cout << GSL_MATRIX_GET(m_working_mat, i, j) << " ";
    }
    // cout << endl;
  } 

//   std::time(&t1);
//   MSG_INFO("TriangleCGMDInterpolation::compute", "matrix inversion: copy took: " << difftime(t1, t0) << " seconds");

  // Compute the LU decomposition
  int signum;
  int error = gsl_linalg_LU_decomp(m_working_mat, m_permutation, &signum);
  if(error)
    throw gError("TriangleCGMDInterpolation::compute", "gsl_linalg_LU_decomp failed with return value " + ObjToString(error));

//   std::time(&t1);
//   MSG_INFO("TriangleCGMDInterpolation::compute", "matrix inversion: LU took: " << difftime(t1, t0) << " seconds");
      
  // Compute the inverse
  error = gsl_linalg_LU_invert(m_working_mat, m_permutation, m_inverse);
  if(error)
    throw gError("TriangleCGMDInterpolation::compute", "gsl_linalg_LU_invert failed with return value " + ObjToString(error));
  
//   std::time(&t1);
//   MSG_INFO("TriangleCGMDInterpolation::compute", "matrix inversion took: " << difftime(t1, t0) << " seconds");

//   std::time(&t1);
//   MSG_INFO("TriangleCGMDInterpolation::compute", "all about matrix inversion took: " << difftime(t1, t0) << " seconds");

//   cout << endl << endl << endl << endl << endl << endl 
//        << "#########################################################" 
//        << endl << endl << endl << endl << endl << endl;

//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "Matrix InverseNmunu:START:");
//   for(size_t i = 0; i < m_shapeMatrix.size(); ++i) {
//     for(size_t j = 0; j < m_shapeMatrix[i].size(); ++j)
//       cout << GSL_MATRIX_GET(m_inverse, i, j) << " ";
//     cout << endl;
//   }
//    MSG_DEBUG("TriangleCGMDInterpolation::compute", "Matrix InverseNmunu:END");

  // -------- END: matrix inverse ------------------


//   size_t debug_nOfNodes = m_elementsOfNode.size();
//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "Finished reading element nodes: nOfNodes = " << debug_nOfNodes); 
//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "This is in my array of elements:");
//   for(size_t element = 0; element < m_cornerIndices.size(); ++element) {
//     for(size_t node = 0; node < 3; ++node) {
//       cout << m_cornerIndices.at(element).at(node) << " ";
//       if(m_cornerIndices.at(element).at(node) > debug_nOfNodes)
// 	throw gError("TriangleCGMDInterpolation::compute", "AFTER matrix inverse: That's bad! nOfNodes = " + ObjToString(debug_nOfNodes));
//     }
//     cout << endl;
//   }

//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "element check after matrix inverse finished");
//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "again especially m_cornerIndices.at(element=864).at(node) = " << m_cornerIndices.at(864).at(0) << ", " << m_cornerIndices.at(864).at(1) << ", " << m_cornerIndices.at(864).at(2) );

//   abort();

  // -------- START: compute m_f ------------------

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_pColour,

     // the element where this particle is inside (was determined previously)
     size_t el = __iSLFE->tag.intByOffset(m_elementSlotOffset);
     vector<size_t>& cornerIndices = m_cornerIndices[el];

     point_t testPos = __iSLFE->r;
     
     // check periodic BCs separately for each direction, i.e. do we have to shift the particle into the element? We check by using the distance to node 0 of the element.
     for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
       // Is the distance to the (possibly shifted) corner 0 periodic? 
       // If yes, shift particle position.
       if(testPos[dir] - m_corners[el][0][dir] > m_boxSize[dir]/2)
	 testPos[dir] -= m_boxSize[dir];
       if(testPos[dir] - m_corners[el][0][dir] < -m_boxSize[dir]/2)
	 testPos[dir] += m_boxSize[dir];
     }

     // for each and every node mu
     for(size_t node_m = 0; node_m < nOfNodes; ++node_m) {
       // for each node nu of the element e_i
       for(size_t node_n = 0; node_n < 3; ++node_n) {	 
	 // f_{\mu i} += Ninv{\mu\nu} * N_{\nu i} 
	 // LHS was previously initialised to zero

	 m_f[__iSLFE->mySlot][node_m] += 
	   GSL_MATRIX_GET(m_inverse, node_m, cornerIndices[node_n])
	     * shapeF(el, node_n, testPos); 
       } // end loop over nodes node_n (nu)
     } // end loop over nodes node_m (mu)

     ); // end loop over particles

  // -------- END: compute m_f ------------------

  // write m_f to file
  if(m_write_f) {
    if(m_fStepCounter == 0) {
      double* fCopy = new double[nOfParticles*nOfNodes];
      for(size_t i = 0; i < nOfParticles; ++i) {
	for(size_t j = 0; j < nOfNodes; ++j) {
	  fCopy[j+i*nOfNodes] = m_f[i][j];
	}}
      ofstream o;
      string fileName = string(m_fName + ObjToString(m_fFileCounter) + ".bin");
      o.open(fileName.c_str(), ios::out | ios::binary);
      o.write((char*) fCopy, nOfParticles*nOfNodes*sizeof(double));
      o.close();
    }
    ++m_fStepCounter;
    if(m_fStepCounter == m_write_fEvery) 
      m_fStepCounter = 0;
    ++m_fFileCounter;
  }


  // ----------START: compute desired interpolated quantities ----------------


  // the following is now (2011-03-21) in call()
//   interpolate();


  // DEBUG-output
//     debugPrintM_F(0);
//     debugPrintM_F(1);
//     debugPrintM_F(10);
//     debugPrintM_F(100);
// //     debugPrintM_F(1000);



//   debugSumPerNode();

//   debugSumPerPart();

//    debugTotalSum();


//   MSG_DEBUG("TriangleCGMDInterpolation::compute", "FINISHED for timestep = " << timestep);

}

void TriangleCGMDInterpolation::interpolate() const
{
  Phase* phase = M_PHASE;
  Controller* controller = M_CONTROLLER;

//    MSG_DEBUG("TriangleCGMDInterpolation::interpolate", "CALLED at time " << controller->time());


  ParticleList& particles = phase->particles(m_pColour);
  ParticleList& nodes = phase->particles(m_nodeColour);
  size_t nOfParticles = phase->particles(m_pColour).size();
  size_t nOfNodes = phase->particles(m_nodeColour).size();

  size_t forceIndex = controller->forceIndex();

  // initialise to zero because we did not want this to be done automatically
  for(size_t mu = 0; mu < nOfNodes; ++mu) {
    Particle& node = nodes[mu];
    for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
//       node.tag.pointByOffset(m_momentum2Offset)[dir] = 0;
//       node.tag.pointByOffset(m_force2Offset)[dir] = 0;
      node.tag.pointByOffset(m_displacementOffset)[dir] = 0;
      node.tag.pointByOffset(m_momentum1Offset)[dir] = 0;
      node.tag.pointByOffset(m_force1Offset)[dir] = 0;
    }
  }

  for(size_t i = 0; i < nOfParticles; ++i) {
    // for m_shapeF we do not have to loop over all nodes
    Particle& p = particles[i];

//     MSG_DEBUG("TriangleCGMDInterpolation::interpolate", "i=" << i << endl << " momentum1-contrib=" << particles[i].v);

//     size_t el = p.tag.intByOffset(m_elementSlotOffset);

//     //---START: local interpolation----------------------
//     // only for nodes m of the same element (where i is in) values
//     // are different from zero 
//     for(size_t node_m = 0; node_m < 3; ++node_m) {
//       size_t global_m = m_cornerIndices[el][node_m];
//       Particle& node = nodes[global_m];
//       //----START:might vanish soon (2011-03-16)-------------
//       // compute here the locally interpolated momentum
//       node.tag.pointByOffset(m_momentum2Offset) +=
// 	// currently (2011-03-16) all masses=1!!! 
// 	//       shapeF(el, node_m, testPos) * i->v;
// 	m_shapeF[i][global_m] * p.v;
//       // compute here the locally interpolated force
//       node.tag.pointByOffset(m_force2Offset) +=
// 	// currently (2011-03-16) all masses=1!!! 
// 	//       shapeF(el, node_m, testPos) * i->force[forceIndex];
// 	m_shapeF[i][global_m] * p.force[forceIndex];
//       //----END:might vanish soon (2011-03-16)-------------
      
//     }
//     //---END: local interpolation----------------------




    for(size_t mu = 0; mu < nOfNodes; ++mu) {
      Particle& node = nodes[mu];
      // displacement
      // CHECK PERIODIC BCs!!!	
      point_t testDisp = p.r - m_rInit[i];
      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
	if(testDisp[dir] > m_boxSize[dir]/2)
	  testDisp[dir] -= m_boxSize[dir];
	if(testDisp[dir] < -m_boxSize[dir]/2)
	    testDisp[dir] += m_boxSize[dir];
      }
      node.tag.pointByOffset(m_displacementOffset) += 
	m_f[i][mu]*(testDisp);
      // non-locally interpolated momentum

//        MSG_DEBUG("TriangleCGMDInterpolation::interpolate", "mu=" << mu << ", i=" << i << endl << " momentum1-contrib=" << particles[i].v << endl << "weight=" << m_f[i][mu] << "momentum1 before:" << nodes[mu].tag.pointByOffset(m_momentum1Offset));

      node.tag.pointByOffset(m_momentum1Offset) +=
	// currently (2011-03-16) all masses=1!!! 
	m_f[i][mu]*p.v;

//        MSG_DEBUG("TriangleCGMDInterpolation::interpolate", "momentum1 after:" << nodes[mu].tag.pointByOffset(m_momentum1Offset));

      // locally interpolated momentum IS COMPUTED ABOVE
      // non-locally interpolated force

//       if(timestep == 9)       
// 	MSG_DEBUG("TriangleCGMDInterpolation::interpolate", "mu=" << mu << ", i=" << i << endl << " force1-contrib=" << particles[i].force[forceIndex] << endl << "weight=" << m_f[i][mu] << "force1 before:" << nodes[mu].tag.pointByOffset(m_force1Offset));
      node.tag.pointByOffset(m_force1Offset) +=
	// currently (2011-03-16) all masses=1!!! 
	m_f[i][mu]*p.force[forceIndex];

//       if(timestep == 9)
// 	MSG_DEBUG("TriangleCGMDInterpolation::interpolate", "force1 after:" << nodes[mu].tag.pointByOffset(m_force1Offset));


      // locally interpolated force IS COMPUTED ABOVE

    }
  }

}

void TriangleCGMDInterpolation::debugSumPerPart() const
{
  MSG_DEBUG("TriangleCGMDInterpolation::debugSumPerPart", "Sum of contributions for");
  for(size_t j = 0; j < m_f.size(); ++j) {
    double sum = 0;
    cout << "atom=" << j << ": ";
    for(size_t node = 0; node < m_f[j].size(); ++node) {
      sum += m_f[j][node];
    }
    cout << sum << endl;
  }
}

void TriangleCGMDInterpolation::debugSumPerNode() const
{
  MSG_DEBUG("TriangleCGMDInterpolation::debugSumPerNode", "Sum of contributions for");
   for(size_t node = 0; node < m_f[0].size(); ++node) {
     double sum = 0;
     cout << "mu=" << node << ": ";
     for(size_t j = 0; j < m_f.size(); ++j) {
       sum += m_f[j][node];
     }
     cout << sum << endl;
   }

}

void TriangleCGMDInterpolation::debugTotalSum() const
{
  MSG_DEBUG("TriangleCGMDInterpolation::debugTotalSum", "Total sum:");
  MSG_DEBUG("TriangleCGMDInterpolation::debugTotalSum", "sum_nodes->sum_particles:");
  
  double sum = 0;
  
  for(size_t j = 0; j < m_f.size(); ++j) {
    for(size_t node = 0; node < m_f[j].size(); ++node) {
    //      cout << "mu=" << node << ": ";
      sum += m_f[j][node];
    }
    //      cout << sum << endl;
  }
  
  cout << sum << endl;
}


void TriangleCGMDInterpolation::debugPrintM_F(size_t node) const
{
//   size_t pColour = M_MANAGER->getColour(m_pSpecies);
//   size_t nodeColour = M_MANAGER->getColour(m_nodeSpecies);

  point_t nodePos = M_PHASE->particles(m_nodeColour)[node].r;


  cout << endl << endl << endl << endl << endl << endl 
       << "#########################################################" 
       << endl << endl << endl << endl << endl << endl;

  MSG_DEBUG("TriangleCGMDInterpolation::debugPrintM_F", "f_mui(" << node << "):START:");
  // loop over all particles in the row of the node
  for(size_t j = 0; j < m_f.size(); ++j) {



    // compute cartesian distance
    point_t distVec = M_PHASE->particles(m_pColour)[j].r - nodePos;
  
    // periodic BCs
    if(m_periodic) {
      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
	if(distVec[dir] > m_boxSize[dir]/2)
	  distVec[dir] -= m_boxSize[dir];
	else if(distVec[dir] < -m_boxSize[dir]/2)
	  distVec[dir] += m_boxSize[dir];
      }
    }
    cout << distVec.abs()  << " " << m_f[j][node] << endl;
  }

  MSG_DEBUG("TriangleCGMDInterpolation::debugPrintM_F", "f_mui(" << node << ")::END");
}
