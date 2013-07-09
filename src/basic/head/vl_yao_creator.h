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




#ifndef __VL_YAO_CREATOR_H
#define __VL_YAO_CREATOR_H


using namespace std;


#include "pair_creator.h"
#include "cell.h"
#include "gen_f.h"
#include "pairdist.h"
#include "pair_list.h"
#include "wall_container.h"
#include "node_many_children.h"
// #include "output.h"
// #include "output_file.h"

#include <list>

#include <fstream>

#define M_PHASE  ((Phase*) m_parent)
#define M_MANAGER  M_PHASE->manager()

/*!
 * Implementation of the modified algorithm for Verlet-list creation based 
 * on partial updating in "dirty" cells as 
 * presented by Yao, Wang, Cheng. It is currently (2013-05-17) deactivated
 * since it is not clear whether it works correctly. For example it is not 
 * clear whether the recalculation of distances of old pairs in non-dirty 
 * cells is done at all. It does not seem so. In addition, the 
 * implementation is very inefficient and can not even compete with the 
 * linked-cell method without Verlet list.
 */
class VLYaoCreator: public PairCreator
{
private:

  void createDistances1stDirty
   (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist);


  void createDistances1stDirtyFrozen
   (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist);


  void createDistances2ndDirty
   (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist);


  void createDistances2ndDirtyFrozen
   (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist);


  void updateAndSwitch(Particle* i, int& size, /*int& sizeP,*/ int size_loc, int size_locP, int capacityP, PairList* pl);

protected:

  /*!
   * The skin size that is added to the cut-off
   */
  double m_skin_size;

  ofstream m_s;

  /*!
   * The columns from the parent data to be written
   */
  list<int> m_columns;

  DataFormat *m_input_format;

  /*!
   * Pointer to the parent phase
   */
  Phase *m_phase;

  //  /*!
  //   * The cut-off radius for pair creation. This is set to the largest cut-off radius
  //   * that has been requested by any force, thermostat, etc. + the skin size set by the verlet list creator
  //   */
  //  double m_cutoff;

  /*!
   * The name of the displacement
   */
  string m_displacement_name;

  /*!
   * Tag offset of the displacement
   */
  size_t m_displacement_o;

  /*!
   * tag offset of helper displacement stored at that time the \a Particle's
   * \a Cell was identified to require an update
   */
  size_t m_displacementOld_o;

  /*!
   * The tag offset to the first element of the neighbour list of a particle
   */
  vector<size_t> m_offset;

  /*!
   * The tag offset to the second element of the neighbour list of a particle
   */
  size_t m_int_offset;

  /*!
   * The number of slots that are needed for saving the neighbour of a particle
   */
  int m_n_slots;

  /*!
   * vector for the offsets of the colour pairs and the free pairs list
   */
  vector< vector<int> > m_offset_free;

  /*!
   * vector for the offsets of the colour pairs and the frozen pairs list
   */
  vector< vector<int> > m_offset_frozen;

  /*!
   * the capacities for each type of colour pair for free particles
   */
  vector< vector<int> > m_capacity_free;

  /*!
   * the capacities for each type of colour pair for frozen particles
   */
  vector< vector<int> > m_capacity_frozen;

  /*!
   * The first time the Verlet list has to be created from scratch. Later the pairs are only partially updated.
   * Initialise to "true" to create the Verlet list at the start of the simulation
   */
  bool m_create_now;

  /*!
   * The distance between the offset to the first neighbour pair index of a particle and the second neighbour pair index
   */
  int offset_dist;

  /*!
   * The "dirty" cells list (cells with dirty flag)
   */
  list<Cell*> m_dtf_cells;

  /*!
   * The "dirty" cells list (cells with dirty flag)
   */
//   list<CellLink*> m_dtf_links;

  /*!
   * Initialize
   */
  void init();

  /*!
   * Calculate the maximum displacement of a particle and set the cell dirty flags
   */
  void setDtf();

  /*!
   * Write the file header
   */
  virtual void writeHeader();


public:

  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   */
  VLYaoCreator(Phase* p);

  /*!
   * Destructor
   */
  virtual ~VLYaoCreator();

  /*!
   * Setup the linked list creator
   */
  virtual void setup();

  /*!
   * Inform the pairCreator that a new free pair has been added
   */
  virtual void pairAddedFree(Pairdist* pd);

  /*!
   * Inform the pairCreator that a new frozen pair has been added
   */
  virtual void pairAddedFrozen(Pairdist* pd);

  /*!
   * Create the pair fast verlet list
   */
  virtual void createDistances();

  /*!
   * The pair information in this pairCreator shouldn't be deleted
   */
  virtual void invalidatePositions();

  /*!
   * Needed to initialise the old displacement stored in the tag
   */
  virtual void setupAfterParticleCreation();

};


#endif
