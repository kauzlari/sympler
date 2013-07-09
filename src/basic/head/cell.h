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



#ifndef __CELL_H
#define __CELL_H

#include <map>
#include <list>
#include <utility>

#include "misc.h"
#include "wall.h"
#include "pair_creator.h"
#include "pairdist.h"
#include "pair_list.h"
#include "integrator.h"
#include "wall_container.h"
#include "integrator_position.h"
#include "integrator_velocity_verlet.h"
#include "integrator_velocity_verlet_disp.h"
#include "integrator_static.h"


#include "wall_triangle.h"

//---- Constants ----

#define NUM_NEIGHBORS 26
#define NUM_HALF_NEIGHBORS 13

#define INV_NEIGHBOR(n)   NUM_NEIGHBORS-n-1

#define OFFSET2NEIGHBOR(off, n) \
    n = (off[0]+1)*9 + (off[1]+1)*3 + (off[2]+1); \
    if (n > NUM_NEIGHBORS/2) \
    n--; \
    while(0)



/*---- CellLink ----
  A cell link describes a link between two cells and stores the
  distances of mutual particles.
*/



class Cell;
/* class PairCreator; */


/*!
 * This class connects two cells for cell subdivision. It handles
 * the creation of the pair list.
 */
class CellLink
{
protected:

  /*!
   * Holds a pointer to the first of the two cells connected by this link.
   */
  Cell *m_first;

  /*!
   * Holds a pointer to the second of the two cells connected by this links.
   */
  Cell *m_second;

  /*!
   * Alignment information for the two cells.
   * E.g. m_alignment = 4 means the second lies to the left in x-direction
   * of first.
   */
  int m_alignment;

  /*!
   * Stores the number of
   * cells that are active. 0 or 1 means the link
   * itself is deactivated.
   */
  int m_n_active_cells;

#ifdef ENABLE_PTHREADS
  /*!
   * Mutex for concurrent access to m_n_active_cells
   */
  pthread_mutex_t m_activation__mutex;
#endif

  /*!
   * Stores the relative distance of the
   * origin of the two cells.
   * Note: This is NOT first->corner1 - second->corner1!
   */
  point_t m_cell_dist;

  /*!
   *  Flag of the CellLink, so that it won't be used two times in the CellLink array
   */
//   bool m_linkUsed;

  /*!
   * Determines whether forces act on particles in the first
   * and/or in the second cell.
   */
  pair<bool, bool> m_acts_on;


#ifdef _OPENMP

  /*!
   * Every CellLink belongs to a thread: m_thread.
   */
  int m_thread;

#endif


public:

  /*!
   * Default constructor. Should not be used.
   */
  CellLink();

  /*!
   * Constructor. Calls set for in order to set the parameters
   * @param first First of the two connected cells
   * @param second Second of the two connected cells
   * @param alignment Where does the second cell lie with respect to the first cell?
   * @param acts_on_first Determines whether forces act on particles in the first cell
   * @param acts_on_second Determines whether forces act on particles in the second cell
   * @see set()
   */
  CellLink(Cell *first, Cell *second, int alignment,
           bool acts_on_first = true, bool acts_on_second = true);

  /*!
   * Copy constructor.
   * @param copy Link to be copied.
   */
  CellLink(const CellLink &copy);

  /*!
   * Destructor
   */
  virtual ~CellLink();

  /*!
   * Set the connected cells and the acts-on information.
   * @param first First of the two connected cells
   * @param second Second of the two connected cells
   * @param alignment Where does the second cell lie with respect to the first cell?
   * @param acts_on_first Determines whether forces act on particles in the first cell
   * @param acts_on_second Determines whether forces act on particles in the second cell
   */
  virtual void set(Cell *first, Cell *second, int alignment,
                   bool acts_on_first = true, bool acts_on_second = true);

  /*!
   * Notifies the link about the activation of one of the two connected cells.
   * If both cells are active, the link is active.
   */
  virtual void cellActivated();

  /*!
   * Notifies the link about the deactivation of one of the two connected cells.
   * If both cells are active, the link is active.
   */
  virtual void cellDeactivated();

  /*!
   * Determine the pointer to the other cell. I.e. usually
   * a Cell will call link->other(this) in order to determine
   * its neighbor.
   * @param c current cell
   */
  inline Cell* other(Cell* c) {
    if (c == m_first)
      return m_second;
    else
      return m_first;
  }

  /*!
   * Return a pointer to the first cell.
   */
  inline Cell *first() {
    return m_first;
  }

  /*!
   * Return a pointer to the second cell.
   */
  inline Cell *second() {
    return m_second;
  }

  /*!
   * Create the pair list.
   * @param t Thread number of this call to createDistances. For multithreaded compilation.
   */
#ifdef _OPENMP   
  virtual void createDistances(int thread_no);
#else
  virtual void createDistances();
#endif  

  /*!
   * Create the verlet pair list.
   * @param t Thread number of this call to createDistances. For multithreaded compilation.
   */
//   virtual void createVerletDistances(size_t t);

  /*!
   * Return actsOn information.
   */
  const pair<bool, bool> &actsOn() const {
    return m_acts_on;
  }

  /*!
   * Return the distance between the two cells, connected by this link in real coordinates
   */
  virtual point_t &mCellDist() {
    return m_cell_dist;
  }

  /*!
   * Next active link.
   */
  CellLink *next;

  /*!
   * Previous active link.
   */
  CellLink *prev;

#ifdef _OPENMP
  /*!
   * Which thread does this CellLink belong to?
   */
  virtual int &mThread() {
  	return m_thread;
  }
#endif

};



class ManagerCell;

/*!
 * A cell represents one simulation cell for cell subdivision.
 */
class Cell: public cuboid_t
{
protected:
  /*!
   * Neighbors of this cell. Can be more than one,
   * because of all this inlet, outlet stuff.
   */
  list<CellLink*> m_neighbors[NUM_NEIGHBORS];

  /*!
   * Outlets tell the cell where the particles go
   * when they leave this cell.
   */
  list<Cell*> m_outlets[NUM_NEIGHBORS];

  /*!
   * Indirect outlets are outlets that do not share a side
   * with the cell.
   */
  list< pair<int, Cell*> > m_indirect_outlets;

  /*!
   * Pairs that exist within this cell
   */
  CellLink* m_local_link;

  /*!
   * Stores all walls that intersect with this cell.
   */
  list<Wall*> m_walls;

  /*!
   * Stores all walls of this and all neighboring cells.
   */
  list<Wall*> m_all_walls;

  /*!
   * List of all particles currently in a cell.
   */
  vector<list<Particle*> > m_particles;

  /*!
   * List of all frozen particles currently in a cell.
   */
  vector<list<Particle*> > m_frozen_particles;

#ifdef ENABLE_PTHREADS
  /*!
   * Mutex for access to m_injected_particles
   */
  vector<pthread_mutex_t> m_injected_particles__mutex;
#endif

  /*!
   * List of all FREE particles to be newly injected.
   */
  vector<list<Particle*> > m_injected_particles;

  /*!
   * Number of particles in this cell.
   */
  size_t m_n_particles;

  /*!
   * Cell's dirty flag
   */
  bool m_dtf;

  /*!
   * A pointer to the manager.
   */
  ManagerCell *m_manager;

  /*!
   * The group all particles in this cell belong to.
   */
  size_t m_group;

  /*!
   *  Flag of the Cell, so that it won't be used two times in the Cell array
   */
  bool m_cellUsed;

  /*!
   * Notification if a particle leaves this cell.
   * Does nothing in default implementation.
   * @param offset In which direction did the particle leave
   * @param n Neigbor index of that direction
   */
  virtual void particleLeftCell(const int_point_t &offset, int n) {
  }

  /*!
   * Connect two neighboring cells.
   * "Where" is the position parameter according to Cell::c_offsets
   * @param neighbor The neighboring cell
   * @param where Neighbor index of the cell alignment
   * @param first Determines whether forces act on particles in the first cell
   * @param second Determines whether forcee act on particle in the second cell (the one given by \a neighbor)
   */
  virtual void establishLink(Cell *neighbor, int where, bool first, bool second);


public:

  /*!
   * Check if the particle collides with a wall.
   * This solves the quadratic equation exactly.
   * @param p Check collision for this particle
 * @param r make possible changes in this position (may be != p->r)
 * @param v make possible changes in this position (may be != p->v)
 * @param force Force on this particle
   */
  virtual void doCollision(Particle *p, point_t& r, point_t& v, const point_t &force, IntegratorPosition *integratorP);

  /*!
   * Constructor.
   * @param mgr Pointer to the manager
   * @param group Group all particles in this cell belong to
   */
  Cell(ManagerCell *mgr, int group = 0);

  /*!
   * Constructor.
   * @param mgr Pointer to the manager
   * @param c1 Top left corner of the cell
   * @param c2 Bottom right corner of the cell
   * @param group Group all particles in this cell belong to
   */
  Cell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int group = 0);

  /*!
   * Destructor.
   */
  virtual ~Cell();

  /*!
   * Basic initialization
   */
  virtual void init();

  /*!
   * Handles cell activation,
   * i.e. the cell is added to the
   * linked list (\a prev/\a next fields)
   */
  virtual void activate();

  /*!
   * Handles cell deactivation,
   * i.e. the cell is removed from the
   * linked list (\a prev/\a next fields)
   */
  virtual void deactivate();

  /*!
   * Setup walls looks for walls in neighboring cells and adds
   * them to the \a m_all_walls list. Furthmore, it looks for indirect
   * outlets, i.e. peridicities. In these cases, additional care
   * has to be taken when calculating collisions with walls.
   */
  virtual void setupWalls();

  /*!
   * \a neighbor stores a pointer to the neighboring cell whereas
   * \a where gives the index of the neighbor's position according
   * to c_offsets.
   * @param neighbor Neighboring cell
   * @param where Alignment of the neighboring cell
   */
  virtual void addNeighbor(Cell *neighbor, int where);

  /*!
   * Adds an unidirectional neighbor, that means
   * particles can fly into \a neighbor and particles in
   * \a neighbor feels forces from particles in this cell,
   * but not the other way round.
   * @param neighbor Neighboring cell
   * @param where Alignment of the neighboring cell
   */
  virtual void addOutlet(Cell *neighbor, int where);

  /*!
   * Adds an outlet without force interaction,
   * i.e. the particle do not "see" each other
   * @param neighbor Neighboring cell
   * @param where Alignment of the neighboring cell
   */
  virtual void addPeriodic(Cell *neighbor, int where);

  /* clearOutlets deletes all outlets in a specific direction
  */
  virtual void clearOutlets(int where);

  /*!
   * Deletes all neighbors in a specific direction
   * @param where Direction in which to delete neighbors
   */
//  virtual void clearNeighbors(int where);

  /*!
   * Interface between Bondaries and Cells. Add a wall to this
   * cell.
   * @param p Pointer to the wall object
   **/
  virtual void addWall(Wall *p) {
    m_walls.push_back(p);
  }

  /*!
   * Return the flag that tells whether the CellLink has been used in the dirty cellLinks array
   */
  virtual bool& cellUsed() {
    return m_cellUsed;
  }

  /*!
   * Set the container
   * @param container Pointer to the container holding the walls
   */
  virtual void assignContainer(WallContainer *container);

  /*!
   * Advance positions of the particle within this cell using
   * the Velocity-Verlet algorithm
   * @param integrator Integrator to use for the update. Fixme!!! Only Velocity-Verlet supported.
   */
  virtual void updatePositions(IntegratorPosition *integrator);

  /*!
   * Add a particle to this cell.
   * @param colour Colour of the particle
   * @param p Pointer to the particle
   */
  void injectFree(size_t colour, Particle *p);

  /*!
   * Add a frozen particle to this cell.
   * @param colour Colour of the particle
   * @param p Pointer to the particle
   */
  void injectFrozen(size_t colour, Particle *p);

  /*!
   * Commit the injection of free particles.
   */
  virtual void commitInjections();

  /*!
   * Clear all tags which contain the force factors.
   */
  virtual void clearTags();

  /*!
   * Return a list of all neighbors in direction \a i
   * @param i Direction
   */
  list<CellLink*> &neighbors(int i) {
    return m_neighbors[i];
  }

  /*!
   * Return a list of all outlets in direction \a i
   * @param i Direction
   */
  list<Cell*> &outlets(int i) {
    return m_outlets[i];
  }

  /*!
   * Return a list of all particle of color \a c
   * @param c Color
   */
  list<Particle*> &particles(size_t c) {
    return m_particles[c];
  }

  /*!
   * Return the CellLink to the cell itself
   */
  virtual CellLink* localLink() {
    return m_local_link;
  }

  /*!
   * Return a list of all frozen particle of color \a c
   * @param c Color
   */
  list<Particle*> &frozenParticles(size_t c) {
    return m_frozen_particles[c];
  }

  /*!
   * Erase a free particle
   * @param e Iterator to the particle
   */
  void eraseParticle(list<Particle*>::iterator e);

  /*!
   * Returns true if the cell has particles
   */
  bool hasParticles() const;

 /*!
  * Returns true if createDistances has to be called for this cell
  */
  virtual bool& dirtyFlag()
  {
    return m_dtf;
  }

  /*!
   * Returns true if the cell has frozen particles
   */
  bool hasFrozenParticles() const;

  /*!
   * Return the manager
   */
  virtual ManagerCell* manager() {
    return m_manager;
  }

  /*!
   * Find the neighbor corresponding to position \a r
   * @param r Position
   */
  int whichNeighbor(const point_t &r) {
    int n;
    int_point_t off;

    for (int j = 0; j < SPACE_DIMS; j++)
      if (r[j] < corner1[j])
        off[j] = -1;
      else if (r[j] >= corner2[j])
        off[j] = 1;

    OFFSET2NEIGHBOR(off, n);

    return n;
  }

  /*!
   * Check if the particle hit a wall
   * @param p Particle of interest
   * @param force Force on the particle
   * @param hit True on return if a hit did occur
   * @param t_travelled Time travelled to that hit
   * @param hit_pos Impact position on the wall
   * @param wall Pointer to the wall that has been hit
   */
  void checkForHit(const Particle *p, const point_t &force,
                   bool &hit, double &t_travelled, point_t &hit_pos, Wall *&wall, IntegratorPosition* integratorP) {
    for (list<Wall*>::iterator i = m_all_walls.begin(); i != m_all_walls.end(); i++) {
      double t;
      point_t h;

      //!!!!
//      ((WallTriangle*) (*i))->initHelpers();

      if ((*i)->hit(p, force, t, h, integratorP)) {
#ifdef TRACK_PARTICLE
        if (p->mySlot == TRACK_PARTICLE) {
          cout << "??? hit @ " << h << " after time " << t << " on wall:" << endl;
          cout << (*i)->toString() << endl;
        }

#endif

        if (t < t_travelled) {
          hit = true;
          t_travelled = t;
          hit_pos = h;
          wall = *i;
        }
      }
    }
  }

  /*!
   * Return the group of this cell
   */
  virtual int group() const {
    return m_group;
  }

  int intersections(const line_t &l, int dir) {
    return WallContainer::intersections(m_walls, l, dir);
  }

  /*!
  * Return all the walls of this cell and its neighbouring cells
  * NOTE: indirect outlet-cells (e.g. periodic cells) are NOT condidered!!!
  */
  list<Wall*>* allWalls()
  {
    return &m_all_walls;
  }
  
  list< pair<int, Cell*> >* indirectOutlets() {
    return &m_indirect_outlets;
  }

  /* Constants */
  static const int_point_t c_offsets[NUM_NEIGHBORS];

  /*!
   * For external use, e.g. when cell sizes are
   * being calculated.
   */
  int_point_t tag;

  /*!
   * Pointer to the next active cell.
   */
  Cell *next;

  /*!
   * Pointer to the previous active cell.
   */
  Cell *prev;
};


#ifdef _OPENMP
inline void addPair(vector<PairList> &distances, double cutoff_sq,
                    int dir, Cell *first_c, Cell *second_c,
                    Particle *first_p, Particle *second_p,
                    bool ao_f, bool ao_s,
                    point_t &cell_dist,
                    int thread_no)
{
#else
inline void addPair(vector<PairList> &distances, double cutoff_sq,
                    int dir, Cell *first_c, Cell *second_c,
                    Particle *first_p, Particle *second_p,
                    bool ao_f, bool ao_s,
                    point_t &cell_dist)
{
#endif
  dist_t d;
  d.abs_square = 0;
  for (int _i = 0; _i < SPACE_DIMS; _i++) {
    d.cartesian[_i] = -dir*cell_dist[_i]
      + first_p -> r[_i] - first_c -> corner1[_i]
      - second_p -> r[_i] + second_c -> corner1[_i];
    d.abs_square += d.cartesian[_i]*d.cartesian[_i];
    }
  /* Take care: The order of *j, *i defines the direction d \
     is pointing to. */
  if (d.abs_square < cutoff_sq) {
    d.abs = sqrt(d.abs_square);
#ifdef _OPENMP    
    distances[thread_no].newPair().set(d, first_p, second_p, ao_f, ao_s);
#else    
    distances[PairCreator::counterTN].newPair().set(d, first_p, second_p, ao_f, ao_s);
#endif

  }
}

#endif
