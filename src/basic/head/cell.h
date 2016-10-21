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
#include "colour_pair.h"
#include "pairdist.h"
#include "pair_list.h"
#include "integrator.h"
#include "wall_container.h"
#include "integrator_position.h"
#include "integrator_velocity_verlet.h"
#include "integrator_velocity_verlet_disp.h"
#include "integrator_static.h"


#include "geometric_primitives.h"
#include "manager_cell.h"
#include "wall_triangle.h"


/*---- AbstractCellLink ----
  A cell link describes a link between two cells and stores the
  distances of mutual particles.
*/



class Cell;
/* class PairCreator; */


/*!
 * This class connects two cells for cell subdivision. It handles
 * the creation of the pair list.
 */
class AbstractCellLink
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
   *
   */
  bool m_cross_regions;

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
   *  Flag of the AbstractCellLink, so that it won't be used two times in the AbstractCellLink array
   */
//   bool m_linkUsed;

  /*!
   * Determines whether forces act on particles in the first
   * and/or in the second cell.
   */
  pair<bool, bool> m_acts_on;

#ifdef _OPENMP

  /*!
   * Every AbstractCellLink belongs to a thread: m_thread.
   */
  int m_thread;

#endif


public:

  /*!
   * Default constructor. Should not be used.
   */
  AbstractCellLink();

  /*!
   * Constructor. Calls set for in order to set the parameters
   * @param first First of the two connected cells
   * @param second Second of the two connected cells
   * @param alignment Where does the second cell lie with respect to the first cell?
   * @param acts_on_first Determines whether forces act on particles in the first cell
   * @param acts_on_second Determines whether forces act on particles in the second cell
   * @see set()
   */
  AbstractCellLink(Cell *first, Cell *second, int alignment,
           bool acts_on_first = true, bool acts_on_second = true, bool cross_regions = false);

  /*!
   * Copy constructor.
   * @param copy Link to be copied.
   */
  AbstractCellLink(const AbstractCellLink &copy);

  /*!
   * A pointer to the manager.
   */
  ManagerCell *m_manager;

  /*!
   * Destructor
   */
  virtual ~AbstractCellLink();

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
  Cell* other(Cell* c) {
    if (c == m_first)
      return m_second;
    else
      return m_first;
  }

  /*!
   * Return a pointer to the first cell.
   */
  Cell *first() {
    return m_first;
  }

  /*!
   * Return a pointer to the second cell.
   */
  Cell *second() {
    return m_second;
  }

  /*!
   * Create the pair list.
   * @param t Thread number of this call to createDistances. For multithreaded compilation.
   */
#ifdef _OPENMP   
  virtual void doForSameColors(ColourPair* cp, int colour, double cutoff_sq, int thread_no) = 0;
#else
  virtual void doForSameColors(ColourPair* cp, int colour, double cutoff_sq) = 0;
#endif

#ifdef _OPENMP   
  virtual void doForDiffColors(ColourPair* cp, int c1, int c2, double cutoff_sq, int thread_no) = 0; 
#else
  virtual void doForDiffColors(ColourPair* cp, int c1, int c2, double cutoff_sq) = 0;
#endif  

#ifdef _OPENMP   
  virtual void createDistances(int thread_no) = 0;
#else
  virtual void createDistances() = 0;
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
   * Next active link.
   */
  AbstractCellLink *next;

  /*!
   * Previous active link.
   */
  AbstractCellLink *prev;

#ifdef _OPENMP
  /*!
   * Which thread does this AbstractCellLink belong to?
   */
  virtual int &mThread() {
  	return m_thread;
  }
#endif

};

/*!
 * This class connects two cells for cell subdivision. It handles
 * the creation of the pair list.
 */

template <typename AddPairCheck_x>
AbstractCellLink *establishLinkTemplate(ManagerCell *a_manager, Cell *first, Cell *second, int where, bool acts_on_first, bool acts_on_second, bool cross_regions =  false);
template <typename AddPairCheck_x>
class CellSelfLink;
template <typename AddPairCheck_x>
class CellLink;
/*!
 * A cell represents one simulation cell for cell subdivision.
 */
class Cell: public cuboid_t
{
friend class ManagerCell;
friend class ManagerCell2x;
template <typename AddPairCheck_x> friend class CellLink;
template <typename AddPairCheck_x> friend class CellSelfLink;
friend class BoundaryCell; //for accessing m_walls in setupwalls
protected:
  /*!
   * Neighbors of this cell. Can be more than one,
   * because of all this inlet, outlet stuff.
   */
  vector<vector<AbstractCellLink*> > m_neighbors;

  /*!
   * Outlets tell the cell where the particles go
   * when they leave this cell.
   */
  vector<vector<Cell*> > m_outlets;

  /*!
   * Indirect outlets are outlets that do not share a side
   * with the cell.
   */
  list<pair <point_t,  Cell*> > m_indirect_outlets;

  /*!
   * Pairs that exist within this cell
   */
  AbstractCellLink* m_local_link;

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
   * Stores the relative distance of the
   * origin of the two cells.
   * Note: This is NOT first->corner1 - second->corner1!
   */
  vector<vector<point_t> > m_cell_dist;

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
  template<typename AddPairCheck_x>
  AbstractCellLink *establishLink(Cell *neighbor, int where, bool first, bool second, bool cross_regions = false)
  {
    return establishLinkTemplate<AddPairCheck_x>(m_manager, this, neighbor, where, first, second, cross_regions);
  }

  template<typename AddPairCheck_x>
  AbstractCellLink *establishLink(Cell *neighbor, int where, bool cross_regions = false)
  {
    return establishLinkTemplate<AddPairCheck_x>(m_manager, this, neighbor, where, true, true, cross_regions);
  }



  virtual bool doCollisionIndirectOutlets(Particle *p, const point_t &force, point_t &r, double &t_travelled, point_t &hit_pos, Wall *&wall, IntegratorPosition *integratorP) {return false;}; 

  inline void commonConstructor(int sameDirOutlets, int_point_t a_tag = {0,0,0});
public:

  /*!
   * Check if the particle collides with a wall.
   * This solves the quadratic equation exactly.
   * @param p Check collision for this particle
   * @param r make possible changes in this position (may be != p->r)
   * @param v make possible changes in this position (may be != p->v)
   * @param force Force on this particle
   */
  void doCollision(Particle *p, point_t& r, point_t& v, const point_t &force, IntegratorPosition *integratorP);

  /*!
   * Constructor.
   * @param mgr Pointer to the manager
   * @param group Group all particles in this cell belong to
   */
  Cell(ManagerCell *mgr, int group = 0, int sameDirOutlets = 1);

  /*!
   * Constructor.
   * @param mgr Pointer to the manager
   * @param c1 Top left corner of the cell
   * @param c2 Bottom right corner of the cell
   * @param group Group all particles in this cell belong to
   */
  Cell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int group = 0, int sameDirOutlets = 1);

  /*!
   * Constructor.
   * @param mgr Pointer to the manager
   * @param c1 Top left corner of the cell
   * @param c2 Bottom right corner of the cell
   * @param a_tag Position of the cell
   * @param group Group all particles in this cell belong to
   */
  Cell(ManagerCell *mgr, cuboid_t cuboid, int_point_t a_tag, int group = 0, int sameDirOutlets = 1);


  /*!
   * Constructor.
   * @param mgr Pointer to the manager
   * @param c1 Top left corner of the cell
   * @param c2 Bottom right corner of the cell
   * @param a_tag Position of the cell
   * @param group Group all particles in this cell belong to
   */
  Cell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int_point_t a_tag, int group = 0, int sameDirOutlets = 1);

  /*!
   * Destructor.
   */
  virtual ~Cell();

  /*!
   * Basic initialization
   */
  template<typename AddPairCheck_x>
  void init()
  {
    m_local_link = new CellSelfLink<AddPairCheck_x>(this, this, -1);

    m_manager->m_links.push_back(m_local_link);
  }

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

  inline void addToMOutlets(Cell *a, Cell *neighbor, int where, bool cross_regions = false)
  {
        a->m_outlets[where][cross_regions] = neighbor;
  //      MSG_DEBUG("addToMOutlets", where << " " << cross_regions << " " <<  neighbor);
  };

  /*!
   * \a neighbor stores a pointer to the neighboring cell whereas
   * \a where gives the index of the neighbor's position according
   * to c_offsets.
   * @param neighbor Neighboring cell
   * @param where Alignment of the neighboring cell
   */
  template<typename AddPairCheck_x>
  AbstractCellLink *addNeighbor(Cell *neighbor, int where, bool cross_regions = false)
  {
      
      AbstractCellLink *link = this->establishLink<AddPairCheck_x>(neighbor, where, true, true, cross_regions);
      if (m_manager->m_divby > 1) where = m_manager->c_2x_1x(where);

    //    MSG_DEBUG("addNeighbor", where << " " << " " <<  neighbor);
      addToMOutlets(this, neighbor, where);
    //    MSG_DEBUG("addNeighbor", INV_DIRECT_NEIGHBOR(where) << " " << " " <<  this);
      addToMOutlets(neighbor, this, INV_DIRECT_NEIGHBOR(where));
      
      return link;
  }

  /*!
   * Adds a unidirectional neighbor, that means
   * particles can fly into \a neighbor and particles in
   * \a neighbor feel forces from particles in this cell,
   * but not the other way round.
   * @param neighbor Neighboring cell
   * @param where Alignment of the neighboring cell
   */
  template<typename AddPairCheck_x>
  AbstractCellLink *addOutlet(Cell *neighbor, int where, bool cross_regions = false)
  {
      AbstractCellLink *link = establishLink<AddPairCheck_x>(neighbor, where, false, true, cross_regions);
       
      if (m_manager->m_divby > 1) where = m_manager->c_2x_1x(where);
      addToMOutlets(this, neighbor, where, cross_regions);
      
      return link;
  }

  /*!
   * Adds an outlet without force interaction,
   * i.e. the particle do not "see" each other
   * @param neighbor Neighboring cell
   * @param where Alignment of the neighboring cell
   */
  virtual void addPeriodic(Cell *neighbor, int where, bool cross_regions = false);

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
   * Return the flag that tells whether the AbstractCellLink has been used in the dirty cellLinks array
   */
  virtual bool& cellUsed() {
    return m_cellUsed;
  }

  /*!
   * Set the container
   * @param container Pointer to the container holding the walls
   */
  virtual void assignContainer(WallContainer *container);

  virtual bool emitIntoOutlets(Particle* p, size_t colour, int_point_t off, int n, IntegratorPosition *integrator, Phase *phase = NULL);

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
  inline virtual void commitInjections();

  /*!
   * Clear all tags which contain the force factors.
   */
  inline virtual void clearTags();

  /*!
   * Return a list of all neighbors in direction \a i
   * @param i Direction
   */
  vector<AbstractCellLink*> &neighbors(int i) {
    return m_neighbors[i];
  }

  /*!
   * Return a list of all outlets in direction \a i
   * @param i Direction
   */
  vector<Cell*> &outlets(int i) {
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
   * Return the AbstractCellLink to the cell itself
   */
  virtual AbstractCellLink* localLink() {
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
  int whichNeighbor(const point_t &r);
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

#ifdef TRACK_PARTICLE
        if (p->mySlot == TRACK_PARTICLE)
          cout << (*i)->toString() << endl;
#endif
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
  
  list< pair< point_t, Cell*> >* indirectOutlets() {
    return &m_indirect_outlets;
  }

  /* Constants */

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

class BoundaryCell : public Cell
{
  public:
  const struct region_t *m_region;

	  BoundaryCell(ManagerCell *mgr, region_t *r, int group = 0);
	  BoundaryCell(ManagerCell *mgr, cuboid_t cuboid, int_point_t a_tag, region_t *r, int group = 0);
	  BoundaryCell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int_point_t a_tag, region_t *r, int group = 0);
  
  virtual bool doCollisionIndirectOutlets(Particle *p, const point_t &force, point_t &r, double &t_travelled, point_t &hit_pos, Wall *&wall, IntegratorPosition *integratorP) override; 
  virtual void setupWalls();
  virtual bool emitIntoOutlets(Particle* p, size_t colour, int_point_t off, int n, IntegratorPosition *integrator, Phase *phase);
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
};

template<typename AddPairCheck_x>
inline static void addPair(vector<PairList> &distances, double cutoff_sq,
  Cell *first_c, Cell *second_c,
  Particle *first_p, Particle *second_p,
  bool ao_f, bool ao_s,
  point_t &cell_dist,
  int thread_no);

//#define ADDPAIR(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist) do { dist_t d; \
  d.abs_square = 0; \
  for (int _i = 0; _i < SPACE_DIMS; _i++) { \
      d.cartesian[_i] = cell_dist[_i] \
        + (first_p) -> r[_i] - (first_c) -> corner1[_i] \
        - (second_p) -> r[_i] + (second_c) -> corner1[_i]; \
      AddPairCheck_x::check(first_c, _i, d.cartesian[_i], first_c->corner2[_i] - first_c->corner1[_i]); \
    d.abs_square += d.cartesian[_i]*d.cartesian[_i]; \
  } \
  if (d.abs_square < cutoff_sq) { \
    d.abs = sqrt(d.abs_square); \
      distances[thread_no].newPair().set(d, first_p, second_p, ao_f, ao_s); \
  }} while(0)



//#ifdef _OPENMP
  #define ADDPAIR(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist) do addPair<AddPairCheck_x>(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist, thread_no); while(0)
//#else
//  #define ADDPAIR(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist) do addPair<AddPairCheck_x>(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist, thread_no); while(0) //PairCreator::counterTN); while(0)
//#endif

//#define CREATE_DISTANCES_FOR_SAME_PARTICLE_LIST(distances, cutoff_sq, c, p, thread_no) do {
#define CREATE_DISTANCES_FOR_SAME_PARTICLE_LIST(distances, cutoff_sq, c, p) do { \
  list<Particle*>::iterator p_end = p.end(); \
  for (list<Particle*>::iterator i = p.begin(); i != p_end; ++i) { point_t precomputed = (*i) -> r - c -> corner1; \
    list<Particle*>::iterator j = i; \
    for (++j; j != p_end; ++j) { \
      ADDPAIR \
         (distances, cutoff_sq, \
         c, c, \
         *i, *j, \
         true, true, \
         precomputed); \
    } \
  }}while(0) 

//#ifdef _OPENMP
//  #define CREATE_DISTANCES_FOR_SAME_PARTICLE_LIST(distances, cutoff_sq, c, p) do createDistancesForSameParticleList_<AddPairCheck_x>(distances, cutoff_sq, c, p, thread_no); while (0)
//#else
//  #define CREATE_DISTANCES_FOR_SAME_PARTICLE_LIST(distances, cutoff_sq, c, p) do createDistancesForSameParticleList_<AddPairCheck_x>(distances, cutoff_sq, c, p); while (0)
//#endif
template<typename AddPairCheck_x>
inline static void createDistancesForSameParticleList_
  (vector<PairList> &distances,
   double cutoff_sq,
   Cell *c,
   list<Particle*> &p
#ifdef _OPENMP
   ,int thread_no
#endif
   )
{
#ifndef _OPENMP   
  int thread_no = PairCreator::counterTN;
#endif
  list<Particle*>::iterator p_end = p.end();
  point_t zero_cell_dist = {0,0,0};
  for (list<Particle*>::iterator i = p.begin(); i != p_end; ++i) {
    list<Particle*>::iterator j = i;
    for (++j; j != p_end; ++j) {
      ADDPAIR//<AddPairCheck_x>
         (distances, cutoff_sq,
         c, c,
         *i, *j,
         true, true,
//         zero_cell_dist, thread_no);
         zero_cell_dist);
    }
  }
}

#define CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist) { \
              list<Particle*>::iterator p1_end = first_p.end(); \
              list<Particle*>::iterator p2_end = second_p.end(); \
              for (list<Particle*>::iterator i = first_p.begin(); i != p1_end; ++i) {point_t precomputed = cell_dist + (*i) -> r - first_c -> corner1; \
                for (list<Particle*>::iterator j = second_p.begin(); j != p2_end; ++j) { \
                  ADDPAIR \
                    (distances, cutoff_sq, \
                     first_c, second_c, \
                     *i, *j, \
                     ao_f, ao_s, \
                     precomputed); \
                } \
              }}while(0)

//#ifdef _OPENMP
//  #define CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist) do createDistancesForDifferentParticleLists_<AddPairCheck_x>(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist, thread_no); while (0)
//#else
//  #define CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist) do createDistancesForDifferentParticleLists_<AddPairCheck_x>(distances, cutoff_sq, first_c, second_c, first_p, second_p, ao_f, ao_s, cell_dist); while (0)
//#endif
template<typename AddPairCheck_x>
inline static void createDistancesForDifferentParticleLists_
  (vector<PairList> &distances,
   double cutoff_sq,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t cell_dist 
#ifdef _OPENMP
   ,int thread_no
#endif
   )
{
#ifndef _OPENMP   
  int thread_no = PairCreator::counterTN;
#endif


  list<Particle*>::iterator p1_end = first_p.end();
  list<Particle*>::iterator p2_end = second_p.end();

  for (list<Particle*>::iterator i = first_p.begin(); i != p1_end; ++i) {
    for (list<Particle*>::iterator j = second_p.begin(); j != p2_end; ++j) {
      ADDPAIR//<AddPairCheck_x
        (distances, cutoff_sq,
         first_c, second_c,
         *i, *j,
         ao_f, ao_s,
//         cell_dist, thread_no);
         cell_dist);        
    }
  }
}

template <typename AddPairCheck_x>
class CellSelfLink: public AbstractCellLink
{
  private:
    point_t m_zero_cell_dist;
  public:
  CellSelfLink(Cell *first, Cell *second, int alignment, \
    bool acts_on_first = true, bool acts_on_second = true/*, bool cross_regions = false*/): \
      AbstractCellLink(first, second, alignment, acts_on_first, acts_on_second, false){ m_zero_cell_dist = {0,0,0};};
  /*!
   * Create the pair list.
   * @param t Thread number of this call to createDistances. For multithreaded compilation.
   */

#ifdef _OPENMP   
  virtual void doForSameColors(ColourPair* cp, int colour, double cutoff_sq, int thread_no) final override {
#else
  virtual void doForSameColors(ColourPair* cp, int colour, double cutoff_sq) final override {
    int thread_no = PairCreator::counterTN;
#endif
    CREATE_DISTANCES_FOR_SAME_PARTICLE_LIST//<AddPairCheck_x>
      (cp->freePairs(),
       cutoff_sq,
       m_first,
       m_first->particles(colour));
//	   m_first->particles(c1), thread_no);

    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->frozenPairs(),
       cutoff_sq,
       m_first, m_first,
       m_first->particles(colour), m_first->frozenParticles(colour),
       true, false,
//	   ((point_t) {0,0,0}));
//	   ((point_t) {0,0,0}), thread_no);
//	   zero_cell_dist, thread_no);
       m_zero_cell_dist);
  }

#ifdef _OPENMP   
  virtual void doForDiffColors(ColourPair* cp, int c1, int c2, double cutoff_sq, int thread_no) final override {
#else
  virtual void doForDiffColors(ColourPair* cp, int c1, int c2, double cutoff_sq) final override {
    int thread_no = PairCreator::counterTN;
#endif
    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->freePairs(),
       cutoff_sq,
       m_first, m_first,
       m_first->particles(c1), m_first->particles(c2),
       true, true,
  //             ((point_t) {0,0,0}));
  //	     ((point_t) {0,0,0}), thread_no);
  //             zero_cell_dist, thread_no);
       m_zero_cell_dist);

    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->frozenPairs(),
       cutoff_sq,
       m_first, m_first,
       m_first->particles(c1), m_first->frozenParticles(c2),
       true, false,
  //             ((point_t) {0,0,0}));
  //	     ((point_t) {0,0,0}), thread_no);
  //             zero_cell_dist, thread_no);
       m_zero_cell_dist);

    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->frozenPairs(),
       cutoff_sq,
       m_first, m_first,
       m_first->frozenParticles(c1), m_first->particles(c2),
       false, true,
  //             ((point_t) {0,0,0}));
  //	     ((point_t) {0,0,0}), thread_no);
  //             zero_cell_dist, thread_no);
       m_zero_cell_dist);
  }

#ifdef _OPENMP   
  virtual void createDistances(int thread_no) final override {
#else
  virtual void createDistances() final override {
    int thread_no = PairCreator::counterTN;
#endif
    ManagerCell *manager = m_first->manager();
    size_t n_colours = manager->nColours();
    /* We are calculating pairs within the same cell */
    /* Loop over the first colour */
    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      /* --- Pairs of the same colour -------------------------------------------------------------- */
      ColourPair *cp = manager->cp(c1, c1);

      if (cp->needPairs()) {
	double cutoff_sq = cp->cutoff() * cp->cutoff();
	CREATE_DISTANCES_FOR_SAME_PARTICLE_LIST//<AddPairCheck_x>
	  (cp->freePairs(),
	   cutoff_sq,
	   m_first,
           m_first->particles(c1));
//	   m_first->particles(c1), thread_no);

        CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
	  (cp->frozenPairs(),
	   cutoff_sq,
	   m_first, m_first,
	   m_first->particles(c1), m_first->frozenParticles(c1),
	   true, false,
//	   ((point_t) {0,0,0}));
//	   ((point_t) {0,0,0}), thread_no);
//	   zero_cell_dist, thread_no);
	   m_zero_cell_dist);
      } 

      for (size_t c2 = c1+1; c2 < n_colours; ++c2) {
	/* --- Pairs of different colour -------------------------------------------------------------- */
	cp = manager->cp(c1, c2);

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();
          CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
	    (cp->freePairs(),
	     cutoff_sq,
	     m_first, m_first,
	     m_first->particles(c1), m_first->particles(c2),
	     true, true,
//             ((point_t) {0,0,0}));
//	     ((point_t) {0,0,0}), thread_no);
//             zero_cell_dist, thread_no);
             m_zero_cell_dist);

          CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
	    (cp->frozenPairs(),
	     cutoff_sq,
	     m_first, m_first,
	     m_first->particles(c1), m_first->frozenParticles(c2),
	     true, false,
//             ((point_t) {0,0,0}));
//	     ((point_t) {0,0,0}), thread_no);
//             zero_cell_dist, thread_no);
             m_zero_cell_dist);

          CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
	    (cp->frozenPairs(),
	     cutoff_sq,
	     m_first, m_first,
	     m_first->frozenParticles(c1), m_first->particles(c2),
	     false, true,
//             ((point_t) {0,0,0}));
//	     ((point_t) {0,0,0}), thread_no);
//             zero_cell_dist, thread_no);
             m_zero_cell_dist);
	}
      } /* Loop over c2 */
    } /* Loop over c1 */
  }

};

/*!
 * This class connects two cells for cell subdivision. It handles
 * the creation of the pair list.
 */
template <typename AddPairCheck_x>
class CellLink: public AbstractCellLink
{ 
  private:
    point_t m_cell_dist, m_minus_cell_dist; 
  public:
    CellLink(Cell *first, Cell *second, int alignment, \
      bool acts_on_first = true, bool acts_on_second = true, bool cross_regions = false): \
        AbstractCellLink(first, second, alignment, acts_on_first, acts_on_second, cross_regions){m_cell_dist = first->m_cell_dist[alignment][cross_regions]; m_minus_cell_dist = -m_cell_dist;};
  /*!
   * Create the pair list.
   * @param t Thread number of this call to createDistances. For multithreaded compilation.
   */
#ifdef _OPENMP   
  virtual void doForSameColors(ColourPair* cp, int colour, double cutoff_sq, int thread_no) final override {
#else
  virtual void doForSameColors(ColourPair* cp, int colour, double cutoff_sq) final override {
    int thread_no = PairCreator::counterTN;
#endif
    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->freePairs(),
       cutoff_sq,
       m_second, m_first,
       m_second->particles(colour), m_first->particles(colour),
       m_acts_on.second, m_acts_on.first,
//               m_cell_dist, thread_no);
       m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);
    if (m_acts_on.first)
      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
        (cp->frozenPairs(),
         cutoff_sq,
         m_second, m_first,
         m_second->frozenParticles(colour), m_first->particles(colour),
         false, true,
//                 m_cell_dist, thread_no);
         m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);
    if (m_acts_on.second)
      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
        (cp->frozenPairs(),
         cutoff_sq,
         m_second, m_first,
         m_second->particles(colour), m_first->frozenParticles(colour),
         true, false,
//                 m_cell_dist, thread_no);
         m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);
  }

#ifdef _OPENMP   
  virtual void doForDiffColors(ColourPair* cp, int c1, int c2, double cutoff_sq, int thread_no) final override {
#else
  virtual void doForDiffColors(ColourPair* cp, int c1, int c2, double cutoff_sq) final override {
    int thread_no = PairCreator::counterTN;
#endif
    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->freePairs(),
       cutoff_sq,
       m_first, m_second,
       m_first->particles(c1), m_second->particles(c2),
       m_acts_on.first, m_acts_on.second,
//               minus_m_cell_dist, thread_no);
//	       -m_cell_dist, thread_no);
       m_minus_cell_dist);
//             -m_first->m_cell_dist[m_alignment][m_cross_regions]);
    CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
      (cp->freePairs(),
       cutoff_sq,
       m_second, m_first,
       m_second->particles(c1), m_first->particles(c2),
       m_acts_on.second, m_acts_on.first,
//               m_cell_dist, thread_no);
       m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);

    if (m_acts_on.first){
      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
        (cp->frozenPairs(),
         cutoff_sq,
         m_first, m_second,
         m_first->particles(c1), m_second->frozenParticles(c2),
         true, false,
         m_minus_cell_dist);
//                 m_minus_cell_dist, thread_no);
//                -m_cell_dist, thread_no);
//                -m_cell_dist);
//                -m_first->m_cell_dist[m_alignment][m_cross_regions]);
      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
        (cp->frozenPairs(),
         cutoff_sq,
         m_second, m_first,
         m_second->frozenParticles(c1), m_first->particles(c2),
         false, true,
//                 m_cell_dist, thread_no);
         m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);
    }

    if (m_acts_on.second){
      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
        (cp->frozenPairs(),
         cutoff_sq,
         m_first, m_second,
         m_first->frozenParticles(c1), m_second->particles(c2),
         false, true,
         m_minus_cell_dist);
//                 m_minus_cell_dist, thread_no);
//                 -m_cell_dist, thread_no);
//                 -m_cell_dist);
//               -m_first->m_cell_dist[m_alignment][m_cross_regions]);
      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
        (cp->frozenPairs(),
         cutoff_sq,
         m_second, m_first,
         m_second->particles(c1), m_first->frozenParticles(c2),
         true, false,
//                 m_cell_dist, thread_no);
         m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);
    }
  }

#ifdef _OPENMP   
  virtual void createDistances(int thread_no) final override {
#else
  virtual void createDistances() final override {
  int thread_no = PairCreator::counterTN;
#endif
    ManagerCell *manager = m_first->manager();
    size_t n_colours = manager->nColours();
    // We have different cells.
    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      for (size_t c2 = 0; c2 < n_colours; ++c2) {
	ColourPair *cp = manager->cp(c1, c2);

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();

	  if (c1 < c2) {
	      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
	      (cp->freePairs(),
	       cutoff_sq,
	       m_first, m_second,
	       m_first->particles(c1), m_second->particles(c2),
	       m_acts_on.first, m_acts_on.second,
//               minus_m_cell_dist, thread_no);
//	       -m_cell_dist, thread_no);
	       m_minus_cell_dist);
//             -m_first->m_cell_dist[m_alignment][m_cross_regions]);

	    if (m_acts_on.first)
	      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
		(cp->frozenPairs(),
		 cutoff_sq,
		 m_first, m_second,
		 m_first->particles(c1), m_second->frozenParticles(c2),
		 true, false,
                 m_minus_cell_dist);
//                 m_minus_cell_dist, thread_no);
//                -m_cell_dist, thread_no);
//                -m_cell_dist);
//                -m_first->m_cell_dist[m_alignment][m_cross_regions]);

	    if (m_acts_on.second)
	      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
		(cp->frozenPairs(),
		 cutoff_sq,
		 m_first, m_second,
		 m_first->frozenParticles(c1), m_second->particles(c2),
		 false, true,
                 m_minus_cell_dist);
//                 m_minus_cell_dist, thread_no);
//                 -m_cell_dist, thread_no);
//                 -m_cell_dist);
//               -m_first->m_cell_dist[m_alignment][m_cross_regions]);
          } else {
	      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
	      (cp->freePairs(),
	       cutoff_sq,
	       m_second, m_first,
	       m_second->particles(c2), m_first->particles(c1),
	       m_acts_on.second, m_acts_on.first,
//               m_cell_dist, thread_no);
               m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);

	    if (m_acts_on.first)
	      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
		(cp->frozenPairs(),
		 cutoff_sq,
		 m_second, m_first,
		 m_second->frozenParticles(c2), m_first->particles(c1),
		 false, true,
//                 m_cell_dist, thread_no);
                 m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);

	    if (m_acts_on.second)
	      CREATE_DISTANCES_FOR_DIFFERENT_PARTICLE_LISTS//<AddPairCheck_x>
		(cp->frozenPairs(),
		 cutoff_sq,
		 m_second, m_first,
		 m_second->particles(c2), m_first->frozenParticles(c1),
		 true, false,
//                 m_cell_dist, thread_no);
                 m_cell_dist);
//               m_first->m_cell_dist[m_alignment][m_cross_regions]);
          }
	}
      }
    }
  }

};

template <typename AddPairCheck_x>
AbstractCellLink *establishLinkTemplate(ManagerCell *a_manager, Cell *first, Cell *second, int where, bool acts_on_first, bool acts_on_second, bool cross_regions)
{

    auto *link = new CellLink<AddPairCheck_x>(first, second, where, acts_on_first, acts_on_second, cross_regions);

    first->neighbors(where).push_back(link);
    second->neighbors(a_manager->num_neighbors()-where-1).push_back(link);
    return link;
}

struct AddPairCheck_Regular {
  static void check(Cell *c, int i, double &cartesian, double width){};
};

struct AddPairCheck_OneCellDims {
  static void check(Cell *c, int i, double &cartesian, double width)
  {
//    int signC = sign(cartesian);
//    if (dynamic_cast<BoundaryCell*> (c)->m_region->m_oneCellPeriodicDims[i] && signC*cartesian > width/2){
//    if (((BoundaryCell*) c)->m_region->m_oneCellPeriodicDims[i] && signC*cartesian > width/2){
//      cartesian -= signC * width;
    if (((BoundaryCell*) c)->m_region->m_oneCellPeriodicDims[i])
      if(cartesian > 0.5*width) cartesian -= width; 
      else if(cartesian < -0.5*width) cartesian += width; 
//      return true;
    //else return false;
  };
};

template<typename AddPairCheck_x>
void addPair(vector<PairList> &distances, double cutoff_sq,
  Cell *first_c, Cell *second_c,
  Particle *first_p, Particle *second_p,
  bool ao_f, bool ao_s,
//  point_t &cell_dist,
  point_t &precomputed,
  int thread_no)
{
  dist_t d;
  d.abs_square = 0;
  for (int _i = 0; _i < SPACE_DIMS; _i++) {
//      d.cartesian[_i] = cell_dist[_i]
//        + first_p -> r[_i] - first_c -> corner1[_i]
      d.cartesian[_i] = precomputed[_i]
        - second_p -> r[_i] + second_c -> corner1[_i];
      AddPairCheck_x::check(first_c, _i, d.cartesian[_i], first_c->corner2[_i] - first_c->corner1[_i]);
    d.abs_square += d.cartesian[_i]*d.cartesian[_i];
  }
  /* Take care: The order of *j, *i defines the direction d \
     is pointing to. */
  if (d.abs_square < cutoff_sq) {
    d.abs = sqrt(d.abs_square);
      distances[thread_no].newPair().set(d, first_p, second_p, ao_f, ao_s);
  }
};




#endif
