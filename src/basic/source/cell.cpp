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



using namespace std;

#include <algorithm>

#include "cell.h"
#include "phase.h"
#include "pairdist.h"
#include "pair_list.h"
#include "controller.h"
#include "manager_cell.h"
#include "wall_triangle.h"
#include "colour_pair.h"
#include <iomanip>
// #include "pair_creator.h"


// const int_point_t Cell::c_offsets[NUM_NEIGHBORS] = {
//   {{-1},{-1},{-1}}, {{-1},{-1}, {0}}, {{-1},{-1}, {1}}, {{-1}, {0},{-1}}, {{-1}, {0}, {0}}, {{-1}, {0}, {1}},
//   {{-1}, {1},{-1}}, {{-1}, {1}, {0}}, {{-1}, {1}, {1}},
//   { {0},{-1},{-1}}, { {0},{-1}, {0}}, { {0},{-1}, {1}}, { {0}, {0},{-1}}, /*{{0}, {0}, {0}},*/ { {0}, {0}, {1}},
//   { {0}, {1},{-1}}, { {0}, {1}, {0}}, { {0}, {1}, {1}},
//   { {1},{-1},{-1}}, { {1},{-1}, {0}}, { {1},{-1}, {1}}, { {1}, {0},{-1}}, { {1}, {0}, {0}}, { {1}, {0}, {1}},
//   { {1}, {1},{-1}}, { {1}, {1}, {0}}, { {1}, {1}, {1}}
// };


// const int_point_t Cell::c_half_offsets[NUM_HALF_NEIGHBORS] = {
//     {-1,-1, 1}, {-1, 0, 1},
//     {-1, 1, 0}, {-1, 1, 1},
//     {0,-1, 1}, {0, 0, 1},
//     {0, 1, 0}, {0, 1, 1},
//     {1,-1, 1}, {1, 0, 0}, {1, 0, 1},
//     {1, 1, 0}, {1, 1, 1}
// };


/*---- Class AbstractCellLink ----*/

AbstractCellLink::AbstractCellLink()
{
  throw gError("AbstractCellLink::AbstractCellLink(default)", "Should not be called. Contact the programmer.");

  set(NULL, NULL, -1, true, true);
}


AbstractCellLink::AbstractCellLink(Cell *first, Cell *second, bool acts_on_first, bool acts_on_second, bool cross_regions)
  : m_n_active_cells(0), m_cross_regions(cross_regions), next(NULL), prev(NULL)/*, m_linkUsed(false)*/
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_init(&m_activation__mutex, &g_mutex_attr);
#endif
  set(first, second, acts_on_first, acts_on_second);
}

AbstractCellLink::AbstractCellLink(const AbstractCellLink &copy)
  : m_n_active_cells(copy.m_n_active_cells), next(NULL), prev(NULL)/*, m_linkUsed(false)*/
{
  throw gError("AbstractCellLink::AbstractCellLink(copy)", "Should not be called. Contact the programmer.");

  set(copy.m_first, copy.m_second, copy.m_acts_on.first, copy.m_acts_on.second);
}


AbstractCellLink::~AbstractCellLink()
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_destroy(&m_activation__mutex);
#endif
}



//---- Methods ----

/* cellDist determines the distance between two cells in real coordinates.
   i.e.
   if first is left of second it takes
            second.corner1 - first.corner2
   and otherwise
            first.corner1 - second.corner2
*/

void AbstractCellLink::set(Cell *first, Cell *second, int alignment, bool acts_on_first, bool acts_on_second)
{
  m_first = first;
  m_second = second;
  m_acts_on.first = acts_on_first;
  m_acts_on.second = acts_on_second;

  assert(first && second && ((alignment != -1) || (first == second)) );
}


/* One of the two cells connected by a cell link is activated.
   If both links are activated, activate this link
*/
void AbstractCellLink::cellActivated() {
#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_activation__mutex);
#endif

  ++m_n_active_cells;
//    MSG_DEBUG ("AbstractCellLink::cellActivated","start");

//    assert(m_n_active_cells >= 0 && m_n_active_cells <= 2);
  if (!(m_n_active_cells >= 0 && m_n_active_cells <= 2)) {
    MSG_DEBUG
      ("AbstractCellLink::cellActivated",
       "m_n_active_cells = " << m_n_active_cells);
    abort();
  }

// Assigning the activated AbstractCellLink in the proper thread list
  if (m_n_active_cells == 2) 
    m_first->manager()->activateCellLink(this);


#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_activation__mutex);
#endif
}


/* One of the two cells connected by a cell link is deactivated.
   If link is active, deactivate link.
*/
void AbstractCellLink::cellDeactivated() {
#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_activation__mutex);
#endif

  --m_n_active_cells;

//    assert(m_n_active_cells >= 0 && m_n_active_cells <= 2);
  if (!(m_n_active_cells >= 0 && m_n_active_cells <= 2)) {
    MSG_DEBUG
      ("AbstractCellLink::cellDeactivated",
       "m_n_active_cells = " << m_n_active_cells << " ! Aborting.");
    abort();
  }

  if (m_n_active_cells == 1)
    m_first->manager()->deactivateCellLink(this);

#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_activation__mutex);
#endif
}


/* ---- Inline helper functions ---- */


/* add a Pairdist to the appropriate distances
   list
*/

// #ifdef _OPENMP
// inline void addPair(vector<PairList> &distances, double cutoff_sq,
//                     int dir, Cell *first_c, Cell *second_c,
//                     Particle *first_p, Particle *second_p,
//                     bool ao_f, bool ao_s,
//                     point_t cell_dist,
//                     int thread_no)
// {
// #else
// inline void addPair(vector<PairList> &distances, double cutoff_sq,
//                     int dir, Cell *first_c, Cell *second_c,
//                     Particle *first_p, Particle *second_p,
//                     bool ao_f, bool ao_s,
//                     point_t cell_dist)
// {
// #endif
//   dist_t d;
//   d.abs_square = 0;
//   for (int _i = 0; _i < SPACE_DIMS; _i++) {
//     d.cartesian[_i] = -dir*cell_dist[_i]
//       + first_p -> r[_i] - first_c -> corner1[_i]
//       - second_p -> r[_i] + second_c -> corner1[_i];
//     d.abs_square += d.cartesian[_i]*d.cartesian[_i];
//     }
//   /* Take care: The order of *j, *i defines the direction d \
//      is pointing to. */
//   if (d.abs_square < cutoff_sq) {
//     d.abs = sqrt(d.abs_square);
// #ifdef _OPENMP    
//     distances[thread_no].newPair().set(d, first_p, second_p, ao_f, ao_s);
// #else    
//     distances[PairCreator::counterTN].newPair().set(d, first_p, second_p, ao_f, ao_s);
// #endif

//   }
// }




/*
void CellLink::addPair_(vector<PairList> &distances, double cutoff_sq,
                    Cell *first_c, Cell *second_c,
                    Particle *first_p, Particle *second_p,
                    bool ao_f, bool ao_s,
                    point_t &cell_dist,
                    int thread_no)
{
  dist_t d;
  d.abs_square = 0;
  for (int _i = 0; _i < SPACE_DIMS; _i++) {
      d.cartesian[_i] = cell_dist[_i]
        + first_p -> r[_i] - first_c -> corner1[_i]
        - second_p -> r[_i] + second_c -> corner1[_i];
      d.abs_square += d.cartesian[_i]*d.cartesian[_i];
    }*/
  /* Take care: The order of *j, *i defines the direction d \
     is pointing to. */
//    else if 
//    if (first_p->mySlot == 46 || second_p->mySlot == 46)// || first_p->mySlot == 554 || second_p->mySlot == 554)
//    MSG_DEBUG("addPair", first_p->mySlot << " " << first_p->r << " " << first_c->tag << " " << second_p->mySlot << " " << second_p->r << " " << second_c->tag << " " << cell_dist << " " << d.cartesian << " " << first_c->corner1 << " " << second_c->corner1 << " " << d.abs_square << " " << cutoff_sq);
/*  if (d.abs_square < cutoff_sq) {
//    if (first_p->mySlot == 46 || second_p->mySlot == 46)// || first_p->mySlot == 554 || second_p->mySlot == 554)
//            MSG_DEBUG("addPair", first_p->mySlot << " " << first_p->r << " " << first_c->tag << " " << second_p->mySlot << " " << second_p->r << " " << second_c->tag << " " << cell_dist << " " << d.cartesian << " " << first_c->corner1 << " " << second_c->corner1 << " " << d.abs_square << " " << cutoff_sq);
    int _i = 2; 
//      d.cartesian[_i] = cell_dist[_i] + first_p->r[_i] - first_c->corner1[_i] - second_p->r[_i] + second_c->corner1[_i];
//    if (first_p->mySlot == 567){
//      MSG_DEBUG("addPair", cell_dist[_i] << " + " << first_p->r[_i] << " = " << cell_dist[_i] + first_p->r[_i]);
//      MSG_DEBUG("addPair", cell_dist[_i] + first_p->r[_i] << " - " << first_c->corner1[_i] << " - " << second_p->r[_i] << " = " << cell_dist[_i] + first_p->r[_i] -first_c->corner1[_i] - second_p->r[_i]);
//      MSG_DEBUG("addPair", cell_dist[_i] + first_p->r[_i] -first_c->corner1[_i] - second_p->r[_i] <<  " + " << second_c->corner1[_i] << " = " <<  cell_dist[_i] + first_p->r[_i] - first_c->corner1[_i] - second_p->r[_i] + second_c->corner1[_i]);
//      MSG_DEBUG("addPair",  d.cartesian[2] << " = " <<  cell_dist[_i] + first_p->r[_i] - first_c->corner1[_i] - second_p->r[_i] + second_c->corner1[_i]);
//}
    //if (second_p->mySlot == 564)
    //        MSG_DEBUG("addPair(", &distances<<",...) " << second_p->mySlot << " " << second_p->r  << " " << first_p->mySlot << " " << first_p->r << " " << -cell_dist << " " << -d.cartesian << " " << second_c->corner1 << " " << first_c->corner1 << " " << d.abs_square << " " << cutoff_sq);
    //else if (first_p->mySlot == 564)
    //        MSG_DEBUG("addPair(", &distances<<",...) " << first_p->mySlot << " " << first_p->r  << " " << second_p->mySlot << " " << second_p->r << " " << cell_dist << " " << d.cartesian << " " << first_c->corner1 << " " << second_c->corner1 << " " << d.abs_square << " " << cutoff_sq);
    d.abs = sqrt(d.abs_square);
//#ifdef _OPENMP    
    distances[thread_no].newPair().set(d, first_p, second_p, ao_f, ao_s);*/
/*#else    
    distances[PairCreator::counterTN].newPair().set(d, first_p, second_p, ao_f, ao_s);
#endif*/
//  }
//}


/*---- Class Cell ----*/

//---- Constructors/Destructor ----

Cell::Cell(ManagerCell *mgr, int group, int sameDirOutlets)
  : m_particles(mgr->nColours()), m_frozen_particles(mgr->nColours()),
    m_injected_particles(mgr->nColours()), m_n_particles(0), m_manager(mgr), m_group(group), m_cellUsed(false),
#ifdef ENABLE_PTHREADS
    m_injected_particles__mutex(mgr->nColours()),
#endif
    next(NULL), prev(NULL)
{
  commonConstructor(sameDirOutlets);
}


Cell::Cell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int group, int sameDirOutlets)
  : cuboid_t(c1, c2), m_particles(mgr->nColours()), m_frozen_particles(mgr->nColours()),
    m_injected_particles(mgr->nColours()),
    m_n_particles(0),
    m_manager(mgr), m_group(group), m_cellUsed(false),
#ifdef ENABLE_PTHREADS
    m_injected_particles__mutex(mgr->nColours()),
#endif
    next(NULL), prev(NULL)
{
  commonConstructor(sameDirOutlets);
}

Cell::Cell(ManagerCell *mgr, cuboid_t cuboid, const int_point_t &a_tag, int group, int sameDirOutlets)
  : cuboid_t(cuboid), tag(a_tag), m_particles(mgr->nColours()), m_frozen_particles(mgr->nColours()),
    m_injected_particles(mgr->nColours()),
    m_n_particles(0),
    m_manager(mgr), m_group(group), m_cellUsed(false),
#ifdef ENABLE_PTHREADS
    m_injected_particles__mutex(mgr->nColours()),
#endif
    next(NULL), prev(NULL)
{
  commonConstructor(sameDirOutlets, a_tag);
}


Cell::Cell(ManagerCell *mgr, const point_t &c1, const point_t &c2, const int_point_t &a_tag, int group, int sameDirOutlets)
  : cuboid_t(c1, c2), tag(a_tag), m_particles(mgr->nColours()), m_frozen_particles(mgr->nColours()),
    m_injected_particles(mgr->nColours()),
    m_n_particles(0),
    m_manager(mgr), m_group(group), m_cellUsed(false),
#ifdef ENABLE_PTHREADS
    m_injected_particles__mutex(mgr->nColours()),
#endif
    next(NULL), prev(NULL)
{
  commonConstructor(sameDirOutlets, a_tag);
}

BoundaryCell::BoundaryCell(ManagerCell *mgr, region_t *r, int group) : Cell(mgr, group, 2) {m_region = r;}

BoundaryCell::BoundaryCell(ManagerCell *mgr, cuboid_t cuboid, int_point_t &a_tag, region_t *r, int group) : Cell(mgr, cuboid, a_tag, group, 2) {m_region = r;}

BoundaryCell::BoundaryCell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int_point_t &a_tag, region_t *r, int group): Cell(mgr, c1, c2, a_tag, group, 2) {m_region = r;}


void Cell::commonConstructor(int sameDirOutlets, const int_point_t &a_tag)//TODO: change sameDirOutlets to sameDir
{
  const int num_neighbors = m_manager->num_neighbors();
  m_neighbors.reserve(num_neighbors);
  m_neighbors.resize(num_neighbors);

  m_cell_dist.reserve(num_neighbors);
  m_cell_dist.resize(num_neighbors);
  for (int i=0; i < num_neighbors; i++)
  {
    m_cell_dist[i].reserve(sameDirOutlets);
    m_cell_dist[i].resize(sameDirOutlets);
  }
  m_outlets.reserve(NUM_DIRECT_NEIGHBORS);
  m_outlets.resize(NUM_DIRECT_NEIGHBORS);
  for (int i=0; i < NUM_DIRECT_NEIGHBORS; i++)
  {
    m_outlets[i].reserve(sameDirOutlets);
    m_outlets[i].resize(sameDirOutlets);
  }
  tag = a_tag;
#ifdef ENABLE_PTHREADS
  for (size_t i = 0; i < m_injected_particles__mutex.size(); i++)
    pthread_mutex_init(&m_injected_particles__mutex[i], &g_mutex_attr);
#endif
}
  

Cell::~Cell()
{
#ifdef ENABLE_PTHREADS
  for (int i = 0; i < m_injected_particles__mutex.size(); i++)
    pthread_mutex_destroy(&m_injected_particles__mutex[i]);
#endif
}

int Cell::whichNeighbor(const point_t &r)
{
  int_point_t off;

  for (int j = 0; j < SPACE_DIMS; j++)
    if (r[j] < corner1[j])
      off[j] = -1;
    else if (r[j] >= corner2[j])
      off[j] = 1;
  return m_manager->offset2neighbor(off);
}

//--- Methods ---




void Cell::addPeriodic(Cell *neighbor, int where, bool cross_regions)
{
    if (m_manager->m_divby > 1) where = m_manager->c_2x_1x(where);
    addToMOutlets(this, neighbor, where, cross_regions);
}

//#define TRACK_PARTICLE 825 
#define TRACK_PARTICLE_COLOUR 0

void Cell::doCollision(Particle *p, point_t& r, point_t& v, const point_t &force, IntegratorPosition *integratorP)
{
	double t_travelled;
	point_t hit_pos, old_r;
	Wall *wall;
	bool hit = true;

#ifdef TRACK_PARTICLE

  if (p->mySlot == TRACK_PARTICLE)
    cout << "Cell::doCollision: tracking particle " << TRACK_PARTICLE << " !!! " << m_all_walls.size() << " walls in cell." << endl;

#endif
  double iterations = 0;
  while (hit) {
    ++iterations;
    if(iterations > 100)
      throw gError("Cell::doCollision", "More than 100 wall collisions for particle " + ObjToString(p->mySlot) + ", " + ObjToString(p->r) + ", colour = " + ObjToString(p->c) + "!!!\nCheck your settings (forces, size of timestep, ...)!");
    hit = false;
    t_travelled = HUGE_VAL;
    old_r = /*p->*/r;

    /* First, find hit that occurs first. */
    // this cell knows about all walls, intersecting this and its neighbouring cells
    checkForHit(p, force, hit, t_travelled, hit_pos, wall, integratorP);

#ifdef TRACK_PARTICLE

    if (p->mySlot == TRACK_PARTICLE && p->c == TRACK_PARTICLE_COLOUR && hit) {
      cout << "--> hit @ " << hit_pos << " after time " << t_travelled << " on wall:" << endl;
      cout << wall->toString() << endl;
    }

#endif
  hit = hit || doCollisionIndirectOutlets(p, force, r, t_travelled, hit_pos, wall, integratorP);

    if (hit) {

//       MSG_DEBUG("Cell::doCollision", "hit = true for " << p->r << ", " << p->mySlot << ", " << p->c << "wall = " << wall->toString());
#ifdef TRACK_PARTICLE

      if (p->mySlot == TRACK_PARTICLE && p->c == TRACK_PARTICLE_COLOUR) {
        cout << "--> reflecting" << endl;
        cout << "before: p->r = " << p->r << ", p->v = " << /*p->*/v << ", p->dt = " << p->dt << endl;
      }

#endif

// MSG_DEBUG("Cell::doCollision", "HITCASE: wall normal = " << wall->normal());
      /* Hit happened inside. Reflect particle. */
      wall->reflector()->reflect(p, /*p->*/r, /*p->*/v, hit_pos, wall->normal(), wall->inPlane());
      p->dt -= t_travelled;
      // next is for the case that the reflector has aborted the collisions
      // by setting p->dt = 0
      if(p->dt < 0) p->dt = 0;
//       assert(p->dt > 0);

#ifdef TRACK_PARTICLE

      if (p->mySlot == TRACK_PARTICLE && p->c == TRACK_PARTICLE_COLOUR) {
        cout << "after: p->r = " << p->r << ", p->v = " << /*p->*/v << ", p->dt = " << p->dt << endl;
      }

#endif

    } else
      /*p->*/r = old_r;

	}
}

bool Cell::emitIntoOutlets(Particle *p, size_t colour, int_point_t &off, int n, IntegratorPosition *integrator, Phase *phase)
{
  if (Cell *c = m_outlets[n][0]) {
    if (c->isInside(p->r)) {
      c->injectFree(colour, p);

      /* Notify any subclass that a particle has left the cell in a certain direction. */
      particleLeftCell(off, n);
#ifdef TRACK_PARTICLE
      if (p->mySlot == TRACK_PARTICLE){
        int force_index = ((Controller*) integrator->parent())->forceIndex();
        cout << "particle id (slot)               = " << p->mySlot << endl;
        cout << "current position                 = " << p->r << endl;
        cout << "current velocity                 = " << p->v << endl;
        cout << "offset                           = " << off << endl;
        cout << "emission direction               = " << n << endl;
        cout << "force                            = " << p->force[force_index] << endl;
        cout << "=== tag ===" << endl;
        cout << p->tag.toString();
        MSG_DEBUG("Cell::emitIntoOutlets", "Particle entered neighboring cell:");
        cout << "=== DEBUGING INFORMATION: CELLS ===" << endl;
        cout << "current cell: tag                        = " << tag << endl;
        cout << "current cell: corner 1                   = " << corner1 << endl;
        cout << "current cell: corner 2                   = " << corner2 << endl;
        cout << "entering cell: tag                       = " << c->tag << endl;
        cout << "entering cell: corner 1                  = " << c->corner1 << endl;
        cout << "entering cell: corner 2                  = " << c->corner2 << endl;
        cout << "offset between cells (current->entering) = " << off << endl;
      }
#endif
    } else {
        int force_index = ((Controller*) integrator->parent())->forceIndex();
        cout << "particle id (slot)               = " << p->mySlot << endl;
        cout << "current position                 = " << p->r << endl;
        cout << "current velocity                 = " << p->v << endl;
        cout << "offset                           = " << off << endl;
        cout << "emission direction               = " << n << endl;
        cout << "force                            = " << p->force[force_index] << endl;
        cout << "=== tag ===" << endl;
        cout << p->tag.toString();
        MSG_DEBUG("Cell::emitIntoOutlets", "Particle did not enter neighboring cell:");
        cout << "=== DEBUGING INFORMATION: CELLS ===" << endl;
        cout << "current cell: tag                        = " << tag << endl;
        cout << "current cell: corner 1                   = " << corner1 << endl;
        cout << "current cell: corner 2                   = " << corner2 << endl;
        cout << "entering cell: tag                       = " << c->tag << endl;
        cout << "entering cell: corner 1                  = " << c->corner1 << endl;
        cout << "entering cell: corner 2                  = " << c->corner2 << endl;
        cout << "offset between cells (current->entering) = " << off << endl;
     throw gError("Cell::emitIntoOutlets", "Particle flew farther than "
                   "a cell in one timestep. Please check the correctness of your forces, decrease dt or "
                   "increase the cut-off radius.", gError::PARTICLEFLEWTOOFAR);
    }
    return false;
  }else
    return true;
}

bool BoundaryCell::emitIntoOutlets(Particle *p, size_t colour, int_point_t &off, int n, IntegratorPosition *integrator, Phase *phase)
{
  point_t old_r = p->r;
  Particle *new_p = p;
  int i=0;
  bool total_erase = true;
  for (vector<Cell*>::iterator c = m_outlets[n].begin(); c != m_outlets[n].end(); c++, i++) if (*c){ total_erase = false;

    if (!new_p) new_p = phase->addParticle(*p);

    /* Update position in case this is periodic, etc. */
//old style    new_p->r = old_r + (*c)->corner1 - corner1 - m_cell_dist[((m_manager->m_divby == 2)?c_2x_direct_neighbors[n]:n)][i];
    new_p->r = old_r + (*c)->corner1 - corner1 - m_cell_dist[m_manager->m_direct_neighbors[n]][i];
    
    if ((*c)->isInside(new_p->r)) {
      (*c)->injectFree(colour, new_p);

      /* Notify any subclass that a particle has left the cell in a certain direction. */
      particleLeftCell(off, n);
#ifdef TRACK_PARTICLE
      if (p->mySlot == TRACK_PARTICLE && p->c == TRACK_PARTICLE_COLOUR) {
        int force_index = ((Controller*) integrator->parent())->forceIndex();
        MSG_DEBUG("BoundaryCell::emitIntoOutlets", "Particle entered neighboring cell:");
        cout << "particle id (slot)               = " << p->mySlot << endl;
        cout << "position before update/collision = " << old_r << endl;
        cout << "current position                 = " << new_p->r << endl;
        cout << "current velocity                 = " << p->v << endl;
        cout << "force                            = " << p->force[force_index] << endl;
        cout << "=== tag ===" << endl;
        cout << p->tag.toString();

        cout << "=== DEBUGING INFORMATION: CELLS ===" << endl;
        cout << "current cell: tag                        = " << tag << endl;
        cout << "current cell: corner 1                   = " << corner1 << endl;
        cout << "current cell: corner 2                   = " << corner2 << endl;
        cout << "entering cell: tag                       = " << (*c)->tag << endl;
        cout << "entering cell: corner 1                  = " << (*c)->corner1 << endl;
        cout << "entering cell: corner 2                  = " << (*c)->corner2 << endl;
        cout << "offset between cells (current->entering) = " << off << " " << n << " " << m_cell_dist[m_manager->m_direct_neighbors[n]][i] << endl;
      }
#endif

    } else {
      MSG_DEBUG("BoundaryCell::emitIntoOutlets", "Particle did not enter neighboring cell:");
      int force_index = ((Controller*) integrator->parent())->forceIndex();

      cout << "=== DEBUGING INFORMATION: PARTICLE ===" << endl;
      cout << "particle id (slot)               = " << p->mySlot << endl;
      cout << "position before update/collision = " << old_r << endl;
      cout << "current position                 = " << new_p->r << endl;
      cout << "current velocity                 = " << p->v << endl;
      cout << "force                            = " << p->force[force_index] << endl;
      cout << "region n_cells                   = " << m_region->n_cells << endl;
      cout << "=== tag ===" << endl;
      cout << p->tag.toString();

      cout << "=== DEBUGING INFORMATION: CELLS ===" << endl;
      cout << "current cell: tag                        = " << tag << endl;
      cout << "current cell: corner 1                   = " << corner1 << endl;
      cout << "current cell: corner 2                   = " << corner2 << endl;
      cout << "entering cell: tag                       = " << (*c)->tag << endl;
      cout << "entering cell: corner 1                  = " << (*c)->corner1 << endl;
      cout << "entering cell: corner 2                  = " << (*c)->corner2 << endl;
      cout << "offset between cells (current->entering) = " << off << endl;
      cout << "neigbor index                            = " << n << endl;
      cout << "region n_cells                   = " << m_region->n_cells << endl;

      throw gError("BoundaryCell::emitIntoOutlets", "Particle flew farther than "
                   "a cell in one timestep. Please check the correctness of your forces, decrease dt or "
                   "increase the cut-off radius.", gError::PARTICLEFLEWTOOFAR);
    }
      new_p = NULL;
  }
  return total_erase;
}

/* Main integration and collision detection logic. */
void Cell::updatePositions(IntegratorPosition *integrator)
{
  Phase *phase = m_manager->phase();

  double dt, dt_div_mass, dt_div2_mass;
  size_t colour;
  //	forces_t *forces;

  dt = integrator->dt();
//   lambda = integrator->lambda();
  dt_div_mass = integrator->dtDivMass();
  dt_div2_mass = integrator->dtDiv2Mass();
  colour = integrator->colour();
  //	forces = integrator->forces();

  list<Particle*>::iterator m_particles_end = m_particles[colour].end();
  for (list<Particle*>::iterator i = m_particles[colour].begin(); i != m_particles_end; ) {

//     MSG_DEBUG("Cell::updatePositions", "now at p =  " << (*i)->mySlot << ", c = " << (*i)->c);
/*    if((*i)->mySlot == 52)
    {
      if(phase->boundary()->isInside((*i)->r))
        MSG_DEBUG("Cell::updatePositionsSTART", "52 INSIDE");
      else MSG_DEBUG("Cell::updatePositionsSTART", "52 NOT INSIDE");
    }*/

    bool erase = false;
    bool total_erase = false;
    Particle *p = *i;
    point_t startr = p->r;

    p->dt = dt;

    integrator->integratePosition(p, this);
    integrator->integrateVelocity(p);


#ifdef TRACK_PARTICLE
      if (p->mySlot == TRACK_PARTICLE && p->c == TRACK_PARTICLE_COLOUR) {
        int force_index = ((Controller*) integrator->parent())->forceIndex();
        MSG_DEBUG("Cell::updatePositions", "");
        cout << "particle id (slot)               = " << p->mySlot << endl;
        cout << "position before update/collision = " << startr << endl;
        cout << "current position                 = " << p->r << endl;
        cout << "current velocity                 = " << p->v << endl;
        cout << "force                            = " << p->force[force_index] << endl;
        cout << "=== tag ===" << endl;
        cout << p->tag.toString();

        cout << "=== DEBUGING INFORMATION: CELLS ===" << endl;
        cout << "current cell: tag                        = " << tag << endl;
        cout << "current cell: corner 1                   = " << corner1 << endl;
        cout << "current cell: corner 2                   = " << corner2 << endl;
      }
#endif

    // Last modified: 2007-12-27: changed from isInside to isInsideEps due to problem with geometrical epsilons. There was an inconsistency due to a small epsilon between the isInside of the cell the particle is leaving and isInside of the cell the particle should enter (see below)
    if (!isInsideEps(p->r, g_geom_eps)) {
//        MSG_DEBUG("Cell::updatePositionsSTART", "52 INSIDE");
      int_point_t off = {0, 0, 0};
      int n;

       for (int j = 0; j < SPACE_DIMS; j++)
        if (p->r[j] < corner1[j])
          off[j] = -1;
        else if (p->r[j] >= corner2[j])
          off[j] = 1;

      OFFSET2DIRECTNEIGHBOR(off, n);

      erase = true;
      total_erase = emitIntoOutlets(p, colour, off, n, integrator, phase);
    } 
//     if((*i)->mySlot == 52)
//     {
//       if(phase->boundary()->isInside((*i)->r))
//         MSG_DEBUG("Cell::updatePositionsEND", "52 INSIDE");
//       else MSG_DEBUG("Cell::updatePositionsEND", "52 NOT INSIDE");
//     }


    if (erase) {
/*      if((*i)->mySlot == 52) MSG_DEBUG("Cell::updatePositions", "52: erase = TRUE");*/

      Particle *p = *i;
      list<Particle*>::iterator e = i;
      i++;
      m_particles[colour].erase(e);

      assert(m_n_particles != 0);

      m_n_particles--;

      if (m_n_particles == 0)
        deactivate();

      if (total_erase) {
//                 cout << ">>> DELETING PARTICLE " << index << " >>>" << endl;
//                 cout << "off = " << off << endl;
//                 cout << "corner1 = " << corner1 << endl;
//                 cout << "corner2 = " << corner2 << endl;
//                 cout << "p->r = " << p->r << endl;
//                 cout << "p->v = " << p->v << endl;
//                 cout << "p->dt = " << p->dt << endl;
//                 cout << "force = " << (*forces)[*i] << endl;

        m_manager->phase()->removeParticle(p);
      }
    } else
      i++;
  }
}

/*
void Cell::createDistances(int t)
{
  m_local_link->createDistances(t);

  map<Cell*, AbstractCellLink*>::iterator m_links_end = m_links.end();
  for (map<Cell*, AbstractCellLink*>::iterator i = m_links.begin(); i != m_links_end; i++)
    i->second->createDistances(t);
}
*/


/*
void Cell::addColour(size_t colour) {
  assert(colour == m_particles.size());

  m_particles.resize(colour+1);
  m_frozen_particles.resize(colour+1);
  m_injected_particles__mutex.resize(colour+1);
  pthread_mutex_init(&m_injected_particles__mutex[colour], NULL);
  m_injected_particles.resize(colour+1);
}
*/

void Cell::injectFree(size_t colour, Particle *p)
{
  if (p->g != m_group)  {
    m_manager->phase()->groupChanged(p, m_group);
    p->g = m_group;
  }

#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_injected_particles__mutex[colour]);
#endif
  m_injected_particles[colour].push_back(p);
#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_injected_particles__mutex[colour]);
#endif
}


void Cell::injectFrozen(size_t colour, Particle *p)
{
  p->g = m_group;
  p->isFrozen = 1;
  m_frozen_particles[colour].push_back(p);

  if (!m_n_particles) {
    activate();
  }
  ++m_n_particles;
}


void Cell::commitInjections()
{
  size_t np = 0;

  for(size_t colour = 0; colour < m_particles.size(); ++colour)
  {
    np += m_injected_particles[colour].size();

    FOR_EACH
      (list<Particle*>,
       m_injected_particles[colour],
       m_particles[colour].push_back(*__iFE);
      );

    m_injected_particles[colour].clear();

  }

  if (np && !m_n_particles)
    activate();
  m_n_particles += np;
}


void Cell::clearTags()
{
//   Phase *phase = m_manager->phase();
  for (size_t c = 0; c < m_particles.size(); ++c) {
    FOR_EACH
      (list<Particle*>,
       m_particles[c],
       (*__iFE)->tag.clear();
      );
// if(c == 1) MSG_DEBUG("Cell::clearTags", "clearing 1");

  }
}


void Cell::clearOutlets(int where)
{
  m_outlets[where].clear();
}

#if 0
void Cell::clearNeighbors(int where)
{
  for (list<AbstractCellLink*>::iterator i = m_neighbors[where].begin();
       i != m_neighbors[where].end(); i++) {
    vector<AbstractCellLink*>::iterator help;

    // Fixme!!! Slow as crap.
    help = find(m_manager->m_links.begin(), m_manager->m_links.end(), *i);

    if (help != m_manager->m_links.end()) {
      m_manager->m_links.erase(help);
    }

    //    m_links.erase(i->first);
  }

  m_neighbors[where].clear();
}
#endif

void Cell::eraseParticle(list<Particle*>::iterator e)
{
  Particle *p = *e;

  m_particles[p->c].erase(e);

  assert(m_n_particles != 0);

  m_n_particles--;

  if (m_n_particles == 0)
    deactivate();

  m_manager->phase()->removeParticle(p);
}


void Cell::assignContainer(WallContainer *container)
{
  for (list<Wall*>::iterator i = container->walls().begin(); i != container->walls().end(); i++) {
    if ((*i)->intersects(*this)) {
      addWall(*i);
    }
    /*else {
      WallTriangle *wt = (WallTriangle*) *i;

      cout << "rejected:" << endl;
      cout << "wt: " << wt->corner(0) << " " << wt->corner(1) << " " << wt->corner(2) << endl;
      cout << "c:  " << corner1 << " " << corner2 << endl;
      }*/
  }
}


void Cell::setupWalls()
{
  copy(m_walls.begin(), m_walls.end(), inserter(m_all_walls, m_all_walls.begin()));

  // size_t counter = ;
  // for (list<Wall*>::iterator ii = m_all_walls.begin(); ii != m_all_walls.end(); ii++) ++counter;
//   if(
//      corner1.x > 19.5 && corner1.x < 19.52 &&
//      corner1.y > -0.1 && corner1.y < 0.1 &&
//     corner1.z > 6.59 && corner1.z < 6.6 &&
//     corner2.x > 20.5 && corner2.x < 20.6 &&
//     corner2.y > 1.2 && corner2.y < 1.3 &&
//     corner2.z > 7.6 && corner2.z < 7.8
// 		   )
 //    MSG_DEBUG("Cell::setupWalls", "elements in m_all_walls before outlet-loop = " << ", NUM_NEIGHBORS= ");

  /* Loop over all outlets. */
  for (int n = 0; n < NUM_DIRECT_NEIGHBORS; n++) {
  //  for (list<Cell*>::iterator i = m_outlets[n].begin(); i != m_outlets[n].end(); i++) {
    if (Cell *i = m_outlets[n][0]){
  //       if(
  // 	 corner1.x > 19.5 && corner1.x < 19.52 &&
  // 	 corner1.y > -0.1 && corner1.y < 0.1 &&
  // 	corner1.z > 6.59 && corner1.z < 6.6 &&
  // 	corner2.x > 20.5 && corner2.x < 20.6 &&
  // 	corner2.y > 1.2 && corner2.y < 1.3 &&
  // 	corner2.z > 7.6 && corner2.z < 7.8
  // 		       )
   //   MSG_DEBUG("Cell::setupWalls", "direct outlet " << n << ": " << i->corner1 << " " << i->corner2);

      for (list<Wall*>::iterator j = i->m_walls.begin(); j != i->m_walls.end(); j++) {
        if (find(m_all_walls.begin(), m_all_walls.end(), *j) == m_all_walls.end())
          m_all_walls.push_back(*j);
      }
    }
  }
}


bool Cell::hasParticles() const
{
  bool has = false;

  FOR_EACH_CONST
    (vector< list<Particle*> >,
     m_particles,
     if (!i->empty())
       has = true;
    );

  return has;
}


bool Cell::hasFrozenParticles() const
{
  bool has = false;

  FOR_EACH_CONST
    (vector<list<Particle*> >,
     m_frozen_particles,
     if (!i->empty())
       has = true;
    );

  return has;
}



typedef map<Cell*, AbstractCellLink*> my_map_t;
void Cell::activate()
{
  /* Update cell list. */
  m_manager->activateCell(this);

  /* Notify local cell link. */
  m_local_link->cellActivated();
  m_local_link->cellActivated();

  for (int j = 0; j < m_manager->num_neighbors(); ++j) {
     m_neighbors[j][0]->cellActivated();
     //MSG_DEBUG("Cell::activate", "neighbor activation" << (m_neighbors[j][1]==NULL));
  }
}
void BoundaryCell::activate()
{
  /* Update cell list. */
  m_manager->activateCell(this);

  /* Notify local cell link. */
  m_local_link->cellActivated();
  m_local_link->cellActivated();//why two?

  for (int j = 0; j < m_manager->num_neighbors(); ++j) {
    FOR_EACH
      (vector<AbstractCellLink*>,
       m_neighbors[j],
//       MSG_DEBUG("BoundaryCell::activate", "neighbor activation " << j << " " << (m_neighbors[j][1]));
       (*__iFE)->cellActivated();
       );
  }
}


void Cell::deactivate()
{
  /* Update cell list. */
  m_manager->deactivateCell(this);

  /* Notify cell links. */
  m_local_link->cellDeactivated();
  m_local_link->cellDeactivated();

  for (int j = 0; j < m_manager->num_neighbors(); ++j) {
    FOR_EACH
      (vector<AbstractCellLink*>,
       m_neighbors[j],
       (*__iFE)->cellDeactivated();
       );
  }
}

void BoundaryCell::deactivate()
{
  /* Update cell list. */
  m_manager->deactivateCell(this);

  /* Notify cell links. */
  m_local_link->cellDeactivated();
  m_local_link->cellDeactivated();

  for (int j = 0; j < m_manager->num_neighbors(); ++j) {
    FOR_EACH
      (vector<AbstractCellLink*>,
       m_neighbors[j],
       (*__iFE)->cellDeactivated();
       );
  }
}

void BoundaryCell::setupWalls()
{
  copy(m_walls.begin(), m_walls.end(), inserter(m_all_walls, m_all_walls.begin()));

  // size_t counter = ;
  // for (list<Wall*>::iterator ii = m_all_walls.begin(); ii != m_all_walls.end(); ii++) ++counter;
//   if(
//      corner1.x > 19.5 && corner1.x < 19.52 &&
//      corner1.y > -0.1 && corner1.y < 0.1 &&
//     corner1.z > 6.59 && corner1.z < 6.6 &&
//     corner2.x > 20.5 && corner2.x < 20.6 &&
//     corner2.y > 1.2 && corner2.y < 1.3 &&
//     corner2.z > 7.6 && corner2.z < 7.8
// 		   )
//     MSG_DEBUG("Cell::setupWalls", "elements in m_all_walls before outlet-loop = " << ", NUM_NEIGHBORS= ");

  /* Loop over all outlets. */
  for (int n = 0; n < NUM_DIRECT_NEIGHBORS; n++) {
    for (int i=0; i < m_outlets[n].size(); i++){ if (Cell *c = m_outlets[n][i]){int n_1x_2x = m_manager->m_direct_neighbors[n]; 
	      /* Check whether cell is a direct or an indirect outlet. */
	      if ((corner1 + m_cell_dist[n_1x_2x][i] - c->corner1).abs() < g_geom_eps) {

	    //       if(
	    // 	 corner1.x > 19.5 && corner1.x < 19.52 &&
	    // 	 corner1.y > -0.1 && corner1.y < 0.1 &&
	    // 	corner1.z > 6.59 && corner1.z < 6.6 &&
	    // 	corner2.x > 20.5 && corner2.x < 20.6 &&
	    // 	corner2.y > 1.2 && corner2.y < 1.3 &&
	    // 	corner2.z > 7.6 && corner2.z < 7.8
	    // 		       )
//		   MSG_DEBUG("Cell::setupWalls", "direct outlet: " << c->corner1 << ", " << c->corner2);

		for (list<Wall*>::iterator j = c->m_walls.begin(); j != c->m_walls.end(); j++) {
		  if (find(m_all_walls.begin(), m_all_walls.end(), *j) == m_all_walls.end())
		    m_all_walls.push_back(*j);
		  }
	      }else {
	    // 	if(
	    // 	   corner1.x > 19.5 && corner1.x < 19.52 &&
	    // 	   corner1.y > -0.1 && corner1.y < 0.1 &&
	    // 	  corner1.z > 6.59 && corner1.z < 6.6 &&
	    // 	  corner2.x > 20.5 && corner2.x < 20.6 &&
	    // 	  corner2.y > 1.2 && corner2.y < 1.3 &&
	    // 	  corner2.z > 7.6 && corner2.z < 7.8
	    // 			 )

//		   MSG_DEBUG("Cell::setupWalls", "this cell: " << corner1 << ", " << corner2);
//		   MSG_DEBUG("Cell::setupWalls", "indirect outlet: " << c->corner1 << ", " << c->corner2);
		m_indirect_outlets.push_back(pair< point_t, Cell*> (m_cell_dist[n_1x_2x][i], c));
	      }
      }
    }
  }
}


// this checks for the indirect outlets = "far neighbours" due to, e.g., periodic BCs
bool BoundaryCell::doCollisionIndirectOutlets(Particle *p, const point_t &force, point_t &r, double &t_travelled, point_t &hit_pos, Wall *&wall, IntegratorPosition *integratorP)
{
  bool hit = false;
  point_t old_r = r;
//  int i = 0;
  for (list<pair< point_t, Cell*> >::iterator c = m_indirect_outlets.begin();
       c != m_indirect_outlets.end(); c++) {
    bool new_hit = false;
//    MSG_DEBUG("BoundaryCell::doCollisionIndirectOutlets", i++);

    /* Update position in case this is periodic, etc. */
    /*p->*/r = old_r + c->second->corner1 - corner1 - c->first;

    c->second->checkForHit(p, force, new_hit, t_travelled, hit_pos, wall, integratorP);

#ifdef TRACK_PARTICLE

    if (p->mySlot == TRACK_PARTICLE && p->c == TRACK_PARTICLE_COLOUR && hit) {
      cout << "--> new hit @ " << hit_pos << " after time " << t_travelled << " on wall:" << endl;
      cout << wall->toString() << endl;
    }

#endif


    if (new_hit) {

      hit = true;

      hit_pos = hit_pos - c->second->corner1 + corner1 + c->first;
    }
  }
  return hit;
}

