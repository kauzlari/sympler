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



#ifndef __INLET_CELL_H
#define __INLET_CELL_H

#include "cell.h"
#include "pc_inlet.h"


/*!
 * This specialization of a \a BoundaryCell does create particles in the inlet whenever
 * a particle leaves this cell through a certain side.
 */
class InletCell_Create: public BoundaryCell
{
protected:
  /*!
   * Request creation of new particles from this \a ParticleCreator
   */
  ParticleCreator *m_pc;

  /*!
   * New particles are created if a particles leaves in direction \a m_dir
   * (+1, -1) of the axis \a m_axis (0, 1, 2)
   */
  int m_axis;
  
  /*!
   * New particles are created if a particles leaves in direction \a m_dir
   * (+1, -1) of the axis \a m_axis (0, 1, 2)
   */
  int m_dir;

  /*!
   * Function called whenever a particle leaves a cell. Check if it leaves into the
   * right direction and call the particle creator.
   * @param off Offset of the cell the particle is entering
   * @param n Index of the direction in which the particle is leaving (supplemental to \a off)
   */
  virtual void particleLeftCell(const int_point_t &off, int n) {
    if (off[m_axis] == m_dir)
      m_pc->createParticle();
  }

public:
  /*!
   * Constructor
   * @param mgr Pointer to the cell subdivision manager
   * @param pc Pointer to the particle creator to use for creation of new particles
   * @param axis Axis
   * @param dir Direction
   * @param group Group of the cell
   */
  InletCell_Create(ManagerCell *mgr, ParticleCreator *pc, int axis, int dir, region_t *r, int group = 0);

  /*!
   * Constructor
   * @param mgr Pointer to the cell subdivision manager
   * @param pc Pointer to the particle creator to use for creation of new particles
   * @param axis Axis
   * @param dir Direction
   * @param c1 Bottom left corner of the cell
   * @param c2 Top right corner of the cell
   * @param group Group of the cell
   */
  InletCell_Create(ManagerCell *mgr, ParticleCreator *pc, int axis, int dir,
            cuboid_t cuboid, int_point_t cell_pos, region_t *r, int group = 0);

  /*!
   * Destructor
   */
  virtual ~InletCell_Create();
};


/*!
 * This specialization of a \a Cell does delete particles in the outlet whenever
 * a particle leaves this cell through a certain side (which means enters the oulet).
 */
class InletCell_Delete: public BoundaryCell
{
protected:
  /*!
   * Request creation of new particles from this \a ParticleCreator
   */
  ParticleCreator *m_pc;

  /*!
   * New particles are created if a particles leaves in direction \a m_dir
   * (+1, -1) of the axis \a m_axis (0, 1, 2)
   */
  int m_axis;
  
  /*!
   * New particles are created if a particles leaves in direction \a m_dir
   * (+1, -1) of the axis \a m_axis (0, 1, 2)
   */
  int m_dir;

  /*!
   * Function called whenever a particle leaves a cell. Check if it leaves into the
   * right direction and call the particle creator.
   * @param off Offset of the cell the particle is entering
   * @param n Index of the direction in which the particle is leaving (supplemental to \a off)
   */
  virtual void particleLeftCell(const int_point_t &off, int n) {
    if (off[m_axis] == m_dir)
      m_pc->deleteParticle();
  }

public:
  /*!
   * Constructor
   * @param mgr Pointer to the cell subdivision manager
   * @param pc Pointer to the particle creator to use for creation of new particles
   * @param axis Axis
   * @param dir Direction
   * @param group Group of the cell
   */ 
 InletCell_Delete(ManagerCell *mgr, ParticleCreator *pc, int axis, int dir, region_t *r, int group = 0);

  /*!
   * Constructor
   * @param mgr Pointer to the cell subdivision manager
   * @param pc Pointer to the particle creator to use for creation of new particles
   * @param axis Axis
   * @param dir Direction
   * @param c1 Bottom left corner of the cell
   * @param c2 Top right corner of the cell
   * @param group Group of the cell
   */
  InletCell_Delete(ManagerCell *mgr, ParticleCreator *pc, int axis, int dir,
            cuboid_t cuboid, int_point_t cell_pos, region_t *r, int group = 0);

  /*!
   * Destructor
   */
  virtual ~InletCell_Delete();
};


#endif
