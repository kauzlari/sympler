/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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



#ifndef __WEIGHTING_FUNCTION_H
#define __WEIGHTING_FUNCTION_H

#include "node.h"
#include "pairdist.h"
#include "smart_enum.h"

using namespace std;

/* --- WeightingFunction --- */

/*!
 * Base class for the definition of weighting functions
 */
class WeightingFunction: public Node
{
 protected:
  /*!
   * The name of this weighting function (for later reference)
   */
  string m_name;

  /*!
   * The cut-off radius of this weighting function
   */
  double m_cutoff;

  /*!
   * Register all basic properties with the property list.
   */
  void init();

 public:
  /*!
  * Default constructor
   */
   WeightingFunction();

  /*!
    * Constructor that should be used
    * @param parent The parent node
   */
  WeightingFunction(Node *parent);

  /*!
   * Destructor
   */
  virtual ~WeightingFunction();

  /*!
   * Return the name of this weighting function
   */
  const string &name() const {
    return m_name; 
  }

  /*!
   * Return the cut-off radius of this weighting function
   */
  double cutoff() const {
    return m_cutoff;
  }


  /* The difference is the following:
     - interpolate   gives W(r, r0) where W is the interpolation function
     - weight        gives w(r, r0) with w defined a -r w = grad_r W
     - localGradient gives u(r, r0) with u = grad_r0 W
  */

  /*!
   * Return the value of the interpolation function.
   * @param r The distance of the two participating particles
   * @param p The position of the second particle (for which the property is not calculated)
   * @param dist_from_wall If not NULL, then the distance to the wall will be returned
   */
  virtual double interpolate(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const = 0;

  /*!
   * Return the value of the weighting function, which is the first derivative of the
   * interpolation function with respect to r. It is given by -rw(r,r0) = grad_r W(r, r0)
   * where W is the interpolation function and w is the weighting function.
   * @param r The distance of the two participating particles
   * @param p The position of the second particle (for which the property is not calculated)
   * @param dist_from_wall If not NULL, then the distance to the wall will be returned
   */
  virtual double weight(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const = 0;

  /*!
   * Return the value of the local gradient function, which is the first derivative of the
   * interpolation function with respect to p. It is given by u(r,r0) = -grad_r0 W(r, r0)
   * where W is the interpolation function and u is the local gradient function.
   * @param r The distance of the two participating particles
   * @param p The position of the second particle (for which the property is not calculated)
   * @param dist_from_wall If not NULL, then the distance to the wall will be returned
   */
  virtual point_t localGradient(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const = 0;

  /*!
   * For a constant surface entropy factor return the gradient
   * grad(\int d^2 r_s W(r_s-r_i, r_i))
   * @param p The position of the particle r_i
   */
  virtual point_t gradientSurfaceWeight(const point_t &p) const {
    throw gError
      ("WeightingFunction::gradientSurfaceWeight",
       "Cannot be used. Please specifiy a different weighting function.");
  }
};



/* --- WeightingFunctionWithWall --- */

/*!
 * Base class for simple weighting functions in the presence of two parallel walls.
 * Note: This is a hack. Write this more general.
 */
class WeightingFunctionWithWall: public WeightingFunction
{
 protected:
  /*!
   * Where is the wall located?
   */
  size_t m_wall_dir;
  
  /*!
   * Position of wall to the left in \a m_wall_dir direction
   */
  double m_left_wall;

  /*!
   * Position of wall to the right in \a m_wall_dir direction
   */
  double m_right_wall;

  /*!
   * Initialize the property list for \a m_wall_dir, \a m_left_wall and \a m_right_wall
   */
  void init();

  /*!
   * Return the value of the interpolation function.
   * @param r The distance of the two participating particles
   * @param normal The normal vector to the wall
   * @param dist The distance from the wall
   */
  virtual double interpolateWithDist(const Pairdist *r, const point_t &normal, double dist) const = 0;

  /*!
   * Return the value of the weighting function, which is the first derivative of the
   * interpolation function with respect to r. It is given by -rw(r,r0) = grad_r W(r, r0)
   * where W is the interpolation function and w is the weighting function.
   * @param r The distance of the two participating particles
   * @param normal The normal vector to the wall
   * @param dist The distance from the wall
   */
  virtual double weightWithDist(const Pairdist *r, const point_t &normal, double dist) const = 0;
 
   /*!
   * Return the value of the local gradient function, which is the first derivative of the
   * interpolation function with respect to p. It is given by u(r,r0) = -grad_r0 W(r, r0)
   * where W is the interpolation function and u is the local gradient function.
   * @param r The distance of the two participating particles
   * @param normal The normal vector to the wall
   * @param dist The distance from the wall
   */
  virtual point_t localGradientWithDist(const Pairdist *r, const point_t &normal, double dist) const = 0;

  /*!
   * For a constant surface entropy factor return the gradient
   * grad(\int d^2 r_s W(r_s-r_i, r_i))
   * @param normal The normal vector to the wall
   * @param dist The distance from the wall
   */
  virtual point_t gradientSurfaceWeightWithDist(const point_t &normal, double dist) const {
    throw gError
      ("WeightingFunctionWithWall::gradientSurfaceWeightWithDist",
       "Cannot be used. Please specifiy a different weighting function.");
  }

 public:
  /*!
   * Constructor
   * @param parent The parent node
   */
  WeightingFunctionWithWall(Node *parent);

  /*!
   * Destructor
   */
  virtual ~WeightingFunctionWithWall();

  /*!
   * Return the value of the interpolation function.
   * @param r The distance of the two participating particles
   * @param p The position of the second particle (for which the property is not calculated)
   * @param dist_from_wall If not NULL, then the distance to the wall will be returned
   */
  virtual double interpolate(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    point_t normal = {{{ 0, 0, 0 }}};

    if(!r) {
      throw gError("WeightingFunctionWithWall::interpolate for " + className(), "r=NULL detected. This is currently not supported! Please use a different WeightingFunction");
      return 0.;
    }
      
    if (r && r->abs() > m_cutoff)
      return 0.;

    if (p[m_wall_dir] < m_left_wall+m_cutoff) {
      normal[m_wall_dir] = 1;

      if (dist_from_wall)
	*dist_from_wall = p[m_wall_dir]-m_left_wall;

      return interpolateWithDist(r, normal, (p[m_wall_dir]-m_left_wall)/m_cutoff);
    } else if (p[m_wall_dir] > m_right_wall-m_cutoff) {
      normal[m_wall_dir] = -1;

      if (dist_from_wall)
	*dist_from_wall = m_right_wall-p[m_wall_dir];

      return interpolateWithDist(r, normal, (m_right_wall-p[m_wall_dir])/m_cutoff);
    } else
      return interpolateWithDist(r, normal, 1);
  }

  /*!
   * Return the value of the weighting function, which is the first derivative of the
   * interpolation function with respect to r. It is given by -rw(r,r0) = grad_r W(r, r0)
   * where W is the interpolation function and w is the weighting function.
   * @param r The distance of the two participating particles
   * @param p The position of the second particle (for which the property is not calculated)
   * @param dist_from_wall If not NULL, then the distance to the wall will be returned
   */
  virtual double weight(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    point_t normal = {{{ 0, 0, 0 }}};

    if (r && r->abs() > m_cutoff)
      return 0.;

    if (p[m_wall_dir] < m_left_wall+m_cutoff) {
      normal[m_wall_dir] = 1;
      
      if (dist_from_wall)
	*dist_from_wall = p[m_wall_dir]-m_left_wall;

      return weightWithDist(r, normal, (p[m_wall_dir]-m_left_wall)/m_cutoff);
    } else if (p[m_wall_dir] > m_right_wall-m_cutoff) {
      normal[m_wall_dir] = -1;

      if (dist_from_wall)
	*dist_from_wall = m_right_wall-p[m_wall_dir];

      return weightWithDist(r, normal, (m_right_wall-p[m_wall_dir])/m_cutoff);
    } else
      return weightWithDist(r, normal, 1);
  }

  /*!
   * Return the value of the local gradient function, which is the first derivative of the
   * interpolation function with respect to p. It is given by u(r,r0) = -grad_r0 W(r, r0)
   * where W is the interpolation function and u is the local gradient function.
   * @param r The distance of the two participating particles
   * @param p The position of the second particle (for which the property is not calculated)
   * @param dist_from_wall If not NULL, then the distance to the wall will be returned
   */
  virtual point_t localGradient(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    point_t normal = {{{ 0, 0, 0 }}};

    if (r && r->abs() > m_cutoff) {
      point_t g = {{{ 0, 0, 0 }}};

      return g;
    }

    if (p[m_wall_dir] < m_left_wall+m_cutoff) {
      normal[m_wall_dir] = 1;

      if (dist_from_wall)
	*dist_from_wall = p[m_wall_dir]-m_left_wall;

      return localGradientWithDist(r, normal, (p[m_wall_dir]-m_left_wall)/m_cutoff);
    } else if (p[m_wall_dir] > m_right_wall-m_cutoff) {
      normal[m_wall_dir] = -1;

      if (dist_from_wall)
	*dist_from_wall = m_right_wall-p[m_wall_dir];

      return localGradientWithDist(r, normal, (m_right_wall-p[m_wall_dir])/m_cutoff);
    } else
      return localGradientWithDist(r, normal, 1);
  }

  /*!
   * For a constant surface entropy factor return the gradient
   * grad(\int d^2 r_s W(r_s-r_i, r_i))
   * @param p The position of the particle r_i
   */
  virtual point_t gradientSurfaceWeight(const point_t &p) const {
    point_t normal = {{{ 0, 0, 0 }}};

    if (p[m_wall_dir] < m_left_wall+m_cutoff) {
      normal[m_wall_dir] = 1;
      return gradientSurfaceWeightWithDist(normal, (p[m_wall_dir]-m_left_wall)/m_cutoff);
    } else if (p[m_wall_dir] > m_right_wall-m_cutoff) {
      normal[m_wall_dir] = -1;
      return gradientSurfaceWeightWithDist(normal, (m_right_wall-p[m_wall_dir])/m_cutoff);
    } else
      return gradientSurfaceWeightWithDist(normal, 1);
  }


  /*!
   * Returns the direction in which to find the wall
   */
  size_t wallDir() const {
    return m_wall_dir;
  }
};


/*--- Factory --- */

class WeightingFunction_Factory: public SmartEnum<WeightingFunction_Factory>
{
 public:
  virtual WeightingFunction *instantiate(Node *parent) const = 0;

 protected:
  WeightingFunction_Factory(const string &name)
    : SmartEnum<WeightingFunction_Factory>(name) { }
};


template <class T>
class WeightingFunction_Register: public WeightingFunction_Factory
{
 public:
  WeightingFunction_Register(const string &name)
    : WeightingFunction_Factory(name) { }

  virtual WeightingFunction *instantiate(Node *parent) const {
    return new T(parent);
  }
};


#endif
