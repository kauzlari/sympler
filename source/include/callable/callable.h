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



#ifndef __CALLABLE_H
#define __CALLABLE_H 

#include "smart_enum.h"
#include "node_many_children.h"

using namespace std;

class Simulation;

/*!
 * A \a Callable is the base class for modules that need to
 * perform tasks during every time step. For example, a specialization
 * of a \a Callable is the \a Thermostat, which modifies the properties
 * of the particles within every time step to ensure the system has 
 * a predefined temperature.
 *
 * A \a Callable is called in the \a Controller directly after the integration
 * is finished by a call to the \a call method.
 */
class Callable : public NodeManyChildren
{
 public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  Callable(Simulation* simulation);

  /*!
   * Destructor
   */
  virtual ~Callable();

  /*!
   * Just throw an error. The default \c Callable cannot have children,
   * subclasses might.
   */
  virtual Node* instantiateChild(const string& name);

  /*!
   * This method is called when the \a Callable should perform its tasks.
   */
  virtual void call(size_t timestep) = 0;
};



//---- Factories ----

class Callable_Factory: public SmartEnum<Callable_Factory>
{
public:
    virtual Callable *instantiate(Simulation *simulation) const = 0;

protected:
    Callable_Factory(const string &name)
    : SmartEnum<Callable_Factory>(name) { }
};


template <class T>
class Callable_Register: public Callable_Factory
{
public:
    Callable_Register(const string &name)
    : Callable_Factory(name) { }

    virtual Callable *instantiate(Simulation *simulation) const;
};



//---- Inline functions ----

template <class T>
inline Callable *Callable_Register<T>::instantiate(Simulation *simulation) const
{
    return new T(simulation);
}


#endif
