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



#ifndef __SMART_POINTER_H
#define __SMART_POINTER_H

#include <stdlib.h>
#include "general.h"

struct UntypedSmartPointer
{
protected:
  void *m_value;
  int *m_ref_count;

  inline void assign(void* value, int *ref_count) {
    release();
    if (value) {
      m_value = value;
      m_ref_count = ref_count;
      (*m_ref_count)++;
    }
  }

public:
  inline UntypedSmartPointer(): m_value(NULL), m_ref_count(NULL) {
  }
  inline UntypedSmartPointer(void *value): m_value(NULL), m_ref_count(NULL) {
    set(value);
  }
  inline UntypedSmartPointer(const UntypedSmartPointer &sp): m_value(NULL), m_ref_count(NULL) {
    assign(sp.m_value, sp.m_ref_count);
  }
  inline ~UntypedSmartPointer() {
    release();
  }

  inline void alloc(size_t size) {
    release();
    m_value = malloc(size);
    m_ref_count = new int;
    (*m_ref_count) = 1;
  }

  inline void set(void *value) {
    release();
    m_value = value;
    m_ref_count = new int;
    (*m_ref_count) = 1;
  }

  inline void release() {
    if (m_value) {
      (*m_ref_count)--;
      if (!*m_ref_count) {
        free(m_value);
        delete m_ref_count;
      }

      /* Make the pointer NULL again in case release is called by hand */
      m_value = NULL;
    }
  }

  inline UntypedSmartPointer &operator=(const UntypedSmartPointer &sp) {
    assign(sp.m_value, sp.m_ref_count);
    return *this;
  }

  inline void *value() const {
    return m_value;
  }

  inline bool isNull() const {
    return m_value == NULL;
  }    
};


/*!
 * A SmartPointer is an auto-releasing pointer, that means when the last reference to
 * the object disappears, the object will be freed. When instantiating a SmartPointer
 * it will always be NULL. To get an actual object use the alloc or set methods.
 * 
 * To generate a new reference to an object, use the copy constructor or the '='
 * operator.
 */
template<typename T>
struct SmartPointer
{
  protected:
  /*!
   * Pointer to the encapsulated object, shared with other SPs
   */
  T *m_value;

  /*!
   * Pointer to the reference counter, shared with other SPs
   */
  int *m_ref_count;

  /*!
   * Assign an object and an existing reference counter to this smart,
   * will also increase the reference counter by one.
   */
  inline void assign(T* value, int *ref_count) {
    release();
    if (value) {
      m_value = value;
      m_ref_count = ref_count;
      (*m_ref_count)++;
    }
  }
    
  public:
  /*!
   * Constructor
   */
  inline SmartPointer(): m_value(NULL), m_ref_count(NULL) {
  }

  /*!
   * Construct a new \a SmartPointer for existing object \a value
   * @param value Pointer to the object to be encapsulated
   */
  inline SmartPointer(T *value): m_value(NULL), m_ref_count(NULL) {
    set(value);
  }

  /*!
   * Shallow copy the \a SmartPointer \a sp
   * @param sp Copy from this \a SmartPointer
   */
  inline SmartPointer(const SmartPointer<T> &sp): m_value(NULL), m_ref_count(NULL) {
    assign(sp.m_value, sp.m_ref_count);
  }

  /*!
   * Destructor
   */
  inline ~SmartPointer() {
    release();
  }

  /*!
   * Allocate a new object and a new reference counter
   */
  inline void alloc() {
    release();
    m_value = new T();
    m_ref_count = new int;
    (*m_ref_count) = 1;
  }

  /*!
   * Make a deep copy of this object and return a new \a SmartPointer holding this
   * copy (with initial reference count = 1, of course)
   */
  inline SmartPointer<T> deepCopy() {
	  return SmartPointer<T>(new T(*m_value));
  }

  /*!
   * Set the encapsulated object to \a value. The reference counter will be
   * reset.
   */
  inline void set(T *value) {
    release();
    m_value = value;
    m_ref_count = new int;
    (*m_ref_count) = 1;
  }

  /*!
   * Release this \a SmartPointer. If the reference count drops to zero,
   * the encapsulated object will be released as well.
   */
 /* inline */void release() {
    if ((m_value!=0) && (m_ref_count!=0)) { 
      (*m_ref_count)--;
      if (!*m_ref_count) {
        delete m_value;
        delete m_ref_count; 
      }
    
      /* Make the pointer NULL again in case release is called by hand */
      m_value = NULL;
    }
  }    

  /*!
   * Shallow copy the \a SmartPointer \a sp
   * @param sp Copy from this \a SmartPointer
   */
  inline SmartPointer<T> &operator=(const SmartPointer<T> &sp) {
    assign(sp.m_value, sp.m_ref_count);
    return *this;
  }

  /*!
   * Compare whether two \a SmartPointer s point to the same object
   * @param sp \a SmartPointer to compare with
   */
  inline bool operator==(const SmartPointer<T> &sp) {
    return (*m_value) == (*sp.m_value);
  }

  /* fixme!!! <, >, <=, etc. to be implemented */

  /*!
   * Return a pointer to the encapsulated object.
   */
  inline T *value() const {
    return m_value;
  }

  /*!
   * Return a pointer to the encapsulated object. This makes ubiquious
   * object->method() calls possible where \a object is a \a SmartPointer
   * \a method() is a method of the encapsulated object.
   */
  inline T *operator->() {
    return m_value;
  }

  /*!
   * Return a pointer to the encapsulated object. This makes ubiquious
   * object->method() calls possible where \a object is a \a SmartPointer
   * \a method() is a method of the encapsulated object.
   */
  inline const T *operator->() const {
    return m_value;
  }

  /*!
   * Return a pointer to the encapsulated object.
   */
  inline T &operator*() {
    return *m_value;
  }

  /*!
   * Return a pointer to the encapsulated object.
   */
  inline const T &operator*() const {
    return *m_value;
  }

  /*!
   * Does this SmartPointer point to NULL?
   */
  inline bool isNull() {
    return m_value == NULL;
  }
  
  /*!
   * Increase reference count by n.
   * @param n increment
   */
  // dirty hack to make memcpy of data storage possible
  inline void incRefCount(int n = 1) {
    if ((m_ref_count!=0)) (*m_ref_count)+=n;
  }
  
  /*!
   * Get reference count.
   */
 
  inline int getRefCount() {
    if (m_ref_count)
    	return (*m_ref_count);
    	else return 0;
  }
};

#endif
