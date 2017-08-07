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



#ifndef __SMART_LIST_H
#define __SMART_LIST_H

#include <list>
#include <vector>

#include "general.h"

using namespace std;

#define CHUNK_LEN 65536
#define CHUNK_SH  16


/* Macros */

#define SLOT2CHUNKID(slot)   (slot >> CHUNK_SH)
#define SLOT2INDEX(slot)     (slot & (CHUNK_LEN-1))


/* Never ever modify slot, prev or next. */
#define SMARTLISTENTRY(T)                       \
size_t mySlot;                                  \
T *prev, *next;


/* Loop over all elements. */
#define SL_FOR_EACH(type, sl, code)                     \
for (type *__iSLFE = sl.first(); __iSLFE != NULL; __iSLFE = __iSLFE->next) {    \
  code                                                  \
} while(0)


/* Loop over all elements that are being hold in memory. */
#define SL_FOR_EACH_IN_MEMORY(type, sl, code)                           \
{                                                                       \
  vector<type*> chunks = sl.chunks();                                   \
  for (vector<type*>::iterator cur_chunk = chunks.begin(); cur_chunk != chunks.end(); \
    ++cur_chunk) { \
    for (int j = 0; j < CHUNK_LEN; ++j) {                               \
      type *i = &(*cur_chunk)[j];                                       \
      code                                                              \
    }                                                                   \
  }                                                                     \
} while(0)


/* SmartList */

/*!
 * A template for a doubly linked list class, which partitions the list
 * into a certain number of chunks to assign unique indices to list entries
 * and ensures that pointers to list entries stay valid all the time.
 */
template <typename T>
class SmartList
{
protected:
  /*
   * Empty list INDEX counter
   */
  size_t m_emptyIndex;
  /*!
   * List capacity
   */
  size_t m_capacity;

  /*!
   * Current list size
   */
  size_t m_size;

  /*!
   * Pointer to the first list entry
   */
  T *m_first;

  /*!
   * Pointer to the last list entry
   */
  T *m_last;

  /*!
   * List of chunks
   */
  vector<T*> m_chunks;

  /*!
   * List of free slots
   */
  list<size_t> m_free_slots;

  /* Thread safety */
#ifdef ENABLE_PTHREAD
  pthread_mutex_t m_mutex;
#endif
    
  /*!
   * Expand the list capacity, i.e., add another chunk and
   * update the free slot information.
   */
  void expandCapacity() {
    /* size_t first_slot, last_slot; */
    T* chunk;

//      MSG_DEBUG("SmartList::expandCapacity", "# chunks now = " << m_chunks.size()+1);

    /* first_slot = m_chunks.size()*CHUNK_LEN; */
    chunk = new T[CHUNK_LEN];
    //        chunk = (T*) malloc(sizeof(T)*CHUNK_LEN);
    /* last_slot = (m_chunks.size()+1)*CHUNK_LEN; */
    
    /*
    for (size_t slot = first_slot; slot < last_slot; slot++) {
      m_free_slots.push_back(slot);
      }*/

    m_chunks.push_back(chunk);

    m_capacity += CHUNK_LEN;
  }

  /*!
   * Initialize the list => create first chunk
   */
  void init() {
#ifdef ENABLE_PTHREAD
    pthread_mutex_init(&m_mutex, &g_mutex_attr);
#endif
    expandCapacity();
  }

public:
  /*!
   * Type of list entries
   */
  typedef T type;

list <size_t> freeSlots()
{
  return m_free_slots;
} 

  /* For compatibility with STL */
  class iterator
  {
  protected:
    T *m_i;
  public:
    iterator(T *i): m_i(i) { }
    ~iterator() { }

    iterator &operator++() {
      m_i = m_i->next;
    }
    bool operator!=(const iterator &i2) {
      return m_i != i2.m_i;
    }

    T& operator->() {
      return &m_i;
    }
    const T& operator->() const {
      return &m_i;
    }

    T& operator*() {
      return &m_i;
    }
    const T& operator*() const {
      return &m_i;
    }
  };
  /* -------------------------- */

  /*!
   * Constructor
   */
  SmartList(): m_emptyIndex(0), m_capacity(0), m_size(0), m_first(NULL), m_last(NULL) {
    init();
  }

  /*!
   * Copy constructor
   * @param copy List to copy from
   */
  SmartList(const SmartList<T> &copy): m_emptyIndex(0), m_capacity(0), m_size(0), m_first(NULL), m_last(NULL) {
    //        MSG_DEBUG("SmartList", "Copy constructor called.");
    if (copy.m_size == 0) {
      /* Everything okay, we just initialize as usual. */
      init();
    } else {
      /* We can't copy a SmartList. */
      throw gError("SmartList::SmartList", "Copying a SmartList is not yet implemented.");
    }
  }

  /*!
   * Destructor
   */
  virtual ~SmartList() {
    typename vector<T*>::iterator chunks_end = m_chunks.end();
    for (typename vector<T*>::iterator i = m_chunks.begin(); i != chunks_end; i++)
      delete [] *i;
#ifdef ENABLE_PTHREAD
    pthread_mutex_destroy(&m_mutex);
#endif
  }

  SmartList<T> &operator=(const SmartList<T> &copy) {
    //        MSG_DEBUG("SmartList", "operator= called.");

    throw gError("SmartList::operator=", "Copying a SmartList is not yet implemented.");
  }

  /*!
   * Create a new entry.
   */
  T &newEntry() {
    size_t slot;

#ifdef ENABLE_PTHREAD
    pthread_mutex_lock(&m_mutex);
#endif

    //printf("m_size %d, m_capacity %d, m_emptyIndex %d\n",m_size, m_capacity,m_emptyIndex);
    if (m_size == m_capacity) {
      //printf("b\n");
      assert(m_free_slots.empty());

      expandCapacity();
    }

    if (m_emptyIndex < m_capacity)
    {
      //printf("c\n");
      slot = m_emptyIndex++;
    }
    else
    {
      //printf("d\n");
      assert(!m_free_slots.empty());
      slot = m_free_slots.front();
      m_free_slots.pop_front();
    }
    ++m_size;

    /* Establish prev and next links. */
    /* Fixme!!! Currently not cache coherent. */
    //        MSG_DEBUG("SmartList::newEntry", "slot = " << slot << ", chunk id = " << SLOT2CHUNKID(slot) << ", index = " << SLOT2INDEX(slot));

    T &entry = m_chunks[SLOT2CHUNKID(slot)][SLOT2INDEX(slot)];

    /* If m_first == NULL this is the first element added to the list. */
    if (!m_first)
      m_first = &entry;

    /* Add element to the end of the list. */
    entry.prev = m_last;
    /* If this is not the first element, tell the last one which one is
       the new last one. */
    if (m_last)
      m_last->next = &entry;
    /* The new one is now last element. */
    entry.next = NULL;
    m_last = &entry;

#ifdef ENABLE_PTHREAD
    pthread_mutex_unlock(&m_mutex);
#endif

    entry.mySlot = slot;

    return entry;
  }

  /*!
   * Delete an entry
   * @param entry Entry to delete
   */
  void deleteEntry(T &entry) {
#ifdef ENABLE_PTHREAD
    pthread_mutex_lock(&m_mutex);
#endif

    assert(m_size > 0);
    assert(m_first != NULL);
    assert(m_last != NULL);

    if (m_first->mySlot == entry.mySlot)
      m_first = m_first->next;
    if (m_last->mySlot == entry.mySlot)
      m_last = m_last->prev;

    m_free_slots.push_back(entry.mySlot);
    --m_size;

    if (entry.prev)
      entry.prev->next = entry.next;
    if (entry.next)
      entry.next->prev = entry.prev;


// if(entry.mySlot == 54)
// MSG_DEBUG("SmartList::basics", "id of this pair = " << &entry << "  previous pair id = " << entry.prev << "  next pair id = " << entry.next << " prev of next = " << entry.next->prev << " next of prev = " << entry.prev->next << " prev of prev of next = " << entry.next->prev->prev << " next of next of prev" << entry.prev->next->next);

#ifdef ENABLE_PTHREAD
    pthread_mutex_unlock(&m_mutex);
#endif
  }

  /*!
   * Delete an entry
   * @param slot Entry slot to delete
   */
  void deleteEntry(int slot) {
    deleteEntry(operator[](slot));
  }

  /*!
   * Delete the whole list
   */
  void clear() {
    m_size = 0;
    m_emptyIndex = 0;
    m_free_slots.clear();

    /*
    for(size_t i = 0; i < m_capacity; ++i)
      {
	m_free_slots.push_back(i);
      }
    */

    m_first = NULL;
    m_last = NULL;
  }

  /*!
   * Return the size (number of entries) of the list
   */
  size_t size() const {
    return m_size;
  }

  /*!
   * Return the current capacity of the list
   */
  size_t capacity() const {
    return m_capacity;
  }

  /*!
   * Return entry with slot \a slot
   * @param slot Slot of the entry to return
   */
  T &operator[](size_t slot) {
    return m_chunks[SLOT2CHUNKID(slot)][SLOT2INDEX(slot)];
  }

  /*!
   * Return entry with slot \a slot
   * @param slot Slot of the entry to return
   */
  const T &operator[](size_t slot) const {
    return m_chunks[SLOT2CHUNKID(slot)][SLOT2INDEX(slot)];
  }

  /*!
   * Return the first element
   */
  T *first() {
    return m_first;
  }

  /*!
   * Return the last element
   */
  T *last() {
    return m_last;
  }

  iterator begin() {
    return iterator(m_first);
  }
  iterator end() {
    return iterator(NULL);
  }

  /*!
   * Return all chunks
   */
  vector<T*> &chunks() {
    return m_chunks;
  }
};

/*!
* This struct is intended to be used for \a SmartList s of int, double, ...
*/
template <typename T>
struct PrimitiveSLEntry
{
  /*!
  * the value we would like to save
  */
  T m_val;
  
  SMARTLISTENTRY(PrimitiveSLEntry<T>)
      
  /*!
   *Copy everything *except* mySlot, prev and next. 
   * @param entry The entry to copy from
   */
  void copyFrom(const PrimitiveSLEntry<T> &entry) {
    m_val = entry.m_val;
  }
  
  /*!
   * constructor
   */
  PrimitiveSLEntry() {
#ifdef ENABLE_PTHREADS
    pthread_mutex_init(&m_mutex, &g_mutex_attr);
#endif
  } 

  /*!
   * Copy constructor
   * @param entry The entry to copy from
   */
  PrimitiveSLEntry(const PrimitiveSLEntry<T>& entry) {
#ifdef ENABLE_PTHREADS
    pthread_mutex_init(&m_mutex, &g_mutex_attr);
#endif
    copyFrom(entry);
  } 

  /*!
   * Assignment operator
   * @param entry The entry to copy from
   */
  PrimitiveSLEntry &operator=(const PrimitiveSLEntry<T> &entry) {
    copyFrom(entry);

    return *this;
  }
  
};

#endif
