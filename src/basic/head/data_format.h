/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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



#ifndef __DATAFORMAT_H
#define __DATAFORMAT_H

#include <stdlib.h>

#include <map>
#include <string> 
#include <vector>

#include "misc.h"
#include "general.h"
#include "smart_pointer.h"

/* Matrix datatypes with TNT library */
#ifdef HAVE_TNT_TNT_H

#define WITH_ARRAY_TYPES
#endif
#ifdef WITH_ARRAY_TYPES

#include <tnt/tnt.h>
#include "MArray2D.h"
#endif

using namespace std;

/* DataFormat stores a description of a data set and Data stores the actual data */

#define MAX_STR_LEN 80

typedef SmartPointer< vector<int> > vector_int_sp;
typedef SmartPointer< vector<double> > vector_double_sp;
typedef SmartPointer< vector<point_t> > vector_point_sp;
typedef SmartPointer< vector<tensor_t> > vector_tensor_sp;

/* Datatypes for use with TNT library *andern*/
#ifdef WITH_ARRAY_TYPES
typedef SmartPointer< MArray2D > array2d_double_sp;
#endif

#define IDDF_DATA_FORMAT "__data_format"

#define DF_UNKNOWN 0           /* Unknown */
#define DF_POINTS 1            /* Collection of points. */
#define DF_STRUCTURED_GRID 2   /* Points on a grid. Number of cells to be taken from 'n_cells' */
#define IDDF_N_CELLS "n_cells"

/* That's it for now. New format should be described here. */

/* Macro to perform operation on all SmartPointer datatypes */
#define NOOP while(0)
#ifdef WITH_ARRAY_TYPES

#define	DATAFORMAT_ARRAY_SWITCH(attr,m_data,oper) \
	case DataFormat::MArray2D:\
		oper(attr,m_data,array2d_double_sp); break
#else
#define	DATAFORMAT_ARRAY_SWITCH(attr,m_data,oper) NOOP
#endif

#define DATAFORMAT_CONTAINER_SWITCH(attr,m_data,oper) switch((attr)->datatype) { \
  case DataFormat::VECTOR_DOUBLE:				       \
    oper((attr),(m_data),vector_double_sp); break;		       \
  case DataFormat::VECTOR_INT:					       \
    oper((attr),(m_data),vector_int_sp); break;			       \
  case DataFormat::VECTOR_POINT:				       \
    oper((attr),(m_data),vector_point_sp); break;		       \
  case DataFormat::VECTOR_TENSOR:				       \
    oper((attr),(m_data),vector_tensor_sp); break;		       \
    DATAFORMAT_ARRAY_SWITCH((attr),(m_data),oper);			\
  default: /*Yes, it will be called for other formats and should then do NOTHING!*/ break; \
}


//---- Classes ----

class Data;
typedef SmartPointer<Data> data_sp;

/*!
 * A \a DataFormat stores the structure of the information stored in a \a Data.
 * A \a Data can be regarded as a run-time extensible structure.
 *
 * DataFormat conventions:
 * For certain output modules the data must contain a field named '__data_format' of
 * type int. This describes the kind of data passed on to postprocessors or
 * output modules.
 * The 'real' data begins AFTER the field '__data_format', additional information such 
 * as 'n_cells' or 'time' have to come BEFORE '__data_format'.
 * Currently, the following identifiers are possible:
 */
class DataFormat
{
public:
  /*!
   * The different data types a \a Data can store.
   */
  enum datatype_t {
      INT = 0, DOUBLE = 1, INT_POINT = 2, POINT = 3, TENSOR = 4, STRING = 5,
      VECTOR_INT = 6, VECTOR_DOUBLE = 7, VECTOR_POINT = 8, VECTOR_TENSOR = 9,
#ifdef WITH_ARRAY_TYPES
       MArray2D= 10, EO_DATATYPE = 11
#else
      EO_DATATYPE = 10
#endif
  };

  /*!
   * Contains the size (in bytes) for the different type a \a Data can store.
   */
  static size_t c_size_of_datatype[EO_DATATYPE];

  /*!
   * Description of an attribute, or entry, in the structure
   */
  struct attribute_t {
    /*!
     * A unique string identifier.
     */
    string name;

    /*!
     * Index of the attribute.
     */
    size_t index;

    /*!
     * Start of the attribute's data from the beginning of the overall data
     * stored in a \a Data.
     */
    size_t offset;

    /*!
     * Type of data.
     */
    datatype_t datatype;

    /*!
     * Can this attribute be cleared using clear?
     */
    bool persistent;

    /*!
     * Short description for use in user defined functions, i.e. n for density
     */
    string symbol;

    /*!
     * Return a string identifier of the datatype
     * Member function version
     */
    string datatypeAsString() const;

    /*!
     * Return a string identifier of the datatype
     * Static member function version
     */
    static string datatypeAsString(const datatype_t& attr);
    
  };

  /*!
   * Setup c_size_of_datatype array in such a way, that the data
   * is always aligned to a 2^align bit boundary. I hope this ensures cache
   * coherency, but someone with more experience in this stuff should have
   * a look at this.
   */
  static void alignDataFor(size_t align);

#ifdef _OPENMP
  /*!
   * Calls a function from Data which counts the number of doubles
   * that each datatype uses.
   */
  static int getNumOfDoubles(datatype_t datatype);
#endif

protected:
  /*!
   * A vector for access to attributes by index
   */
  vector<attribute_t> m_attr_by_index;

  /*!
   * A map for access to attributes by name
   */
  map<string, attribute_t> m_attr_by_name;

  /*!
   * Size of the whole data structure when allocated
   */
  size_t m_size;

public:
  /*!
   * Constructor
   */
  DataFormat();

  /*!
   * Copy constructor
   */
  DataFormat(const DataFormat &copy_dataFormat);

  /*!
   * Destructor
   */
  virtual ~DataFormat();

  /*!
   * Add an additional attribute
   * @param name Attribute's name
   * @param datatype Attribute's type
   * @param persistent Can this attribut be cleared using \a clear? Default = false
   * @param symbol Short name. Default = name
   */
  virtual attribute_t addAttribute
    (const string &name, datatype_t datatype, bool persistent = false, string symbol = "");

  /*!
   * Allocate and initialize a new data block for this structure (of size \a m_size)
   * @param alloc_sp Whether smartpointers should be initialized (default: true)
   */
  void *alloc(bool alloc_sp = true) const;

  /*!
   * Release a data block allocated by this \a DataFormat
   * @param data Pointer to the data block to release
   */
  void release(void *&data) const;

  /*!
   * Clear all field of the data block \a data, except those where persistent = true
   * @param data Pointer to the data block to clear
   */
  void clear(void *data) const;

  /*!
   * Clear all field of the data block \a data, including those where persistent = false
   * @param data Pointer to the data block to clear  
   */
  void clearAll(void *data) const;

  /*!
   * Return a string containing the names of all attributes. For passing
   * error messages to the user.
   */
  string toString() const;

  /*!
   * Check if the attribute with name \a name does already exist.
   * @param name Name of the attribute
   */
  bool attrExists(const string &name) const;

  /*!
   * Return an attribute by specifying a name
   * @param name Name of the attribute to return
   */
  inline const attribute_t &attrByName(const string &name) {
    if (!attrExists(name))
      throw gError
        ("DataFormat::attrByName",
         string("Attribute '") + name + "' does not exists. Possibilities: " + toString());

    return m_attr_by_name[name];
  }

  /*!
   * Return an attribute by specifying the index
   * @param i Index of the attribute to return
   */
  inline const attribute_t &attrByIndex(int i) const {
    return m_attr_by_index[i];
  }
    
  /*!
   * Return an attribute by specifying the index
   * @param i Index of the attribute to return
   */
  inline attribute_t &attrByIndex(int i) {
    return m_attr_by_index[i];
  }

  /*!
   * \a indexOf returns the index of the attribute of name \a name
   * ONLY IF the attribute exists and has datatype \a datatype
   * @param name Name of the attribute to find
   * @param datatype Type of the attribute to find
   */
  inline int indexOf(const string &name, datatype_t datatype) {
    if (attrExists(name)) {
      attribute_t attr;

      attr = attrByName(name);
      if (attr.datatype == datatype)
        return attr.index;
      else
        throw gError
          ("DataFormat::indexOf",
           "Attribute '" + name + "' has wrong data type.");
    } else
      throw gError
        ("DataFormat::indexOf",
         "Attribute '" + name + "' does not exist. Possibilities: " + toString());
  }

  /*!
   * Return the offset of the attribute with name \a name
   * @param name Name of the attribute
   */
  inline ptrdiff_t offsetByName(const string &name) {
    return m_attr_by_name[name].offset;
  }

  /*!
   * Return the of offset of the attribute with index \a i
   * @param i Index of the attribute
   */
  inline ptrdiff_t offsetByIndex(int i) {
    return m_attr_by_index[i].offset;
  }

  /*!
   * Return the pointer to the attribute with name \a name in data block \a data
   * @param name Name of the attribute
   * @param data Start of the data block where to find the actual data
   */
  inline void *ptrByName(const string &name, void *data) {
    return (void*) ((size_t) data + offsetByName(name));
  }

  /*!
   * Return the pointer to the attribute with index \a i in data block \a data
   * @param i Index of the attribute
   * @param data Start of the data block where to find the actual data
   */
  inline void *ptrByIndex(int i, void *data) {
    return (void*) ((size_t) data + offsetByIndex(i));
  }

  /*!
   * Return the pointer to the attribute \a attr in data block \a data
   * @param attr Attribute information
   * @param data Start of the data block where to find the actual data
   */
  inline static void *ptrByAttr(const attribute_t &attr, void *data) {
    return (void*) ((size_t) data + attr.offset);
  }

  /*!
   * Return the pointer to the attribute with offset \a offset in data block \a data
   * @param offset Offset of the attribute
   * @param data Start of the data block where to find the actual data
   */
  inline static void *ptrByOffset(ptrdiff_t offset, void *data) {
    return (void*) ((size_t) data + offset);
  }

  /*!
   * Return the number of attributes managed by this \a DataFormat
   */
  inline size_t rows() {
    return m_attr_by_index.size();
  }

  /*!
   * Return the size (in bytes) of structures allocated by this \a DataFormat
   */
  inline size_t size() {
    return m_size;
  }

  /*!
   * Return a new \a Data object
   */
  data_sp newData();
};


/*!
 * The \a Data object manages the actual data. It therefore stores a pointer
 * to the appropriate \a DataFormat and to the actual block of data.
 */
class Data
{
protected:
  /*!
   * Pointer to the format description
   */
  DataFormat *m_format;

  /*!
   * Pointer to the actual data
   */
  void *m_data;

public:
  /*!
   * Default constructor. 
   */
  Data();

  /*!
   * Create a \a Data for format description \a format
   * @param format Format description for this \a Data
   */
  Data(DataFormat *format);

  /*!
   * Copy constructor. Copy over all data stored in \a copy_data
   * @param copy_data Copy from this object
   */
  Data(const Data &copy_data);

  /*!
   * Destructor
   */
  virtual ~Data();

  /*!
   * Assigment operator
   * @param copy_data Copy from this object
   */
  Data &operator=(const Data &copy_data);

#ifdef _OPENMP
  /*!
   * Define how many doubles we need for the ValCalculators,
   * based on the m_datatype they use. The function is called by Symbol.
   */
  static int countDoubles(DataFormat::datatype_t dtype);
#endif

  /*!
   * Add new attribute to the Data field 
   */
  virtual DataFormat::attribute_t addAttribute
  (const string &name, DataFormat::datatype_t datatype, bool persistent = false, string symbol = "");

  /*!
   * Check if the attribute with name \a name does already exist.
   * @param name Name of the attribute
   */
  bool attrExists(const string &name) const;
 
  /*!
   * Return the offset of the attribute with name \a name
   * @param name Name of the attribute
   */
  inline ptrdiff_t offsetByName(const string &name) {
    return m_format->offsetByName(name);
  }

  /*!
   * Return the of offset of the attribute with index \a i
   * @param i Index of the attribute
   */
  inline ptrdiff_t offsetByIndex(int i) {
    return m_format->offsetByIndex(i);
  }

  /*!
   * Return a string containing the names of all attributes. For passing
   * error messages to the user.
   */
  string toStr() const;

  /*!
   * \a indexOf returns the index of the attribute of name \a name
   * ONLY IF the attribute exists and has datatype \a datatype
   * @param name Name of the attribute to find
   * @param datatype Type of the attribute to find
   */
  inline int indexOf(const string &name, DataFormat::datatype_t datatype) {
    return m_format->indexOf(name, datatype);
  }

  /*!
   * Return an attribute by specifying a name
   * @param name Name of the attribute to return
   */
  inline const DataFormat::attribute_t &attrByName(const string &name) {
    return m_format->attrByName(name);
  }

  /*!
   * Return an attribute by specifying the index
   * @param i Index of the attribute to return
   */
  inline const DataFormat::attribute_t &attrByIndex(int i) const {
    return m_format->attrByIndex(i);
  }

  /*!
   * Is this data allocated?
   */
  bool isNull() const {
    return m_data == NULL;
  }

  /*!
   * Make this Data persistent
   */
  inline void protect(size_t i) {
    m_format->attrByIndex(i).persistent = true;
  }
  
  /*!
   * Remove the persistency of this Data
   */
  inline void unprotect(size_t i) {
    m_format->attrByIndex(i).persistent = false;
  }

  /*!
   * Return the number of attributes managed by this \a DataFormat
   */
  inline size_t rows() {
    return m_format->rows();
  }
  
  /*!
   * Set the format description to \a format and allocate
   * and initialize a new memory block
   * @param format Format description to use from now on
   */
  void setFormatAndAlloc(DataFormat *format);

  /*!
   * Analog to setFormatAndAlloc
   */
//   void setFormatAndCopy(Data *defaultData);

  /*!
   * Allocate the memory block
   */
  inline void alloc() {
    m_data = m_format -> alloc();
  }

  /*!
   * Release and allocate the memory block
   */
  void reAlloc() {
    m_format -> release(m_data);
    alloc();
  }

  /*!
   * Release the memory block
   */
  void release() {
    m_format -> release(m_data);
    m_data = NULL;
  }

  /*!
   * Clear all attributes, except those with persistent = true
   */
  void clear() {
    m_format->clear(m_data);
  }

  /*!
   * Clear all attributes, including those with persistent = true
   */
  void clearAll() {
    m_format->clearAll(m_data);
  }
   
  /*!
   * Return the integer with name \a name
   * @param name Name of the attribute to return
   */
  inline int &intByName(const string &name) const {
    return *((int*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the double with name \a name
   * @param name Name of the attribute to return
   */
  inline double &doubleByName(const string &name) const {
    return *((double*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the integer point with name \a name
   * @param name Name of the attribute to return
   */
  inline int_point_t &intPointByName(const string &name) const {
    return *((int_point_t*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the point with name \a name
   * @param name Name of the attribute to return
   */
  inline point_t &pointByName(const string &name) const {
    return *((point_t*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the tensor with name \a name
   * @param name Name of the attribute to return
   */
  inline tensor_t &tensorByName(const string &name) const {
    return *((tensor_t*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the string with name \a name
   * @param name Name of the attribute to return
   */ 
  inline string &stringByName(const string &name) const {
    return *((string*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the double vector with name \a name
   * @param name Name of the attribute to return
   */
  inline vector_double_sp &vectorDoubleByName(const string &name) const {
    return *((vector_double_sp*) m_format->ptrByName(name, m_data));
  }    

  /*!
   * Return the point vector with name \a name
   * @param name Name of the attribute to return
   */
  inline vector_point_sp &vectorPointByName(const string &name) const {
    return *((vector_point_sp*) m_format->ptrByName(name, m_data));
  }

  /*!
   * Return the tensor vector with name \a name
   * @param name Name of the attribute to return
   */
  inline vector_tensor_sp &vectorTensorByName(const string &name) const {
    return *((vector_tensor_sp*) m_format->ptrByName(name, m_data));
  }

#ifdef WITH_ARRAY_TYPES


  /*!
   * Return the array2d_double matrix with name \a name
   * @param name Name of the attribute to return
   */
  inline array2d_double_sp &array2dDoubleByName(const string &name) const {
	  return *((array2d_double_sp*) m_format->ptrByName(name, m_data));
  }
#endif

  /*!
   * Return a stringified version of the attribute with name \a name
   * @param name Name of the attribute
   */
  virtual string toStringByName(const string &name, const std::string &eol = "\n") const;

  /*!
   * Return a mathematicaized version of the attribute with name \a name
   * @param name Name of the attribute
   */
  virtual string toMathematicaByName(const string &name, const std::string &eol = "\n") const;

  /*!
   * Set the attribute with name \a name from the string information in \a value
   * @param name Name of the attribute to set
   * @param value String to convert to the appropriate format
   */
  virtual void fromStringByName(const string &name, const string &value);

  /*!
   * Return the integer with index \a i
   * @param i Index of the attribute to return
   */
  inline int &intByIndex(int i) const {
    return *((int*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the double with index \a i
   * @param i Index of the attribute to return
   */
  inline double &doubleByIndex(int i) const {
    return *((double*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the integer point with index \a i
   * @param i Index of the attribute to return
   */
  inline int_point_t &intPointByIndex(int i) const {
    return *((int_point_t*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the point with index \a i
   * @param i Index of the attribute to return
   */
  inline point_t &pointByIndex(int i) const {
    return *((point_t*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the tensor with index \a i
   * @param i Index of the attribute to return
   */
  inline tensor_t &tensorByIndex(int i) const {
    return *((tensor_t*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the string with index \a i
   * @param i Index of the attribute to return
   */
  inline string &stringByIndex(int i) const {
    return *((string*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the double vector with index \a i
   * @param i Index of the attribute to return
   */
  inline vector_double_sp &vectorDoubleByIndex(int i) const {
    return *((vector_double_sp*) m_format->ptrByIndex(i, m_data));
  }    

  /*!
   * Return the integer vector with index \a i
   * @param i Index of the attribute to return
   */
  inline vector_int_sp &vectorIntByIndex(int i) const {
    return *((vector_int_sp*) m_format->ptrByIndex(i, m_data));
  }    

  /*!
   * Return the point vector with index \a i
   * @param i Index of the attribute to return
   */
  inline vector_point_sp &vectorPointByIndex(int i) const {
    return *((vector_point_sp*) m_format->ptrByIndex(i, m_data));
  }

  /*!
   * Return the tensor vector with index \a i
   * @param i Index of the attribute to return
   */
  inline vector_tensor_sp &vectorTensorByIndex(int i) const {
    return *((vector_tensor_sp*) m_format->ptrByIndex(i, m_data));
  }

#ifdef WITH_ARRAY_TYPES
  
  /*!
   * Return the array2d_double vector with index \a i
   * @param i Index of the attribute to return
   */
  inline array2d_double_sp &array2dDoubleByIndex(int i) const {
	  return *((array2d_double_sp*) m_format->ptrByIndex(i, m_data));
  }
#endif
    
  /*!
   * Return a stringified version of the attribute with index \a i
   * @param i Index of the attribute
   */
  virtual string toStringByIndex(int i, const std::string &eol = "\n") const;

  /*!
   * Return a mathematicaized version of the attribute with index \a i
   * @param i Index of the attribute
   */
  virtual string toMathematicaByIndex(int i, const std::string &eol = "\n") const;

  /*!
   * Set the attribute with index \a i from the string information in \a value
   * @param i Index of the attribute to set
   * @param value String to convert to the appropriate format
   */
  virtual void fromStringByIndex(int i, const string &value);

  /*!
   * Return the integer with offset \a o
   * @param o Offset of the attribute to return
   */
  inline int &intByOffset(int o) const {
    return *((int*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the double with offset \a o
   * @param o Offset of the attribute to return
   */
  inline double &doubleByOffset(int o) const {
    return *((double*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the integer point with offset \a o
   * @param o Offset of the attribute to return
   */
  inline int_point_t &intPointByOffset(int o) const {
    return *((int_point_t*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the point with offset \a o
   * @param o Offset of the attribute to return
   */
  inline point_t &pointByOffset(int o) const {
    return *((point_t*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the tensor with offset \a o
   * @param o Offset of the attribute to return
   */
  inline tensor_t &tensorByOffset(int o) const {
    return *((tensor_t*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the string with offset \a o
   * @param o Offset of the attribute to return
   */
  inline string &stringByOffset(int o) const {
    return *((string*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the double vector with offset \a o
   * @param o Offset of the attribute to return
   */
  inline vector_double_sp &vectorDoubleByOffset(int o) const {
    return *((vector_double_sp*) DataFormat::ptrByOffset(o, m_data));
  }    

  /*!
   * Return the point vector with offset \a o
   * @param o Offset of the attribute to return
   */
  inline vector_point_sp &vectorPointByOffset(int o) const {
    return *((vector_point_sp*) DataFormat::ptrByOffset(o, m_data));
  }

  /*!
   * Return the tensor vector with offset \a o
   * @param o Offset of the attribute to return
   */
  inline vector_tensor_sp &vectorTensorByOffset(int o) const {
    return *((vector_tensor_sp*) DataFormat::ptrByOffset(o, m_data));
  }

#ifdef WITH_ARRAY_TYPES

  /*!
   * Return the array2d_double vector with offset \a o
   * @param o Offset of the attribute to return
   */
  inline array2d_double_sp &array2dDoubleByOffset(int o) const {
	  return *((array2d_double_sp*) DataFormat::ptrByOffset(o, m_data));
  }
#endif

  /*!
   * Return a list of all attributes. For passing information to the user.
   */
  virtual string toString() const;

  /*!
   * Return the pointer to the actual data
   */
  inline void* data() {
    return m_data;
  }
  
  /*!
   * Return the format description
   */
  inline DataFormat *format() {
    return m_format;
  }
};



inline data_sp DataFormat::newData()
{
  data_sp h(new Data(this));

  return h;
}

#endif
