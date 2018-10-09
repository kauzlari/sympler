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



#ifndef __PROPERTYLIST_H
#define __PROPERTYLIST_H

#include <climits>
#include <map>
#include <string>
#include <vector>

#include <libxml/tree.h>

#include "geometric_primitives.h"
#include "general.h"

using namespace std;

/* PropertyList manages the properties that can be given in an XML file */

//---- Classes ----

class FunctionFixed;
class FunctionPair;
class HelpNode;

class PropertyListConstraint;

/*!
 * The property ties properties of an XML file to real variables in the memory.
 * It automatically reads an XML file and sets the appropriate variable.
 */
class PropertyList
{
public:
  /*!
   * Data type of a \a PropertyList entry
   */
  enum datatype_t {
      INT = 0, DOUBLE = 1, STRING = 2, BOOLEAN = 3, FUNCTIONPAIR = 4, FUNCTIONFIXED = 5, POINT = 6, EO_DATATYPE = 7
  };


  /*!
   * List of string identifiers corresponding to the data type
   */
  static string datatype_names[EO_DATATYPE];
    
  /*!
   * An entry in the property list
   */
  struct property_t {
    /*!
     * Name of this property
     */
    string name;

    /*!
     * Index of this property within the \a PropertyList
     */
    size_t index;

    /*!
     * Type of this property
     */
    datatype_t datatype;
    
    /*!
     * Pointer to where the data read from the XML file should be written to
     */
    void *ptr;

    /*!
     * Any constraints on the possible values of this property?
     */
    PropertyListConstraint *constraint;

    /*!
     * Description string for automatic help file generation
     */
    string description;

    /*!
     * Set the help string for this \a property_t
     * @param description New help string
     */
    void setDescription(string newDescription) {
    	description = newDescription;
    }

    /*!
     * Read/Write an integer value
     */
    int &asInt() {
          return *((int*) ptr);
      }



    /*!
     * Read an integer value
     */
    const int &asInt() const {
      return *((int*) ptr);
      }

    
        
    /*!
     * Read/Write a double value
     */
    double &asDouble() {
       return *((double*) ptr);
      }


   
    /*!
     * Read a double value
     */
    const double &asDouble() const {
      return *((double*) ptr);
      }


        
    /*!
     * Read/Write a string value
     */
    string &asString() {
      return *((string*) ptr);
    }

    /*!
     * Read a string value
     */
    const string &asString() const {
      return *((string*) ptr);
    }

    /*!
     * Read/Write a bool value
     */        
    bool &asBool() {
      return *((bool*) ptr);
    }

    /*!
     * Read a bool value
     */        
    const bool &asBool() const {
      return *((bool*) ptr);
    }        

    /*!
     * Read/Write a \a FunctionPair
     */
    FunctionPair &asFunctionPair() {
	    return *((FunctionPair*) ptr);
    }

    /*!
     * Read a \a FunctionPair
     */
    const FunctionPair &asFunctionPair() const {
	    return *((FunctionPair*) ptr);
    }
    
    /*!
     * Read/Write a \a FunctionFixed
     */
    FunctionFixed &asFunctionFixed() {
	    return *((FunctionFixed*) ptr);
    }

    /*!
     * Read a \a FunctionFixed
     */
    const FunctionFixed &asFunctionFixed() const {
	    return *((FunctionFixed*) ptr);
    }

    /*!
     * Read/Write a point
     */
    point_t &asPoint() {
      return *((point_t*) ptr);
    }

    /*!
     * Read a point
     */
    const point_t &asPoint() const {
      return *((point_t*) ptr);
    }

    /*!
     * Return a stringified version of the data (for XML file generation)
     */
    string toString() const;
/*    {
      stringstream str;

      switch (datatype) {
      case INT:
          str << asInt();
            break;
      case DOUBLE:
          str << asDouble();
          break;
      case STRING:
          str << "'" << asString() << "'";
          break;
      case BOOLEAN:
          str << (asBool() ? "yes" : "no");
          break;
	    case FUNCTIONPAIR:
          str << "'" << asFunctionPair().expression() << "'";
          if(asFunctionPair().isNull()) str << "UNDEFINED";
          break;
	    case FUNCTIONFIXED:
          str << "'" << asFunctionFixed().expression() << "'";
          if(asFunctionFixed().isNull()) str << "UNDEFINED";
          break;
      case POINT:
          str << asPoint();
          break;
      default:
          str << "undefined.";
          break;
      }

      return str.str();
    }*/
  }; // end of definition struct property_t...
    
protected:
  /*!
   * The name of this object (initially identical to the class name)
   */
  string m_name;

  /*!
   * The name of the class this property belongs to.
   */
  string m_class_name;

  /*!
   * Description of the class.
   */
  string m_description;
         
  /*!
   * A vector for access by index
   */
  vector<property_t> m_prop_by_index;
    
  /*!
   * A map for access by name
   * FIXME: together with \a m_prop_by_index this generates redundant
   * and to be synchronised information! Remove redundancy!
   */
  map<string, property_t> m_prop_by_name;

  /*!
   * Allow properties that are undefined?
   */
  bool m_allow_unknown;

  /*!
   * Name/value pairs of all undefined properties for delayed evaluation
   */
  map<string, string> m_unknown_props;
    
public:
  /*!
   * Constructor
   */
  PropertyList();

  /*!
   * Destructor
   */
  virtual ~PropertyList();
    
  /*!
   * Add a new property
   * @param name Name of the property
   * @param datatype Type of the property
   * @param ptr Pointer to where the property should be written to in memory
   * @param constraint Constraints for this property
   * @param description Description of this property
   */
  virtual property_t addProperty
    (const string &name, datatype_t datatype,
     void *ptr, PropertyListConstraint *constraint = NULL, string description = "");
    
  /*!
   * Return property of name \a name
   * @param name Name of the property
   */
  property_t &propByName(const string &name) {
    return m_prop_by_name[name];
  }

  /*!
   * Return property of name \a name. Const version.
   * @param name Name of the property
   */
  const property_t &propByName(const string &name) const {
  	// map-operator [] can not return a const reference, so we do not use it in
  	// this const version of \a propByName
    return m_prop_by_name.find(name) -> second;
  }

  /*!
   * Return property with index \a i
   * @param i Index of the property
   */
  property_t &propByIndex(int i) {
    return m_prop_by_index[i];
  }
    
  /*!
   * Does the property \a name exist?
   * @param name Name of the property to look for
   */
  bool exists(const string &name) const {
    return m_prop_by_name.find(name) != m_prop_by_name.end();
  }

  /*!
   * Return the number of entries in this \a PropertyList
   */
  int size() const {
    return m_prop_by_index.size();
  }

  /*!
   * Set the name of the object this property list belongs to
   * @param name The new name
   */
  void setName(string name) {
    m_name = name;
  }

  /*!
   * Return the name of the object this property list belongs to
   */
  const string &name() const {
    return m_name;
  }

  /*!
   * Set the name of the class this property list belongs to
   * @param class_name The new class name
   */
  void setClassName(string class_name) {
    m_class_name = class_name;
    m_name = class_name;
  }

  /*!
   * Return the name of the class this property list belongs to
   */
  const string &className() const {
    return m_class_name;
  }

  /*!
   * Set the help string for this \a PropertyList
   * @param description New help string
   */
  void setDescription(string description) {
    m_description = description;
  }

  /*!
   * Return the help string of this \a PropertyList
   */
  const string &description() const {
    return m_description;
  }

  /*!
   * Set the help string for specified \a property_t
   * @param[in] name Name of \a property_t in \a m_prop_by_name and 
   * \a m_prop_by_index (FIXME: redundancy!)
   * @param[in] description New help string
   */
  void setPropDescription(string name, string description);

  /*!
   * Call this to prevent unknown properties to raise an exception
   */
  void allowUnknown() {
    m_allow_unknown = true;
  }

  /*!
   * Return the name/value pairs of unknown propertie
   */
  const map<string, string> &unknown() const {
    return m_unknown_props;
  }

  /*!
   * Create a help tree
   */
  void help(HelpNode *node) const;
    
  /*!
   * Read XML file
   */
  virtual void fromXML(const xmlNode *xmln);

  /*!
   * Write XML header
   */
  virtual ostream &toXML_begin(ostream &s, int shift = 0);

  /*!
   * Write XML footer
   */
  virtual ostream &toXML_end(ostream &s, int shift = 0);

  /*!
   * Produce an error message if an unknown attribute was found while parsing
   * @param s String forming the first part of the error message
   */
  void throwListIfUnknown(string s);
};



/*!
 * Class for constraint checking within the \a PropertyList
 */
class PropertyListConstraint
{
 public:
  /*!
   * Constructor
   */
  PropertyListConstraint();

  /*!
   * Destructor
   */
  virtual ~PropertyListConstraint();
    
  /*!
   * Check if property \a property matches this constraint
   * @param property Property to check
   */
  virtual bool check(const PropertyList::property_t &property) = 0;
};


/* Different generic constraints */

#define PLC_COMPARISON(name, datatype, access_function, operator) \
    class name: public PropertyListConstraint \
    { \
    protected: \
        datatype m_than; \
    \
    public: \
        name(datatype than): m_than(than) { } \
        virtual bool check(const PropertyList::property_t &property) { \
            return property.access_function() operator m_than; \
        } \
    };
    
#define PLC_BINARY(name, operator) \
    class name: public PropertyListConstraint \
    { \
    protected: \
        PropertyListConstraint *m_first, *m_second; \
        \
    public: \
        name(PropertyListConstraint *first, PropertyListConstraint *second)\
            : m_first(first), m_second(second) { } \
        ~name() { \
            delete m_first; delete m_second; \
        } \
        virtual bool check(const PropertyList::property_t &property) { \
            return m_first->check(property) operator m_second->check(property); \
        } \
    };


PLC_COMPARISON(PLCIntEqual, int, asInt, ==)
//PLC_COMPARISON(PLCDoubleEqual, double, asDouble, ==)
PLC_COMPARISON(PLCStringEqual, string, asString, ==)

PLC_COMPARISON(PLCIntGreater, int, asInt, >)
PLC_COMPARISON(PLCDoubleGreater, double, asDouble, >)
PLC_COMPARISON(PLCIntGreaterEqual, int, asInt, >=)
PLC_COMPARISON(PLCDoubleGreaterEqual, double, asDouble, >=)

PLC_COMPARISON(PLCIntLess, int, asInt, <)
PLC_COMPARISON(PLCDoubleLess, double, asDouble, <)
PLC_COMPARISON(PLCIntLessEqual, int, asInt, <=)
PLC_COMPARISON(PLCDoubleLessEqual, double, asDouble, <=)

PLC_BINARY(PLCAnd, &&)
PLC_BINARY(PLCOr, ||)

#define PLCIntRange(a, b) PLCAnd(new PLCIntGreaterEqual(a), new PLCIntLessEqual(b))
#define PLCDoubleRange(a, b) PLCAnd(new PLCDoubleGreaterEqual(a), new PLCDoubleLessEqual(b))

#define PLCStringList2(a, b) PLCOr(new PLCStringEqual(a), new PLCStringEqual(b))
#define PLCStringList3(a, b, c) PLCOr(new PLCStringEqual(a), new PLCStringList2(b, c))
#define PLCStringList4(a, b, c, d) PLCOr(new PLCStringEqual(a), new PLCStringList3(b, c, d))
#define PLCStringList5(a, b, c, d, e) PLCOr(new PLCStringEqual(a), new PLCStringList4(b, c, d, e))

#define PLCIntList2(a, b) PLCOr(new PLCIntEqual(a), new PLCIntEqual(b))
#define PLCIntList3(a, b, c) PLCOr(new PLCIntEqual(a), new PLCIntList2(b, c))

#endif
