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



#include <algorithm>

#include "data_format.h"

//---- DataFormat ----

size_t DataFormat::c_size_of_datatype[EO_DATATYPE] = {
    sizeof(int), sizeof(double), sizeof(int_point_t), sizeof(point_t), sizeof(tensor_t), sizeof(string),
    sizeof(vector_int_sp), sizeof(vector_double_sp), sizeof(vector_point_sp), sizeof(vector_tensor_sp) 
#ifdef WITH_ARRAY_TYPES
    ,  sizeof(array2d_double_sp)
#endif
};

//---- Constructors/Destructor ----

DataFormat::DataFormat(): m_size(0)
{
}


DataFormat::DataFormat(const DataFormat &copy_dataFormat): m_size(0)
{
  m_attr_by_index = vector<attribute_t>(copy_dataFormat.m_attr_by_index);
  m_attr_by_name = map<string, attribute_t>(copy_dataFormat.m_attr_by_name);
  m_size = copy_dataFormat.m_size;
}


DataFormat::~DataFormat()
{
}


void DataFormat::alignDataFor(size_t align)
{

  for (int i = 0; i < EO_DATATYPE; ++i) {
    c_size_of_datatype[i] = (((c_size_of_datatype[i] - 1) >> align) + 1) << align;
  }
}


//---- Methods ----

DataFormat::attribute_t DataFormat::addAttribute
  (const string &name, datatype_t datatype, bool persistent, string symbol)
{
  map<string, attribute_t>::iterator it = m_attr_by_name.find(name);
  if(it == m_attr_by_name.end())
  {
    string s = symbol;
  
    if (s == "")
      s = name;
    attribute_t attr = {
      name, m_attr_by_index.size(), m_size, datatype, persistent, s
    };
      
    m_size += c_size_of_datatype[datatype];
    
    m_attr_by_index.push_back(attr);
    m_attr_by_name[name] = attr;
//     if(s == "D")
//     MSG_DEBUG("DataFormat::addAttribute", "symbol = " << s << ", datatype = " << attr.datatype);   
    return attr;
  }
  else
  {
    MSG_DEBUG("DataFormat::addAttribute", "Will do nothing because attribute already existing: " << endl << "name = " << name << endl << "symbol = " << symbol);
    if(it->second.datatype != datatype) 
	    throw gError("DataFormat::addAttribute", "Type mismatch during request of attribute. Requested type: " + ObjToString(datatype) + ". Registered type: " + ObjToString(it->second.datatype) + ". With: INT = 0, DOUBLE = 1, INT_POINT = 2, POINT = 3, TENSOR = 4, STRING = 5, VECTOR_INT = 6, VECTOR_DOUBLE = 7, VECTOR_POINT = 8, VECTOR_TENSOR = 9, "+
#ifdef WITH_ARRAY_TYPES
	     "MArray2D = 10,  EO_DATATYPE = 11"
#else
	    "EO_DATATYPE = 10"
#endif
	    );
    return it->second;
  }
}

string DataFormat::attribute_t::datatypeAsString() const {

  /*static*/ DataFormat::attribute_t::datatypeAsString(this->datatype);

}

/*static*/ string DataFormat::attribute_t::datatypeAsString(const datatype_t& datatype) {
  switch (datatype) {
  case DataFormat::INT:
    return string("INT");
  case DataFormat::DOUBLE:
    return string("DOUBLE");
  case DataFormat::POINT:
    return string("POINT");
  case DataFormat::TENSOR:
    return string("TENSOR");
  case DataFormat::STRING:
    return string("STRING");
  case DataFormat::VECTOR_DOUBLE:
    return string("VECTOR_DOUBLE");
  case DataFormat::VECTOR_INT:
    return string("VECTOR_INT");
  case DataFormat::VECTOR_POINT:
    return string("VECTOR_POINT");
#ifdef WITH_ARRAY_TYPES
  case DataFormat::MArray2D:
    return string("MArray2D");
#endif
  default:
    throw gError
      ("DataFormat::attribute_t::datatypeAsString",
       "Don't know how to handle data type: "
       + ObjToString(datatype));
  }
}

struct alloc_smart_pointer: public unary_function<DataFormat::attribute_t, void> 
{
    void *m_data;
    alloc_smart_pointer(void *data): m_data(data) { }
    void operator()(const DataFormat::attribute_t &attr) {
      if (attr.datatype == DataFormat::VECTOR_DOUBLE) {
        ((vector_double_sp*) ((size_t) m_data + attr.offset))->alloc();
      } 
      else if (attr.datatype == DataFormat::VECTOR_INT) {
        ((vector_int_sp*) ((size_t) m_data + attr.offset))->alloc();
      } 
      else if (attr.datatype == DataFormat::VECTOR_POINT) {
            ((vector_point_sp*) ((size_t) m_data + attr.offset))->alloc();
        } 
        else if (attr.datatype == DataFormat::VECTOR_TENSOR) {
            ((vector_tensor_sp*) ((size_t) m_data + attr.offset))->alloc();
        }
#ifdef WITH_ARRAY_TYPES      

        else if (attr.datatype == DataFormat::MArray2D) {
        	reinterpret_cast<array2d_double_sp*>
        		((size_t) m_data + attr.offset)->alloc(); // = new array2d_double_sp();
        }
#endif
    }
};


void *DataFormat::alloc(bool alloc_sp) const
{
    void *data = NULL;

    if (m_size) {
        data = malloc(m_size);
        memset((void*) data, (char) 0, m_size);

	    if(alloc_sp) {
           for_each(m_attr_by_index.begin(), m_attr_by_index.end(), alloc_smart_pointer(data));
	    }
    }

    return data;
}


struct release_smart_pointer: public unary_function<DataFormat::attribute_t, void>
{
    void *m_data;
    release_smart_pointer(void *data): m_data(data) { }
    void operator()(const DataFormat::attribute_t &attr) {
      if (attr.datatype == DataFormat::VECTOR_DOUBLE) {
        ((vector_double_sp*) ((size_t) m_data + attr.offset))->release();
      } 
      else if (attr.datatype == DataFormat::VECTOR_INT) {
        ((vector_int_sp*) ((size_t) m_data + attr.offset))->release();
      } 
      else if (attr.datatype == DataFormat::VECTOR_POINT) {
            ((vector_point_sp*) ((size_t) m_data + attr.offset))->release();
      } 
      else if (attr.datatype == DataFormat::VECTOR_TENSOR) {
            ((vector_tensor_sp*) ((size_t) m_data + attr.offset))->release();
      }
#ifdef WITH_ARRAY_TYPES      

      else if (attr.datatype == DataFormat::MArray2D) {
            ((array2d_double_sp*) ((size_t) m_data + attr.offset))->release();
      } 
#endif
    }
};

void DataFormat::release(void *&data) const
{
    if (data) {

        for_each(m_attr_by_index.begin(), m_attr_by_index.end(), release_smart_pointer(data));

        free(data);
        data = NULL;
    }
}

void DataFormat::clear(void *data) const
{
    vector<attribute_t>::const_iterator attr_end = m_attr_by_index.end();
    for (vector<attribute_t>::const_iterator i = m_attr_by_index.begin();
         i != attr_end; i++) {
        if (!i->persistent) {
#define 	DATAFORMATCLEAR_RELEASE(attr,m_data,finaltype) \
				((finaltype*)((size_t) m_data+attr->offset))->release()
        	DATAFORMAT_CONTAINER_SWITCH(i,data,DATAFORMATCLEAR_RELEASE);
            memset((void*) ((size_t) data + i->offset), (char) 0, c_size_of_datatype[i->datatype]);
        }
    }

}


void DataFormat::clearAll(void *data) const
{
    vector<attribute_t>::const_iterator attr_end = m_attr_by_index.end();
    for (vector<attribute_t>::const_iterator i = m_attr_by_index.begin();
         i != attr_end; i++) {
         	DATAFORMAT_CONTAINER_SWITCH(i,data,DATAFORMATCLEAR_RELEASE);
            memset((void*) ((size_t) data + i->offset), (char) 0, c_size_of_datatype[i->datatype]);
    }

}



bool DataFormat::attrExists(const string &name) const 
{
  return m_attr_by_name.find(name) != m_attr_by_name.end();
}
 
string DataFormat::toString() const
{
  string s("");
  map<string, attribute_t>::const_iterator attr_end = m_attr_by_name.end();

  for (map<string, attribute_t>::const_iterator i = m_attr_by_name.begin();
       i != attr_end; i++) {
    if (s != "")
      s += ", '" + i->first + "'";
    else
      s += "'" + i->first + "'";
  }

  return s;
}


#ifdef _OPENMP
int DataFormat::getNumOfDoubles(datatype_t datatype)
{
  return Data::countDoubles(datatype);
}
#endif



//---- Data ----

//---- Constructors/Destructor ----

Data::Data(): m_format(NULL), m_data(NULL)
{
}


Data::Data(DataFormat *format): m_format(format)
{
    m_data = m_format->alloc();
}


Data::Data(const Data &copy_data): m_format(NULL), m_data(NULL)
{

    if (copy_data.m_format) {
    	// Copy pointer to data format, DataFormat instance is shared
        m_format = copy_data.m_format;

        // Create new data storage space and copy old data,
        // but do not assign smartpointers
        m_data = m_format->alloc(false);
        
        if(m_format->size()) {
	  memcpy((void*) m_data, (void*) copy_data.m_data, m_format->size());
	  
	  // Iterate over all attributes and deep copy smartpointers where necessary 
	  size_t iend = (m_format->rows());
	  for(size_t i = 0; i < iend; i++) {
	    DataFormat::attribute_t a = m_format->attrByIndex(i);
	    
	    // incRefCount: dirty hack to make memcpy possible
#define	DATAFORMATCOPY_DEEP(a,m_data,finaltype)				\
            {reinterpret_cast<finaltype*>((size_t) m_data + a->offset)->incRefCount(); \
	      *reinterpret_cast<finaltype*>((size_t) m_data + a->offset) = \
		reinterpret_cast<finaltype*>((size_t) m_data + a->offset)->deepCopy();} NOOP 
	    DATAFORMAT_CONTAINER_SWITCH(&a, m_data, DATAFORMATCOPY_DEEP);
	  } // for
	} // if
    } // if
}


Data::~Data()
{
    
    if (m_format)
        m_format->release(m_data);
}


Data &Data::operator=(const Data &copy_data)
{

    if (m_format != copy_data.m_format) {
        if (m_format)
            m_format->release(m_data);

        m_format = copy_data.m_format;

        if (m_format)
            m_data = m_format->alloc(false);
    }

    if (m_format) {
        memcpy((void*) m_data, (void*) copy_data.m_data, m_format->size());

        // Iterate over all attributes and deep copy smartpointers where necessary 
        size_t iend = (m_format->rows());
        for(size_t i = 0; i < iend; i++) {
        	DataFormat::attribute_t a = m_format->attrByIndex(i);
		DATAFORMAT_CONTAINER_SWITCH(&a, m_data, DATAFORMATCOPY_DEEP);
          ;} // for         
    } //if
        
    return *this;
}



//---- Methods ----

void Data::setFormatAndAlloc(DataFormat *format)
{

    if (m_format)
        m_format->release(m_data);

    m_format = format;
    alloc();

}


DataFormat::attribute_t Data::addAttribute
  (const string &name, DataFormat::datatype_t datatype, bool persistent, string symbol)
{
  size_t old_size = m_format->size();
  DataFormat::attribute_t tempAttr = m_format->addAttribute(name, datatype, persistent, symbol);
  size_t new_size = m_format->size();

  if (old_size != new_size) {
    void* oldData = m_data;
    m_data = malloc(new_size);
    memset((char*) m_data, (char) 0, new_size);
    memcpy(m_data, oldData, old_size);
    free(oldData);
    (alloc_smart_pointer(m_data ))(tempAttr);
  }
  return tempAttr;
}


bool Data::attrExists(const string &name) const 
{
  return m_format->attrExists(name);
}


string Data::toStr() const
{
  return m_format->toString();
}


struct accumulate_double: public unary_function<DataFormat::attribute_t, double>
{
  string m_str, m_before, m_after, m_sep;
  bool m_comma;
  accumulate_double(string before, string after, string sep)
  : m_str(""), m_before(before), m_after(after), m_sep(sep), m_comma(false) { }
  void operator()(double d) {
    char s[80];
        
    if (m_comma)
      m_str += m_sep;
    m_comma = true;

    sprintf(s, "%f", d);
    m_str += s;
  }
  inline string str() {
    return m_before + m_str + m_after;
  }    
};

struct accumulate_double_scientific: public unary_function<DataFormat::attribute_t, double>
{
  string m_str, m_before, m_after, m_sep;
  bool m_comma;
  accumulate_double_scientific(string before, string after, string sep)
  : m_str(""), m_before(before), m_after(after), m_sep(sep), m_comma(false) { }
  void operator()(double d) {
    char s[80];
        
    if (m_comma)
      m_str += m_sep;
    m_comma = true;

    sprintf(s, "%g", d);
    m_str += s;
  }
  inline string str() {
    return m_before + m_str + m_after;
  }    
};

struct accumulate_int: public unary_function<DataFormat::attribute_t, int>
{
  string m_str, m_before, m_after, m_sep;
  bool m_comma;
  accumulate_int(string before, string after, string sep)
  : m_str(""), m_before(before), m_after(after), m_sep(sep), m_comma(false) { }
  void operator()(int num) {
    char s[80];
        
    if (m_comma)
      m_str += m_sep;
    m_comma = true;

    sprintf(s, "%i", num);
    m_str += s;
  }
  inline string str() {
    return m_before + m_str + m_after;
  }    
};

struct accumulate_point: public unary_function<DataFormat::attribute_t, point_t>
{
  string m_str, m_before, m_after, m_sep;
  bool m_comma;
  accumulate_point(string before, string after, string sep)
    : m_str(""), m_before(before), m_after(after), m_sep(sep), m_comma(false) { }
  void operator()(const point_t &p) {
    if (m_comma)
      m_str += m_sep;
    m_comma = true;

    m_str += m_before;
    for (int i = 0; i < SPACE_DIMS; i++) {
      char s[80];
      sprintf(s, "%f", p[i]);
      m_str += s;
      if (i != SPACE_DIMS-1)
        m_str += m_sep;
    }
    m_str += m_after;
  }
  inline string str() {
    return m_before + m_str + m_after;
  }
};

struct accumulate_point_scientific: public unary_function<DataFormat::attribute_t, point_t>
{
  string m_str, m_before, m_after, m_sep;
  bool m_comma;
  accumulate_point_scientific(string before, string after, string sep)
    : m_str(""), m_before(before), m_after(after), m_sep(sep), m_comma(false) { }
  void operator()(const point_t &p) {
    if (m_comma)
      m_str += m_sep;
    m_comma = true;

    m_str += m_before;
    for (int i = 0; i < SPACE_DIMS; i++) {
      char s[80];
      sprintf(s, "%g", p[i]);
      m_str += s;
      if (i != SPACE_DIMS-1)
        m_str += m_sep;
    }
    m_str += m_after;
  }
  inline string str() {
    return m_before + m_str + m_after;
  }
};

string Data::toStringByName(const string &name, const std::string &eol) const
{
    int index = m_format->attrByName(name).index;

    return toStringByIndex(index, eol);
}

string Data::toMathematicaByName(const string &name, const std::string &eol) const
{
    int index = m_format->attrByName(name).index;

    return toMathematicaByIndex(index, eol);
}


void Data::fromStringByName(const string &name, const string &value)
{
    int index = m_format->attrByName(name).index;

    fromStringByIndex(index, value);
}


string Data::toStringByIndex(int i, const std::string &eol) const
{
  char s[80];
  stringstream strstr;
  vector<double> *vd;
  vector<int> *vi;
  vector<point_t> *vp;
#ifdef WITH_ARRAY_TYPES
  MArray2D *ma2d;
#endif
    
  switch (m_format->attrByIndex(i).datatype) {
  case DataFormat::INT:
      sprintf(s, "%i", *((int*) m_format->ptrByIndex(i, m_data)));
      return string(s);
  case DataFormat::DOUBLE:
      sprintf(s, "%g", *((double*) m_format->ptrByIndex(i, m_data)));
      return string(s);
  case DataFormat::POINT:
      strstr << *((point_t*) m_format->ptrByIndex(i, m_data));
      return strstr.str();
  case DataFormat::TENSOR:
      strstr << *((tensor_t*) m_format->ptrByIndex(i, m_data));
      return strstr.str();
  case DataFormat::STRING:
      return *((string*) m_format->ptrByIndex(i, m_data));
    case DataFormat::VECTOR_DOUBLE:
      vd = ((vector_double_sp*) m_format->ptrByIndex(i, m_data))->value();
      return for_each(vd->begin(), vd->end(), accumulate_double_scientific(" ", " ", " ")).str();
    case DataFormat::VECTOR_INT:
      vi = ((vector_int_sp*) m_format->ptrByIndex(i, m_data))->value();
      return for_each(vi->begin(), vi->end(), accumulate_int(" ", " ", " ")).str();
    case DataFormat::VECTOR_POINT:
      vp = ((vector_point_sp*) m_format->ptrByIndex(i, m_data))->value();
      return for_each(vp->begin(), vp->end(), accumulate_point_scientific(" ", " ", " ")).str();
#ifdef WITH_ARRAY_TYPES
    case DataFormat::MArray2D:
    	ma2d = ((array2d_double_sp*) m_format->ptrByIndex(i, m_data))->value();
    	return ma2d->toMTXString(eol);
#endif
  default:
    throw gError
      ("Data::toStringByIndex",
       "Don't know how to handle data type: "
       + ObjToString(m_format->attrByIndex(i).datatype));
  }
}


#ifdef _OPENMP
int Data::countDoubles(DataFormat::datatype_t dtype)
{
  switch (dtype) {
    case DataFormat::DOUBLE:
        return 1;
    case DataFormat::POINT:
        return SPACE_DIMS;
    case DataFormat::TENSOR:
        return SPACE_DIMS_SQUARED;

    default:
      throw gError("Data::countDoubles", "This datatype shouldn't be used. Please contact the programmer.");
  }
}
#endif


string Data::toMathematicaByIndex(int i, const std::string &eol) const
{
  char s[80];
  stringstream strstr;
  vector<double> *vd;
  vector<int> *vi;
  vector<point_t> *vp;
#ifdef WITH_ARRAY_TYPES
  MArray2D *ma2d;
#endif

  int d1, d2;
    
  switch (m_format->attrByIndex(i).datatype) {
  case DataFormat::INT:
      sprintf(s, "%i", *((int*) m_format->ptrByIndex(i, m_data)));
      return string(s);
  case DataFormat::DOUBLE:
      sprintf(s, "%f", *((double*) m_format->ptrByIndex(i, m_data)));
      return string(s);
  case DataFormat::TENSOR: {
      tensor_t &t = *((tensor_t*) m_format->ptrByIndex(i, m_data));
      strstr << "{";
      for (int a = 0; a < SPACE_DIMS; a++) {
        strstr << "{";
        for (int b = 0; b < SPACE_DIMS; b++) {
          strstr << t(a, b);
          if (b != SPACE_DIMS-1)
            strstr << ", ";
        }
        strstr << "}";
        if (a != SPACE_DIMS-1)
          strstr << ", ";
      }
      strstr << "}";
      return strstr.str();
  }
  case DataFormat::STRING:
      return *((string*) m_format->ptrByIndex(i, m_data));
    case DataFormat::VECTOR_DOUBLE:
      vd = ((vector_double_sp*) m_format->ptrByIndex(i, m_data))->value();
      return for_each(vd->begin(), vd->end(), accumulate_double("{", "}", ", ")).str();
    case DataFormat::VECTOR_INT:
      vi = ((vector_int_sp*) m_format->ptrByIndex(i, m_data))->value();
      return for_each(vi->begin(), vi->end(), accumulate_int("{", "}", ", ")).str();
    case DataFormat::VECTOR_POINT:
      vp = ((vector_point_sp*) m_format->ptrByIndex(i, m_data))->value();
      return for_each(vp->begin(), vp->end(), accumulate_point("{", "}", ", ")).str();            
#ifdef WITH_ARRAY_TYPES
    case DataFormat::MArray2D:
		ma2d = ((array2d_double_sp*) m_format->ptrByIndex(i, m_data))->value();
		d1 = ma2d->dim1();
		d2 = ma2d->dim2();
		strstr << "{";
		for (int k = 0; k < d1; ++k) {
			strstr << "{";
			for (int j = 0; j < d2; ++j) {
				strstr << ((*ma2d)[k])[j];
				if (j != d2-1)
					strstr << ", ";
			}
			strstr << "}";
			if (k != d1-1)
				strstr << ", ";
			strstr << eol;
		}
	    strstr << "}";
		return strstr.str();
#endif
  default:
    throw gError
      ("Data::toMathematicaByIndex",
       "Don't know how to handle data type: "
       + ObjToString(m_format->attrByIndex(i).datatype));
  }
}


void Data::fromStringByIndex(int i, const string &value)
{
  switch (m_format->attrByIndex(i).datatype) {
  case DataFormat::INT:
    intByIndex(i) = atoi(value.c_str());
    break;
  case DataFormat::DOUBLE:
    doubleByIndex(i) = atof(value.c_str());
    break;
  case DataFormat::STRING:
    stringByIndex(i) = value;
  case DataFormat::POINT: {
    point_t p;
    int pos;

    pos = value.find('(')+1;
    for (int j = 0; j < SPACE_DIMS; j++) {
      int endpos;

      if (j != SPACE_DIMS-1)
	endpos = value.find(',', pos);
      else
	endpos = value.find(')', pos);

      p[j] = atof(string(value, pos, endpos-pos).c_str());
      pos = endpos+1;

    }
    pointByIndex(i) = p;

  }
    break;
  case DataFormat::TENSOR: {

    tensor_t t;
    int pos;

    pos = value.find('(')+1;

    for (int k = 0; k < SPACE_DIMS; k++) {
      pos = value.find('(', pos)+1;

      for (int j = 0; j < SPACE_DIMS; j++) {
	int endpos;

	if (j != SPACE_DIMS-1)
	  endpos = value.find(',', pos);
	else
	  endpos = value.find(')', pos);

	t(k, j) = atof(string(value, pos, endpos-pos).c_str());
	pos = endpos+1;
      }

      if (k != SPACE_DIMS-1)
	pos = value.find(',', pos)+1;
      else
	pos = value.find(')', pos)+1;
    }
    tensorByIndex(i) = t;
  }
    break;
  default:
    throw gError
      ("Data::fromStringByIndex",
       "Unsupported data format.");
  }
}


/* For debugging purposes */
string Data::toString() const 
{
  stringstream strstr;

  for (size_t i = 0; i < m_format->rows(); i++)
    strstr << m_format->attrByIndex(i).name << ": " << toStringByIndex(i) << endl;

  return strstr.str();
}
