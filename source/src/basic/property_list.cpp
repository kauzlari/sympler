/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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



#include <climits>
#include <string>
#include "property_list.h"
#include "function_pair.h"
#include "function_fixed.h"

#include "help.h"

string PropertyList::datatype_names[EO_DATATYPE] = 
    {"Integer", "Double", "String", "Boolean", "FunctionPair", "FunctionFixed", "Point"};


/* ---- Helper functions ---- */

point_t string2point(const string &str)
{
  stringstream s(str);
  point_t p;
  char c;

  s >> skipws;

  s.get(c);
  if (c != '(')
    throw gError("string2point", "This is not a point: '" + str + "'");

  s >> skipws;
  s >> p.x;
  s >> skipws;
  
  s.get(c);
  if (c != ',')
    throw gError("string2point", "This is not a point: '" + str + "'");

  s >> skipws;
  s >> p.y;
  s >> skipws;

  s.get(c);
  if (c != ',')
    throw gError("string2point", "This is not a point: '" + str + "'");

  s >> skipws;
  s >> p.z;
  s >> skipws;

  s.get(c);
  if (c != ')')
    throw gError("string2point", "This is not a point: '" + str + "'");

  return p;
}


/*---- PropertyList ---- */
//---- Constructors/Destructor ----

PropertyList::PropertyList(): m_allow_unknown(false)
{
}


PropertyList::~PropertyList()
{
}



//---- Methods ----

PropertyList::property_t PropertyList::addProperty
  (const string &name, datatype_t datatype,
   void *ptr, PropertyListConstraint *constraint, string description)
{
  property_t prop = {
      name, m_prop_by_index.size(), datatype, ptr, constraint, description
  };

  if (exists(name)) {
    MSG_DEBUG("PropertyList::addProperty", "gError");
    throw gError
      ("PropertyList::addProperty",
       "class name = " + m_class_name + ", name = " + m_name + ". " +
       "Property '" + name + "' has been defined for the second time. " + 
       "Please contact the programmer.");
  }

  m_prop_by_index.push_back(prop);
  m_prop_by_name[name] = prop;

  return prop;
}


void PropertyList::help(HelpNode *node) const
{
  HelpNode *n;

  n = new HelpNode(node, m_name);
  new HelpNode(n, m_description);

  for (vector<property_t>::const_iterator i = m_prop_by_index.begin();
       i != m_prop_by_index.end(); i++) {
    HelpNode *m;

    m = new HelpNode(n, i->name + " (" + datatype_names[i->datatype] + ")");

    // changed on 2013-06-17 because transformation from "_" to " " was removed
    new HelpNode(m, i->description + " (default: " + i->toString() + ")");
//     new HelpNode(m, i->description + " (default:_" + i->toString() + ")");
  }
}

void PropertyList::setPropDescription(string name, string description) {

 	propByName(name).setDescription(description);

	// FIXME: the following is necessary because we store the properties
	// redundantly in \a m_prop_by_name and \a m_prop_by_index. Remove the
	// redundancy!
	vector<property_t>::iterator i;
  for (i = m_prop_by_index.begin(); i != m_prop_by_index.end(); i++) {
  	if (i -> name == name) {
  		i -> description = description;
  		// Multiple occurrences of name should have been prevented by
  		// function PropertyList::addProperty
  		break;
  	}
  }

  if (i == m_prop_by_index.end()) {
  	string s = string("PropertyList::setPropDescription") + "Requested property "
  			"name '" + name + "' not found in indexed properties. Aborting. "
				"Possible properties in the current context are: ";

  	throwListIfUnknown(s);
  }
}


string PropertyList::property_t::toString() const {
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
}

/*
void PropertyList::fromSGML(const SGML &e)
{
  for (map_string_string_ci i = e.attr.begin(); i != e.attr.end(); i++) {
    if (exists(i->first)) {
      property_t p = propByName(i->first);

      switch (p.datatype) {
      case INT:
          p.asInt() = atoi(i->second.c_str());
          break;
      case DOUBLE:
          p.asDouble() = atof(i->second.c_str());
          break;
      case STRING:
          p.asString() = i->second;
          break;
      case BOOLEAN:
          if (i->second == "true" || i->second == "yes" || i->second == "1") {
            p.asBool() = true;
            break;
          } else if (i->second == "false" || i->second == "no" || i->second == "0") {
            p.asBool() = false;
            break;
          } else
            throw gError("PropertyList::fromSGML", "Property '" + i->first + 
                         "' of class '" + m_name + "' has an invalid value. "
                         "Can be 'true'|'yes'|1 or 'false'|'no'|0.");
      case FUNCTION:
          p.asFunction().setExpression(i->second);
          break;
      case EO_DATATYPE:
          throw gError("It is forbidden to use the datatype 'EO_DATATYPE'.");
      }
            
      if (p.constraint)
        if (!p.constraint->check(p))
          throw gError("PropertyList::fromSGML", "Property '" + i->first + "' has an invalid value. Constraint on value: " + p.description);
    } else {
      string s = "Unknown property '" + i->first + "' of class '" + m_name + "'. Possibilities are";
            
      if(m_prop_by_index.empty()) s += " none.";
            
      for (vector<property_t>::iterator pi = m_prop_by_index.begin(); pi != m_prop_by_index.end(); pi++)
        s += "\n  " + pi->name + "\n    " + pi->description;;
        
      throw gError("PropertyList::fromSGML", s);
    }
  }
}
*/

void PropertyList::fromXML(const xmlNode *xmln)
{
   string bigString = ObjToString(LONG_MAX);
   string smallString = ObjToString(LONG_MIN);
   int comparison;
  for (xmlAttr *xmlp = xmln->properties; xmlp; xmlp = xmlp->next) {
    if (exists((const char*) xmlp->name)) {
      property_t p = propByName((const char*) xmlp->name);
      string s = string((const char*) xmlp->children->content);

      switch (p.datatype) {
      case INT:
   if (s.size()>bigString.size()) {
      throw gError("Integer Value \"" + s + "\" out of admissible range! Please correct your input data.");
   }
   else {
        if (s.size()<bigString.size()) {
        p.asInt() = atoi(s.c_str());
        }
        else {           //s.size()==bigString.size() || s.size()==smallString.size()
             comparison = s.compare(bigString);
                if (comparison<=0){
               p.asInt() = atoi(s.c_str());
                }
                else {
                throw gError("Integer Value \"" + s + "\" out of admissible range! Please correct your input data.");
                }
             }
         }
          break;
      case DOUBLE:
          p.asDouble() = atof(s.c_str());
          break;
      case STRING:
          p.asString() = s;
          break;
      case BOOLEAN:
          if (s == "true" || s == "yes" || s == "1") {
            p.asBool() = true;
            break;
          } else if (s == "false" || s == "no" || s == "0") {
            p.asBool() = false;
            break;
          } else
            throw gError("PropertyList::fromXML", "(module " + m_name + ") Property '" + string((const char*) xmlp->name) + 
                         "' of class '" + m_name + "' has an invalid value. "
                         "Can be 'true'|'yes'|1 or 'false'|'no'|0.");
      case FUNCTIONPAIR:
          p.asFunctionPair().setExpression(s);
          break;
      case FUNCTIONFIXED:
          p.asFunctionFixed().setExpression(s);
          break;
      case POINT:
          p.asPoint() = string2point(s);
          break;
      case EO_DATATYPE:
          throw gError("It is forbidden to use the datatype 'EO_DATATYPE'.");
      }
            
      if (p.constraint)
        if (!p.constraint->check(p))
          throw gError("PropertyList::fromXML", "(module " + m_name + ") Property '" + string((const char*) xmlp->name) + "' has an invalid value. Constraint on value: " + p.description);
    } else {
      if (m_allow_unknown) {
        m_unknown_props[string((const char*) xmlp->name)] = string((const char*) xmlp->children->content);
      } else {
	
        string s = "(module " + m_name + ") Unknown property '" + 
          string((const char*) xmlp->name) + "'. Possibilities are:";
            
// 	MSG_DEBUG("PropertyList::fromXML",  "(" + m_name + ") Unknown property '" + string((const char*) xmlp->name) + "'");
	throwListIfUnknown(s);

      }
    }
  }


  /* Now check if all properties have been set correctly. */

  for (vector<property_t>::iterator i = m_prop_by_index.begin(); i != m_prop_by_index.end(); ++i) {
     switch (i->datatype) {
     case FUNCTIONPAIR:
         if (i->asFunctionPair().isNull()) {
           throw gError
             ("PropertyList::fromXML",
              "(module " + m_name + ") Please provide a value for property '" + i->name + "'.");
         }
         break;
     case FUNCTIONFIXED:
         if (i->asFunctionFixed().isNull()) {
           throw gError
             ("PropertyList::fromXML",
              "(module " + m_name + ") Please provide a value for property '" + i->name + "'.");
         }
         break;
     default: break;
#if 0
     default:
         /* We just assume that if a) the value is not set in the input file and b)
            the constraint is not met the programmer wants the value to be provided
            always (no default value). */
         if (i->constraint)
           if (!i->constraint->check(*i)) {
             throw gError
               ("PropertyList::fromXML",
                "(module " + m_name + ") Please provide a value for property '" + i->name + "'.");
           }
         break;
#endif
     }
  }
}

void PropertyList::throwListIfUnknown(string s)
{
  
  if(m_prop_by_index.empty()) s += " none.";
  
  for (vector<property_t>::iterator pi = m_prop_by_index.begin(); pi != m_prop_by_index.end(); pi++)
    s += "\n  " + pi->name + "\n    " + pi->description;
  
  throw gError("PropertyList::throwListIfUnknown", s);
//   MSG_DEBUG("PropertyList::fromXML", s);
//   abort();
}

ostream &PropertyList::toXML_begin(ostream &s, int shift)
{
  s << PutTab(shift) << "<" << m_class_name;
    
  for (vector<property_t>::iterator i = m_prop_by_index.begin(); i != m_prop_by_index.end(); i++) {
    s << endl << PutTab(shift+1) << i->name << " = ";
        
    switch (i->datatype) {
    case INT:
        s << "\"" << i->asInt() << "\"";
        break;
    case DOUBLE:
        s << "\"" << i->asDouble() << "\"";
        break;
    case STRING:
        s << "\"" << i->asString() << "\"";
        break;
    case BOOLEAN:
        if (i->asBool())
          s << "\"yes\"";
        else
          s << "\"no\"";
        break;
    case FUNCTIONPAIR:
        s << "\"" << i->asFunctionPair().expression() << "\"";
        break;
    case FUNCTIONFIXED:
        s << "\"" << i->asFunctionFixed().expression() << "\"";
        break;
    case POINT:
        s << "\"" << i->asPoint() << "\"";
        break;
    case EO_DATATYPE:
        throw gError("It is forbidden to use the datatype 'EO_DATATYPE'.");
    }
  }
    
  s << endl << PutTab(shift) << ">" << endl;
  return s;
}


ostream &PropertyList::toXML_end(ostream &s, int shift)
{
    s << PutTab(shift) << "</" << m_class_name << ">" << endl;
    return s;
}



/*---- PropertyListConstraint ----*/
//---- Constructors/Destructor ----

PropertyListConstraint::PropertyListConstraint()
{
}


PropertyListConstraint::~PropertyListConstraint()
{
}


