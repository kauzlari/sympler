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


#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "function_compiler.h"

// #include "valgrind/memcheck.h"

int FunctionCompiler::s_function_counter = 0;


// bool fexists(const char *filename)
// {
//   FILE *file = fopen(filename,"r");
// 
//   if(file) {
//     fclose(file);
//     return 1;
//   }
//   else {
//     return 0;
//   }
// }
// 


/* FunctionCompiler */

FunctionCompiler::FunctionCompiler(): m_handle(NULL), m_fn(NULL)
{
}


FunctionCompiler::~FunctionCompiler()
{
  close();
}


void FunctionCompiler::close()
{
  if (m_handle) {
    dlclose(m_handle);
  }
}


/* Fixme!!! We have dynamically determine
   - which compiler to use
   - where to compile (i.e. /tmp directory changes)
*/
bool FunctionCompiler::compile()
{
  ofstream f;
  stringstream s;
  int r;
  string ext;
 // srand ( time(NULL) );
 // int numberForCounter = rand() % 1000000 + 1;
  
   	char *pTemp;
	if(!(pTemp = std::getenv("TMP"))){
//	if (strlen(pTemp)==0){
		pTemp=(char *)"/tmp";
		}
  /* Now we _hard_ compile the function. Hopefully gcc won't boil out. */

#ifdef __APPLE__
  ext = ".dylib";
#elif WIN32
  /* Fixme!!! This correct? Check! */
  ext = ".dll";
#else
  ext = ".so";
#endif

//  s << "/tmp/" C_FC_PREFIX << s_function_counter << "." << numberForCounter;

s << pTemp << "/" << C_FC_PREFIX << s_function_counter;

  s_function_counter++;
  struct stat dummy_stat;
  while /*(fexists((s.str()+ext).c_str()) || fexists((s.str()+".c").c_str()))*/ 
  (stat((s.str()+ext).c_str(), &dummy_stat) == 0 || stat((s.str()+".c").c_str(), &dummy_stat) == 0) {
    s.seekp(0);
    s << pTemp << "/" << C_FC_PREFIX << s_function_counter;
    s_function_counter++;
    MSG_DEBUG("FunctionCompiler::compile", "counter for existing source files in " << pTemp <<  ":" << s_function_counter);
  }

  m_so_filename = s.str()+ext;

  f.open((s.str()+".c").c_str());
  f << "#include <math.h>" << endl;
  f << "#include <stdlib.h>" << endl;
  f << "void " << C_FC_FN_NAME << "(" << m_header << ")" << endl;

//   MSG_DEBUG("FunctionCompiler::compile", "m_parser->expression = " << m_parser->expression());

  Variant v = m_parser->toC();
  
//   MSG_DEBUG("FunctionCompiler::compile", "1st C-expression: " << v.strings()[0]);
  
  assert(v.typeId() == Variant::SCALAR_STRING ||
	 v.typeId() == Variant::VECTOR_STRING ||
	 v.typeId() == Variant::TENSOR_STRING);

//   assert(v.typeId() == Variant::TENSOR_STRING);
  
  if (v.strings().size() != m_result_vars.size()) {
    f.close();


    throw gError
      ("FunctionCompiler::compile",
       "Type mismatch in return value. Expression: '" + m_parser->expression() + "'. "
       "I expected a return type of size (entries) " + ObjToString(m_result_vars.size()) + ", "
       "but got " + ObjToString(v.strings().size()) + " entrie(s).");
  }

  f << "{" << endl;
  
/*  f << "printf(\"RHS0= %f\\n \", " << v.strings()[0] << ");" << endl;
  f << "printf(\"RHS0: " << v.strings()[0] << "\\n \");" << endl;
  f << "printf(\"RHS0= %f\\n \", " << v.strings()[0] << ");" << endl;
  f << "printf(\"SPECIAL Before: %f\\n \", (*((double*) ((char*) particle_tag + 8))));" << endl;
  f << "printf(\"Result Before: %f\\n \", (*((double*) ((char*) result + 0))));" << endl;*/

  //   f << "printf(\"Before: %f\\n \", (*((double*) ((char*) result + 8))));" << endl;
//   f << "printf(\"Before: %f\\n \", (*((double*) ((char*) result + 16))));" << endl;
//   f << "printf(\"Before ni =  %f\\n \", (*((double*) ((char*) first_tag + 16))));" << endl;
  
  for (size_t i = 0; i < m_result_vars.size(); ++i) {
    FUNCTION_PARSER_LOG("FunctionCompiler::compile","exporting \"" << m_result_vars[i] << " = " << v.strings()[i] << ";\"");
//     cout << "FunctionCompiler::compile exporting \"" << m_result_vars[i] << " = " << v.strings()[i] << ";\"" << endl;

    f << "   " << m_result_vars[i] << " = " << v.strings()[i] << ";" << endl;
  }
  
/*  f << "printf(\"SPECIAL After: %f\\n \", (*((double*) ((char*) particle_tag + 8))));" << endl;
  f << "printf(\"Results After: %f\\n \", (*((double*) ((char*) result + 0))));" << endl;*/
  
//   f << "printf(\"After: %f\\n \", (*((double*) ((char*) result + 8))));" << endl;
//   f << "printf(\"After: %f\\n \", (*((double*) ((char*) result + 16))));" << endl;
//   f << "printf(\"After ni =  %f\\n \", (*((double*) ((char*) first_tag + 16))));" << endl;
  f << "}" << endl;

  f.close();

  /* Fixme!!! Check for installed compilers, temp directory, etc... */
#ifdef __APPLE__
  r = system
    (("gcc -O3 -dynamiclib -fPIC -nostartfiles -o " + m_so_filename + " " + s.str() + ".c -lm").c_str());
#else
  r = system
    (("gcc -O3 -shared -fPIC -nostartfiles -o " + m_so_filename + " " + s.str() + ".c -lm").c_str());
#endif

  remove((s.str() + ".c").c_str());

  return r == 0;
}


void FunctionCompiler::setParserAndCompile(FunctionParser *parser)
{
//   MSG_DEBUG("FunctionCompiler::setParserAndCompile", "START");
  close();

  m_parser = parser;
  
  if (compile()) {
//     MSG_DEBUG("FunctionCompiler::setParserAndCompile", "IF: true-case");
    char *error_str;

//     error_str = "";
    
    // experience with Ubuntu hardy (2009-08-06):
    // when using valgrind --tool=memcheck, we get a message "Invalid read of size 4" 
    // for the following line of code. This is known when searching with google. To 
    // get rid of this message we have to install libc6-dbg; Since this is suring 
    // initialisation, we can probably afford the small performance loss (Of course 
    // this is only true if this lib isn't used somewhere else, which I currently 
    // don't know). On the other hand, the error does not seem to have any side effect. 
    m_handle = dlopen(m_so_filename.c_str(), RTLD_NOW);
    
//    MSG_DEBUG("FunctionCompiler::setParserAndCompile", "BEFORE handle-if, m_handle = " << m_handle);
    if (!m_handle)
      throw gError
          ("FunctionCompiler::setParserAndCompile", "error opening " + m_so_filename + ": " + dlerror());

    // this call is for discarding the last error message
    dlerror();
    m_fn = (gfnptr_t) dlsym(m_handle, C_FC_FN_NAME);
    if (!m_fn)
	    throw gError
          ("FunctionCompiler::setParserAndCompile", "Didn't get a function from " + m_so_filename + ".");
    error_str = dlerror();
//     MSG_DEBUG("FunctionCompiler::setParserAndCompile", "BEFORE error_str-if: error_str = " /*<< error_str*/);
    if (error_str)
	    throw gError
          ("FunctionCompiler::setParserAndCompile", "error during linking of " + m_so_filename + ": " + error_str);
    
//     MSG_DEBUG("FunctionCompiler::setParserAndCompile", "BEFORE remove .so");

    remove(m_so_filename.c_str());
//     MSG_DEBUG("FunctionCompiler::setParserAndCompile", "AFTER remove .so");
  
  } else
      throw gError
        ("FunctionCompiler::setParserAndCompile",
        "Failed to compile the function using gcc.");
//   MSG_DEBUG("FunctionCompiler::setParserAndCompile", "END");
}
