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



/*
   This file contains geomerty entities, like points, lines, etc.
   and certains operations that can be performed.
*/

#ifndef __GEOMETRIC_PRIMITIVES_H
#define __GEOMETRIC_PRIMITIVES_H

#include <ostream>
#include <iostream>
// for gcc4.3
#include "string.h"

#include <math.h>

using namespace std;

#include "misc.h"
//---- Defines ----

#define SPACE_DIMS 3

#define TOCELLINDEX(p, n)  ((p.z * n.y) + p.y) * n.x + p.x
#define TOCELLINDEXD(x, y, z, n)  ((z * n.y) + y) * n.x + x

/* Epsilon value for geometry calculations. */
extern double g_geom_eps;

template<typename T> struct math_tensor_t;
typedef math_tensor_t<double> tensor_t;

template<typename T> struct math_tensor3_t;
typedef math_tensor3_t<double> tensor3_t;




//---- Types ----

/*!
 * A vector type which implements scalar multiplication, cross products, etc.
 */
template<typename T>
struct math_vector_t {
/* Datatypes */
  union {
    T coords[SPACE_DIMS];
    struct {
      T x, y, z;
    };
  };

        
  /*!
   * Access to coordinates
   * @param i Return coordinate with this index
   */
  inline T &operator[](int i) {
    return coords[i];
  }

  /*!
   * Access to coordinates
   * @param i Return coordinate with this index
   */
  inline const T &operator[](int i) const {
    return coords[i];
  }

  /*!
   * Length of the vector squared
   */
  inline T absSquare() const {
    T h = 0;
    
    for (int i = 0; i < SPACE_DIMS; ++i)
      h+= coords[i]*coords[i];
    
    return h;
  }

  /*!
   * Length of the vector
   */	
  inline T abs() const {
    return sqrt(absSquare());
  }
	
  /*!
   * Comparison operator
   * @param p Compare to this vector
   */
  inline bool operator==(const math_vector_t<T> &p) const {
    bool id = true;
    
    for (int i = 0; id && (i < SPACE_DIMS); ++i)
      id = (coords[i] == p[i]);
    
    return id;
  }    

  /*!
   * Addition
   * @param p Add this vector to the current vector
   */
  inline math_vector_t<T> operator+(const math_vector_t<T> &p) const {
    math_vector_t<T> n;
    
    for (int i = 0; i < SPACE_DIMS; ++i)
      n[i] = coords[i] + p[i];
		
    return n;
  }
	
  inline math_vector_t<T> &operator+=(const math_vector_t<T> &p) {
    for (int i = 0; i < SPACE_DIMS; ++i)
      coords[i] += p[i];
    
    return *this;
  }    	

    inline math_vector_t<T> operator-(const math_vector_t<T> &p) const {
        math_vector_t<T> n;

        for (int i = 0; i < SPACE_DIMS; ++i)
            n[i] = coords[i] - p[i];

        return n;
    }

    inline math_vector_t<T> operator-() const {
        math_vector_t<T> n;

        for (int i = 0; i < SPACE_DIMS; ++i)
            n[i] = -coords[i];

        return n;
    }

    inline math_vector_t<T> operator-=(const math_vector_t<T> &p) {
        for (int i = 0; i < SPACE_DIMS; ++i)
            coords[i] -= p[i];
			
		return *this;
    }    

/* Scalar product */
    inline T operator*(const math_vector_t<T> &p) const {
        T result = 0;

        for (int i = 0; i < SPACE_DIMS; ++i)
            result += coords[i]*p[i];

        return result;
    }    
	
/* Multiplication with a scalar */
    inline math_vector_t<T> operator*(T d) const {
        math_vector_t<T> result;
		
        for (int i = 0; i < SPACE_DIMS; ++i)
            result[i] = coords[i] * d;
		
        return result;
    }
    inline math_vector_t<T> operator*=(T d) {
        for (int i = 0; i < SPACE_DIMS; ++i)
            coords[i] *= d;
		
		return *this;
    }
	
    inline math_vector_t<T> operator/(T d) const {
        math_vector_t<T> n;

        for (int i = 0; i < SPACE_DIMS; ++i)
            n[i] = coords[i] / d;

        return n;
    }
    inline math_vector_t<T> operator/=(T d) {
        for (int i = 0; i < SPACE_DIMS; ++i)
            coords[i] /= d;
			
		return *this;
    }

    inline void /*math_vector_t<T>*/ assign(T v) {
        for (int i = 0; i < SPACE_DIMS; ++i)
            coords[i] = v;
// return *this;
    }

/* Cross product */
	inline math_vector_t<T> cross(const math_vector_t<T> &p) const {
		math_vector_t<T> result;
		
		for (int i = 0; i < SPACE_DIMS; ++i) {
			int n, nn;
			n = (i + 1) % SPACE_DIMS;
			nn = (i + 2) % SPACE_DIMS;
			result[i] = coords[n]*p[nn] - coords[nn]*p[n];
		}

		return result;
	}

/* Divide element by element */
	inline math_vector_t<T> listDivide(const math_vector_t<T> &p) const {
		math_vector_t<T> result;
		for (int i = 0; i < SPACE_DIMS; ++i) {
			result[i] = coords[i]/p[i];
		}
		return result;
	}

/* Multiply element by element */
	inline math_vector_t<T> listMultiply(const math_vector_t<T> &p) const {
		math_vector_t<T> result;
		for (int i = 0; i < SPACE_DIMS; ++i) {
			result[i] = coords[i]*p[i];
		}
		return result;
	}

  /* outer product */
	inline math_tensor_t<T> outer(const math_vector_t<T> &p) const {
		math_tensor_t<T> result;
		
		for (int i = 0; i < SPACE_DIMS; ++i) {
		  for (int j = 0; j < SPACE_DIMS; ++j) {
		    
		    result.tensor[j+i*SPACE_DIMS] = coords[i]*p[j];
		  }
		}
		return result;
	}

};

typedef math_vector_t<int> int_point_t;
typedef math_vector_t<bool> bool_point_t;
// GLOBAL FUNCTIONS

/* << operator */

template<typename T>
ostream& operator<<(ostream &out, const math_vector_t<T> v) {
	out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	
	return out;
}

template<typename T>
  inline math_vector_t<T> operator*(T d, const math_vector_t<T> &p) {
    return p*d;
  }
      


/* Specific geometric quantities, i.e. the above with T = double */

typedef math_vector_t<double> point_t;

inline point_t operator*(int n, const point_t &p) {
	return p*n;
}
inline point_t operator/(double d, const point_t &p) {
	point_t tmp_p;
	tmp_p.x = p.x / d;
	tmp_p.y = p.y /d;
	tmp_p.z = p.z / d;
	return tmp_p;
}

// GLOBAL FUNCTIONS

inline void g_pointTimesEqInt(point_t& p, const int& i) {
  p.x *= i;
  p.y *= i;
  p.z *= i;
}

inline void g_pointTimesEqDbl(point_t& p, const double& d) {
  p.x *= d;
  p.y *= d;
  p.z *= d;
}

inline void g_pointTimesInt(point_t& p, const point_t& p2, const int& i) {
  p.x = i*p2.x;
  p.y = i*p2.y;
  p.z = i*p2.z;
}

inline void g_pointTimesDbl(point_t& p, const point_t& p2, const double& d) {
  p.x = d*p2.x;
  p.y = d*p2.y;
  p.z = d*p2.z;
}

inline void g_pointPlusEqPoint(point_t& p, const point_t& p2) {
  p.x += p2.x;
  p.y += p2.y;
  p.z += p2.z;
}

inline void g_pointPlusPoint(point_t& p, const point_t& p2, const point_t& p3) {
  p.x = p2.x + p3.x;
  p.y = p2.y + p3.y;
  p.z = p2.z + p3.z;
}


/* Geometric objects */

struct line_t {
	union {
		point_t points[2];
		struct {
			point_t from, to;
		};
	};
	
	/* [] access */
    inline point_t &operator[](int i) {
        return points[i];
    }
	
    inline const point_t &operator[](int i) const {
        return points[i];
    }	
};


struct cuboid_t {
    /* All coordinates of corner1 have to be smaller than
    those of corner2. */
	union {
		point_t corners[2];
		struct {
			point_t corner1, corner2;
		};
	};

    cuboid_t(){}

    cuboid_t(point_t a_c1, point_t a_c2)
    {
        corner1 = a_c1;
        corner2 = a_c2;
    }
	
	/* [] access */
    inline point_t &operator[](int i) {
        return corners[i];
    }
	
    inline const point_t &operator[](int i) const {
        return corners[i];
    }		
    
    /* For sorting! fixme!!! That's total bullshit */
  /*
    inline bool operator<(const cuboid_t &c) const {
        return corner1 < c.corner1;
    }
  */

    inline point_t size() const {
        return corner2-corner1;
    }

    inline bool isInside(const point_t &pos) const {
        for (int i = 0; (i < SPACE_DIMS); i++) {
            if ((pos[i] < corner1[i]) || (pos[i] >= corner2[i]))
	    {		
               	return false;
            }
        }
        
		return true;
    }

    inline bool isInsideEps(const point_t &pos, double eps) const {
        for (int i = 0; (i < SPACE_DIMS); i++) {
            if ((pos[i] < corner1[i]-eps) || (pos[i] >= corner2[i]+eps))
                return false;
        }
        
		return true;
    }


    /* Fixme!!! Intersect was designed for WallTriangle
       the >, >=, <, <= or not the way it should be. */
 #define CHECK_I(c1, c2)                                                 \
    if (x >= 0 && x <= 1) {                                             \
        point_t h = a + d*x;                                            \
                                                                      \
        if (h[c1] >= corner1[c1] && h[c2] >= corner1[c2] && h[c1] <= corner2[c1] && h[c2] <= corner2[c2]) \
	  { \
\
	    /* if			     \
  (corner1[0] > 16.1 && corner1[0] < 17.9 && \
   corner1[1] > 1.6 && corner1[1] < 3.4 && \
	    corner1[2] > 0.5 && corner1[2] < 1.625) \
     cout << "g_primitives:intersect: " << "TRUE: rl=" << rl << "(corner1-a="  << corner1-a << ", corner2-a=" << corner2-a << ")" << endl  << "corner1="  << corner1 << ", corner2="  << corner2 << ", c1=" << c1 << ", c2=" << c2 << endl; \
	    */								\
return true;				\
	  }								\
} while(0)								\

    inline bool intersects(const point_t &a, const point_t &b) const {
        /* Test all six sides for intersection. */
        double x;
        point_t rl, d;

        if (isInside(a) || isInside(b))
            return true;

        d = b-a;
        rl = corner1-a;

	if(d.x != 0)
	  {
	    x = rl.x/d.x;
	    CHECK_I(1, 2);
	  }
	if(d.y != 0)
	  {
	    x = rl.y/d.y;
	    CHECK_I(0, 2);
	  }
	if(d.z != 0)
	  {
	    x = rl.z/d.z;
	    CHECK_I(0, 1);
	  }
        rl = corner2-a;

	if(d.x != 0)
	  {
	    x = rl.x/d.x;
	    CHECK_I(1, 2);
	  }
	if(d.y != 0)
	  {
	    x = rl.y/d.y;
	    CHECK_I(0, 2);
	  }
	if(d.z != 0)
	  {
	    x = rl.z/d.z;
	    CHECK_I(0, 1);
	  }
        return false;
    }
#undef CHECK_I
};

#define SPACE_DIMS_SQUARED SPACE_DIMS*SPACE_DIMS

       
/*!
* Tensor definition with single C-array
*/
template<typename T>
struct math_tensor_t {
/* Datatypes */
	union {
		T tensor[SPACE_DIMS_SQUARED]/*[SPACE_DIMS]*/;
		struct {
			T xx, xy, xz, yx, yy, yz, zx, zy, zz;
		};
	};

 // if used in a tag, the following would have to be set
//   size_t size;
 
  math_tensor_t<T>() {
    memset(tensor, 0, SPACE_DIMS_SQUARED*sizeof(T));
//     size = SPACE_DIMS_SQUARED;
  
//     cout << "{math_tensor_t<T>::math_tensor_t<T>}" << "size = " << size << endl;
  }

  math_tensor_t<T>(const math_tensor_t<T> &copy) {
    memcpy(tensor, copy.tensor, SPACE_DIMS_SQUARED*sizeof(T));
//     size = SPACE_DIMS_SQUARED;
//     cout << "{math_tensor_t<T>::copy-constructor}" << "size = " << size << endl;
  }

  math_tensor_t<T> operator=(const math_tensor_t<T> &copy) {
    memcpy(tensor, copy.tensor, SPACE_DIMS_SQUARED*sizeof(T));
//     size = SPACE_DIMS_SQUARED;
    return *this;
//     cout << "{math_tensor_t<T>::operator=}" << "size = " << size << endl;
  }

  T &operator()(size_t i, size_t j) {
//     return tensor[i][j];
    return tensor[j+i*SPACE_DIMS];
  }
  const T &operator()(size_t i, size_t j) const {
//     return tensor[i][j];
    return tensor[j+i*SPACE_DIMS];
  }

  inline math_tensor_t<T> operator+(math_tensor_t<T> d) const {
    math_tensor_t<T> n;
//     size_t sdq = SPACE_DIMS_SQUARED;
    
//     for (int i = 0; i < SPACE_DIMS; ++i)
//       for (int j = 0; j < SPACE_DIMS; ++j)
//         n.tensor[i][j] = tensor[i][j] + d.tensor[i][j];

    for(size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
      n.tensor[i] = tensor[i] + d.tensor[i];
    return n;
  }
  
  inline math_tensor_t<T> operator+=(math_tensor_t<T> d) {
//     size_t sdq = SPACE_DIMS_SQUARED;

    /*    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        tensor[i][j] += d.tensor[i][j];*/
    for(size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
      tensor[i] += d.tensor[i];
//     cout << "{math_tensor_t::operator +=}" << "size = " << size << endl << "tensor = " << tensor << endl << "d.tensor = " << d.tensor;    			
		return *this;
  }

  inline math_tensor_t<T> operator-(math_tensor_t<T> d) const {
    math_tensor_t<T> n;
//     size_t sdq = SPACE_DIMS_SQUARED;

/*    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        n.tensor[i][j] = tensor[i][j] - d.tensor[i][j];*/
    for(size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
      n.tensor[i] = tensor[i] - d.tensor[i];

    return n;
  }
  
  inline math_tensor_t<T> operator-=(math_tensor_t<T> d) {
//     size_t sdq = SPACE_DIMS_SQUARED;
/*    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        tensor[i][j] -= d.tensor[i][j];*/
    for (int i = 0; i < SPACE_DIMS_SQUARED; ++i)
        tensor[i] -= d.tensor[i];
			
		return *this;
  }

  inline math_tensor_t<T> operator-(T d) const {
    math_tensor_t<T> n;
    size_t incr = SPACE_DIMS+1;
/*    for (int i = 0; i < SPACE_DIMS; ++i)
      n.tensor[i][i] = tensor[i][i] - d;*/
    for (int i = 0; i < SPACE_DIMS_SQUARED; i+=incr)
      n.tensor[i] = tensor[i] - d;

    return n;
  }
  inline math_tensor_t<T> operator-=(T d) {
    size_t incr = SPACE_DIMS+1;
/*    for (int i = 0; i < SPACE_DIMS; ++i)
      tensor[i][i] -= d;*/
    for (size_t i = 0; i < SPACE_DIMS_SQUARED; i+=incr)
      tensor[i] -= d;
			
		return *this;
  }

  inline math_tensor_t<T> operator*(T d) const {
    math_tensor_t<T> n;

/*    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        n.tensor[i][j] = tensor[i][j] * d;*/
    for (size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
        n.tensor[i] = tensor[i] * d;

    return n;
  }
  inline math_tensor_t<T> operator*=(T d) {
/*    for (int i = 0; i < SPACE_DIMS; ++i)
    for (int j = 0; j < SPACE_DIMS; ++j)
    tensor[i][j] *= d;*/
    for (size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
      tensor[i] *= d;
			
    return *this;
  }

  inline math_tensor_t<T> operator*=(math_tensor_t<T> t) {
/*    for (int i = 0; i < SPACE_DIMS; ++i)
    for (int j = 0; j < SPACE_DIMS; ++j)
    tensor[i][j] *= d;*/
    for (size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
      tensor[i] *= t.tensor[i];
			
    return *this;
  }

  inline math_tensor_t<T> operator/(T d) const {
    math_tensor_t<T> n;

    for (size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
        n.tensor[i] = tensor[i] / d;
//     for (int i = 0; i < SPACE_DIMS; ++i)
//       for (int j = 0; j < SPACE_DIMS; ++j)
//         n.tensor[i][j] = tensor[i][j] / d;

    return n;
  }
  inline math_tensor_t<T> operator/=(T d) {
    for (size_t i = 0; i < SPACE_DIMS_SQUARED; ++i)
        tensor[i] /= d;
/*    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        tensor[i][j] /= d;*/
			
		return *this;
  }

  inline T trace() const {
    T __t = 0;
//     cout << "{math_tensor_t::trace()}" << "CALLED" << endl;
    size_t __incr = SPACE_DIMS+1;
    
    for (size_t __i = 0; __i < SPACE_DIMS_SQUARED; __i+=__incr)
      __t += tensor[__i];
//     for (int __i = 0; __i < SPACE_DIMS; ++__i)
//       __t += tensor[__i][__i];

    return __t;
  }

  /* matrix vector operations*/
  
  inline math_vector_t<T> operator*(math_vector_t<T> vec) {
    math_vector_t<T> tempVec = {{{0, 0, 0}}};
    size_t k = 0;
    for (int i = 0; i < SPACE_DIMS; ++i)
    {
      for (int j = 0; j < SPACE_DIMS; ++j)
      {
        tempVec[i] += vec[j]*tensor[k/*j+i*SPACE_DIMS*/]/*[j]*/;
        ++k;
      }
    }
    return tempVec;
  }
  /* end: matrix vector operations*/
   
  inline math_vector_t<T> contract() {
    math_vector_t<T> t = { { { 0, 0, 0 } } };
    size_t k = 0;

/*    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        t[i] += tensor[i][j];*/
    for (int i = 0; i < SPACE_DIMS; ++i)
    {
      for (int j = 0; j < SPACE_DIMS; ++j)
      {
        t[i] += tensor[k/*j+i*SPACE_DIMS*/]/*[j]*/;
        ++k;
      }
    }

    return t;
  }

  inline double Q() {
    double q = 0;

    for (size_t i = 0; i < SPACE_DIMS_SQUARED/*SPACE_DIMS*/; ++i)
//       for (int j = 0; j < SPACE_DIMS; ++j)
        q += tensor[i]/*[j]*/*tensor[i]/*[j]*/;

    return q;
  }

  void assign(T v) {
//     memset(&tensor, sizeof(T)*SPACE_DIMS*SPACE_DIMS, 0);
    for (size_t i = 0; i < SPACE_DIMS_SQUARED/*SPACE_DIMS*/; ++i)
//       for (int j = 0; j < SPACE_DIMS; ++j)
        tensor[i]/*[j]*/ = v;

  }
};

#define SPACE_DIMS_CUBED SPACE_DIMS*SPACE_DIMS*SPACE_DIMS
       
/*!
* 3rd order tensor definition with single C-array
*/
template<typename T>
struct math_tensor3_t {
/* Datatypes */
/* 	union { */
		T tensor[SPACE_DIMS_CUBED];
/* 		struct { */
/* 			T xxx, xxy, xxz, xyx, xyy, xyz, xzx, xzy, xzz;*/
/* 			T yxx, yxy, yxz, yyx, yyy, yyz, yzx, yzy, yzz;*/
/* 			T zxx, zxy, zxz, zyx, zyy, zyz, zzx, zzy, zzz;*/
/* 		}; */
/* 	}; */
};



template<typename T>
ostream& operator<<(ostream &out, const math_tensor_t<T> t)
{
  out << "tensor(";
  for (int i = 0; i < SPACE_DIMS; i++) {
    out << "(";
      for (int j = 0; j < SPACE_DIMS; j++) {
        out << t(i, j);
        if (j != SPACE_DIMS-1)
          out << ", ";
      }
    out << ")";
    if (i != SPACE_DIMS-1)
      out << ", ";
  }
  out << ")";
	
	return out;
}



/*
template<typename T>
inline math_vector_t<T> operator*(T d, const math_tensor_t<T> &p) {
	return p*d;
}
*/

// was previously HERE
/* typedef math_tensor_t<double> tensor_t; */

inline tensor_t operator*(double d, const tensor_t &t) {
	return t*d;
}


inline void g_outerPointPoint(tensor_t& t, const point_t& p2, const point_t& p3) {
  for (int i = 0; i < SPACE_DIMS; ++i) {
    int idx = i*SPACE_DIMS;
    for (int j = 0; j < SPACE_DIMS; ++j) {
      t.tensor[j+idx] = p2[i]*p3[j];
    }
  }
}

inline void g_outerTensorPoint(tensor3_t& t, const tensor_t& t2, const point_t& p3) {
  for (int i = 0; i < SPACE_DIMS; ++i) {
    int idx = i*SPACE_DIMS;
    for (int j = 0; j < SPACE_DIMS; ++j) {
      int jdx = j+idx;
      int jdxx = SPACE_DIMS*jdx;
      for (int k = 0; k < SPACE_DIMS; ++k) {
	t.tensor[k+jdxx] = t2.tensor[jdx]*p3[k];
      }
    }
  }
}

/*!
 * Here the point_t contributes to the second index of the resulting tensor3_t
 */
inline void g_outer2TensorPoint(tensor3_t& t, const tensor_t& t2, const point_t& p3) {
  /* commented out versions are from g_outerTensorPoint */
  for (int i = 0; i < SPACE_DIMS; ++i) {
    int idx = i*SPACE_DIMS;
    for (int j = 0; j < SPACE_DIMS; ++j) {
/*       int jdx = j+idx; */
/*       int jdxx = SPACE_DIMS*jdx; */
      int jdxx = SPACE_DIMS*(j+idx);
      for (int k = 0; k < SPACE_DIMS; ++k) {
/* 	t.tensor[k+jdxx] = t2.tensor[jdx]*p3[k]; */
	t.tensor[k+jdxx] = t2.tensor[k+idx]*p3[j];
      }
    }
  }
}

inline void g_outerPointTensor(tensor3_t& t, const point_t& p2, const tensor_t& t3) {
  for (int i = 0; i < SPACE_DIMS; ++i) {
    int idx = i*SPACE_DIMS;
    for (int j = 0; j < SPACE_DIMS; ++j) {
      int jdx = j*SPACE_DIMS;
      int jdxx = SPACE_DIMS*idx+jdx;
      for (int k = 0; k < SPACE_DIMS; ++k) {
	t.tensor[k+jdxx] = p2[i]* t3.tensor[k+jdx];
      }
    }
  }
}

inline void g_tensorTimesEqInt(tensor_t& t, const int& i) {
  t.xx *= i;
  t.xy *= i;
  t.xz *= i;
  t.yx *= i;
  t.yy *= i;
  t.yz *= i;
  t.zx *= i;
  t.zy *= i;
  t.zz *= i;
}

inline void g_tensorTimesEqDbl(tensor_t& t, const double& d) {
  t.xx *= d;
  t.xy *= d;
  t.xz *= d;
  t.yx *= d;
  t.yy *= d;
  t.yz *= d;
  t.zx *= d;
  t.zy *= d;
  t.zz *= d;
}

inline void g_tensorTimesDbl(tensor_t& t, const tensor_t& t2, const double& d) {
  t.xx = t2.xx * d;
  t.xy = t2.xy * d;
  t.xz = t2.xz * d;
  t.yx = t2.yx * d;
  t.yy = t2.yy * d;
  t.yz = t2.yz * d;
  t.zx = t2.zx * d;
  t.zy = t2.zy * d;
  t.zz = t2.zz * d;
}

inline void g_tensorPlusEqTensor(tensor_t& t, const tensor_t& t2) {
  t.xx += t2.xx;
  t.xy += t2.xy;
  t.xz += t2.xz;
  t.yx += t2.yx;
  t.yy += t2.yy;
  t.yz += t2.yz;
  t.zx += t2.zx;
  t.zy += t2.zy;
  t.zz += t2.zz;
}

inline void g_tensorPlusTensor(tensor_t& t, const tensor_t& t2, const tensor_t& t3) {
  t.xx = t2.xx + t3.xx;
  t.xy = t2.xy + t3.xy;
  t.xz = t2.xz + t3.xz;
  t.yx = t2.yx + t3.yx;
  t.yy = t2.yy + t3.yy;
  t.yz = t2.yz + t3.yz;
  t.zx = t2.zx + t3.zx;
  t.zy = t2.zy + t3.zy;
  t.zz = t2.zz + t3.zz;
}

inline void g_idMinusEqTensor(tensor_t& t) {
  t.xx = 1-t.xx;
  t.xy = -t.xy;
  t.xz = -t.xz;
  t.yx = -t.yx;
  t.yy = 1-t.yy;
  t.yz = -t.yz;
  t.zx = -t.zx;
  t.zy = -t.zy;
  t.zz = 1-t.zz;
}

inline void g_idMinusEqTensor2D(tensor_t& t) {
  t.xx = 1-t.xx;
  t.xy = -t.xy;
  t.xz = 0;
  t.yx = -t.yx;
  t.yy = 1-t.yy;
  t.yz = 0;
  t.zx = 0;
  t.zy = 0;
  t.zz = 0;
}

inline void g_tensorDotPoint(point_t& p, const tensor_t& t2, const point_t& p2) {
  size_t k = 0;
  for (int i = 0; i < SPACE_DIMS; ++i)
    {
      p[i] = 0;
      for (int j = 0; j < SPACE_DIMS; ++j)
	{
	  p[i] += p2[j]*t2.tensor[k/*j+i*SPACE_DIMS*/];
	  ++k;
	}
    }
}


inline void g_transpose(tensor_t& t, const tensor_t& t2) {
  t.xx = t2.xx;
  t.xy = t2.yx;
  t.xz = t2.zx;
  t.yx = t2.xy;
  t.yy = t2.yy;
  t.yz = t2.zy;
  t.zx = t2.xz;
  t.zy = t2.yz;
  t.zz = t2.zz;
}


inline void g_tensor3TimesEqInt(tensor3_t& t, const int& i) {
  size_t max = SPACE_DIMS_CUBED;
  for(size_t __i = 0; __i < max; ++__i)
    t.tensor[__i] *= i;
}

inline void g_tensor3TimesEqDbl(tensor3_t& t, const double& d) {
  size_t max = SPACE_DIMS_CUBED;
  for(size_t __i = 0; __i < max; ++__i)
    t.tensor[__i] *= d;
}

inline void g_tensor3TimesDbl(tensor3_t& t, tensor3_t& t2, const double& d) {
  size_t max = SPACE_DIMS_CUBED;
  for(size_t __i = 0; __i < max; ++__i)
    t.tensor[__i] = d*t2.tensor[__i];
}

inline void g_tensor3PlusEqTensor3(tensor3_t& t, tensor3_t& t2) {
  size_t max = SPACE_DIMS_CUBED;
  for(size_t __i = 0; __i < max; ++__i)
    t.tensor[__i] += t2.tensor[__i];
}

inline void g_tensor3PlusTensor3(tensor3_t& t, tensor3_t& t2, tensor3_t& t3) {
  size_t max = SPACE_DIMS_CUBED;
  for(size_t __i = 0; __i < max; ++__i)
    t.tensor[__i] = t2.tensor[__i] + t3.tensor[__i];
}

inline void g_tensor3DotPoint(tensor_t& t, tensor3_t& t2, point_t& p3) {
  size_t m = 0;
  size_t n = 0;
  for (size_t i = 0; i < SPACE_DIMS; ++i) {
    for (size_t j = 0; j < SPACE_DIMS; ++j) {
      t.tensor[j+i*SPACE_DIMS] = 0;
      for (size_t k = 0; j < SPACE_DIMS; ++j) {
	t.tensor[n/*j+i*SPACE_DIMS */] 
	  += p3[k]*t2.tensor[m/*k+j*SPACE_DIMS+i*SPACE_DIMS_SQUARED*/];
	++m;
      }
      ++n;
    }
  }
}

template<typename T>
ostream& operator<<(ostream &out, const math_tensor3_t<T> t)
{
  out << "tensor3_t: ";
  for (int h = 0; h < SPACE_DIMS; h++) {
    out << "(";
    for (int i = 0; i < SPACE_DIMS; i++) {
      out << "(";
      for (int j = 0; j < SPACE_DIMS; j++) {
        out << t.tensor[j+i*SPACE_DIMS+h*SPACE_DIMS*SPACE_DIMS];
        if (j != SPACE_DIMS-1)
          out << ", ";
      }
      out << ")";
      if (i != SPACE_DIMS-1)
	out << ", ";
    }
  out << ")" << endl;
  }
	
	return out;
}


struct triangle_t {
	union {
		point_t corners[3];
		struct {
			point_t first, second, third;
		};
	};
	
	/* [] access */
    inline point_t &operator[](int i) {
        return corners[i];
    }
	
    inline const point_t &operator[](int i) const {
        return corners[i];
    }		
};

#endif
