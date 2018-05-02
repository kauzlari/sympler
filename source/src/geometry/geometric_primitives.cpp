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



#include "general.h"
#include "geometric_primitives.h"

double g_geom_eps = 1e-10;


//     /*inline*/ bool cuboid_t::isInside(const point_t &pos) const {
//       for (int i = 0; (i < SPACE_DIMS); i++) {
//         if ((pos[i] < corner1[i]) || (pos[i] >= corner2[i]))
//         {		
//           if(pos[i] < corner1[i])
//             MSG_DEBUG("cuboid_t::isInside", "<: pos = " << pos[i] << ", corner = " << corner1[i] << ", diff = " << pos[i]-corner1[i]);
//           if(pos[i] >= corner2[i])
//             MSG_DEBUG("cuboid_t::isInside", ">=: pos = " << pos[i] << ", corner = " << corner2[i] << ", diff = " << pos[i]-corner2[i]);
//           return false;
//         }
//       }
//         
//       return true;
//     }
// 
//     /*inline*/ bool cuboid_t::isInsideEps(const point_t &pos, double eps) const {
//       for (int i = 0; (i < SPACE_DIMS); i++) {
//         if ((pos[i] < corner1[i]-eps) || (pos[i] >= corner2[i]+eps))
//           return false;
//       }
//         
//       return true;
//     }

// template<typename T>
// ostream& operator<<(ostream &out, const math_tensor_t<T> t) 
// {
//   out << "(";
//   for (int i = 0; i < SPACE_DIMS; i++) {
//     out << "(";
//       for (int j = 0; j < SPACE_DIMS; j++) {
//         out << t(i, j);
//         if (j != SPACE_DIMS-1)
//           out << ", ";
//       }
//     out << ")";
//     if (i != SPACE_DIMS-1)
//       out << ", ";
//   }
//   out << ")";
	
// 	return out;
// }


// template<typename T>
// ostream& operator<<(ostream &out, const math_tensor3_t<T> t) {
// // ostream& operator<<(ostream &out, const tensor3_t t) {
// //   out << "tensor3_t: "; 
//   for (int h = 0; h < SPACE_DIMS; h++) {
//     out << "(";
//     for (int i = 0; i < SPACE_DIMS; i++) {
//       out << "(";
//       for (int j = 0; j < SPACE_DIMS; j++) {
//         out << t.tensor[j+i*SPACE_DIMS+h*SPACE_DIMS*SPACE_DIMS];
//         if (j != SPACE_DIMS-1)
//           out << ", ";
//       }
//       out << ")";
//       if (i != SPACE_DIMS-1)
// 	out << ", ";
//     }
//     out << ")" ;
//     if (h != SPACE_DIMS-1)
//       out << ", ";
//     out << endl;
//   }
	
// 	return out;
// }
