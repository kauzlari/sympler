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
#ifdef WITH_ARRAY_TYPES
#include "MArray2D.h"
#ifdef HAVE_JAMA_JAMA_LU_H
#include "jama/jama_lu.h"
#endif
//empty consturctor
MArray2D::MArray2D() {
	p_array2d=new TNT::Array2D<double>;
	if (!p_array2d) {
		throw gError("MArray2D","Allocation error!");
	}
}

//normal constructor, init importing a file
MArray2D::MArray2D(string filename) {
	p_array2d=new TNT::Array2D<double>;
	if (!p_array2d) {
		throw gError("MArray2D","Allocation error!");
	}
	path.assign(filename);
	p_array2d=read_MMheader(path);

}

//copy constructor
MArray2D::MArray2D(const MArray2D &o) {
	if (!o.p_array2d) {
		p_array2d=new TNT::Array2D<double>();
	} else
		p_array2d=new TNT::Array2D<double>(*o.p_array2d);
	if (!p_array2d) {
		throw gError("MArray2D","Allocation error!");
	}
}

MArray2D &MArray2D::operator=(const MArray2D &copy_marray2d) {
	if (p_array2d)
		delete p_array2d;
	if (!copy_marray2d.p_array2d) {
		p_array2d=new TNT::Array2D<double>();
	} else
		p_array2d=new TNT::Array2D<double>(*copy_marray2d.p_array2d);
	if (!p_array2d) {
		throw gError("MArray2D","Allocation error!");
	}
	return *this;

}

//init with MArray2D (row,column)
MArray2D::MArray2D(int r, int c) {
	p_array2d=new TNT::Array2D<double>(r,c);
	if (!p_array2d) {
		throw gError("MArray2D","Allocation error!");
	}
}

//init with MArray2D (row,column,constant the elements are initialized to)
MArray2D::MArray2D(int r, int c, int z) {
	p_array2d=new TNT::Array2D<double>(r,c,z);
	if (!p_array2d) {
		throw gError("MArray2D","Allocation error!");
	}
}

//destructor
MArray2D::~MArray2D() {
	if(p_array2d)
		delete p_array2d;
}

TNT::Array2D< double>* read_MMheader(string path) {
	TNT::Array2D<double>* mat;
	char banner[MM_MAX_TOKEN_LENGTH];
	char mtx[MM_MAX_TOKEN_LENGTH];
	char crd[MM_MAX_TOKEN_LENGTH];
	char data_type[MM_MAX_TOKEN_LENGTH];
	char storage_scheme[MM_MAX_TOKEN_LENGTH];

	string comment;
	int column, row;
	ifstream mm_file(path.c_str()); //file handle for reading

	//mm_file.exceptions(std::ifstream::eofbit| std::ifstream::failbit| std::ifstream::badbit);

	if (!mm_file) //File not found?!
	{
		throw gError("MArray2D","The input Matrix Market file was not found");

	}

	//try {	

	mm_file>> banner>> mtx>> crd >> data_type >> storage_scheme;

	//cheking for the right file
	if (strcmp(banner, MMBanner)!=0)
		throw gError("read_MMheader","The matrix input file is not in Market Martix format!");
	if (strcmp(mtx, MM_MTX_STR)!=0)
		throw gError("read_MMheader","Not a matrix defined in the input file!Can not read!");
	if (strcmp(crd, MM_DENSE_STR)!=0)
		throw gError("read_MMheader","The matrix is not a dense matrix!Can not be read!");
	if (strcmp(data_type, MM_REAL_STR)!=0)
		throw gError("read_MMheader","The matrix is not a real matrix!Can not be read!");
	if (strcmp(storage_scheme, MM_GENERAL_STR)!=0)
		throw gError("read_MMheader","The matrix is not a general one!Can not be read!");

	//continue reading of the file
	getline(mm_file, comment);
	getline(mm_file, comment);

	//get rows and columns
	mm_file>>row;
	mm_file>>column;

	mat=new TNT::Array2D<double>(row,column);

	while (!mm_file.eof()) {

		for (int i=0; i<column; i++) {
			for (int k=0; k<row; k++) {

				mm_file>>(*mat)[k][i];
			}
		}
	}

	//}
	//	catch (std::ifstream::failure e) {
	//	throw gError("read_MMheader","something wrong");
	//		}

	mm_file.close();
	return mat;
}

//Functions for exporting a row or a column from the matrix as a MArray2D "vector"
MArray2D &MArray2D::row(int i) {

	int row_num=(p_array2d)->dim1();//dim1 is the number of rows in the given matrix
	int column_num=(p_array2d)->dim2();//dim2 is the number of columns in the given matrix
	if (i>row_num)
		throw gError("row","Row number out of matrix boundaries!");

	MArray2D* temp=new MArray2D(1, column_num);

	for (int j=0; j<column_num; j++) {
		(*temp->p_array2d)[0][j]=(*p_array2d)[i][j];

	}
	return *temp;
}



MArray2D operator-(const MArray2D &m1, const MArray2D &m2) {
	int d1_m1=(m1.p_array2d)->dim1();
	int d2_m1=(m1.p_array2d)->dim2();

	if ((m2.p_array2d)->dim1() != d1_m1 || (m2.p_array2d)->dim2() != d2_m1)
		throw gError("operator+","Matrix dimensions not consistent!");

	else {
		MArray2D temp(d1_m1, d2_m1);

		for (int i=0; i<d1_m1; i++) {
			for (int j=0; j<d2_m1; j++)
				(*temp.p_array2d)[i][j] = (*m1.p_array2d)[i][j] - (*m2.p_array2d)[i][j];
		}
		return temp;
	}
}
MArray2D MArray2D::scalarmult(double s) {

	int row_num=(p_array2d)->dim1();//dim1-number of rows
	int column_num=(p_array2d)->dim2();//dim2- number of columns

	MArray2D temp(row_num, column_num);
	for (int j=0; j<row_num; j++) {
		for (int i=0; i<column_num; i++)
			(*temp.p_array2d)[j][i]=(*p_array2d)[j][i]*s;
	}
	return temp;
}


MArray2D matmult(const MArray2D &m1, const MArray2D &m2) {
	{
		if ((m1.p_array2d)->dim2() != (m2.p_array2d)->dim1())
			throw gError("matmult","Dimensions not consistent!dim1 mat1 differs from dim2 mat2!");

		int d1_m1 = (m1.p_array2d)->dim1();
		int d2_m1 = (m1.p_array2d)->dim2();
		int d2_m2 = (m2.p_array2d)->dim2();

		MArray2D temp(d1_m1, d2_m2);

		for (int i=0; i<d1_m1; i++)
			for (int j=0; j<d2_m2; j++) {
				double sum = 0;

				for (int k=0; k<d2_m1; k++)
					sum += (*m1.p_array2d)[i][k] * (*m2.p_array2d)[k][j];

				(*temp.p_array2d)[i][j] = sum;
			}

		return temp;
	}
}

///*matrix mal vector = vector*/
//MArray2D matvecmult(const MArray2D &m1, const MArray2D &m2) {
//	{
//		if ((m1.p_array2d)->dim2() != (m2.p_array2d)->dim1())
//			throw gError("matmult","Dimensions not consistent!dim1 mat1 differs from dim2 mat2!");
//
//		int d1_m1 = (m1.p_array2d)->dim1();
//		int d2_m1 = (m1.p_array2d)->dim2();
//		int d2_m2 = (m2.p_array2d)->dim2();
//
//		MArray2D temp(d1_m1, d2_m2);
//
//		for (int i=0; i<d1_m1; i++)
//			for (int j=0; j<d2_m2; j++) {
//				double sum = 0;
//
//				for (int k=0; k<d2_m1; k++)
//					sum += (*m1.p_array2d)[i][k] * (*m2.p_array2d)[k][j];
//
//				(*temp.p_array2d)[i][j] = sum;
//			}
//
//		return temp;
//	}
//}

std::ostream& operator<<(std::ostream &s, MArray2D &m) {
	int d1=(m.p_array2d)->dim1();
	int d2=(m.p_array2d)->dim2();

	for (int i=0; i<d1; i++) {
		for (int j=0; j<d2; j++) {
			s << m[i][j] << " ";
		}
		s << "\n";
	}
	return s;
}

MArray2D operator+(const MArray2D &m1, const MArray2D &m2) {
	int d1_m1=(m1.p_array2d)->dim1();
	int d2_m1=(m1.p_array2d)->dim2();

	if ((m2.p_array2d)->dim1() != d1_m1 || (m2.p_array2d)->dim2() != d2_m1)
		throw gError("operator+","Matrix dimensions not consistent!");

	else {
		MArray2D temp(d1_m1, d2_m1);

		for (int i=0; i<d1_m1; i++) {
			for (int j=0; j<d2_m1; j++)
				(*temp.p_array2d)[i][j] = (*m1.p_array2d)[i][j] + (*m2.p_array2d)[i][j];
		}
		return temp;
	}
}

MArray2D &MArray2D::column(int i) {

	int row_num=(p_array2d)->dim1();//dim1 is the number of rows
	int col_num=(p_array2d)->dim2();//dim2 is the number of columns

	if (i>col_num)
		throw gError("column","Column number out of matrix boundaries!");

	MArray2D* temp=new MArray2D(row_num, 1);

	for (int j=0; j<row_num; j++) {
		(*temp->p_array2d)[j][0]=(*p_array2d)[j][i];

	}
	return *temp;
}

#ifdef HAVE_JAMA_JAMA_LU_H
//function for inverting a matrix
//return a TNT Array, so (MArray2D.p_array2d)=&(invert(m_matrixZ));
MArray2D &invert(const MArray2D &m1) {

	//create identity matrix
	size_t d1 = (m1.p_array2d)->dim1();
	size_t d2 = (m1.p_array2d)->dim2();
	if (d1 != d2)
		throw gError("MArray2D::"+PRETTYFUNC_INFO,
				"This function can only invert square matrices, matrix"
				"dimensions are: " + ObjToString(d1) + ", " + ObjToString(d2) + ".");
	MArray2D matrixU(d1, d1, 0);
	for (size_t k=0; k < d1; k++) {
		matrixU[k][k]=1;
	}
	JAMA::LU<double> lu(*m1.p_array2d);
	// solves A * A_inv = Identity
	MArray2D* result = new MArray2D(d1,d1);
	result->p_array2d = new TNT::Array2D<double>(lu.solve(*(matrixU.p_array2d)));
	return *result;
}

std::string MArray2D::toMTXString(const std::string &eol) {
	int d1 = dim1(), d2 = dim2();
	std::ostringstream outstr;
	outstr << string(MMBanner) + " " + MM_MTX_STR + " " + 
			MM_DENSE_STR + " " + MM_REAL_STR + " " + MM_GENERAL_STR + eol;
	outstr << d1 << " " << d2 << eol;
	for (int j = 0; j < d2; ++j)
		for (int i = 0; i < d1; ++i)
			outstr << (*p_array2d)[i][j] << eol;
	return outstr.str();
}

#endif /*HAVE_JAMA_JAMA_LU_H*/
#endif /*WITH_ARRAY_TYPES*/
