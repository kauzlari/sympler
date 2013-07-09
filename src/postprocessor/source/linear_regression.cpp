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



#include <math.h>

#include "simulation.h"
#include "linear_regression.h"

#define IDX_SLOPE 0
//#define IDX_SLOPE_ERROR 1
#define IDX_INTERCEPT 1
//#define IDX_INTERCEPT_ERROR 3
#define IDX_VALUE 2

/* Register this Postprocessor with the factory. */
const Postprocessor_Register<LinearRegression> linear_regression("LinearRegression");



//---- Constructors/Destructor ----

LinearRegression::LinearRegression(Node *parent, Simulation *simulation)
	: Postprocessor(parent, simulation), m_idx_x(0), m_idx_y(1)
{
    init();
}


LinearRegression::~LinearRegression()
{
}



//---- Methods ----

void LinearRegression::init()
{
  /* Properties */
  m_properties.setClassName("LinearRegression");
    
  m_properties.setDescription(
    "Performes a linear regression. Output columns are "
    "'slope', 'intercept' and 'value'. 'value' is 'slope' with *postfn* applied."
  );

  STRINGPC
    (x, m_xcol,
     "Name of the column containing the x-value.");
  STRINGPC
    (y, m_ycol,
     "Name of the column containing the y-values.");
  STRINGPC
    (prefn, m_str_prefn,
     "Function applied before the regression.");
  STRINGPC
    (postfn, m_str_postfn,
     "Function applied after the regression (x-value is the slope the the regression).");
    
  m_xcol = m_ycol = "---";
  m_str_prefn = m_str_postfn = "x";

  /* Format for the output pipe */
  m_output.addAttribute("slope", DataFormat::DOUBLE);
//    m_output.addAttribute("slope_error", DataFormat::DOUBLE);
  m_output.addAttribute("intercept", DataFormat::DOUBLE);
//    m_output.addAttribute("intercept_error", DataFormat::DOUBLE);
  m_output.addAttribute("value", DataFormat::DOUBLE);
}


void LinearRegression::setup()
{
  Postprocessor::setup();

  if (m_input_format->attrExists(m_xcol))
    m_idx_x = m_input_format->indexOf(m_xcol, DataFormat::DOUBLE);

  if (m_input_format->attrExists(m_ycol))
    m_idx_y = m_input_format->indexOf(m_ycol, DataFormat::DOUBLE);
}


void LinearRegression::describeInput(DataFormat *input_format)
{
    m_input_format = input_format;
}


void LinearRegression::push(data_sp data)
{
    m_data.push_back(data);
}


void LinearRegression::flush()
{
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0, sxx, syy, sxy;
    double slope, slope_error, intercept;
    int n = 0;

	/* Initialize our post- and preprocessor functions */
	for (int i = 0; i < SPACE_DIMS; i++) {
		string s = "box" + string(1, 'X'+i);
		m_prefn.constant(s) = m_simulation->phase()->boundary()->boundingBox().size()[i];
		m_postfn.constant(s) = m_simulation->phase()->boundary()->boundingBox().size()[i];
	}

	m_prefn.FromString(m_str_prefn);
	m_postfn.FromString(m_str_postfn);

	/* Apply m_prefn and do the linear regression */
	for (list<data_sp>::iterator i = m_data.begin(); i != m_data.end(); i++) {
        double x, y;

        x = (*i)->doubleByIndex(m_idx_x);
        y = m_prefn.value((*i)->doubleByIndex(m_idx_y));
        
        sum_x += x;
        sum_y += y;
        sum_xy += x*y;
        sum_x2 += x*x;
        sum_y2 += y*y;
        n++;
    }

    if(n < 3) 
       throw gError
         ("LinearRegression::flush", "I need at least 3 values, but have only " 
           + ObjToString(n));
    sxy = sum_xy - sum_x*sum_y/n;
    sxx = sum_x2 - sum_x*sum_x/n;
    syy = sum_y2 - sum_y*sum_y/n;

    slope = sxy/sxx;
    slope_error = sqrt((syy - (slope*slope*sxx)) / (n - 2));
    intercept = sum_y/n - slope*sum_x/n;

    data_sp data = m_output.newData();
    data->doubleByIndex(IDX_SLOPE) = slope;
//    data->doubleByIndex(IDX_SLOPE_ERROR) = slope_error;
    data->doubleByIndex(IDX_INTERCEPT) = intercept;
//    data->doubleByIndex(IDX_INTERCEPT_ERROR) = 0;
    data->doubleByIndex(IDX_VALUE) = m_postfn.value(slope);
    distribute(data);
}


