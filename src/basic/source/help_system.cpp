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



#include <string>

#include "help_system.h"

#include "help.h"
#include "simulation.h"
#include "grid_meter.h"
#include "reflector.h"
#include "particle_creator.h"
#include "pair_creator.h"

#define HELP_CLASS(name, factory_class, cmd) \
    class name: public Help<factory_class> { \
    protected: \
        virtual Node *instantiate(const factory_class &factory) const { \
            return factory.cmd; \
        } \
    }


HELP_CLASS(Force_Help, GenFType, makeGenF(NULL));
HELP_CLASS(Boundary_Help, Boundary_Factory, instantiate(NULL));
HELP_CLASS(Reflector_Help, Reflector_Factory, instantiate(NULL));
HELP_CLASS(ParticleCreator_Help, ParticleCreator_Factory, instantiate(NULL));
HELP_CLASS(Meter_Help, Meter_Factory, instantiate(NULL));
HELP_CLASS(Postprocessor_Help, Postprocessor_Factory, instantiate(NULL, NULL));
HELP_CLASS(Callable_Help, Callable_Factory, instantiate(NULL));
HELP_CLASS(GridMeter_Help, GridMeter_Factory, instantiate(NULL));
HELP_CLASS(Integrator_Help, Integrator_Factory, instantiate(NULL));
HELP_CLASS(Symbol_Help, SymbolFactory, instantiate(NULL));
HELP_CLASS(WeightingFunction_Help, WeightingFunction_Factory, instantiate(NULL));
HELP_CLASS(PairCreator_Help, PairCreator_Factory, instantiate(NULL));



void syntax()
{
    cout << "Usage: sympler input-file\n\nOr: sympler --help\nOr: sympler --help all\nOr: sympler --help workflow\nOr: sympler --help expressions\n";
}


string main_help_str
    ("Usage: sympler INPUTFILE\n"
     "\n"
     "This screen:                  sympler --help\n"
     "Module description:           sympler --help [module]\n"
     "              e.g.:           sympler --help Meters\n"
     "Help for all modules:         sympler --help all\n"
     "Workflow description:         sympler --help workflow\n"
     "Help on symbolic expressions: sympler --help expressions\n"
     "\n"
     "The file given as INPUTFILE has to be of the following format: Objects are to be defined "
     "between the paranthesis '<' and '>' followed by attribute definitions.\n\n"
     "Example:\n"
     "  <Object\n"
     "      attribute1 = \"4*sin(x)\"\n"
     "      attribute2 = \"12.3\"\n"
     "  >\n"
     "      ... child objects ...\n"
     "  </Object>\n\n"
     "The body of any object can contain child objects which usually take and process data from "
     "the parent object or provide additional information for the parent object. \n"
     "The objects have the following parent-child structure:\n\n"
     " ------------  \n"
     "| Simulation | \n"
     " ------------  \n"
     "  | \n"
     "  |   --------------------  \n"
     "  |--| WeightingFunctions | \n"
     "  |   --------------------  \n"
     "  | \n"
     "  |   ---------  \n"
     "  |--| Symbols | \n"
     "  |   ---------  \n"
     "  | \n"
     "  |   ------------  \n"
     "  |--| Controller | \n"
     "  |   ------------  \n"
     "  |    |\n"
     "  |    |   -------------  \n"
     "  |    |--| Integrators | \n"
     "  |        -------------  \n"
     "  | \n"
     "  |   --------  \n"
     "  |--| Forces | \n"
     "  |   --------  \n"
     "  | \n"
     "  |   -----------  \n"
     "  |--| Callables | \n"
     "  |   -----------  \n"
     "  | \n"
     "  |   -------  \n"
     "  |--| Phase | \n"
     "  |   -------  \n"
     "  |    | \n"
     "  |    |   --------------  \n"
     "  |    |--| PairCreators | \n"
     "  |    |   --------------  \n"
     "  |    | \n"
     "  |    |   ------------  \n"
     "  |    |--| Boundaries | \n"
     "  |        ------------  \n"
     "  |         | \n"
     "  |         |   ------------------  \n"
     "  |         |--| ParticleCreators | \n"
     "  |         |   ------------------  \n"
     "  |         | \n"
     "  |         |   ------------  \n"
     "  |         |--| Reflectors | \n"
     "  |             ------------  \n"
     "  | \n"
     "  |   --------  \n"
     "  |--| Meters | \n"
     "      -------   \n"
     "       | \n"
     "       |   ------------ \n"
     "       |--| GridMeters | \n"
     "       |   ------------  \n"
     "       | \n"
     "       |   ----------------  \n"
     "       |--| Postprocessors | \n"
     "           ----------------  \n\n"
     "This is the general object structure. For certain special classes there are additional "
     "child objects. Those are mentioned in the help screen for the specific object.\n"
     "Please note that forces and callables have to be declared before and meters "
     "after the phase.\n"
     "Type 'sympler --help [module]', where [module] is of the ones given in the tree above, to get "
     "further information.\n"
     "Type 'sympler --help all', to get information on all modules at once.\n"
     "Type 'sympler --help expressions', for learning about accepted algebraic expressions "
     "to be runtime compiled.\n"
     "Type 'sympler --help workflow', for more information about the order in which the main operations are performed.\n\n"
     );

string flow_help_str
    ("It follows a list of basic computations that are performed, in the order given, from initialisation to end. This may help for example to figure out the reason for unexpected behaviour, for example if particle-values are unexpectedly zero:\n\n"
     " - set Forces and unprotected Symbols = 0 (protected Symbols are those declared by an Integrator.)\n"
     " - calculate neighbour-list\n"
     " - calculate Symbols of stage 1 (see help option --help Symbols for setting stages)\n"
     " - compute Forces\n\n"
     " - --- START OF MAIN LOOP OVER TIMESTEPS ------\n\n"
     " - calculate Symbols of stage 0 (see help option --help Symbols for setting stages)\n"
     " - call Meters (see help option --help Meters)\n"
     " - step1 of time integration, for example particle positions, predictor-step for velocity or other quantity\n"
     " - set out-dated Forces and unprotected Symbols = 0 (protected Symbols are those declared by an Integrator.)\n"
     " - calculate neighbour-list\n"
     " - calculate Symbols of stage 1 (see help option --help Symbols for setting stages)\n"
     " - compute Forces\n"
     " - step2 of time integration, for example corrector step for final value of velocity or other quantity\n"
     " - run Callables (see help option --help Callables)\n\n"
     " - --- END OF MAIN LOOP OVER TIMESTEPS ------\n\n"
     " - calculate Symbols of stage 0 (see help option --help Symbols for setting stages)\n"
     " - final call of Meters (see help option --help Meters)\n"
     " - END of program\n"
     "\n\n");


string expr_help_str
    (
     "Many modules understand algebraic expressions which are then runtime compiled.\n "
     "\n"
     " * \"Vector\" is a vector with three components\n"
     " * \"Matrix\" is a 3x3 matrix\n"
     "\n"
     " * Notation in the input file:\n"
     "  * \"aij\" is a scalar. The indices \"i\", \"j\" indicate that it is stored for each PAIR of particles \"i\" and \"j\". If the variable is to be introduced by the user, s-he will just write 'symbol = \"a\"' in the corresponding module, e.g., in a <PairScalar/>\n"
     "  * \"[ai]\" is a \"Vector\" as defined above. The index \"i\", indicates that it appears in an expression for a pair of particles \"i\" and \"j\" and belongs to particle \"i\". It is stored for each PARTICLE. If the variable is to be introduced by the user, s-he will just write 'symbol = \"a\"' in the corresponding module, e.g., in a <ParticleVector/>\n"
"  * \"{aj}\" is a \"Matrix\" as defined above. The index \"j\", indicates that it appears in an expression for a pair of particles \"i\" and \"j\" and belongs to particle \"j\". It is stored for each PARTICLE. If the variable is to be introduced by the user, s-he will just write 'symbol = \"a\"' in the corresponding module, e.g., in a <ParticleTensor/>\n"
     "  * If the symbol appears in a particle-based expression of modules such as <ParticleScalar/>, <ParticleVector/>, <ParticleTensor/>, then, indices \"i\" or \"j\" are not used, because it is clear that the symbol in the expression belongs to the very same particle as the output symbol of the module. So, we simply write \"a\", \"[a]\", \"{a}\".\n"
     "\n"
     " * Symbols which are always known to\n"
     "  * particle-expressions:\n"
     "   * i: particle index\n"
     "   * [r]: particle position\n"
     "   * [v]: particle velocity\n"
     "  * pair-expressions:\n"
     "   * i: index of first particle\n"
     "   * [ri]: position of first particle\n"
     "   * [vi]: velocity of first particle\n"
     "   * j: index of second particle\n"
     "   * [rj]: position of second particle\n"
     "   * [vj]: velocity of second particle\n"
     "   * rij: scalar absolute pair-distance\n"
     "   * [rij]: distance vector, computed as [ri]-[rj]\n"
     "\n"
     " * Operators:\n"
     "   * \":\"\n"
     "        * Vector:Vector = Scalar (scalar product), e.g. \"[ai]:[bj]\"\n"
     "        * Matrix:Matrix = Scalar (contraction on tensors), e.g. \"{ai}:{bj}\"\n"
     "        * Tensor:Vector = Vector (Matrix-Vector product), e.g. \"{ai}:[bj]\"\n"
     "   * \"@\"\n"
     "        * Vector@Vector = Tensor (outer product), e.g. \"[ai]@[bj]\"\n"
     "   * \"^\" = exponentiation, e.g., \"a^b\"\n"
     "   * \"*\" = multiplication, component wise on Matrix and Vector, e.g., \"a*b\" or also \"{ai}*{bj}\"\n"
     "   * \"�\"\n"
     "       * Matrix�Matrix = Matrix (standard matrix product), e.g. \"{a}�{b}\"\n"
     "         Clarification in case the operator does not display correctly\n"
     "         in your editor or console: The operator corresponds do the\n"
     "         \"degrees\" symbol, i.e., 90� = \"90 degrees\".\n"
     "\n"
     " * Functions of a scalar (or scalar-components):\n"
     "   * abs(Scalar): self-explanatory, component-wise for non-scalars\n"
     "   * atan(Scalar): self-explanatory, component-wise for non-scalars\n"
     "   * cos(Scalar): self-explanatory, component-wise for non-scalars, e.g., \"cos({aj})\"\n"
     "   * hcos(Scalar): (cosh) self-explanatory, component-wise for non-scalars\n"
     "   * exp(Scalar): self-explanatory, component-wise for non-scalars\n"
     "   * idMat(Scalar) = Scalar times identity-matrix, e.g., \"idMat(ai)\"\n"
     "   * idVec(Scalar) = Vector: [Scalar, Scalar, Scalar], e.g., \"idVec(a)\"\n"
     "   * round(Scalar): rounds to nearest integer, component-wise for non-scalars\n"
     "   * sin(Scalar): self-explanatory, component-wise for non-scalars, e.g., \"sin(a)\"\n"
     "   * hsin(Scalar): (sinh) self-explanatory, component-wise for non-scalars\n"
     "   * step(Scalar): unit step function with the following definition: 1 for x>0 and 0 otherwise, component-wise for non-scalars\n"
     "   * stpVal(Scalar): step function with the following definition: x for x>0 and 0 otherwise, component-wise for non-scalars\n"
     "   * sqrt(Scalar): square root of a positive real number, component-wise for non-scalars\n"
     "   * tan(Scalar): self-explanatory, component-wise for non-scalars\n"
     "   * htan(Scalar): (tanh) self-explanatory, component-wise for non-scalars\n"   
     "   * unitMat(Scalar) = Matrix: Scalar*([1,1,1]@[1,1,1]), e.g., \"unitMat(a)\"\n"
     "   * uran(Scalar) = Uniform random number in the interval [0,1]; the argument is not used at the moment, e.g., \"uran(1)\"\n"
     "   * uVecX|Y|Z(Scalar) = [Scalar, 0, 0]|[0,Scalar,0]|[0,0,Scalar], e.g., \"uVecY(a)\"\n"
     "\n"
     " * Functions of a vector:\n"
     "   * diagMat(Vector) = Matrix: Matrix with vector components in diagonal, e.g., \"diagMat([a])\"\n"
     "   * x|y|zCoord(Vector) = Scalar: self-explanatory, e.g., \"zCoord([a])\"\n"
     "\n"
     " * Functions of a matrix:\n"
     "   * det(Matrix) = Scalar: determinant of Matrix, e.g., \"det({ai})\"\n"
     "   * Q(Matrix) = Scalar: Matrix:Matrix, e.g., \"Q({ai})\"\n"
     "   * T(Matrix) = Matrix: transpose of Matrix, e.g., \"T({ai})\"\n"
     "   * trace(Matrix) = Scalar: self-explanatory, e.g., \"trace({ai})\"\n"
     "   * xyMat(Matrix) = Matrix: zeroes the third row and column, e.g. \"xyMat({ai})\"\n" 
     "\n\n");


#define FOUND_HELP(Factory, Helper) \
    node = new HelpNode(parent, Factory::description()); \
    Helper().asHelpNode(node)


HelpNode *find_help(string keyword, HelpNode *parent = NULL)
{
    HelpNode *node = NULL;

    if (keyword == "Simulation") {
        node = new HelpNode(parent, "");
        Simulation().help(node);
    } else if (keyword == "Controller") {
        node = new HelpNode(parent, "");
        Controller(NULL).help(node);
    } else if (keyword == "WeightingFunctions") {
      FOUND_HELP(WeightingFunction_Factory, WeightingFunction_Help);
    } else if (keyword == "Symbols") {
      FOUND_HELP(SymbolFactory, Symbol_Help);
    } else if (keyword == "Integrators") {
      FOUND_HELP(Integrator_Factory, Integrator_Help);
    } else if (keyword == "Phase") {
        node = new HelpNode(parent, "");
        Phase(NULL).help(node);
    } else if (keyword == "Forces") {
        FOUND_HELP(GenFType, Force_Help);
    } else if (keyword == "Boundaries") {
        FOUND_HELP(Boundary_Factory, Boundary_Help);
    } else if (keyword == "Reflectors") {
        FOUND_HELP(Reflector_Factory, Reflector_Help);
    } else if (keyword == "ParticleCreators") {
      FOUND_HELP(ParticleCreator_Factory, ParticleCreator_Help);
    } else if (keyword == "PairCreators") {
      FOUND_HELP(PairCreator_Factory, PairCreator_Help);
    } else if (keyword == "Meters") {
        FOUND_HELP(Meter_Factory, Meter_Help);
    } else if (keyword == "Postprocessors") {
        FOUND_HELP(Postprocessor_Factory, Postprocessor_Help);
    } else if (keyword == "Callables") {
        FOUND_HELP(Callable_Factory, Callable_Help);
    } else if (keyword == "GridMeters") {
        FOUND_HELP(GridMeter_Factory, GridMeter_Help);
    }

    return node;
}


const char *keywords[] =
    { "Simulation", "Controller", "WeightingFunctions", "Symbols", "Integrators", "Phase", 
    "Forces", "Callables",
    "Boundaries", "ParticleCreators", "PairCreators", "Reflectors", "Meters", "GridMeters", 
      "Postprocessors",
      NULL };


HelpNode *all_help()
{
    HelpNode *node = new HelpNode(NULL, main_help_str);

    for (int i = 0; keywords[i]; i++) {
        HelpNode *n = new HelpNode(node, keywords[i]);
        find_help(keywords[i], n);
    }

    return node;
}


void help(string keyword)
{
    HelpNode *node;
    int h;

    if (keyword == "all") {
        node = all_help();
        h = 0;
    } else if (keyword == "workflow") {
      node = new HelpNode(NULL, flow_help_str);
      h = 0;
    } else if (keyword == "expressions") {
      node = new HelpNode(NULL, expr_help_str);
      h = 0;
    } else {
        node = find_help(keyword);
        h = 2;
    }

    if (node) {
        HelpFormatScreen(cout, *node, h);

        delete node;
    } else {
        cout << "help(" << keyword  << "): Unknown keyword: " << keyword << endl
            << "Possibilities are:" << endl;
        size_t i = 0;
        while(keywords[i])
        {
          if(i!=0) cout << ", ";
          cout << keywords[i];
          ++i;
        }
        cout << endl;
	cout << "Or one of \"all\", \"workflow\", \"expressions\".";
        cout << endl;
    }
}


void main_help()
{
    cout << block(main_help_str, 0, 80) << endl;
}

