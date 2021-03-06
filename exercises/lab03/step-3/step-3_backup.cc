/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */



#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>
#include<deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

#include <math.h>

//to define functions
#include <deal.II/base/function_parser.h>

//for dynamical mesh refinement
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/error_estimator.h>


#define CYCLES 8
using namespace dealii;

//I modify this file along the following lines:
//1. the RHS of the EOM becomes sin(2*pi*x)sin(6*pi*y)
//2. estimate the error of the numerical solution w.r.t. to exact one with varios norms

//Moreover, I decided to implement by myself a quadrature-weighted L2_norm
//and I compared my result with the one from the built-in function "integrate_difference"

//all the norms are computed in the function error_norms()
class Step3
{
public:
  Step3 ();

  void run ();


private:
  void refine_grid ();
  void setup_system ();
  void assemble_system (const FunctionParser<2>& Dfp);
  void solve (const FunctionParser<2>& fp);
  void output_results (const unisgned int cycle) const;
  void error_norms (const FunctionParser<2>& fp) const;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;

  AffineConstraints<double> constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
  Vector<double>       exact_sol;
};


Step3::Step3 ()
  :
  fe (1),
  dof_handler (triangulation)
{}



void Step3::refine_grid ()
{
  Vector<float> error_per_cell(triangulation.n_active_cells());

  KellyErrorEstimator<2>::estimate
	(
	dof_handler,QGauss<2>(2),
	std::map<types::boundary_id,
	const Function<2>*>(),solution,
	error_per_cell
	);

  GridRefinement::refine_and_coarsen_fixed_number
	(triangulation,error_per_cell,0.3,0.03);
  triangulation.execute_coarsening_and_refinement();
}




void Step3::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

  //re-enumerate the dofs
  DoFRenumbering::Cuthill_McKee(dof_handler);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  //clear the constraints and reinit them
  constarints.clear();
  DoFTolls::make_hanging_node_constraints(dof_handler,constraimts);
  VectorTools::interpolate_boundary_values(dof_handler,0,
		Function::ZeroFunction<2>(),constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit (sparsity_pattern);
  
}



void Step3::assemble_system (const FunctionParser<2>& Dfp)
{
  QGauss<2>  quadrature_formula(2);
  //update_quadrature_points to use them as function arguments
  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_gradients | update_JxW_values | update_quadrature_points);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        {
	//map the local Q_point
	const auto& xx = fe_values.quadrature_point(q_index);
  
	for (unsigned int i=0; i<dofs_per_cell; ++i){
            cell_rhs(i) += (fe_values.shape_value (i, q_index) *
                            Dfp.value(xx)*
                            fe_values.JxW (q_index));
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                   fe_values.shape_grad (j, q_index) *
                                   fe_values.JxW (q_index));
	}
      }
   }

      cell->get_dof_indices (local_dof_indices);

  constraints.distribute_local_to_global(cell_matrix, cell_rhs,
			local_dof_indices,system_matrix,system_rhs);
}



void Step3::solve (const FunctionParser<2>& fp)
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);

  solver.solve (system_matrix, solution, system_rhs,
                PreconditionIdentity());

  //interpolate exact solution
  exact_sol.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,fp,exact_sol);

  constraints.distribute(solution);

}



void Step3::output_results (const unsigned int cycle) const
{
  if(cycle!=CYCLES) return;

  //plot numerical and exact solutions
  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.add_data_vector (exact_sol, "exact_sol");

  data_out.build_patches ();

  std::ofstream output ("solution.vtk");
  data_out.write_vtk (output);

}

void Step3::error_norms (const FunctionParser<2>& fp) const {
//compute error norms
  Vector<double>  norm_vec(exact_sol);
  norm_vec-= solution; //sol - exact_sol
  std::cout << "L_infty norm: " << (norm_vec).linfty_norm() << std::endl;
  double norm_factor =  1./dof_handler.n_dofs(); //normalization for Lp_norms
  std::cout << "Normalized L_1 norm: " << (norm_vec).l1_norm()*norm_factor << std::endl;
  std::cout << "Normalized L_2 norm: " << (norm_vec).l2_norm()*sqrt(norm_factor) << std::endl;
//--------my L2_norm -------
  //I want to estimate error via a Quadrature version of the L2_norm
  //There is a built in function that does this: integrate_difference
  //However, it is a good exercise to do it manually and compare the two results

  //First, set variables of general utility
  double my_L2_norm = 0.;
  QGauss<2>  quadrature_formula(2);
  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_JxW_values | update_quadrature_points);

  const unsigned int n_q_points = quadrature_formula.size();
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

  //Next, declare a vector to store the numerical solution
  //evaluated at quadrature points of all cells
  //Q_solution is later filled cell-by-cell via the get_function_values
  std::vector<double> Q_numerical(n_q_points);

  for (; cell!=endc; ++cell){

      fe_values.reinit (cell);
      //evaluate Q_solution at the Q_points of the cell
      fe_values.get_function_values(solution,Q_numerical);

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
                const auto& xx = fe_values.quadrature_point(q_index);
                auto Q_elem = Q_numerical[q_index]-fp.value(xx);
                my_L2_norm += Q_elem*Q_elem*fe_values.JxW(q_index);
      }

  }

  //I choose as the normalization factor the "volume" of the system
  norm_factor = 1./16;
  my_L2_norm = sqrt(my_L2_norm*norm_factor);
  std::cout << "My Quadrature_L2_norm: " << my_L2_norm << std::endl;

//-------Deal.II L2 norm ------
  Vector<double> norm_vec_2;
  norm_vec_2.reinit(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,solution,fp,norm_vec_2,quadrature_formula,VectorTools::NormType::L2_norm);
  std::cout << "Deal.II Quadrature_L2_norm: " << (norm_vec_2).l2_norm()*sqrt(norm_factor) << std::endl;
}

void Step3::run ()
{
  //initialize RHS function...
  std::string variables = "x,y";
  std::map<std::string,double> constants;
  constants["pi"] = numbers::PI;
  std::string expression = "sin(2*pi*x)*sin(6*pi*y)";
  FunctionParser<2> fp(1);
  fp.initialize(variables, expression, constants);
  
  //...and its (minus-the-)laplacian
  expression = "40*pi*pi*sin(2*pi*x)*sin(6*pi*y)";
  FunctionParser<2> Dfp(1);
  Dfp.initialize(variables, expression, constants);

  for(unsigned int cycle = 0; cycle < CYCLES; ++cycle){
        if(cycle==0){
                GridGenerator::hyper_cube (triangulation, -1, 1);
                triangulation.refine_global (1);
        }
        else
                refine_grid();

   setup_system ();
   assemble_system (Dfp);
   solve (fp);
   output_results (cycle);
   //error_norms(fp);
  }
  error_norms(fp);
}



int main ()
{
  deallog.depth_console (2);

  Step3 laplace_problem;
  laplace_problem.run ();

  return 0;
}
