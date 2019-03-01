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

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

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


using namespace dealii;



class Step3
{
public:
  Step3 ();

  void run ();


private:
  void make_grid ();
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results () const;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
  // exact solution
  //not here: here is private! unless you flag as mutable...
  //Vector<double>       exact_sol;
};


Step3::Step3 ()
  :
  fe (1),
  dof_handler (triangulation)
{}



void Step3::make_grid ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);

  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
}




void Step3::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
  
  //initialize the exact solution
  //exact_sol.reinit(dof_handler.n_dofs());
}



void Step3::assemble_system ()
{
  QGauss<2>  quadrature_formula(2);
  //update_quadrature_points to use as function arguments
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
	//map the quadrature points
	const auto& xx = fe_values.quadrature_point(q_index);
	//function on the RHS
	const double FF = 40*M_PI*M_PI*sin(2*M_PI*xx[0])*sin(6*M_PI*xx[1]);     
  
	for (unsigned int i=0; i<dofs_per_cell; ++i){
            cell_rhs(i) += (fe_values.shape_value (i, q_index) *
                            FF*
                            fe_values.JxW (q_index));
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                   fe_values.shape_grad (j, q_index) *
                                   fe_values.JxW (q_index));
	}
      }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i){
	system_rhs(local_dof_indices[i]) += cell_rhs(i);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));
      	}

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<2>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
  }
}



void Step3::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);

  solver.solve (system_matrix, solution, system_rhs,
                PreconditionIdentity());
}



void Step3::output_results () const
{

//compute exact solution
//initialize function
  std::string variables = "x,y";
  std::map<std::string,double> constants;
  constants["pi"] = numbers::PI;
  std::string expression = "sin(2*pi*x)*sin(6*pi*y)";
  FunctionParser<2> fp(1);
  fp.initialize(variables, expression, constants);

//plot
  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

//interpolate function
  Vector<double>       exact_sol;
  exact_sol.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,fp,exact_sol);
  data_out.add_data_vector (exact_sol, "exact_sol");
  
  data_out.build_patches ();

  std::ofstream output ("solution.vtk");
  data_out.write_vtk (output);

  Vector<double>  norm_vec(exact_sol);
  norm_vec-= solution;
  std::cout << "L_infty norm: " << (norm_vec).linfty_norm() << std::endl;
  const double norm_factor =  1./dof_handler.n_dofs();
  std::cout << "Normalized L_1 norm: " << (norm_vec).l1_norm()*norm_factor << std::endl;
  std::cout << "Normalized L_2 norm: " << (norm_vec).l2_norm()*sqrt(norm_factor) << std::endl;

//my L2_norm -------
  //I want to estimate error via a Quadrature_L2_norm
  //There is a built in function that does this: integrate_difference
  //However, it is a good exercise to do it manually and compare the two results
  //First of all, I set variables of general utility
  double my_L2_norm = 0.;
  QGauss<2>  quadrature_formula(2);
  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_JxW_values | update_quadrature_points);

  const unsigned int n_q_points = quadrature_formula.size();
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

  //Next, declare a vector to store the value of the numerical solution at cells quadrature points
  std::vector<double> Q_solution(n_q_points);

  for (; cell!=endc; ++cell){

      fe_values.reinit (cell);
      //The following function call maps the solution vector into the cell Q_points
      fe_values.get_function_values(solution,Q_solution);

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index){
         	const auto& xx = fe_values.quadrature_point(q_index);
	 	auto Q_elem = Q_solution[q_index]-sin(2*M_PI*xx[0])*sin(6*M_PI*xx[1]);
		my_L2_norm += Q_elem*Q_elem*fe_values.JxW(q_index);
      }

  }

  //normalize the L2_norm through a form factor
  my_L2_norm = sqrt(my_L2_norm);
  std::cout << "My_L2_norm: " << my_L2_norm << std::endl;

//Her L2 norm ------
  Vector<double> norm_vec_2;
  norm_vec_2.reinit(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,solution,fp,norm_vec_2,quadrature_formula,VectorTools::NormType::L2_norm);
  std::cout << (norm_vec_2).l2_norm() << std::endl;

}



void Step3::run ()
{
  make_grid ();
  setup_system ();
  assemble_system ();
  solve ();
  output_results ();
}



int main ()
{
  deallog.depth_console (2);

  Step3 laplace_problem;
  laplace_problem.run ();



  return 0;
}
