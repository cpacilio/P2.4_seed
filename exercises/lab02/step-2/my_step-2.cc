#include<iostream>
#include<cmath>
#include<fstream>

#include<deal.II/grid/tria.h>
#include<deal.II/grid/tria_accessor.h>
#include<deal.II/grid/tria_iterator.h>

#include<deal.II/grid/grid_generator.h>

#include<deal.II/grid/grid_out.h>

#include<deal.II/dofs/dof_handler.h>
#include<deal.II/dofs/dof_tools.h>
#include<deal.II/dofs/dof_renumbering.h>

#include<deal.II/fe/fe_q.h>

#include<deal.II/lac/sparse_matrix.h>
#include<deal.II/lac/dynamic_sparsity_pattern.h>

using namespace dealii;

void make_grid(Triangulation<2>& t);

void print_grid(Triangulation<2>& t);

void distribute_dofs(DoFHandler<2>& dh);

int main(){
  //define the triangulation
  Triangulation<2> triangulation;
  //make the grid_type
  make_grid(triangulation);
  //print the grid
  //print_grid(triangulation);
  //handle the dofs
  DoFHandler<2> dof_handler(triangulation);
  //distribute the dofs
  distribute_dofs(dof_handler);
  
}

void make_grid(Triangulation<2>& t){
  GridGenerator::hyper_cube(t); 
  t.refine_global(4);
}

void print_grid(Triangulation<2>& t){
  std::ofstream out("my_grid.svg");
  GridOut grid_out;
  grid_out.write_svg(t,out);
}

void distribute_dofs(DoFHandler<2>& dh){
  FE_Q<2> fe(1);
  dh.distribute_dofs(fe);
  //renumber dofs
  DoFRenumbering::Cuthill_McKee(dh);
  //sparsity pattern
  DynamicSparsityPattern ds_pattern(dh.n_dofs(),dh.n_dofs());
  DoFTools::make_sparsity_pattern(dh,ds_pattern);
  SparsityPattern s_pattern;
  s_pattern.copy_from(ds_pattern);
  //print pattern
  std::ofstream out("my_pattern.svg");
  s_pattern.print_svg(out);
}
