#ifndef POISSON_2D_HPP
#define POISSON_2D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
// TODO: include grid_in.h
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

// TODO: include filesystem
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Poisson2D
{
public:
  // Physical dimension (1D, 2D, 3D)
  // TODO: fix problem dimension.

  // TODO: implement a class for the Dirichlet boundary condition.

  // TODO: define a constructor that takes the mesh file name, polynomial degree, diffusion coefficient, forcing term, and Neumann boundary condition as parameters.

  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve();

  // Output.
  void
  output() const;

protected:
  // Name of the mesh.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Diffusion coefficient.
  std::function<double(const Point<dim> &)> mu;

  // Forcing term.
  std::function<double(const Point<dim> &)> f;

  // TODO: add a Neumann boundary condition.

  // Triangulation.
  Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // TODO: add a quadrature formula for boundary integrals.

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // Sparsity pattern.
  SparsityPattern sparsity_pattern;

  // System matrix.
  SparseMatrix<double> system_matrix;

  // System right-hand side.
  Vector<double> system_rhs;

  // System solution.
  Vector<double> solution;
};

#endif