#ifndef POISSON_2D_HPP
#define POISSON_2D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
// TODO: include grid_in.h
#include <deal.II/grid/grid_in.h>
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
#include <filesytem>
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
  static constexpr unsigned int dim = 2;

  // TODO: implement a class for the Dirichlet boundary condition.
  class FunctionG : public Function<dim>
  {
    public:
      // Constructor
      FunctionG()
      {}

      // Evualtion
      virtual double
      value(const Point<dim> &p,
      const unsigned int /*component*/ = 0) const override
      {
        return p[0] + p[1]; // x + y
      }
  };

  // TODO: define a constructor that takes the mesh file name, polynomial degree,
  // diffusion coefficient, forcing term, and Neumann boundary condition as parameters.
  Poisson2D(const std:string &mesh_file_name_,
    const unsigned int &r_,
    const std::function<double(const Point<dim> &)> &mu_,
    const std::function<double(const Point<dim> &)> &f_,
    const std::function<double(const Point<dim> &)> &h_)
  : mesh_file_name(mesh_file_name_)
  , r(r_)
  , mu(mu_)
  , f(f_)
  , h(h_)
  {}

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
  std::function<double(const Point<dim> &)> h;

  // Triangulation.
  Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // TODO: add a quadrature formula for boundary integrals.
  std::unique_ptr<Quadrature<dim-1>> quadrature_boundary;

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