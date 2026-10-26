#ifndef POISSON_1D_HPP
#define POISSON_1D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
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

#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Poisson1D
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 1;

  // Constructor.
  Poisson1D(const unsigned int                              &N_el_,
            const unsigned int                              &r_,
            const std::function<double(const Point<dim> &)> &mu_,
            const std::function<double(const Point<dim> &)> &f_)
    : N_el(N_el_)
    , r(r_)
    , mu(mu_)
    , f(f_)
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
  // Number of elements.
  const unsigned int N_el;

  // Polynomial degree.
  const unsigned int r;

  // Diffusion coefficient.
  std::function<double(const Point<dim> &)> mu;

  // Forcing term.
  std::function<double(const Point<dim> &)> f;
};

#endif