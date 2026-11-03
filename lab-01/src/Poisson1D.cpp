#include "Poisson1D.hpp"

void
Poisson1D::setup()
{
  std::cout << "Initialize the mesh." << std::endl;
  {
    GridGenerator::subdivided_hyper_cube(mesh, N_el, 0.0, 1.0, true);
    std::cout << " Number of elements = " <<
      mesh.n_active_cells() << std::endl;

    // Write the mesh to file.
    const std::string mesh_file_name = "mesh-" + std::to_string(N_el)
      + ".vtk";
    GridOut grid_out;
    std::ofstream grid_out_file(mesh_file_name);
    grid_out.write_vtk(mesh, grid_out_file);
    std::cout << "Mesh saved to " << mesh_file_name << std::endl;
  }

  std::cout << "Initialize the FEM space." << std::endl;
  {
    fe = std::make_unique<FE_SimplexP<dim>>(r);
    std::cout << " Degree = " << fe->degree << std::endl;
    std::cout << " DoFs per cell = " << fe->dofs_per_cell << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r+1);
    std::cout << " Quadrature points per cell = " << quadrature->size()
      << std::endl;
  }

  std::cout << "Initialize the dof handler." << std::endl;
  {
    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);
    std::cout << " Number DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  std::cout << "Initialize the system." << std::endl;
  {
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    rhs.reinit(dof_handler.n_dofs());
    solution.reinit(dof_handler.n_dofs());

    std::cout << "Initialization completed!" << std::endl;
  }
}

void
Poisson1D::assemble()
{}

void
Poisson1D::solve()
{}

void
Poisson1D::output() const
{}