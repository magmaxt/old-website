#include <base/quadrature_lib.h>
#include <base/utilities.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <fe/fe_q.h>
#include <grid/grid_out.h>


#include <lac/constraint_matrix.h>

#include <grid/grid_refinement.h>

#include <numerics/error_estimator.h>

using namespace dealii;



template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();
    
  private:
    void setup_system ();
    void assemble_matrix ();
    void assemble_primal_rhs ();
    void assemble_dual_rhs ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       primal_solution;
    Vector<double>       dual_solution;
    Vector<double>       primal_rhs;
    Vector<double>       dual_rhs;
};





template <int dim>
LaplaceProblem<dim>::LaplaceProblem ()
		:
		dof_handler (triangulation),
                fe (1)
{}



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem () 
{
  dof_handler.clear ();
}



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  primal_solution.reinit (dof_handler.n_dofs());
  dual_solution.reinit (dof_handler.n_dofs());
  primal_rhs.reinit (dof_handler.n_dofs());
  dual_rhs.reinit (dof_handler.n_dofs());

  
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);

  hanging_node_constraints.close ();

  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
}


template <int dim>
void LaplaceProblem<dim>::assemble_matrix () 
{  
  const QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_gradients | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;

      fe_values.reinit (cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
				 fe_values.shape_grad(j,q_point) *
				 fe_values.JxW(q_point));

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add (local_dof_indices[i],
			     local_dof_indices[j],
			     cell_matrix(i,j));
    }

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (primal_rhs);
  hanging_node_constraints.condense (dual_rhs);

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      primal_solution,
				      primal_rhs);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      dual_solution,
				      dual_rhs);
}



template <int dim>
void LaplaceProblem<dim>::assemble_primal_rhs () 
{  
  const QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values    |  
			   update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_rhs = 0;

      fe_values.reinit (cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			  (1./2. * numbers::PI * numbers::PI *
			   std::cos(numbers::PI/2 *
				fe_values.quadrature_point(q_point)[0]) *
			   std::cos(numbers::PI/2 *
				    fe_values.quadrature_point(q_point)[1])) *
			  fe_values.JxW(q_point));

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	primal_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
}



template <int dim>
void LaplaceProblem<dim>::assemble_dual_rhs () 
{
  dual_rhs = 0;
  
/* Option 1: Point value at (1/2,1/2) */
/*        2: Point derivative at (1/2,1/2) */
/*        3: boundary integral */
  const unsigned int functional = 1;
  
  switch (functional)
    {
      case 1:
      case 2:
      {
	MappingQ1<dim> mapping;
	std::pair<typename DoFHandler<dim>::active_cell_iterator,
	  Point<dim> >
	  cell_and_point
	  = GridTools::find_active_cell_around_point (mapping,
						      dof_handler,
						      Point<dim> (0.5, 0.5));
	const Quadrature<dim> quadrature (cell_and_point.second);
	FEValues<dim> fe_values (fe, quadrature,
				 update_values | update_gradients);

	fe_values.reinit (cell_and_point.first);
  
	std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
	cell_and_point.first->get_dof_indices (local_dof_indices);

	if (functional == 1)
	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    dual_rhs(local_dof_indices[i]) = fe_values.shape_value(i,0);
	else
	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    dual_rhs(local_dof_indices[i]) = fe_values.shape_grad(i,0)[0];

	break;
      }

      case 3:
      {
	const QGauss<dim-1>  quadrature_formula(2);

	FEFaceValues<dim> fe_face_values (fe, quadrature_formula, 
					  update_gradients |
					  update_JxW_values);

	const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	const unsigned int   n_q_points    = quadrature_formula.size();

	Vector<double>       cell_rhs (dofs_per_cell);

	std::vector<unsigned int> local_dof_indices (dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	  cell = dof_handler.begin_active(),
	  endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	    if (cell->face(f)->at_boundary()
		&&
		(cell->face(f)->center()[0] == 1))
	      {
		cell_rhs = 0;
		fe_face_values.reinit (cell, f);

		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    cell_rhs(i) += (fe_face_values.shape_grad(i,q_point)[0] *
				    fe_face_values.JxW(q_point));

		cell->get_dof_indices (local_dof_indices);
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  dual_rhs(local_dof_indices[i]) += cell_rhs(i);
	  }
	break;
      }
      
      default:
	    Assert (false, ExcNotImplemented());
    }
}





template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, primal_solution, primal_rhs,
	    preconditioner);
  cg.solve (system_matrix, dual_solution, dual_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (primal_solution);
  hanging_node_constraints.distribute (dual_solution);
}



template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
/* Option 1: Global refinement */
/*        2: Kelly */
/*        3: DWR */
  const int refinement_option = 3;

  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  switch (refinement_option)
    {
      case 1:
      {
	triangulation.refine_global (1);
	return;
      }

      case 2:
      {
	KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss<dim-1>(3),
					    typename FunctionMap<dim>::type(),
					    primal_solution,
					    estimated_error_per_cell);
	break;
      }


      case 3:
      {
	Vector<float> primal_error_per_cell (triangulation.n_active_cells());
	KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss<dim-1>(3),
					    typename FunctionMap<dim>::type(),
					    primal_solution,
					    primal_error_per_cell);
	Vector<float> dual_error_per_cell (triangulation.n_active_cells());
	KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss<dim-1>(3),
					    typename FunctionMap<dim>::type(),
					    dual_solution,
					    dual_error_per_cell);

	unsigned int index=0;
	for (typename Triangulation<dim>::active_cell_iterator
	       cell = triangulation.begin_active();
	     cell != triangulation.end();
	     ++index, ++cell)
	  estimated_error_per_cell(index)
	    = (primal_error_per_cell(index) *
	       dual_error_per_cell(index) *
	       cell->diameter());
	
	break;
      }

      default:
	    Assert (false, ExcNotImplemented());
    }
  

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  {
    std::string filename = "grid-";
    filename += Utilities::int_to_string(cycle,2);
    filename += ".eps";
  
    std::ofstream output (filename.c_str());

    GridOut grid_out;
    grid_out.write_eps (triangulation, output);
  }
  
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss<dim-1>(3),
					    typename FunctionMap<dim>::type(),
					    primal_solution,
					    estimated_error_per_cell);
  unsigned int index=0;
  for (typename Triangulation<dim>::active_cell_iterator 
         cell=triangulation.begin_active(); cell!=triangulation.end(); 
       ++cell, ++index)
   estimated_error_per_cell(index) *= std::pow (cell->diameter(), -2);

  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (primal_solution, "primal_solution");
    data_out.add_data_vector (dual_solution, "dual_solution");
    data_out.add_data_vector (estimated_error_per_cell, "grad2_u");
    data_out.build_patches ();
  
    std::string filename = "solution-";
    filename += Utilities::int_to_string(cycle,2);
    filename += ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
  }

  double point_value, point_derivative, boundary_integral;

  point_value = VectorTools::point_value (dof_handler,
					  primal_solution,
					  Point<dim>(0.5, 0.5));
  {
    const double h = GridTools::minimal_cell_diameter (triangulation);
    point_derivative = (VectorTools::point_value (dof_handler,
						  primal_solution,
						  Point<dim>(0.5+h/100, 0.5))
			-
			VectorTools::point_value (dof_handler,
						  primal_solution,
						  Point<dim>(0.5-h/100, 0.5)))
		       / (2*h/100);
  }

  {
    boundary_integral = 0;

    const QGauss<dim-1>  quadrature_formula(2);
    FEFaceValues<dim> fe_face_values (fe, quadrature_formula, 
				      update_gradients |
				      update_JxW_values);

    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<Tensor<1,dim> > cell_gradients (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	if (cell->face(f)->at_boundary()
	    &&
	    (cell->face(f)->center()[0] == 1))
	  {
	    fe_face_values.reinit (cell, f);
	    fe_face_values.get_function_gradients (primal_solution,
						   cell_gradients);
	    
	    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	      boundary_integral += cell_gradients[q_point][0] *
				   fe_face_values.JxW (q_point);
	  }
    
  }
  

  const double exact_point_value
    = (std::cos(numbers::PI/2 * 0.5) *
       std::cos(numbers::PI/2 * 0.5));
  const double exact_point_derivative
    = (-numbers::PI/2 *
       std::sin(numbers::PI/2 * 0.5) *
       std::cos(numbers::PI/2 * 0.5));
  const double exact_boundary_integral = -2;
  
  std::cout << dof_handler.n_dofs() << " \t"
	    << std::setw(12)
	    << point_value << " \t"
	    << std::setw(12)
    	    << point_value-exact_point_value << " \t"
	    << std::setw(12)
	    << point_derivative << " \t"
	    << std::setw(12)
	    << point_derivative-exact_point_derivative << " \t"
	    << std::setw(12)
    	    << boundary_integral << " \t"
	    << std::setw(12)
    	    << boundary_integral-exact_boundary_integral << " \t"
	    << std::endl;
}




template <int dim>
void LaplaceProblem<dim>::run () 
{
  std::cout << "n_dofs    u(x0)         error            du(x0)/dx       error             boundary int.   error" << std::endl;
  std::cout.precision(8);
  unsigned int cycle=0;
  do
    {
      if (cycle == 0)
	{
	  if (true)
	    {
	      GridGenerator::hyper_cube (triangulation, -1, 1);
	      triangulation.refine_global (2);
	    }
	  else
	    {
	      static const Point<2> vertices_1[]
		= {  Point<2> (-1.,   -1.),
		     Point<2> (-1./2, -1.),
		     Point<2> (0.,    -1.),
		     Point<2> (+1./2, -1.),
		     Point<2> (+1,    -1.),
	     
		     Point<2> (-1.,   -1./2.),
		     Point<2> (-1./2, -1./2.),
		     Point<2> (0.,    -1./2.),
		     Point<2> (+1./2, -1./2.),
		     Point<2> (+1,    -1./2.),
	     
		     Point<2> (-1.,   0.),
		     Point<2> (-1./2, 0.),
		     Point<2> (+1./2, 0.),
		     Point<2> (+1,    0.),
	     
		     Point<2> (-1.,   1./2.),
		     Point<2> (-1./2, 1./2.),
		     Point<2> (0.,    1./2.),
		     Point<2> (+1./2, 1./2.),
		     Point<2> (+1,    1./2.),
	     
		     Point<2> (-1.,   1.),
		     Point<2> (-1./2, 1.),
		     Point<2> (0.,    1.),			  
		     Point<2> (+1./2, 1.),
		     Point<2> (+1,    1.)    };
	      const unsigned int
		n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

	      const std::vector<Point<dim> > vertices (&vertices_1[0],
						       &vertices_1[n_vertices]);

	      static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
		= {{0, 1, 5, 6},
		   {1, 2, 6, 7},
		   {2, 3, 7, 8},
		   {3, 4, 8, 9},
		   {5, 6, 10, 11},
		   {8, 9, 12, 13},
		   {10, 11, 14, 15},
		   {12, 13, 17, 18},
		   {14, 15, 19, 20},
		   {15, 16, 20, 21},
		   {16, 17, 21, 22},
		   {17, 18, 22, 23}};
	      const unsigned int
		n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

	      std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
	      for (unsigned int i=0; i<n_cells; ++i) 
		{
		  for (unsigned int j=0;
		       j<GeometryInfo<dim>::vertices_per_cell;
		       ++j)
		    cells[i].vertices[j] = cell_vertices[i][j];
		  cells[i].material_id = 0;
		}

	      triangulation.create_triangulation (vertices,
						  cells,
						  SubCellData());
	      triangulation.refine_global (1);
	    }
	}
      else
	refine_grid ();
	  
      
      setup_system ();
      assemble_primal_rhs ();
      assemble_dual_rhs ();
      assemble_matrix ();
      solve ();
      output_results (cycle);

      ++cycle;
    }
  while (dof_handler.n_dofs() < 10000);
}



int main () 
{

  try
    {
      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
