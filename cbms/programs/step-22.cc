/* $Id: step-22.cc 18724 2009-04-23 23:06:37Z bangerth $ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2008 */

/*    $Id: step-22.cc 18724 2009-04-23 23:06:37Z bangerth $       */
/*                                                                */
/*    Copyright (C) 2008, 2009 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


                        
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/constraint_matrix.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

#include <lac/sparse_direct.h>

#include <lac/sparse_ilu.h>

#include <fstream>
#include <sstream>

using namespace dealii;

                        
template <int dim>
struct InnerPreconditioner;

template <>
struct InnerPreconditioner<2> 
{
    typedef SparseDirectUMFPACK type;
};

template <>
struct InnerPreconditioner<3> 
{
    typedef SparseILU<double> type;
};


                    
template <int dim>
class StokesProblem 
{
  public:
    StokesProblem (const unsigned int degree);
    void run ();
    
  private:
    void setup_dofs ();
    void assemble_system ();
    void solve (const BlockVector<double> &rhs,
		const ConstraintMatrix &constraints,
		BlockVector<double> &solution);
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();
    
    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     primal_constraints;
    ConstraintMatrix     dual_constraints;
    
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> primal_solution;
    BlockVector<double> dual_solution;
    BlockVector<double> primal_rhs;
    BlockVector<double> dual_rhs;

    std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
};


                    
template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
};


template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>  &p,
			    const unsigned int component) const 
{
  Assert (component < this->n_components,
	  ExcIndexRange (component, 0, this->n_components));
  
  if (component == 1)
    return p[dim-1];
  return 0;
}


template <int dim>
void
BoundaryValues<dim>::vector_value (const Point<dim> &p,
				   Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = BoundaryValues<dim>::value (p, c);
}



template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
    
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                           const unsigned int /*component*/) const 
{
  return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide<dim>::value (p, c);
}


                        
                        
                        
template <class Matrix, class Preconditioner>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix         &m,
                   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const SmartPointer<const Preconditioner> preconditioner;
};


template <class Matrix, class Preconditioner>
InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						     const Preconditioner &preconditioner)
		:
		matrix (&m),
		preconditioner (&preconditioner)
{}


                    
template <class Matrix, class Preconditioner>
void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						  const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
  SolverCG<>    cg (solver_control);

  dst = 0;

  cg.solve (*matrix, dst, src, *preconditioner);
}



template <class Preconditioner>
class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		     const InverseMatrix<SparseMatrix<double>, Preconditioner> &A_inverse);

    void vmult (Vector<double>       &dst,
		const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>, Preconditioner> > A_inverse;
    
    mutable Vector<double> tmp1, tmp2;
};



template <class Preconditioner>
SchurComplement<Preconditioner>::
SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		 const InverseMatrix<SparseMatrix<double>,Preconditioner> &A_inverse)
		:
		system_matrix (&system_matrix),
		A_inverse (&A_inverse),
		tmp1 (system_matrix.block(0,0).m()),
		tmp2 (system_matrix.block(0,0).m())
{}


template <class Preconditioner>
void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
					     const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  A_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}


                        

template <int dim>
StokesProblem<dim>::StokesProblem (const unsigned int degree)
                :
                degree (degree),
                triangulation (Triangulation<dim>::maximum_smoothing),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1),
                dof_handler (triangulation)
{}


                        
				 
template <int dim>
void StokesProblem<dim>::setup_dofs ()
{
  A_preconditioner.reset ();
  system_matrix.clear ();
  
  dof_handler.distribute_dofs (fe);  
  DoFRenumbering::Cuthill_McKee (dof_handler);

  std::vector<unsigned int> block_component (dim+1,0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise (dof_handler, block_component);

  {
    primal_constraints.clear ();
    dual_constraints.clear ();

    ConstraintMatrix     only_hanging_node_constraints;

    std::vector<bool> component_mask (dim+1, true);
    component_mask[dim] = false;

    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      BoundaryValues<dim>(),
					      primal_constraints,
					      component_mask);
    VectorTools::interpolate_boundary_values (dof_handler,
					      2,
					      ZeroFunction<dim>(dim+1),
					      primal_constraints,
					      component_mask);

    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      ZeroFunction<dim>(dim+1),
					      dual_constraints,
					      component_mask);
    VectorTools::interpolate_boundary_values (dof_handler,
					      2,
					      ZeroFunction<dim>(dim+1),
					      dual_constraints,
					      component_mask);

    DoFTools::make_hanging_node_constraints (dof_handler,
					     only_hanging_node_constraints);
    primal_constraints.merge (only_hanging_node_constraints);
    dual_constraints.merge (only_hanging_node_constraints);
  }

  primal_constraints.close ();
  dual_constraints.close ();

  std::vector<unsigned int> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);  
  const unsigned int n_u = dofs_per_block[0],
                     n_p = dofs_per_block[1];

  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
            << std::endl;
      
  {
    BlockCompressedSimpleSparsityPattern csp (2,2);

    csp.block(0,0).reinit (n_u, n_u);
    csp.block(1,0).reinit (n_p, n_u);
    csp.block(0,1).reinit (n_u, n_p);
    csp.block(1,1).reinit (n_p, n_p);
  
    csp.collect_sizes();    
  
    DoFTools::make_sparsity_pattern (dof_handler, csp, primal_constraints, false);
    sparsity_pattern.copy_from (csp);
  }
  
  system_matrix.reinit (sparsity_pattern);
                                   
  primal_solution.reinit (2);
  primal_solution.block(0).reinit (n_u);
  primal_solution.block(1).reinit (n_p);
  primal_solution.collect_sizes ();
  
  dual_solution.reinit (2);
  dual_solution.block(0).reinit (n_u);
  dual_solution.block(1).reinit (n_p);
  dual_solution.collect_sizes ();
  
  primal_rhs.reinit (2);
  primal_rhs.block(0).reinit (n_u);
  primal_rhs.block(1).reinit (n_p);
  primal_rhs.collect_sizes ();

  dual_rhs.reinit (2);
  dual_rhs.block(0).reinit (n_u);
  dual_rhs.block(1).reinit (n_p);
  dual_rhs.collect_sizes ();
}


                        
template <int dim>
void StokesProblem<dim>::assemble_system () 
{
  system_matrix=0;
  primal_rhs=0;
  dual_rhs=0;
  
  QGauss<dim>   quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);
  
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const RightHandSide<dim>          right_hand_side;
  std::vector<Vector<double> >      rhs_values (n_q_points,
                                                Vector<double>(dim+1));

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<SymmetricTensor<2,dim> > phi_grads_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);
				   
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    { 
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;
      
      right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                        rhs_values);
      
      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    {
	      phi_grads_u[k] = fe_values[velocities].symmetric_gradient (k, q);
	      div_phi_u[k]   = fe_values[velocities].divergence (k, q);
	      phi_p[k]       = fe_values[pressure].value (k, q);
	    }

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<=i; ++j)
		{
		  local_matrix(i,j) += (phi_grads_u[i] * phi_grads_u[j]
					- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j]
					+ phi_p[i] * phi_p[j])
				       * fe_values.JxW(q);     

		}

	      const unsigned int component_i =
		fe.system_to_component_index(i).first;
	      local_rhs(i) += fe_values.shape_value(i,q) * 
			      rhs_values[q](component_i) *
			      fe_values.JxW(q);
	    }
	}


      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=i+1; j<dofs_per_cell; ++j)
	  local_matrix(i,j) = local_matrix(j,i);

      cell->get_dof_indices (local_dof_indices);
      primal_constraints.distribute_local_to_global (local_matrix, local_rhs,
						     local_dof_indices, 
						     system_matrix, primal_rhs);
    }

  {
    MappingQ1<dim> mapping;
    std::pair<typename DoFHandler<dim>::active_cell_iterator,
      Point<dim> >
      cell_and_point
      = GridTools::find_active_cell_around_point (mapping,
						  dof_handler,
						  Point<dim> (3.5, 3.5, 0.5));
    const Quadrature<dim> quadrature (cell_and_point.second);
    FEValues<dim> fe_values (fe, quadrature,
			     update_values | update_gradients);

    fe_values.reinit (cell_and_point.first);
  
    std::vector<unsigned int> local_dof_indices (fe.dofs_per_cell);
    cell_and_point.first->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      if (fe.system_to_component_index(i).first == 1)
	dual_rhs(local_dof_indices[i]) = -fe_values.shape_value(i,0);

    dual_constraints.distribute (dual_rhs);
  }
  
  
  std::cout << "   Computing preconditioner..." << std::endl << std::flush;
      
  A_preconditioner
    = std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
  A_preconditioner->initialize (system_matrix.block(0,0),
				typename InnerPreconditioner<dim>::type::AdditionalData());

}



                        
template <int dim>
void StokesProblem<dim>::solve (const BlockVector<double> &rhs,
				const ConstraintMatrix &constraints,
				BlockVector<double> &solution) 
{
  const InverseMatrix<SparseMatrix<double>,
                      typename InnerPreconditioner<dim>::type>
    A_inverse (system_matrix.block(0,0), *A_preconditioner);
  Vector<double> tmp (solution.block(0).size());
  
  {
    Vector<double> schur_rhs (solution.block(1).size());
    A_inverse.vmult (tmp, rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= rhs.block(1);
  
    SchurComplement<typename InnerPreconditioner<dim>::type>
      schur_complement (system_matrix, A_inverse);
    
    SolverControl solver_control (solution.block(1).size(),
				  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);
    
    SparseILU<double> preconditioner;
    preconditioner.initialize (system_matrix.block(1,1), 
      SparseILU<double>::AdditionalData());
  
    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      m_inverse (system_matrix.block(1,1), preconditioner);
    
    cg.solve (schur_complement, solution.block(1), schur_rhs,
	      m_inverse);
  
    constraints.distribute (solution);
  
    std::cout << "   "
	      << solver_control.last_step()
	      << " outer CG Schur complement iterations for pressure"
	      << std::flush
	      << std::endl;    
  }
    
  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += rhs.block(0);
  
    A_inverse.vmult (solution.block(0), tmp);

    constraints.distribute (solution);
  }
}


                        

template <int dim>
void
StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
{
  std::vector<std::string> primal_solution_names (dim, "primal_velocity");
  primal_solution_names.push_back ("primal_pressure");
  std::vector<std::string> dual_solution_names (dim, "dual_velocity");
  dual_solution_names.push_back ("dual_pressure");
  
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
      
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);  
  data_out.add_data_vector (primal_solution, primal_solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.add_data_vector (dual_solution, dual_solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.build_patches ();
  
  std::ostringstream filename;
  filename << "solution-"
           << Utilities::int_to_string (refinement_cycle, 2)
           << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);

  Vector<double> point_values (dim+1);
  VectorTools::point_value (dof_handler,
			    primal_solution,
			    Point<dim>(3.5, 3.5, 0.5),
			    point_values);
  
  std::cout << "   y-velocity at (3.5,3.5,0.5): ["
	    << point_values(0) << ',' << point_values(1) << ',' << point_values(2) << ']'
	    << std::endl;
}


                        
template <int dim>
void
StokesProblem<dim>::refine_mesh () 
{
/* Option 1: Global refinement */
/*        2: Kelly */
/*        3: DWR */
  const int refinement_option = 3;

  std::vector<bool> component_mask (dim+1, false);
  component_mask[dim] = true;
  
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
					    estimated_error_per_cell,
					    component_mask);
	break;
      }


      case 3:
      {
	Vector<float> primal_error_per_cell (triangulation.n_active_cells());
	KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss<dim-1>(3),
					    typename FunctionMap<dim>::type(),
					    primal_solution,
					    primal_error_per_cell,
					    component_mask);
	Vector<float> dual_error_per_cell (triangulation.n_active_cells());
	KellyErrorEstimator<dim>::estimate (dof_handler,
					    QGauss<dim-1>(3),
					    typename FunctionMap<dim>::type(),
					    dual_solution,
					    dual_error_per_cell,
					    component_mask);

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
						   0.15, 0.1);

  triangulation.execute_coarsening_and_refinement ();
}


                        
template <int dim>
void StokesProblem<dim>::run () 
{
  {
    static const Point<dim> vertices_1[]
      = {Point<dim>(0,0,0),
	 Point<dim>(1,0,0),
	 Point<dim>(2,0,0),
	 Point<dim>(3,0,0),
	 Point<dim>(4,0,0),
	 Point<dim>(5,0,0),
	 Point<dim>(6,0,0),
	 Point<dim>(7,0,0),
	 Point<dim>(0,1,0),
	 Point<dim>(1,1,0),
	 Point<dim>(2,1,0),
	 Point<dim>(3,1,0),
	 Point<dim>(4,1,0),
	 Point<dim>(5,1,0),
	 Point<dim>(6,1,0),
	 Point<dim>(7,1,0),
	 Point<dim>(0,2,0),
	 Point<dim>(1,2,0),
	 Point<dim>(2,2,0),
	 Point<dim>(3,2,0),
	 Point<dim>(4,2,0),
	 Point<dim>(5,2,0),
	 Point<dim>(6,2,0),
	 Point<dim>(7,2,0),
	 Point<dim>(0,3,0),
	 Point<dim>(1,3,0),
	 Point<dim>(2,3,0),
	 Point<dim>(3,3,0),
	 Point<dim>(4,3,0),
	 Point<dim>(5,3,0),
	 Point<dim>(6,3,0),
	 Point<dim>(7,3,0),
	 Point<dim>(0,4,0),
	 Point<dim>(1,4,0),
	 Point<dim>(2,4,0),
	 Point<dim>(3,4,0),
	 Point<dim>(4,4,0),
	 Point<dim>(5,4,0),
	 Point<dim>(6,4,0),
	 Point<dim>(7,4,0),
	 Point<dim>(0,5,0),
	 Point<dim>(1,5,0),
	 Point<dim>(2,5,0),
	 Point<dim>(3,5,0),
	 Point<dim>(4,5,0),
	 Point<dim>(5,5,0),
	 Point<dim>(6,5,0),
	 Point<dim>(7,5,0),
	 Point<dim>(0,6,0),
	 Point<dim>(1,6,0),
	 Point<dim>(2,6,0),
	 Point<dim>(3,6,0),
	 Point<dim>(4,6,0),
	 Point<dim>(5,6,0),
	 Point<dim>(6,6,0),
	 Point<dim>(7,6,0),
	 Point<dim>(0,7,0),
	 Point<dim>(1,7,0),
	 Point<dim>(2,7,0),
	 Point<dim>(3,7,0),
	 Point<dim>(4,7,0),
	 Point<dim>(5,7,0),
	 Point<dim>(6,7,0),
	 Point<dim>(7,7,0),
	 Point<dim>(0,0,1),
	 Point<dim>(1,0,1),
	 Point<dim>(2,0,1),
	 Point<dim>(3,0,1),
	 Point<dim>(4,0,1),
	 Point<dim>(5,0,1),
	 Point<dim>(6,0,1),
	 Point<dim>(7,0,1),
	 Point<dim>(0,1,1),
	 Point<dim>(1,1,1),
	 Point<dim>(2,1,1),
	 Point<dim>(3,1,1),
	 Point<dim>(4,1,1),
	 Point<dim>(5,1,1),
	 Point<dim>(6,1,1),
	 Point<dim>(7,1,1),
	 Point<dim>(0,2,1),
	 Point<dim>(1,2,1),
	 Point<dim>(2,2,1),
	 Point<dim>(3,2,1),
	 Point<dim>(4,2,1),
	 Point<dim>(5,2,1),
	 Point<dim>(6,2,1),
	 Point<dim>(7,2,1),
	 Point<dim>(0,3,1),
	 Point<dim>(1,3,1),
	 Point<dim>(2,3,1),
	 Point<dim>(3,3,1),
	 Point<dim>(4,3,1),
	 Point<dim>(5,3,1),
	 Point<dim>(6,3,1),
	 Point<dim>(7,3,1),
	 Point<dim>(0,4,1),
	 Point<dim>(1,4,1),
	 Point<dim>(2,4,1),
	 Point<dim>(3,4,1),
	 Point<dim>(4,4,1),
	 Point<dim>(5,4,1),
	 Point<dim>(6,4,1),
	 Point<dim>(7,4,1),
	 Point<dim>(0,5,1),
	 Point<dim>(1,5,1),
	 Point<dim>(2,5,1),
	 Point<dim>(3,5,1),
	 Point<dim>(4,5,1),
	 Point<dim>(5,5,1),
	 Point<dim>(6,5,1),
	 Point<dim>(7,5,1),
	 Point<dim>(0,6,1),
	 Point<dim>(1,6,1),
	 Point<dim>(2,6,1),
	 Point<dim>(3,6,1),
	 Point<dim>(4,6,1),
	 Point<dim>(5,6,1),
	 Point<dim>(6,6,1),
	 Point<dim>(7,6,1),
	 Point<dim>(0,7,1),
	 Point<dim>(1,7,1),
	 Point<dim>(2,7,1),
	 Point<dim>(3,7,1),
	 Point<dim>(4,7,1),
	 Point<dim>(5,7,1),
	 Point<dim>(6,7,1),
	 Point<dim>(7,7,1),
	 Point<dim>(0,0,2),
	 Point<dim>(1,0,2),
	 Point<dim>(2,0,2),
	 Point<dim>(3,0,2),
	 Point<dim>(4,0,2),
	 Point<dim>(5,0,2),
	 Point<dim>(6,0,2),
	 Point<dim>(7,0,2),
	 Point<dim>(0,1,2),
	 Point<dim>(1,1,2),
	 Point<dim>(2,1,2),
	 Point<dim>(3,1,2),
	 Point<dim>(4,1,2),
	 Point<dim>(5,1,2),
	 Point<dim>(6,1,2),
	 Point<dim>(7,1,2),
	 Point<dim>(0,2,2),
	 Point<dim>(1,2,2),
	 Point<dim>(2,2,2),
	 Point<dim>(3,2,2),
	 Point<dim>(4,2,2),
	 Point<dim>(5,2,2),
	 Point<dim>(6,2,2),
	 Point<dim>(7,2,2),
	 Point<dim>(0,3,2),
	 Point<dim>(1,3,2),
	 Point<dim>(2,3,2),
	 Point<dim>(3,3,2),
	 Point<dim>(4,3,2),
	 Point<dim>(5,3,2),
	 Point<dim>(6,3,2),
	 Point<dim>(7,3,2),
	 Point<dim>(0,4,2),
	 Point<dim>(1,4,2),
	 Point<dim>(2,4,2),
	 Point<dim>(3,4,2),
	 Point<dim>(4,4,2),
	 Point<dim>(5,4,2),
	 Point<dim>(6,4,2),
	 Point<dim>(7,4,2),
	 Point<dim>(0,5,2),
	 Point<dim>(1,5,2),
	 Point<dim>(2,5,2),
	 Point<dim>(3,5,2),
	 Point<dim>(4,5,2),
	 Point<dim>(5,5,2),
	 Point<dim>(6,5,2),
	 Point<dim>(7,5,2),
	 Point<dim>(0,6,2),
	 Point<dim>(1,6,2),
	 Point<dim>(2,6,2),
	 Point<dim>(3,6,2),
	 Point<dim>(4,6,2),
	 Point<dim>(5,6,2),
	 Point<dim>(6,6,2),
	 Point<dim>(7,6,2),
	 Point<dim>(0,7,2),
	 Point<dim>(1,7,2),
	 Point<dim>(2,7,2),
	 Point<dim>(3,7,2),
	 Point<dim>(4,7,2),
	 Point<dim>(5,7,2),
	 Point<dim>(6,7,2),
	 Point<dim>(7,7,2),
	 Point<dim>(0,0,3),
	 Point<dim>(1,0,3),
	 Point<dim>(2,0,3),
	 Point<dim>(3,0,3),
	 Point<dim>(4,0,3),
	 Point<dim>(5,0,3),
	 Point<dim>(6,0,3),
	 Point<dim>(7,0,3),
	 Point<dim>(0,1,3),
	 Point<dim>(1,1,3),
	 Point<dim>(2,1,3),
	 Point<dim>(3,1,3),
	 Point<dim>(4,1,3),
	 Point<dim>(5,1,3),
	 Point<dim>(6,1,3),
	 Point<dim>(7,1,3),
	 Point<dim>(0,2,3),
	 Point<dim>(1,2,3),
	 Point<dim>(2,2,3),
	 Point<dim>(3,2,3),
	 Point<dim>(4,2,3),
	 Point<dim>(5,2,3),
	 Point<dim>(6,2,3),
	 Point<dim>(7,2,3),
	 Point<dim>(0,3,3),
	 Point<dim>(1,3,3),
	 Point<dim>(2,3,3),
	 Point<dim>(3,3,3),
	 Point<dim>(4,3,3),
	 Point<dim>(5,3,3),
	 Point<dim>(6,3,3),
	 Point<dim>(7,3,3),
	 Point<dim>(0,4,3),
	 Point<dim>(1,4,3),
	 Point<dim>(2,4,3),
	 Point<dim>(3,4,3),
	 Point<dim>(4,4,3),
	 Point<dim>(5,4,3),
	 Point<dim>(6,4,3),
	 Point<dim>(7,4,3),
	 Point<dim>(0,5,3),
	 Point<dim>(1,5,3),
	 Point<dim>(2,5,3),
	 Point<dim>(3,5,3),
	 Point<dim>(4,5,3),
	 Point<dim>(5,5,3),
	 Point<dim>(6,5,3),
	 Point<dim>(7,5,3),
	 Point<dim>(0,6,3),
	 Point<dim>(1,6,3),
	 Point<dim>(2,6,3),
	 Point<dim>(3,6,3),
	 Point<dim>(4,6,3),
	 Point<dim>(5,6,3),
	 Point<dim>(6,6,3),
	 Point<dim>(7,6,3),
	 Point<dim>(0,7,3),
	 Point<dim>(1,7,3),
	 Point<dim>(2,7,3),
	 Point<dim>(3,7,3),
	 Point<dim>(4,7,3),
	 Point<dim>(5,7,3),
	 Point<dim>(6,7,3),
	 Point<dim>(7,7,3),
	 Point<dim>(0,0,4),
	 Point<dim>(1,0,4),
	 Point<dim>(2,0,4),
	 Point<dim>(3,0,4),
	 Point<dim>(4,0,4),
	 Point<dim>(5,0,4),
	 Point<dim>(6,0,4),
	 Point<dim>(7,0,4),
	 Point<dim>(0,1,4),
	 Point<dim>(1,1,4),
	 Point<dim>(2,1,4),
	 Point<dim>(3,1,4),
	 Point<dim>(4,1,4),
	 Point<dim>(5,1,4),
	 Point<dim>(6,1,4),
	 Point<dim>(7,1,4),
	 Point<dim>(0,2,4),
	 Point<dim>(1,2,4),
	 Point<dim>(2,2,4),
	 Point<dim>(3,2,4),
	 Point<dim>(4,2,4),
	 Point<dim>(5,2,4),
	 Point<dim>(6,2,4),
	 Point<dim>(7,2,4),
	 Point<dim>(0,3,4),
	 Point<dim>(1,3,4),
	 Point<dim>(2,3,4),
	 Point<dim>(3,3,4),
	 Point<dim>(4,3,4),
	 Point<dim>(5,3,4),
	 Point<dim>(6,3,4),
	 Point<dim>(7,3,4),
	 Point<dim>(0,4,4),
	 Point<dim>(1,4,4),
	 Point<dim>(2,4,4),
	 Point<dim>(3,4,4),
	 Point<dim>(4,4,4),
	 Point<dim>(5,4,4),
	 Point<dim>(6,4,4),
	 Point<dim>(7,4,4),
	 Point<dim>(0,5,4),
	 Point<dim>(1,5,4),
	 Point<dim>(2,5,4),
	 Point<dim>(3,5,4),
	 Point<dim>(4,5,4),
	 Point<dim>(5,5,4),
	 Point<dim>(6,5,4),
	 Point<dim>(7,5,4),
	 Point<dim>(0,6,4),
	 Point<dim>(1,6,4),
	 Point<dim>(2,6,4),
	 Point<dim>(3,6,4),
	 Point<dim>(4,6,4),
	 Point<dim>(5,6,4),
	 Point<dim>(6,6,4),
	 Point<dim>(7,6,4),
	 Point<dim>(0,7,4),
	 Point<dim>(1,7,4),
	 Point<dim>(2,7,4),
	 Point<dim>(3,7,4),
	 Point<dim>(4,7,4),
	 Point<dim>(5,7,4),
	 Point<dim>(6,7,4),
	 Point<dim>(7,7,4)};

    const unsigned int
      n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

    const std::vector<Point<dim> > vertices (&vertices_1[0],
					     &vertices_1[n_vertices]);

    static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
      = {
	  {0,1,8,9,64,65,72,73},
	  {1,2,9,10,65,66,73,74},
	  {2,3,10,11,66,67,74,75},
	  {3,4,11,12,67,68,75,76},
	  {4,5,12,13,68,69,76,77},
	  {5,6,13,14,69,70,77,78},
	  {6,7,14,15,70,71,78,79},
	  {8,9,16,17,72,73,80,81},
	  {9,10,17,18,73,74,81,82},
	  {10,11,18,19,74,75,82,83},
	  {11,12,19,20,75,76,83,84},
	  {12,13,20,21,76,77,84,85},
	  {13,14,21,22,77,78,85,86},
	  {14,15,22,23,78,79,86,87},
	  {16,17,24,25,80,81,88,89},
	  {17,18,25,26,81,82,89,90},
	  {18,19,26,27,82,83,90,91},
	  {19,20,27,28,83,84,91,92},
	  {20,21,28,29,84,85,92,93},
	  {21,22,29,30,85,86,93,94},
	  {22,23,30,31,86,87,94,95},
	  {24,25,32,33,88,89,96,97},
	  {25,26,33,34,89,90,97,98},
//	  {26,27,34,35,90,91,98,99},
	  {27,28,35,36,91,92,99,100},
//	  {28,29,36,37,92,93,100,101},
	  {29,30,37,38,93,94,101,102},
	  {30,31,38,39,94,95,102,103},
	  {32,33,40,41,96,97,104,105},
	  {33,34,41,42,97,98,105,106},
	  {34,35,42,43,98,99,106,107},
	  {35,36,43,44,99,100,107,108},
	  {36,37,44,45,100,101,108,109},
	  {37,38,45,46,101,102,109,110},
	  {38,39,46,47,102,103,110,111},
	  {40,41,48,49,104,105,112,113},
	  {41,42,49,50,105,106,113,114},
	  {42,43,50,51,106,107,114,115},
	  {43,44,51,52,107,108,115,116},
	  {44,45,52,53,108,109,116,117},
	  {45,46,53,54,109,110,117,118},
	  {46,47,54,55,110,111,118,119},
	  {48,49,56,57,112,113,120,121},
	  {49,50,57,58,113,114,121,122},
	  {50,51,58,59,114,115,122,123},
	  {51,52,59,60,115,116,123,124},
	  {52,53,60,61,116,117,124,125},
	  {53,54,61,62,117,118,125,126},
	  {54,55,62,63,118,119,126,127},
	  {64,65,72,73,128,129,136,137},
	  {65,66,73,74,129,130,137,138},
	  {66,67,74,75,130,131,138,139},
	  {67,68,75,76,131,132,139,140},
	  {68,69,76,77,132,133,140,141},
	  {69,70,77,78,133,134,141,142},
	  {70,71,78,79,134,135,142,143},
	  {72,73,80,81,136,137,144,145},
	  {73,74,81,82,137,138,145,146},
	  {74,75,82,83,138,139,146,147},
	  {75,76,83,84,139,140,147,148},
	  {76,77,84,85,140,141,148,149},
	  {77,78,85,86,141,142,149,150},
	  {78,79,86,87,142,143,150,151},
	  {80,81,88,89,144,145,152,153},
	  {81,82,89,90,145,146,153,154},
	  {82,83,90,91,146,147,154,155},
	  {83,84,91,92,147,148,155,156},
	  {84,85,92,93,148,149,156,157},
	  {85,86,93,94,149,150,157,158},
	  {86,87,94,95,150,151,158,159},
	  {88,89,96,97,152,153,160,161},
	  {89,90,97,98,153,154,161,162},
//	  {90,91,98,99,154,155,162,163},
//	  {91,92,99,100,155,156,163,164},
//	  {92,93,100,101,156,157,164,165},
	  {93,94,101,102,157,158,165,166},
	  {94,95,102,103,158,159,166,167},
	  {96,97,104,105,160,161,168,169},
	  {97,98,105,106,161,162,169,170},
	  {98,99,106,107,162,163,170,171},
	  {99,100,107,108,163,164,171,172},
	  {100,101,108,109,164,165,172,173},
	  {101,102,109,110,165,166,173,174},
	  {102,103,110,111,166,167,174,175},
	  {104,105,112,113,168,169,176,177},
	  {105,106,113,114,169,170,177,178},
	  {106,107,114,115,170,171,178,179},
	  {107,108,115,116,171,172,179,180},
	  {108,109,116,117,172,173,180,181},
	  {109,110,117,118,173,174,181,182},
	  {110,111,118,119,174,175,182,183},
	  {112,113,120,121,176,177,184,185},
	  {113,114,121,122,177,178,185,186},
	  {114,115,122,123,178,179,186,187},
	  {115,116,123,124,179,180,187,188},
	  {116,117,124,125,180,181,188,189},
	  {117,118,125,126,181,182,189,190},
	  {118,119,126,127,182,183,190,191},
	  {128,129,136,137,192,193,200,201},
	  {129,130,137,138,193,194,201,202},
	  {130,131,138,139,194,195,202,203},
	  {131,132,139,140,195,196,203,204},
	  {132,133,140,141,196,197,204,205},
	  {133,134,141,142,197,198,205,206},
	  {134,135,142,143,198,199,206,207},
	  {136,137,144,145,200,201,208,209},
	  {137,138,145,146,201,202,209,210},
	  {138,139,146,147,202,203,210,211},
	  {139,140,147,148,203,204,211,212},
	  {140,141,148,149,204,205,212,213},
	  {141,142,149,150,205,206,213,214},
	  {142,143,150,151,206,207,214,215},
	  {144,145,152,153,208,209,216,217},
	  {145,146,153,154,209,210,217,218},
	  {146,147,154,155,210,211,218,219},
	  {147,148,155,156,211,212,219,220},
	  {148,149,156,157,212,213,220,221},
	  {149,150,157,158,213,214,221,222},
	  {150,151,158,159,214,215,222,223},
	  {152,153,160,161,216,217,224,225},
	  {153,154,161,162,217,218,225,226},
	  {154,155,162,163,218,219,226,227},
	  {155,156,163,164,219,220,227,228},
	  {156,157,164,165,220,221,228,229},
	  {157,158,165,166,221,222,229,230},
	  {158,159,166,167,222,223,230,231},
	  {160,161,168,169,224,225,232,233},
	  {161,162,169,170,225,226,233,234},
	  {162,163,170,171,226,227,234,235},
	  {163,164,171,172,227,228,235,236},
	  {164,165,172,173,228,229,236,237},
	  {165,166,173,174,229,230,237,238},
	  {166,167,174,175,230,231,238,239},
	  {168,169,176,177,232,233,240,241},
	  {169,170,177,178,233,234,241,242},
	  {170,171,178,179,234,235,242,243},
	  {171,172,179,180,235,236,243,244},
	  {172,173,180,181,236,237,244,245},
	  {173,174,181,182,237,238,245,246},
	  {174,175,182,183,238,239,246,247},
	  {176,177,184,185,240,241,248,249},
	  {177,178,185,186,241,242,249,250},
	  {178,179,186,187,242,243,250,251},
	  {179,180,187,188,243,244,251,252},
	  {180,181,188,189,244,245,252,253},
	  {181,182,189,190,245,246,253,254},
	  {182,183,190,191,246,247,254,255},
	  {192,193,200,201,256,257,264,265},
	  {193,194,201,202,257,258,265,266},
	  {194,195,202,203,258,259,266,267},
	  {195,196,203,204,259,260,267,268},
	  {196,197,204,205,260,261,268,269},
	  {197,198,205,206,261,262,269,270},
	  {198,199,206,207,262,263,270,271},
	  {200,201,208,209,264,265,272,273},
	  {201,202,209,210,265,266,273,274},
	  {202,203,210,211,266,267,274,275},
	  {203,204,211,212,267,268,275,276},
	  {204,205,212,213,268,269,276,277},
	  {205,206,213,214,269,270,277,278},
	  {206,207,214,215,270,271,278,279},
	  {208,209,216,217,272,273,280,281},
	  {209,210,217,218,273,274,281,282},
	  {210,211,218,219,274,275,282,283},
	  {211,212,219,220,275,276,283,284},
	  {212,213,220,221,276,277,284,285},
	  {213,214,221,222,277,278,285,286},
	  {214,215,222,223,278,279,286,287},
	  {216,217,224,225,280,281,288,289},
	  {217,218,225,226,281,282,289,290},
	  {218,219,226,227,282,283,290,291},
	  {219,220,227,228,283,284,291,292},
	  {220,221,228,229,284,285,292,293},
	  {221,222,229,230,285,286,293,294},
	  {222,223,230,231,286,287,294,295},
	  {224,225,232,233,288,289,296,297},
	  {225,226,233,234,289,290,297,298},
	  {226,227,234,235,290,291,298,299},
	  {227,228,235,236,291,292,299,300},
	  {228,229,236,237,292,293,300,301},
	  {229,230,237,238,293,294,301,302},
	  {230,231,238,239,294,295,302,303},
	  {232,233,240,241,296,297,304,305},
	  {233,234,241,242,297,298,305,306},
	  {234,235,242,243,298,299,306,307},
	  {235,236,243,244,299,300,307,308},
	  {236,237,244,245,300,301,308,309},
	  {237,238,245,246,301,302,309,310},
	  {238,239,246,247,302,303,310,311},
	  {240,241,248,249,304,305,312,313},
	  {241,242,249,250,305,306,313,314},
	  {242,243,250,251,306,307,314,315},
	  {243,244,251,252,307,308,315,316},
	  {244,245,252,253,308,309,316,317},
	  {245,246,253,254,309,310,317,318},
	  {246,247,254,255,310,311,318,319}};

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

  }

  if (false)
    {
      std::vector<unsigned int> subdivisions (dim, 7);
      subdivisions[dim-1] = 4;

      const Point<dim> bottom_left = (dim == 2 ?
				      Point<dim>(0,0) :
				      Point<dim>(0,0,0));
      const Point<dim> top_right   = (dim == 2 ?
				      Point<dim>(7,4) :
				      Point<dim>(7,7,4));
    
      GridGenerator::subdivided_hyper_rectangle (triangulation,
						 subdivisions,
						 bottom_left,
						 top_right);
    
      std::ofstream x ("x");
      for (typename std::vector<Point<dim> >::const_iterator
	     p = triangulation.get_vertices().begin();
	   p != triangulation.get_vertices().end();
	   ++p)
	{
	  x << "Point<dim>(";
	  for (unsigned int i=0; i<3; ++i)
	    x << (*p)[i] << ',';
	  x << ")," << std::endl;
	}
    
    
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end();
	   ++cell)
	{
	  for (unsigned int i=0; i<8; ++i)
	    x << cell->vertex_index(i) << ',';
	  x << std::endl;
	}
    }
  
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->center()[1] == 0)
	cell->face(f)->set_all_boundary_indicators(1);
      else
	if (cell->at_boundary(f)
	    &&
	    (cell->face(f)->center()[0] > 0.1)
	    &&
	    (cell->face(f)->center()[0] < 6.9)
	    &&
	    (cell->face(f)->center()[1] > 0.1)
	    &&
	    (cell->face(f)->center()[1] < 6.9)
	    &&
	    (cell->face(f)->center()[2] > 0.1)
	    &&
	    (cell->face(f)->center()[2] < 3.9))
	  cell->face(f)->set_all_boundary_indicators(2);
	else
	  if (cell->at_boundary(f)
	      &&
	      (cell->face(f)->center()[2] == 0))
	    cell->face(f)->set_all_boundary_indicators(2);
  
  
  triangulation.refine_global (4-dim);

  for (unsigned int refinement_cycle = 0; refinement_cycle<5;
       ++refinement_cycle)
    {
      std::cout << "Refinement cycle " << refinement_cycle << std::endl;
      
      if (refinement_cycle > 0)
        refine_mesh ();
      
      setup_dofs ();

      std::cout << "   Assembling..." << std::endl << std::flush;
      assemble_system ();      

      std::cout << "   Solving..." << std::endl;
      Threads::ThreadGroup<> t;
      t += Threads::spawn (*this, &StokesProblem<dim>::solve)(primal_rhs, primal_constraints, primal_solution);
      t += Threads::spawn (*this, &StokesProblem<dim>::solve)(dual_rhs, dual_constraints, dual_solution);
      t.join_all ();
      
      output_results (refinement_cycle);

      std::cout << std::endl;
    }
}



int main () 
{
  try
    {
      deallog.depth_console (0);

      StokesProblem<3> flow_problem(1);
      flow_problem.run ();
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
