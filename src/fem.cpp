#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const /* renvoie positions globales des vertices en fonction de x et de y*/
    {
        assert( i < nb_points() ) ; //arrête le programme si l'expression est fausse
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = { // 1 point
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {  //3 points
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        if (border) 
	// Get the mapping of two vertices, one for each points of an edge.
        {	
        	for (int nb_ver=0; nb_ver<2; nb_ver++)
        	{
        		vertices_.push_back(M.get_edge_vertex(i, nb_ver));
        	}
        }
        else 
	// Get the mapping of three vertices, one for each points of a triangle.
        {
        	for (int nb_ver=0; nb_ver<3; nb_ver++)
        	{
        		vertices_.push_back(M.get_triangle_vertex(i, nb_ver));
        	}
        }

    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        vertex r ;
        if (border_) 
	/* Transforms a point from the reference element
           edge to the world space. */
        {
	        r.x = (1 - x_r.x )*vertices_[0].x + x_r.x*vertices_[1].x;
	        r.y = (1 - x_r.x )*vertices_[0].y + x_r.x*vertices_[1].y;
        }
        else
	/* Transforms a point from the reference element
           triangle to the world space. */
        {
        	r.x = (1 - x_r.x - x_r.y)*vertices_[0].x + x_r.x*vertices_[1].x + x_r.y*vertices_[2].x;
	        r.y = (1 - x_r.x - x_r.y)*vertices_[0].y + x_r.x*vertices_[1].y + x_r.y*vertices_[2].y;
        }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        DenseMatrix J ;
        if (border_)
	// Computes the jacobian matrix of the mapping of an edge.
        {
        	J.set_size(2,1); 
        	J.set(0, 0, -vertices_[0].x + vertices_[1].x);/* Derivative of the x coordinates of x_r with respect to xie */
        	J.set(1, 0, -vertices_[0].y + vertices_[1].y);/* Derivative of the y coordinates of x_r with respect to xie*/
        }
        else
	// Computes the jacobian matrix of the mapping of a triangle.
        {
        	J.set_size(2,2); 
        	J.set(0, 0, -vertices_[0].x + vertices_[1].x);/* Derivative of the x coordinates of x_r with respect to xie */
        	J.set(1, 0, -vertices_[0].y + vertices_[1].y);/* Derivative of the y coordinates of x_r with respect to xie*/
        	J.set(0, 1, -vertices_[0].x + vertices_[2].x);/* Derivative of the x coordinates of x_r with respect to eta*/
        	J.set(1, 1, -vertices_[0].y + vertices_[2].y);/* Derivative of the y coordinates of x_r with respect to eta*/
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
	// Creation of the jacobian matrix of x_r's vertex.
        DenseMatrix J = jacobian_matrix(x_r); 
	   
	double det = 0.; 
        if (border_)
	// Computes the determinant of the jacobian matrix of the mapping of an edge.
        {
        	det = sqrt( J.get(0,0)*J.get(0,0) + J.get(1,0)*J.get(1,0) );
        }
        else
	// Computes the determinant of the jacobian matrix of the mapping of a triangle.
        {
        	det = J.det_2x2(); //Calcul du déterminant
        }
        return det ;
    }

   /** 
   * \brief Get the coordinates of vertices.
   *        Used to test ElementMapping.
   */
   int ElementMapping::get_vertices()
   {
   	for (int j=0; j< vertices_.size() ; j++)
   	{
   		std::cout << "v" << j << " : " << vertices_[j].x << "	" << vertices_[j].y << std::endl;
   	}
   	return 0;
   }
   

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/

    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {  
	// Verification of values for dim_ and order_
        assert( (dim_ ==1) || (dim_ ==2) );
        assert( order_== 1 );
    }

    int ShapeFunctions::nb_functions() const
    {
        int nb_func = 0;
        if (dim_==1) nb_func = 2;
        else nb_func = 3 ;
        return nb_func ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        double shape_func=0.;
        if (dim_==1)
        {
        	if (i==0) return 1 - x_r.x;
        	else return x_r.x;
        }
        else
        {
		if (i==0) return 1 - x_r.x - x_r.y;
        	else if (i==1) return x_r.x;
        	else return x_r.y; 
        }
        return 0. ; // should not be reached

    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        vec2 g ;
       
        if (dim_==1)
        {
        	if ( i == 0)
        	{
        		g.x = -1;
        	}
        	else 
        	{
        		g.x = 1;
        	}
        	g.y = 0;	
        }
        else 
        {
        	if ( i == 0)
        	{
        		g.x = -1;
        		g.y = -1;
        	}
        	else if ( i == 1 )
        	{
        		g.x = 1;
        		g.y = 0;
        	}
        	else 
        	{
        		g.x = 0;
        		g.y = 1;
        	}	
        }      
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/

    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping, // an element: a segment or a triangle
        const ShapeFunctions& reference_functions, /* shape functions associated with vertices of the element */ 
        const Quadrature& quadrature, 
        double (*coefficient)(vertex), 
        DenseMatrix& Ke ) 
    {
        for (int i=0 ; i < reference_functions.nb_functions(); i++) /*for each shape function*/
        {
        	for (int j=0; j < reference_functions.nb_functions(); j++)
        	{
       			double sum_Ke = 0.;
        		for (int q =0 ; q < quadrature.nb_points() ; q++)
        		{
        			vertex pt_integration = quadrature.point(q); 		
        			double wq = quadrature.weight(q);
				
        			/* Jacobian matrix with respect to Gauss' points (here they are equivalent to reference points)*/
        			DenseMatrix jacob_mat = elt_mapping.jacobian_matrix( pt_integration );
        			DenseMatrix inv_jacob_mat = jacob_mat.invert_2x2();
        			DenseMatrix trans_inv_jacob_mat = inv_jacob_mat.transpose(); 
        		
        			/* Gradient of shape function of index i*/
        			vec2 grad_shape_func_i = reference_functions.evaluate_grad(i, pt_integration );
        		
        			/* Gradient of shape function of index t*/
        			vec2 grad_shape_func_j = reference_functions.evaluate_grad(j, pt_integration );
        			
        			/* Jacobian matrix determinant of elt_mapping*/
        			double det_jacob_mat = elt_mapping.jacobian( pt_integration );
 
 				/***************** Sum for Ke *****************/
 				
 				/* Scalar product */   
 				vec2 vect_je_grad_i = trans_inv_jacob_mat.mult_2x2_2( grad_shape_func_i );  
 				vec2 vect_je_grad_j = trans_inv_jacob_mat.mult_2x2_2( grad_shape_func_j );
 				
 				double prod_scal = dot( vect_je_grad_i , vect_je_grad_j );
 				
 				sum_Ke += wq * coefficient( pt_integration ) * prod_scal * det_jacob_mat ;
        		}	
        		Ke.set( i, j, sum_Ke ); 
        	}
        }
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        for (int i=0; i<3; i++)
        {
		int glob_ind1 = M.get_triangle_vertex_index( t, i);
        	for (int j=0; j<3; j++)
        	{ 
        		int glob_ind2= M.get_triangle_vertex_index( t, j) ;
        		K.add(glob_ind1, glob_ind2, Ke.get(i,j));
        	}
        }
       	//K.print();  // Undo comment to compute test_local_to_global_matrix()
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        for (int i=0; i< reference_functions.nb_functions(); i++)
        {
        	double sum_Fe=0.;
        	for (int q=0; q < quadrature.nb_points(); q++) /* nombre de point d'intégration est nb de pt quadra*/
        	{
        		vertex pt_integration = quadrature.point(q); 
        			
        		double wq = quadrature.weight(q);
        		
        		double val_shape_func = reference_functions.evaluate( i, pt_integration );
        		
        		/* Value of f at the point of Me transformation*/
        		double f = source( pt_integration);
        		
        		/*Jacobian matrix determinant of elt_mapping*/
        		double det_jacob_mat = elt_mapping.jacobian( pt_integration );
        		
        		sum_Fe += wq * val_shape_func * f * det_jacob_mat;
        	}
        	Fe.push_back( sum_Fe );
    	}
    }
    

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        for (int i=0; i< reference_functions_1D.nb_functions(); i++)
        {
        	double sum_Fe=0.;
        	for (int q=0; q < quadrature_1D.nb_points(); q++) 
		/* The number of integration points are the number of quadrature points*/
        	{
        		vertex pt_integration = quadrature_1D.point(q); 
        			
        		double wq = quadrature_1D.weight(q);
        		
        		double val_shape_func = reference_functions_1D.evaluate( i, pt_integration );
        		
        		/* Value of h at the transformation point Me*/
        		double h = neumann( pt_integration);
        		
        		/* Jacobian matrix determinant of elt_mapping_1D*/
        		double racine_jacob = elt_mapping_1D.jacobian( pt_integration);
        		sum_Fe += wq * val_shape_func * h * racine_jacob;
        	}
        	Fe.push_back( sum_Fe );
    	}
    }
    	

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::vector<int> ind_global; 
        std::cout << "Assemblage: " << std::endl;
        if (border)
        {
        	for (int ind_local =0; ind_local < Fe.size(); ind_local++)
        	{
        		F[ M.get_edge_vertex_index( i, ind_local) ] = Fe[ind_local];
        		std::cout << "le numéro global du point " << ind_local << " est " << M.get_edge_vertex_index( i, ind_local) << "." << std::endl;
        	}
        }
        else
        {
        	for (int ind_local =0; ind_local < Fe.size(); ind_local++)
        	{
        		F[ M.get_triangle_vertex_index( i, ind_local) ] = Fe[ind_local];
        		std::cout << "le numéro global du point " << ind_local << " est " << M.get_triangle_vertex_index( i, ind_local) << "." << std::endl;
        	}
        }
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        /*vector of size of values filled with false*/
       	std::vector<bool> processed_vertices(values.size(), false); 
       	/* penalty_coefficient is a big coefficient*/
        double penalty_coefficient =10000.;
        
        for (int i=0; i< M.nb_edges(); i++)
        {
        	/* For each edge, we get the values of the attribute*/
        	int edge_attribute = M.get_edge_attribute(i);
        	
        	/* 
        	In attribut_is_diriclet, there is the information about if it's an edge or not.
	 	True if it's an edge, false otherwise.
        	 */
        	if (attribute_is_dirichlet[edge_attribute])
        	{
        	/* For each vertices associated with the edge*/
        	    for( int v = 0; v < 2; v++ ) 
		    {
                        int vertex_index = M.get_edge_vertex_index(i, v);
                        if( !processed_vertices[vertex_index] ) 
		        {
                            processed_vertices[vertex_index] = true;
                            K.add(vertex_index, vertex_index, penalty_coefficient);
                            F[vertex_index] += penalty_coefficient*values[vertex_index]; /* values correspond à g*/
                        }  
                    }
               }
        }
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
