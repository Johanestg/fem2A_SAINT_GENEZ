#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            std::cout << "TO BE IMPLEMENTED !!!" << std::endl;

	    
	    /*** Initiation des variables pour la création de la matrice K ***/
	    Mesh M;
            M.load(mesh_filename);
            
	    ShapeFunctions fonction_forme(2, 1);
	    Quadrature q= q.get_quadrature(2, false);
	    
	    DenseMatrix Ke;
	    Ke.set_size( 3, 3 );
	    
	    SparseMatrix K(M.nb_vertices());
	    	    
	    /*** Création de la matrice K pour tous les triangles ***/
	    for (int triangle = 0; triangle < M.nb_triangles(); triangle++)
	    {
	    	ElementMapping element_map_triangle(M, false, triangle);
	    	assemble_elementary_matrix(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Ke );
	    	local_to_global_matrix( M, triangle, Ke, K );
	    }
	    
	    /*** Création de g et de F pour ensuite appliquer les conditions de Dirichlet ***/
	    std::vector<double> g;
	    std::vector<double> F;
	    for (int vert_ind =0; vert_ind < M.nb_vertices() ; vert_ind++)
	    {
	    	g.push_back(M.get_vertex( vert_ind ).x + M.get_vertex( vert_ind ).y ); /* u= x + y */
	    	F.push_back(0); /* car pas de conditions de Neumann et f=0 */
	    }
	    
	    /*** Création de attribute_is_Dirichlet ***/
	    M.set_attribute( unit_fct, 0 , false ); //affecte la valeur d'attribut 0 à tous les triangles
	    M.set_attribute( unit_fct, 1 , true ); /*affecte la valeur d'attribut 1 à tous les edge
	    					      1 correspond à la condition de Dirichlet*/
	    
	    std::vector< bool > attribute_is_dirichlet;
	    attribute_is_dirichlet.push_back( false );
	    attribute_is_dirichlet.push_back( true );
	    
	    /*** Application des conditions de Dirichlet ***/
	    apply_dirichlet_boundary_conditions( M, attribute_is_dirichlet, g, K, F );
	    
	    /*** Résolution du système ***/
	    std::vector< double > x(M.nb_vertices());
	    solve( K, F, x );
	    std::string out_name = "pure_dirichlet";
	    M.save(out_name+".mesh");
	    save_solution(x, out_name+".bb");

        }

    }

}
