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
        
        double sinus_bump_fct( vertex v)
        {
            return 2 * pow(M_PI,2) * sin(M_PI * v.x) * sin(M_PI * v.y);
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
	    	F.push_back(0); /* car pas de terme source */
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
	    std::string out_name = "pure_dirichlet_fine";
	    M.save("./simulation/"+out_name+".mesh");  /* créer un dossier simulation où stocker les mesh et la solution des simulation*/
	    save_solution(x, "./simulation/"+out_name+".bb");

        }
        
        
        
        void dirichlet_terme_source( const std::string& mesh_filename, bool verbose , bool border)
        {
        	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
            	      	
           	/*** Initiation des variables pour la création de la matrice K ***/
	    	Mesh M;
            M.load(mesh_filename);
            
	    ShapeFunctions fonction_forme(2, 1);
	    Quadrature q= q.get_quadrature(2, border);
	    
	    DenseMatrix Ke;
	    Ke.set_size( 3, 3 );
	    
	    SparseMatrix K(M.nb_vertices());
	    
	    std::vector< double > F(M.nb_vertices()); /*taille de F est nb de points d'intégration*/
	    	
	    	    
	    /*** Création de la matrice K pour tous les triangles ***/
	    for (int triangle = 0; triangle < M.nb_triangles(); triangle++)
	    {
	    	std::vector< double > Fe;
	    	ElementMapping element_map_triangle(M, false, triangle);
	    	assemble_elementary_matrix(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Ke );	
	    	local_to_global_matrix( M, triangle, Ke, K );

	    	assemble_elementary_vector(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Fe );
		local_to_global_vector(M, border, triangle, Fe, F );
	    }
	    
	    /*** Création de g et de F pour ensuite appliquer les conditions de Dirichlet ***/
	    std::vector<double> g;
	    for (int vert_ind =0; vert_ind < M.nb_vertices() ; vert_ind++)
	    {
	    	g.push_back(0); /* On impose partout u= 0 */
	    }
	    
	    /*** Création de attribute_is_Dirichlet ***/
	    M.set_attribute( unit_fct, 0 , false ); /*affecte la valeur d'attribut 0 à tous les triangles*/
	    M.set_attribute( unit_fct, 1 , true ); /*affecte la valeur d'attribut 1 à tous les edge1 correspond à la condition de Dirichlet*/
	    
	    std::vector< bool > attribute_is_dirichlet;
	    attribute_is_dirichlet.push_back( false );
	    attribute_is_dirichlet.push_back( true );
	    
	    /*** Application des conditions de Dirichlet ***/
	    apply_dirichlet_boundary_conditions( M, attribute_is_dirichlet, g, K, F );
	    
	    /*** Résolution du système ***/
	    std::vector< double > x(M.nb_vertices());
	    solve( K, F, x );
	    std::string out_name = "Dirichlet_terme_source";
	    M.save("./simulation/"+out_name+".mesh");  /* créer un dossier simulation où stocker les mesh et la solution des simulation*/
	    save_solution(x, "./simulation/"+out_name+".bb");
        }
        
        
        void sinus_bump(const std::string& mesh_filename, bool verbose , bool border)
        {
        	std::cout << "Solving a Dirichlet problem with source term" << std::endl;
            	if ( verbose ) 
            	{
                	std::cout << " with lots of printed details..." << std::endl;
            	}
            	      	
           	/*** Initiation des variables pour la création de la matrice K ***/
	    	Mesh M;
            M.load(mesh_filename);
            
	    ShapeFunctions fonction_forme(2, 1);
	    Quadrature q= q.get_quadrature(2, border);
	    
	    DenseMatrix Ke;
	    Ke.set_size( 3, 3 );
	    
	    SparseMatrix K(M.nb_vertices());
	    
	    std::vector< double > F(M.nb_vertices()); /*taille de F est nb de points d'intégration*/
	    	
	    	    
	    /*** Création de la matrice K pour tous les triangles ***/
	    for (int triangle = 0; triangle < M.nb_triangles(); triangle++)
	    {
	    	std::vector< double > Fe;
	    	ElementMapping element_map_triangle(M, false, triangle);
	    	assemble_elementary_matrix(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Ke );	
	    	local_to_global_matrix( M, triangle, Ke, K );

	    	assemble_elementary_vector(element_map_triangle, fonction_forme, q, FEM2A::Simu::sinus_bump_fct, Fe );
		local_to_global_vector(M, border, triangle, Fe, F );
	    }
	    
	    /*** Création de g et de F pour ensuite appliquer les conditions de Dirichlet ***/
	    std::vector<double> g;
	    for (int vert_ind =0; vert_ind < M.nb_vertices() ; vert_ind++)
	    {
	    	g.push_back(0); /* On impose partout u= 0 */
	    }
	    
	    /*** Création de attribute_is_Dirichlet ***/
	    M.set_attribute( unit_fct, 0 , false ); /*affecte la valeur d'attribut 0 à tous les triangles*/
	    M.set_attribute( unit_fct, 1 , true ); /*affecte la valeur d'attribut 1 à tous les edge1 correspond à la condition de Dirichlet*/
	    
	    std::vector< bool > attribute_is_dirichlet;
	    attribute_is_dirichlet.push_back( false );
	    attribute_is_dirichlet.push_back( true );
	    
	    /*** Application des conditions de Dirichlet ***/
	    apply_dirichlet_boundary_conditions( M, attribute_is_dirichlet, g, K, F );
	    
	    /*** Résolution du système ***/
	    std::vector< double > x(M.nb_vertices());
	    solve( K, F, x );
	    std::string out_name = "Sinus_bump";
	    M.save("./simulation/"+out_name+".mesh");  /* créer un dossier simulation où stocker les mesh et la solution des simulation*/
	    save_solution(x, "./simulation/"+out_name+".bb");
        
        	
        
        }
      

    }

}
