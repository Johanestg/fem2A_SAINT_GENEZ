#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

#include "simu.h"

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }


    /****************************************************************/
    /* Tests for Quadrature */
    /****************************************************************/

	bool test_quadrature(int ordre, bool bord)
	{
		/* Declare a quadrature space */
		Quadrature q= q.get_quadrature(ordre, bord);
		std::cout << "nombre de points: " << q.nb_points() << std::endl;
		double somme=0;
		for (int i=0; i< q.nb_points(); i++)
		{
			std::cout << "poid de chaque point: " << q.weight(i) << std::endl;
			somme+=q.weight(i);
		}
		std::cout << "Somme des poids: " << somme << std::endl;
		return true;
 	}


    /****************************************************************/
    /* Tests for ElementMapping */
    /****************************************************************/

	bool test_constructeur_elementmapping()
	{
		Mesh mesh;
		mesh.load("data/square.mesh");
		// Test for triangle
		ElementMapping element_map_triangle(mesh, false, 4);
		element_map_triangle.get_vertices();
		// Test for an edge
		ElementMapping element_map_edge(mesh, true, 4);
		element_map_edge.get_vertices();
		
		return true; 
		
	}
	
	bool test_transform_elementmapping()
	{
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Test for a triangle
		ElementMapping element_map_triangle(mesh, false, 4);
		vertex ref_v;
		ref_v.x=0.2;
		ref_v.y=0.4;
		vertex global_v= element_map_triangle.transform(ref_v);
		std::cout << "Mapping du point:" << std::endl << global_v.x << " " << global_v.y << std::endl;
		// Test for an edge
		ElementMapping element_map_edge(mesh, true, 4);
		vertex global_v2= element_map_edge.transform(ref_v);
		std::cout << "Mapping du point:" << std::endl << global_v2.x << " " << global_v2.y << std::endl;
	
		return true; 
		
	}
	
	bool test_jacobian_matrix_elementmapping()
	{
		// Load of a mesh
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Test for a triangle
		ElementMapping element_map_triangle(mesh, false, 4);
		vertex ref_v;
		ref_v.x=0.2;
		ref_v.y=0.4;
		DenseMatrix jacob_matrix= element_map_triangle.jacobian_matrix(ref_v);
		// Jacobian matrix of the reference point for a triangle
		jacob_matrix.print();
		// Determinant of jacobian matrix
		std::cout << element_map_triangle.jacobian(ref_v) << std::endl;
		
		// Test for an edge
		ElementMapping element_map_edge(mesh, true, 4);
		DenseMatrix jacob_matrix2= element_map_edge.jacobian_matrix(ref_v);
		// Jacobian matrix of the reference point for an edge
		jacob_matrix2.print();
		// Determinant of jacobian matrix
		std::cout << element_map_edge.jacobian(ref_v) << std::endl;
		
		return true; 
		
	}


    /****************************************************************/
    /* Tests for ShapeFunctions */
    /****************************************************************/

	bool test_constructeur_shapefunction(int dim, int order)
	{
		ShapeFunctions fonction_forme(dim, order );
		
		return true; 
	}
	
	bool test_nbfunction_shapefunction(int dim, int order)
	{
		ShapeFunctions fonction_forme(dim, order );
		int nb_func = fonction_forme.nb_functions();
		std::cout << nb_func << std::endl;
		
		return true; 
	}
	
	bool test_evaluate_shapefunction(int dim, int order, int i)
	{
		ShapeFunctions fonction_forme(dim, order );
		vertex ref_v;
		ref_v.x=0.2;
		ref_v.y=0.4;
		
		double eval_shape_func = fonction_forme.evaluate(i, ref_v);
		std::cout << eval_shape_func << std::endl;
		
		return true; 
	}
	
	bool test_evaluategrad_shapefunction(int dim, int order, int i)
	{
		ShapeFunctions fonction_forme(dim, order );
		vertex ref_v;
		ref_v.x=0.2;
		ref_v.y=0.4;
		
		vec2 evalgrad_shape_func = fonction_forme.evaluate_grad(i, ref_v);
		std::cout << evalgrad_shape_func.x << std::endl << evalgrad_shape_func.y << std::endl;
		
		return true; 
	}


    /****************************************************************/
    /* Tests for Finite Element functions */
    /****************************************************************/
	
	bool test_assemble_elementary_matrix()
	{
		// Load a mesh
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Mapping for a triangle
		ElementMapping element_map_triangle(mesh, false, 4);
		
		// Shape function 
		ShapeFunctions fonction_forme(2, 1);
		
		// Quadrature
		Quadrature q= q.get_quadrature(2, false);
		
		// Ke
		DenseMatrix Ke;
		Ke.set_size( 3, 3 );
				
		assemble_elementary_matrix(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Ke ); 
		return true;
	}
	
	bool test_local_to_global_matrix()
	{
		// Load a mesh
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Mapping for a triangle
		ElementMapping element_map_triangle(mesh, false, 4);
		
		// Shape function 
		ShapeFunctions fonction_forme(2, 1);
		
		// Quadrature
		Quadrature q= q.get_quadrature(2, false);
		
		// Crée Ke et calcul ses valeurs
		DenseMatrix Ke;
		Ke.set_size( 3, 3 );
		assemble_elementary_matrix(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Ke );
				
		// Crée K et lui donne ses valeurs pour l'élément t
		int t = 4;
		SparseMatrix K(mesh.nb_vertices());		
		local_to_global_matrix( mesh , t, Ke, K );
		return true;
	}
	
	bool test_assemble_elementary_vector(bool border)
	{
		// Charge un mesh
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Test pour un triangle
		ElementMapping element_map_triangle(mesh, border, 4);
		
		// Shape function 
		int dim =0;
		int& ref_dim= dim;
		if (border) dim=1;
		else dim=2;
		ShapeFunctions fonction_forme(dim, 1);
	
		// Quadrature
		Quadrature q= q.get_quadrature(2, border);
		
		// Fe
		std::vector< double > Fe;	
		assemble_elementary_vector(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Fe ); 
		
		for (int i = 0; i< Fe.size() ; i++)
		{
			std::cout << Fe[i] << std::endl;
		}
		return true;
	}
	
	bool test_local_to_global_vector(bool border)
	{
		// Charge un mesh
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Test pour un triangle
		ElementMapping element_map_triangle(mesh, border, 4);
		
		// Shape function 
		int dim =0;
		int& ref_dim= dim;
		if (border) dim=1;
		else dim=2;
		ShapeFunctions fonction_forme(dim, 1);
		
		// Quadrature
		Quadrature q= q.get_quadrature(2, border);
		
		// Fe
		std::vector< double > Fe;			
		assemble_elementary_vector(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Fe );
				
		// Crée F et lui donne ses valeurs pour l'élément i
		int i = 4;
		std::vector< double > F( q.nb_points()); /*taille de F est nb de points d'intégration*/	
		local_to_global_vector(mesh, border, i, Fe, F );
		return true;
	}
	
	bool test_assemble_elementary_neumann_vector(bool border)
	{
		// Charge un mesh
		Mesh mesh;
		mesh.load("data/square.mesh");
		
		// Test pour un triangle
		ElementMapping element_map_triangle(mesh, border, 4);
		
		// Shape function 
		int dim =0;
		int& ref_dim= dim;
		if (border) dim=1;
		else dim=2;
		ShapeFunctions fonction_forme(dim, 1);
	
		// Quadrature
		Quadrature q= q.get_quadrature(2, border);
		
		// Fe
		std::vector< double > Fe;	
		assemble_elementary_neumann_vector(element_map_triangle, fonction_forme, q, FEM2A::Simu::unit_fct, Fe ); 
		
		for (int i = 0; i< Fe.size() ; i++)
		{
			std::cout << Fe[i] << std::endl;
		}
		return true;
	}
    }
}
