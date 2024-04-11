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

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            /*std::cout << "Vertices <x> <y> <att>" << std::endl;
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
            }*/

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }

	bool test_quadrature(int ordre, bool bord)
	{
		Quadrature q= q.get_quadrature(ordre, bord);/*On déclare un espace quadrature (Quadrature q) et on définie q en lui donnant des valeurs true quand on a un segment, false quand on a un triangle */
		std::cout << "nombre de points: " << q.nb_points() << std::endl;
		double somme=0;
		for (int i=0; i< q.nb_points(); i++)
		{
			std::cout << "poid de chaque point: " << q.weight(i) << std::endl;
			somme+=q.weight(i);
		}
		std::cout << "somme: " << somme << std::endl;
		return true;
	}
	
	bool test_constructeur_elementmapping()
	{
		Mesh mesh;
		mesh.load("data/square.mesh");
		ElementMapping element_map_triangle(mesh, false, 4);
		std:: cout << "vertices= " << element_map_triangle.get_vertices() << std::endl;
		
		
		return true; 
		
	}
    }
}
