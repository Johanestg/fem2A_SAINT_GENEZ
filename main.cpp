#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A; /*pour utiliser les fct, les classes définit dans cet espace.*/


/* Fonction pour lancer les tests*/
void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quadra = false; 
    const bool t_constructeur_element_mapping = false;
    const bool t_transform_element_mapping = false;
    const bool t_jacobian_matrix_elementmapping = true;
    const bool t_constructeur_shapefunction = false;
    const bool t_nbfunction_shapefunction = false;
    const bool t_evaluate_shapefunction = false;
    const bool t_evaluategrad_shapefunction = false;
    const bool t_assemble_elementary_matrix = false;
    const bool t_local_to_global_matrix= false;
    const bool t_assemble_elementary_vector= false;
    const bool t_local_to_global_vector= false;
    const bool t_assemble_elementary_neumann_vector= false;
    
    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if (t_quadra) Tests::test_quadrature(2,true);
    if (t_constructeur_element_mapping) Tests::test_constructeur_elementmapping();
    if (t_transform_element_mapping) Tests::test_transform_elementmapping();
    if (t_jacobian_matrix_elementmapping) Tests::test_jacobian_matrix_elementmapping();
    if (t_constructeur_shapefunction) Tests::test_constructeur_shapefunction(2,1);
    if (t_nbfunction_shapefunction) Tests::test_nbfunction_shapefunction(1,1);
    if (t_evaluate_shapefunction) Tests::test_evaluate_shapefunction(1, 1, 0);
    if (t_evaluategrad_shapefunction) Tests::test_evaluategrad_shapefunction(1, 1, 0);
    if (t_assemble_elementary_matrix) Tests::test_assemble_elementary_matrix();
    if (t_local_to_global_matrix) Tests::test_local_to_global_matrix();
    if (t_assemble_elementary_vector) Tests::test_assemble_elementary_vector(false); /* pour segment ou triangle*/
    if (t_local_to_global_vector) Tests::test_local_to_global_vector(true); /* pour segment ou triangle*/
    if (t_assemble_elementary_neumann_vector) Tests::test_assemble_elementary_neumann_vector(false);
    
} /* comme les trois true, on lance les trois fonctions test */

/* Fonction pour lancer la simulation*/
void run_simu()
{

    const bool simu_pure_dirichlet = false;
    const bool simu_dirichlet = false;
    const bool simu_dirichlet_sinus_bump = true;


    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) Simu::pure_dirichlet_pb("data/square_fine.mesh", verbose);
    if (simu_dirichlet) Simu::dirichlet_terme_source("data/square_fine.mesh", verbose, false);
    if (simu_dirichlet_sinus_bump) Simu::sinus_bump("data/square_fine.mesh", verbose, false);
}


/* !!!!!!!!! DÉBUT DU MAIN !!!!!!!!!*/

int main( int argc, const char * argv[] ) /* argc est le nombre d'argument et argv est un 							pointeur sur le premier argument d'un tableau de char*/
{
    /* Command line parsing : argument est le vecteur de argv */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }
    

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
