/**
 * Main.cc
 *
 *  @date Created on: Dec 19, 2016
 *  @author Keith D Brauss
 */

#include<iostream>
#include<vector>
#include<cmath>

#include "Point.h"
#include "Util.h"
#include "Potential.h"
#include "FmmTree.h"

// bartu edit
#include <math.h>
#include <string>
#include <chrono>

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>

#include <Eigen/Core>

#include "cppoptlib/function.h"
#include "cppoptlib/solver/solver.h"
#include "cppoptlib/solver/bfgs.h"
#include "cppoptlib/solver/gradient_descent.h"
#include "cppoptlib/solver/newton_descent.h"
#include "cppoptlib/solver/conjugated_gradient_descent.h"

using namespace cppoptlib;
using namespace Eigen;

typedef Matrix<double,Dynamic,Dynamic,RowMajor> RowMatrixXd;
using FunctionXd = cppoptlib::function::Function<double>;

struct adjacency{
    std::vector<int> connected_vertex_indices;
    std::vector<double> connected_edge_lengths;
};


MatrixXd V;         // Vertices
MatrixXi F;         // Faces
MatrixXd E;         // Edges
std::vector<adjacency> adjacency_list;

// CONSTANTS (TO BE CHANGED) TODO: convert these to arguments
const int STARTING_ITER = 0;            // Doesn't affect functionality, just to save iteration points (TODO: remove)
const int CHECKPOINT_ITER = 10;        // Save the result every CHECKPOINT_ITER iterations
const int T_STEPS = 300;                // Number of iterations to optimize the objective function
const std::string OBJ_MESH_PATH = "/Users/bartu/Desktop/Coulomb-FFM/test meshes/16_decimated.obj";   // Path to the input mesh .obj

// --------------------------------------------------------------------------------------------------------------

// Get the edge lengths according to adjacency list
std::vector<adjacency> edge_length_adjacency(const MatrixXd& V, const MatrixXd& E){
    std::vector<adjacency> adjacency_list(V.rows());
    
    for(int e_idx = 0; e_idx < E.rows(); e_idx++){
        
        // Get the indices in V, of two points of an edge
        int v1_idx = E(e_idx,0);
        int v2_idx = E(e_idx, 1);
        
        // Update the adjacency lists
        adjacency_list[v1_idx].connected_vertex_indices.insert(adjacency_list[v1_idx].connected_vertex_indices.end(),
                                                               v2_idx);
        
        adjacency_list[v2_idx].connected_vertex_indices.insert(adjacency_list[v2_idx].connected_vertex_indices.end(),
                                                               v1_idx);
        
        // Get edge length
        double edge_length = (V.row(v1_idx) - V.row(v2_idx)).norm();
        
        // Insert edge length
        // TODO: remove?
        adjacency_list[v1_idx].connected_edge_lengths.insert(adjacency_list[v1_idx].connected_edge_lengths.end(),
                                                             edge_length);
        
        adjacency_list[v2_idx].connected_edge_lengths.insert(adjacency_list[v2_idx].connected_edge_lengths.end(),
                                                             edge_length);
        
    }
    return adjacency_list;
}
// --------------------------------------------------------------------------------------------------------------
// Eqn.6
void apply_edge_constraint(MatrixXd &X, std::vector<adjacency> &adjacency_list){
    
    // Loop over the vertices
    for(int i = 0; i < X.rows(); i++){
        
        adjacency adj = adjacency_list[i];
        int v_i = int(adj.connected_vertex_indices.size());
        VectorXd sum_vec(3);
        sum_vec << 0, 0, 0;
        
        VectorXd x_i = X.row(i);
        // Sum over neighbours
        for(int j = 0; j < v_i; j++){
            VectorXd x_j = X.row(adj.connected_vertex_indices[j]);
            VectorXd diff_ij = x_i - x_j;
            double d_ij = adj.connected_edge_lengths[j];
            
            // Eqn. 6, inside the paranthesis
            sum_vec += x_j + d_ij * ( diff_ij / diff_ij.norm() );
        }
        // Take the average of sum
        X.row(i) = sum_vec / v_i;
    }
    return;
}

// --------------------------------------------------------------------------------------------------------------
// X is the flat vector containing vertices (x1, y1, z1, x2, y2, z2,...) where x_i is the x coordinate of vertex_i
// E is the edge vector containing pair of vertices
// edge_lengths is an empty vector with num_edges size, to be filled
// squared is a boolean, to compute squared norm of edges intead of actual norm
// TODO: do we really need get_inverse? doesn't Eigen provide element-wise operations like a = 1 / a where a is a flat vector?
void get_edge_lengths(const VectorXd& X, const MatrixXd& E, VectorXd& edge_lengths, bool squared = false, bool get_inverse = false){
    
    for(int i = 0; i < E.rows(); i++){
        
        int p1_begin = int(E(i,0)) * 3;
        int p2_begin = int(E(i,1)) * 3;
        
        VectorXd p1{{ X(p1_begin), X(p1_begin + 1), X(p1_begin + 2)}}; // {x,y,z} coordinates of p1
        VectorXd p2{{ X(p2_begin), X(p2_begin + 1), X(p2_begin + 2)}}; // {x,y,z} coordinates of p2
        
        if(squared)
            edge_lengths(i) = (p2 - p1).squaredNorm();
        else
            edge_lengths(i) = (p2 - p1).norm();
        
        if(get_inverse)
            edge_lengths(i) = 1.0 / edge_lengths(i);
    }
    
    return;
}

// Convert VectorXd Eigen vector to std::vector<Point>
void eigen2point(const VectorXd& X_vector, std::vector<Point>& X_points){
    
    // Assert the size of vectors
    assert(X_vector.rows() / 3 == X_points.size());
    
    size_t num_points = X_points.size();
    for(size_t i = 0; i < num_points; i++){
        X_points[i].setX(X_vector(i * 3));
        X_points[i].setY(X_vector(i * 3 + 1));
        X_points[i].setZ(X_vector(i * 3 + 2));
    }
}

// end of bartu edit


int main()
{
    // ------------------------------------------------------------------------------------------------------
    //                                        INPUT MESH PRECOMPUTATION
    // ------------------------------------------------------------------------------------------------------
    
    // Read input mesh
    //igl::readOBJ("../../../result.obj", V, F);    // Use this to resume where you left
    igl::readOBJ(OBJ_MESH_PATH, V, F);
    
    // Scale for numerical stabilitiy
    // V = V * 1000;
    
    igl::edges(F, E);
    
    // Precompute inverse edge lengths
    adjacency_list = edge_length_adjacency(V, E);
    
    // TODO: How to get rid of matrix to vector conversions (and stick to a single notation either matrix or vector)?
    // Convert V to X, where V is a matrix and X is a vector, flattened version of V
    RowMatrixXd V_RowMajor(V);                                     // Eigen matrices are Column-Major by default
    Map<RowVectorXd> V_map(V_RowMajor.data(), V_RowMajor.size());  // Flatten Eigen matrix into a vector
    VectorXd x_eigen = V_map;
    
    // Get the maximum edge length over the entire mesh,
    // to be utilized in step size, aka c(t) of eqn.4
    VectorXd edge_lengths(E.rows());
    get_edge_lengths(x_eigen, E, edge_lengths);
    double max_d =edge_lengths.maxCoeff();
    
  // ------------------------------------------------------------------------------------------------------
  //                                                F F M
  // ------------------------------------------------------------------------------------------------------

  // note that the DEFAULT_NUM_LEVEL count starts on l = 1
  // Example:  If the default number of FMM levels is 5
  //           Then the levels are l = 0, l = 1, l = 2, l = 3, and l = 4
  int DEFAULT_NUM_LEVEL = 4;
  // abs_alpha is also the order_of_approximation - highest order of derivatives used in Taylor series approximation
  unsigned int abs_alpha = 4;
  // number of Taylor series terms used to approximate the potential function
  int p = (abs_alpha+1)*(abs_alpha+1);

  std::vector<Point>  x, y;
  std::vector<double> u;

  // Below we use 8 source particles (2 in each coordinate direction)
  // per cell with pow(8,L-1) cells/boxes partitioning the domain
  size_t nSourceParticles = V.rows(); // pow(2, 0);
  x.resize(nSourceParticles);
 
  // One target point
  y.resize(nSourceParticles);

  // u holds the charges of source particles TODO: i forgot how to initialize all elements to same..
  u.resize(nSourceParticles);
  for(size_t i = 0; i < nSourceParticles; i++){
      u[i] = 1.0;
  }

  // ********************************************************************
  // Creating the Source Points *****************************************
  // ********************************************************************
    
    eigen2point(x_eigen, x);
    eigen2point(x_eigen, y);

   // ********************************************************************
  // Creating the Target Points *****************************************
  // ********************************************************************
  // Declaring the Potential Function needed by fmmtree data structure
  std::cout << ">> We are here before constructor for the potential function" << std::endl;
  Potential potential(p,abs_alpha);

  std::cout << ">> We are here before fmmtree construct" << std::endl;
  std::cout << ">> The number of source particles is " << x.size() << std::endl;
  std::cout << ">> The number of target particles is " << y.size() << std::endl;
  std::cout << ">> The potential function has p = " << potential.getP() << std::endl;
  std::cout << ">> The potential function has order_of_approximation = " << potential.getOrderApprox() << std::endl;

  // Declaring the FMM Tree data structure
  auto start = std::chrono::steady_clock::now();
  FmmTree fmmtree(DEFAULT_NUM_LEVEL,x,y,potential);
  auto end = std::chrono::steady_clock::now();
  std::cout << ">>> FFM tree construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    
    
    /*
  start = std::chrono::steady_clock::now();
  std::vector<double> direct = fmmtree.solveDirect(u);
  end = std::chrono::steady_clock::now();
    
  std::cout << ">>> Direct solver took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    */
    
    

    /*
    start = std::chrono::steady_clock::now();
    std::vector<std::vector<double> > direct_grad = fmmtree.solveDirectGrad(u);
    end = std::chrono::steady_clock::now();
    
    std::cout << ">>> Direct gradient took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
   
    */
    Eigen::VectorXd copied_x(nSourceParticles*3);
    Eigen::VectorXd copied_grad(nSourceParticles*3);
    // Iteratively compute Eqn.4 and 6 for T_STEPS
    for(int t = 0; t < T_STEPS; t++){
        std::cout << "Step " << t << "...\n";
        // Evaluate energy function
        std::vector<std::vector<double>> direct_grad = fmmtree.solveDirectGrad(u);
        
        //copied_x = Eigen::Map<RowVectorXd>(&x[nSourceParticles][3], nSourceParticles*3 );
        //copied_grad = Eigen::Map<RowVectorXd>(&direct_grad[nSourceParticles][3], nSourceParticles*3 );
        // TODO: why doesn't mapping of copied_x above work?
        // i combined my previous code on Eigen with this FMM library but now it does unnecessary copying...
        for(size_t i = 0; i < nSourceParticles ; i++){
            copied_x(i * 3) = x[i].getX();
            copied_x(i * 3 + 1) = x[i].getY();
            copied_x(i * 3 + 2) = x[i].getZ();

            copied_grad(i * 3) = direct_grad[0][i];
            copied_grad(i * 3 + 1) = direct_grad[1][i];
            copied_grad(i * 3 + 2) = direct_grad[2][i];
        }
        double c_t = 0.1  * (max_d / copied_grad.maxCoeff());   // step-size heuristic used by the paper
        copied_x = copied_x - c_t * copied_grad;
        
        // Eqn. 6
        // I tried disabling the constraint but then the result is garbage.
        MatrixXd x_mat = RowMatrixXd::Map(copied_x.data(), int(copied_x.rows()/3), 3);
        apply_edge_constraint(x_mat, adjacency_list);
        
        // TODO: adjust apply_edge_constraints to work on flat vectors
        RowMatrixXd x_mat_rowmajor(x_mat);
        Map<RowVectorXd> x_map(x_mat_rowmajor.data(), x_mat_rowmajor.size());
        copied_x = x_map;
        
        // convert them back TODO: get rid of conversions...
        eigen2point(x_eigen, x);
        eigen2point(x_eigen, y);
        
        if((t+1) % CHECKPOINT_ITER == 0){
            MatrixXd V_new = RowMatrixXd::Map(copied_x.data(), int(copied_x.rows()/3), 3);
            igl::writeOBJ("../../../result_" + std::to_string(STARTING_ITER+t+1) + ".obj", V_new, F);
        }
    }
    
    // Unflatten X vector back to V matrix
    // TODO: stick to either vector or matrix representation?
    MatrixXd V_new = RowMatrixXd::Map(copied_x.data(), int(copied_x.rows()/3), 3);
    
    // ------------------------------------------------------------------------------------------------------
    //                                          OUTPUT MESH
    // ------------------------------------------------------------------------------------------------------
    // Write the output mesh into .obj
    // See result.obj created in the project folder
    igl::writeOBJ("../../../result.obj", V_new, F);
    
    
    
  std::cout << ">> Finished." << "\n";

    
    return EXIT_SUCCESS;
}
