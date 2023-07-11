//
//  main.cpp
//  Coulomb-Energy-Minimizer
//
//  Naive implementation of Coulomb Shapes
//  in O(n^2) time.

#include <math.h>
#include <string>
#include <iostream>

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/edges.h>

#include <Eigen/Core>

#include <cppoptlib/function.h>
#include <cppoptlib/solver/solver.h>
#include <cppoptlib/solver/bfgs.h>
#include <cppoptlib/solver/gradient_descent.h>
#include <cppoptlib/solver/newton_descent.h>
#include <cppoptlib/solver/conjugated_gradient_descent.h>

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

// ===============================================================================================================================================
// ===============================================================================================================================================
// CONSTANTS (TO BE CHANGED) TODO: convert these to arguments
const int STARTING_ITER = 0;            // Doesn't affect functionality, just to save iteration points (TODO: remove)
const int CHECKPOINT_ITER = 100;        // Save the result every CHECKPOINT_ITER iterations
const int T_STEPS = 300;                // Number of iterations to optimize the objective function
const std::string OBJ_MESH_PATH = "../../../test meshes/kid01_decimated.obj";   // Path to the input mesh .obj
// ===============================================================================================================================================
// ===============================================================================================================================================


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

// TODO: squared is never needed ig, it can be either removed or utilized.
double get_coulomb_energy(const VectorXd& X, bool squared = false, bool get_inverse = false){
    
    double distance;
    double tot_energy = 0.0;
    size_t num_particles = X.rows() / 3;
    
    for(size_t i = 0; i < num_particles; i++){
        for(int j = 0; j < num_particles; j++){
            if(i == j)
                continue;
            
            VectorXd p1{{ X(i), X(i + 1), X(i + 2)}}; // {x,y,z} coordinates of p1
            VectorXd p2{{ X(j), X(j + 1), X(j + 2)}}; // {x,y,z} coordinates of p2
            
            if(squared)
                distance = (p2 - p1).squaredNorm();
            else
                distance = (p2 - p1).norm();
            
            if(get_inverse){
                distance += __DBL_EPSILON__;    // For numerical stability
                distance = 1.0 / distance;
            }
            tot_energy += distance;
        }
    }
    return tot_energy;
}


// --------------------------------------------------------------------------------------------------------------

class CoulombEnergy : public FunctionXd {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    // Definition the Coulomb Energy in Eqn.1
    double operator()(const Eigen::VectorXd &X) const override{
        
        //VectorXd inv_edge_lengths(E.rows());
        //get_edge_lengths(X, E, inv_edge_lengths, false, true);
        
        //double energy = inv_edge_lengths.sum();
        //return energy;
        return get_coulomb_energy(X, false, true);
    }
    
    // For computing gradients numerically, this function can be commented
    // However, a single iteration of numerical solver takes more than 40 minutes,
    // whereas provided gradients run 10 iterations under a minute.
    void Gradient(const vector_t &X, vector_t *grad) const override {
        
        // assert(X.rows() % 3 == 0);        // X should contain 3 coordinates for every vertex
        // assert(X.rows() == grad->rows()); // ensured by library
        
        size_t N = X.rows() / 3;          // number of vertices
        size_t idx, idx2;                 // used for finding the two points connected by an edge
                
        // Loop over every vertex
        for(size_t i = 0; i < N; i++){
            
            idx = 3*i;          // index of the x coordinate of the vertex_i
            
            int k = 0;          // counter for the loop below
            (*grad)[idx] = 0;   // d_xix
            (*grad)[idx+1] = 0; // d_xiy
            (*grad)[idx+2] = 0; // d_xiz
            double edge_len = 0.0;
            double cache = 0.0; // cache to be utilized in multiple derivatives
            for(size_t j = 0; j < N; j++){
                if(i == j)
                    continue;
                
                idx2 = 3 * j;     // index of the x coordinate of the neighbour vertex_j
                
                // TODO: Just slice the vector using Eigen::block?
                VectorXd p1{{ X(idx), X(idx + 1), X(idx + 2)}};    // {x,y,z} coordinates of p1
                VectorXd p2{{ X(idx2), X(idx2 + 1), X(idx2 + 2)}}; // {x,y,z} coordinates of p2
                
                edge_len = (p2 - p1).squaredNorm(); // I used squared norm in the derivation
                cache = -0.5 * pow(edge_len, -1.5);
                
                // Computing the gradients for each x,y,z coordinate of the vertex_i, over the neighbour vertex_j
                (*grad)[idx]   += cache * (2 * (X(idx)   - X(idx2)));       // gradient for x coordinate
                (*grad)[idx+1] += cache * (2 * (X(idx+1) - X(idx2+1)));     // gradient for y coordinate
                (*grad)[idx+2] += cache * (2 * (X(idx+2) - X(idx2+2)));     // gradient for z coordinate
                
                k++;
            }
        }
    }
};
// --------------------------------------------------------------------------------------------------------------


int main(int argc, char const *argv[]) {
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
    
    // ------------------------------------------------------------------------------------------------------
    //                                      OPTIMIZATION (Eqn.6 and Eqn.4)
    // ------------------------------------------------------------------------------------------------------
    using Solver = cppoptlib::solver::ConjugatedGradientDescent<CoulombEnergy>;
    
    CoulombEnergy coulomb_energy;
    // TODO: How to get rid of matrix to vector conversions (and stick to a single notation either matrix or vector)?
    // Convert V to X, where V is a matrix and X is a vector, flattened version of V
    RowMatrixXd V_RowMajor(V);                                     // Eigen matrices are Column-Major by default
    Map<RowVectorXd> V_map(V_RowMajor.data(), V_RowMajor.size());  // Flatten Eigen matrix into a vector
    VectorXd x = V_map;
    

    // Get the maximum edge length over the entire mesh,
    // to be utilized in step size, aka c(t) of eqn.4
    VectorXd edge_lengths(E.rows());
    get_edge_lengths(x, E, edge_lengths);
    double max_d =edge_lengths.maxCoeff();
    
    // Iteratively compute Eqn.4 and 6 for T_STEPS
    for(int t = 0; t < T_STEPS; t++){
        // Evaluate energy function
        auto state = coulomb_energy.Eval(x, 1);
        
        // Print objective function's value
        if((t+1) % 1 == 0){
            std::cout << "iter " << (t+1) << " - energy :" << state.value << "\n";
        }
        
        // Eqn. 4
        auto grad = state.gradient;
        double c_t = 0.1 * (max_d / grad.maxCoeff());   // step-size heuristic used by the paper
        x = state.x - c_t * grad;
        
        // Eqn. 6
        // I tried disabling the constraint but then the result is garbage.
        MatrixXd x_mat = RowMatrixXd::Map(x.data(), int(x.rows()/3), 3);
        apply_edge_constraint(x_mat, adjacency_list);
        
        // TODO: adjust apply_edge_constraints to work on flat vectors
        RowMatrixXd x_mat_rowmajor(x_mat);
        Map<RowVectorXd> x_map(x_mat_rowmajor.data(), x_mat_rowmajor.size());
        x = x_map;
        
        if((t+1) % CHECKPOINT_ITER == 0){
            MatrixXd V_new = RowMatrixXd::Map(x.data(), int(x.rows()/3), 3);
            igl::writeOBJ("../../../result_" + std::to_string(STARTING_ITER+t+1) + ".obj", V_new, F);
        }
    }
    
    // Unflatten X vector back to V matrix
    // TODO: stick to either vector or matrix representation?
    MatrixXd V_new = RowMatrixXd::Map(x.data(), int(x.rows()/3), 3);
    
    // ------------------------------------------------------------------------------------------------------
    //                                          OUTPUT MESH
    // ------------------------------------------------------------------------------------------------------
    // Write the output mesh into .obj
    // See result.obj created in the project folder
    igl::writeOBJ("../../../result.obj", V_new, F);
    
    // TODO: Run this once to verify your gradients --> Result: correct :) (for decimated mesh)
    // Note: Running the line below takes about ~1.5 hours, depending on the mesh size.
    // This is due to the speed of numerical computation of gradients.
    // std::cout << "Is gradient correct? " << cppoptlib::utils::IsGradientCorrect(coulomb_energy, x) << std::endl;
    
    return EXIT_SUCCESS;
}

