#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>



#include "RotPriorHTypes.h"

namespace RotPriorH
{

    
    /** This struct contains the various parameters that control the algorithm 
        Based on: SESync struct*/
    struct HomographyEstimationOptions {

      /// OPTIMIZATION STOPPING CRITERIA
      /** Stopping tolerance for the norm of the Riemannian gradient */
      double tol_grad_norm = 1.e-9;

       /** Stopping criterion based upon the norm of an accepted update step */
      double preconditioned_grad_norm_tol = 1.e-9;

      /** Stopping criterion based upon the relative decrease in function value */
      double rel_func_decrease_tol = -1e-10;

       /** Stopping criterion based upon the norm of an accepted update step */
      double stepsize_tol = 1e-05;

       /** Gradient tolerance for the truncated preconditioned conjugate gradient
   * solver: stop if ||g|| < kappa * ||g_0||.  This parameter should be in the
   * range (0,1). */
      double STPCG_kappa = 0.7;

      /** Gradient tolerance based upon a fractional-power reduction in the norm of
   * the gradient: stop if ||g|| < ||kappa||^{1+ theta}.  This value should be
   * positive, and controls the asymptotic convergence rate of the
   * truncated-Newton trust-region solver: specifically, for theta > 0, the TNT
   * algorithm converges q-superlinearly with order (1+theta). */
      double STPCG_theta = .5;

      /** Maximum permitted number of (outer) iterations of the RTR algorithm */
      unsigned int max_RTR_iterations = 2000;
      /** Maximum number of inner (truncated conjugate-gradient) iterations to
      * perform per out iteration */
      unsigned int max_tCG_iterations = 5;


      // verbose = {0, 1}
      unsigned int estimation_verbose = 0;
      
      unsigned int verbose = 0;  // for TNT


      /// DEFAULT CONSTRUCTOR with default values

       HomographyEstimationOptions() {};
       

};  // end of struct: EssentialEstimationOptions


    
/** This struct contains the output of the Essential Estimation */
struct HomographyEstimationResult {

  
  // Primal objective value
  double f_hat;

  // The norm of the Riemannian gradient at x_hat
  double gradnorm;

   // Elapsed time for the initialisation
   double elapsed_init_time;
   // Elapsed time for C
   double elapsed_C_time;
   // Elapsed time for the optimization on manifold
   double elapsed_iterative_time;
   // Elapsed time for the whole estimation: initialisation + optimization + verification
   double elapsed_estimation_time;

   
   // Output rotation matrix
   Matrix3 H_opt;

   // Output rotation matrix
   Matrix2 R_opt;
   // Output translation vector
   Vector3 t_opt;
   // Output translation vector
   Vector3 n_opt;

    bool status_riem;  // condition termination
    /* Default constructor */
   HomographyEstimationResult() {};

}; // end of EssentialEstimationResult struct




class RotPriorHClass
{
       
 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
    /* Default constructor */
    RotPriorHClass(void){};

                    
    RotPriorHClass( const  Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                   const  Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                   const  Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations,
                    const HomographyEstimationOptions & options = HomographyEstimationOptions()): p0_(bearingVectors0), 
                    p1_(bearingVectors1), w_(weights_observations), options_(options) 
                    {};        
                           

    ~RotPriorHClass(void){};

    // void estimateEssentialMatrix(void);

    HomographyEstimationResult getResults(Matrix34 & Rt);

    void printResult(HomographyEstimationResult & results);

    private:

        HomographyEstimationOptions options_;
        Eigen::Matrix<double, 3, Eigen::Dynamic> p0_, p1_; 
        Eigen::Matrix<double, 1, Eigen::Dynamic> w_;

};  //end of RotPrior class




}  // end of RotPrior namespace

