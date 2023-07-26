#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "RotPriorHTypes.h"
#include "RotPriorH.h"
#include "RotPriorHUtils.h"


// for the noisy Y-rotation
#include "../utils/generatePointCloud.h"


using namespace std;
using namespace Eigen;
using namespace RotPriorH; 



int main(int argc, char** argv)
{
    
    std::cout << "Essential Matrix Estimation with rotation prior!\n";


    double N_iter = 1;
    //set experiment parameters
    double noise = .5;
    const size_t n_points = 50;
    double FoV = 100;  // in degrees
    double min_depth = 1.;     // in meters
    double max_depth = 8.;       // in meters
    double outlier_fraction = 0;
    double focal_length = 512; 

    std::srand(std::time(nullptr));


    RotPriorH::Vector3 translation;
    RotPriorH::Matrix3 rotation;

    Eigen::MatrixXd points_3D(3, n_points);

    // define struct with params
    UtilsRelPose::StrInputGenProblem str_in = UtilsRelPose::StrInputGenProblem(); 
    str_in.FoV = FoV;  
    str_in.min_depth = min_depth; 
    str_in.max_depth = max_depth; 
    str_in.focal_length = focal_length; 
    str_in.n_points = n_points; 
    str_in.noise = noise; 
    str_in.parallax = 0.0; 
    // str_in.use_custom_t = options.use_custom_t; 
    // str_in.max_rotation = options.max_rotation;                                      
    // str_in.custom_t = options.custom_t; 
                                           
    // generate problem
    UtilsRelPose::StrOutputGenProblem str_out = UtilsRelPose::createSyntheticExperiment(str_in, 
                                               UtilsRelPose::generateRandomTranslationDefault, 
                                               UtilsRelPose::generatePerturbRotationY); 
    // extract data
    translation = str_out.translation; 
    rotation = str_out.rotation; 

    Eigen::Matrix<double, 3, n_points> p0,p1;
    Eigen::Matrix<double, 1, n_points> weights = Eigen::Matrix<double, 1, n_points>::Ones();  
                          
                              
    // normalize observations        
    for (int i=0; i < n_points; i++)
    {
            Vector3 obs1 = str_out.points_correspondences[i].bearing_vector_0; 
            p0.col(i) = obs1.normalized(); 

            Vector3 obs2 = str_out.points_correspondences[i].bearing_vector_1; 
            p1.col(i) = obs2.normalized(); 

    }                                  
                   

    // GT homography
    RotPriorH::Matrix3 Hgt = rotation;
                   
    std::vector<Matrix2> rot_dlt_s; 
    std::vector<Vector3> trans_dlt_s, n_dlt_s; 
    std::vector<Matrix3> H_dlt_s; 
    computeInitPattern(p0, p1, weights, &H_dlt_s, &rot_dlt_s, &trans_dlt_s, &n_dlt_s);   


    double best_cost = 100;
    Matrix2 R_opt; 
    Vector3 t_opt; 
    Vector3 n_opt; 
    Matrix3 H_opt;
          
      
    for (size_t ii=0; ii<H_dlt_s.size(); ii++)
    {

    Eigen::Matrix2d R22 = rot_dlt_s[ii]; 
    RotPriorH::Vector3 t22 = trans_dlt_s[ii], n2 = n_dlt_s[ii];      


    std::cout << "Rotation from init:\n" << R22 << std::endl; 
    std::cout << "Translation from init:\n" << t22.transpose() << std::endl; 
    std::cout << "Normal from init:\n" << n2.transpose() << std::endl;         

    std::cout << "Ground truth rotation:\n" << rotation << std::endl;


    RotPriorH::Matrix34 Rtinit; 
    Rtinit.setZero();
    Rtinit.block<2,2>(0,0) = R22;
    Rtinit.block<3,1>(0,2) = t22;
    Rtinit.block<3,1>(0,3) = n2.normalized();

        
    HomographyEstimationOptions options;
    options.verbose = 1;
    options.estimation_verbose=1;
    options.stepsize_tol = 1e-07;
               

    /* Run essential with prior estimation and DLT + pattern*/       
    RotPriorHClass h_prior_est_pattern(p0,p1,weights, options); 
    // run the actual estimation
    HomographyEstimationResult my_result_prior_pattern = h_prior_est_pattern.getResults(Rtinit); 
     


    std::cout << "Ground truth rotation:\n" << rotation << std::endl;
    std::cout << "Ground truth translation Tgt:\n" << translation.transpose() << std::endl;
     

    std::cout << "+++++++++++++++\n\nResults from H with prior estimation\n"; 
    // ess_prior_est_pattern.printResult(my_result_prior_pattern);
    std::cout << "Rotation:\n" << my_result_prior_pattern.R_opt << std::endl; 
    std::cout << "Traslation:\n" << my_result_prior_pattern.t_opt.transpose() << std::endl;  


    if (my_result_prior_pattern.f_hat < best_cost)
    {
      std::cout << "Update best solution\n"; 
      best_cost = my_result_prior_pattern.f_hat;                     
      // save bets solution                            
      R_opt     = my_result_prior_pattern.R_opt; 
      t_opt     = my_result_prior_pattern.t_opt; 
      n_opt     = my_result_prior_pattern.n_opt; 
      H_opt     = my_result_prior_pattern.H_opt;   
    }  

}
       

    return 0;

}
