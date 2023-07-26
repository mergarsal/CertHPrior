#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "RotPriorHTypes.h"
#include "RotPriorH.h"
#include "RotPriorHUtils.h"


// point cloud 
#include "../utils/generatePointCloudPlanar.h"

using namespace std;
using namespace Eigen;
using namespace RotPriorH; 



RotPriorH::Vector3 returnTranslation(const double max_parallax, 
                                     const RotPriorH::Vector3 & dir_parallax)
{
        return (max_parallax * dir_parallax);

}
        
int main(int argc, char** argv)
{
    
    std::cout << "[EXAMPLE] Homography Estimation with rotation prior!\n";



    double noise = 0.5;        
    double max_parallax = 2.0;  // in meters
    double width_plane = 3;              // this is half the plane so: [-width, +width]
    double height_plane = 3; 
    double focal_length = 512; 
    size_t size_img = 1024; 
    double d_plane = 3;

    const int n_points = 10;

    std::srand(std::time(nullptr));



    RotPriorH::Vector3 translation;
    RotPriorH::Matrix3 rotation;

    Eigen::MatrixXd points_3D(3, n_points);

    // define struct with params
    TwoViewPlanar::PCPlanarParams str_in = TwoViewPlanar::PCPlanarParams();   
    str_in.focal_length = focal_length; 
    str_in.size_img = size_img;
    str_in.noise = noise; 
    str_in.width_plane = width_plane; 
    str_in.height_plane = height_plane;
    str_in.max_parallax = max_parallax; 
    str_in.d_plane = d_plane;
    str_in.N_points = n_points;
    bool use_custom_t = false; 

    Vector3 custom_t; 
    custom_t << 2, 0, 0;
    str_in.dir_parallax = custom_t; 
    double noise_Ry = 0.0;
    str_in.dir_rotation << noise_Ry, 0, 0;
       
                                               
    // generate problem
    TwoViewPlanar::PCPlanarRes str_out = TwoViewPlanar::PCPlanarRes(n_points); 
                                           
    if (use_custom_t)
    {
            str_out = TwoViewPlanar::generatePointCloudPlanar(str_in, 
                            returnTranslation, TwoViewPlanar::generatePerturbRotationY);
    }
    else
    {
            str_out = TwoViewPlanar::generatePointCloudPlanar(str_in, 
                    TwoViewPlanar::generateRandomTranslationDefault,        
                    TwoViewPlanar::generatePerturbRotationY);
    }


    // extract data
    translation = str_out.translation.normalized(); 
    rotation = str_out.rotation; 

    RotPriorH::Matrix3 Tx; 
    Tx << 0, -translation(2), translation(1), translation(2), 0, -translation(0), -translation(1), translation(0), 0; 
    RotPriorH::Matrix3 E = Tx * rotation;

    RotPriorH::Vector3 n; 
    n << 0, 0, 1; 

       
    Eigen::Matrix<double, 3, n_points> p0,p1;
    Eigen::Matrix<double, 1, n_points> weights = Eigen::Matrix<double, 1, n_points>::Ones(); 
                          
    // normalize observations 
    Matrix3 K = Matrix3::Zero(); 
    K(0,0) = focal_length; 
    K(1,1) = focal_length; 
    K(0,2) = size_img / 2; 
    K(1,2) = size_img / 2; 
    K(2,2) = 1; 

    RotPriorH::Matrix3 iK; 
    iK = K.inverse(); 

    for (int i=0; i < n_points; i++)
    {
        Vector3 obs1 = iK * str_out.obs1.col(i); 
        p0.col(i) = obs1.normalized(); 

        Vector3 obs2 = iK * str_out.obs2.col(i); 
        p1.col(i) = obs2.normalized(); 
    }                                  
   


    RotPriorH::Matrix3 R1, R2; 
    RotPriorH::Vector3 t1, t2; 
    R1 = str_out.R1; 
    R2 = str_out.R2; 
    t1 = str_out.t1; 
    t2 = str_out.t2;

    // GT homography
    RotPriorH::Matrix3 Hgt; 
    computeHomography(n, d_plane, R1, R2, t1, t2, Hgt); 
           
           
           
    // init solution from DLT       
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


        RotPriorH::Matrix34 Rtinit; 
        Rtinit.setZero();
        Rtinit.block<2,2>(0,0) = R22;
        Rtinit.block<3,1>(0,2) = t22;
        Rtinit.block<3,1>(0,3) = n2.normalized();

            
        HomographyEstimationOptions options;
        options.verbose = 0;
        options.estimation_verbose=0;
        options.stepsize_tol = 1e-07;
                   

        /* Run essential with prior estimation and DLT + pattern*/       
        RotPriorHClass h_prior_est_pattern(p0, p1, weights, options); 
        // run the actual estimation
        HomographyEstimationResult my_result_prior_pattern = h_prior_est_pattern.getResults(Rtinit); 
         
   
        std::cout << "[EXAMPLE] Ground truth rotation:\n" << rotation << std::endl;
        std::cout << "[EXAMPLE] Ground truth translation Tgt:\n" << translation << std::endl;
        std::cout << "[EXAMPLE] Ground truth normal:\n" << R1 * n << std::endl;
        // std::cout << "Rotation from init:\n" << R22 << std::endl; 
        // std::cout << "Translation from init:\n" << t22 << std::endl; 
        // std::cout << "Normal from init:\n" << n2 << std::endl;  

        std::cout << "[EXAMPLE] \n+++++++++++++++\n\nResults from E with prior estimation PATTERN\n"; 
        // h_prior_est_pattern.printResult(my_result_prior_pattern);
        std::cout << "[EXAMPLE] Rotation:\n" << my_result_prior_pattern.R_opt << std::endl; 
        std::cout << "[EXAMPLE] Traslation:\n" << my_result_prior_pattern.t_opt << std::endl;   
        std::cout << "[EXAMPLE] Normal:\n" << my_result_prior_pattern.n_opt << std::endl;  
        std::cout << "[EXAMPLE] Init from DLT:\n" << Rtinit << std::endl; 

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
