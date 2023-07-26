#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "RotPriorHTypes.h"
#include "RotPriorH.h"
#include "RotPriorHUtils.h"
#include "RotPriorHConstraints.h"


// point cloud 
#include "../utils/generatePointCloudPlanar.h"

#include "IterCertAlg/SymmCert.h"
#include "IterCertAlg/SymmCert.cpp"
#include "IterCertAlg/RankYCert.h"
#include "IterCertAlg/RankYCert.cpp"

using namespace std;
using namespace Eigen;
using namespace RotPriorH; 

typedef Eigen::Matrix<double, 11, 1> Vector11;
typedef Eigen::Matrix<double, 21, 1> Vector21;
typedef Eigen::Matrix<double, 69, 1> Vector69;
typedef Eigen::Matrix<double, 97, 1> Vector97;

typedef Eigen::Matrix<double, 8, 8> Matrix8;
typedef Eigen::Matrix<double, 16, 16> Matrix16;

int main(int argc, char** argv)
{
    
    std::cout << "Homography (plane) Estimation with rotation prior!\n";



    double noise = 5.0;        
    double max_parallax = 1.5;  // in meters
    double width_plane = 3;              // this is half the plane so: [-width, +width]
    double height_plane = 3; 
    double focal_length = 512; 
    size_t size_img = 1024; 
    double d_plane = 4;

    const int n_points = 20;

    std::srand(std::time(nullptr));



    RotPriorH::Vector3 translation;
    Eigen::Matrix3d rotation;

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
                                               
    // generate problem
    TwoViewPlanar::PCPlanarRes str_out = TwoViewPlanar::generatePointCloudPlanar(str_in, 
                        TwoViewPlanar::generateRandomTranslationDefault, 
                        TwoViewPlanar::generatePerturbRotationY);
    // extract data
    translation = str_out.translation.normalized(); 
    rotation = str_out.rotation; 

    Eigen::Matrix3d Tx; 
    Tx << 0, -translation(2), translation(1), translation(2), 0, -translation(0), -translation(1), translation(0), 0; 
    Eigen::Matrix3d E = Tx * rotation;

    Vector3 n; 
    n << 0, 0, 1; 
              

    Eigen::Matrix<double, 3, n_points> p0,p1;
    Eigen::Matrix<double, 1, n_points> weights = Eigen::Matrix<double, 1, n_points>::Ones(); 
                          
    // normalize observations 
    Eigen::Matrix3d K = Matrix3::Zero(); 
    K(0,0) = focal_length; 
    K(1,1) = focal_length; 
    K(0,2) = size_img / 2; 
    K(1,2) = size_img / 2; 
    K(2,2) = 1; 

    Eigen::Matrix3d iK; 
    iK = K.inverse(); 

    for (int i=0; i < n_points; i++)
    {
            Vector3 obs1 = iK * str_out.obs1.col(i); 
            p0.col(i) = obs1.normalized(); 

            Vector3 obs2 = iK * str_out.obs2.col(i); 
            p1.col(i) = obs2.normalized(); 

    }                                  


    Eigen::Matrix2d rot_init_pattern; 
    RotPriorH::Vector3 trans_init_pattern, n_init_pattern; 
    Eigen::Matrix3d H_init_prior; 
      
    // std::cout << "Rinit:\n" << Rtinit << std::endl; 

    // std::cout << "R*t:\n" << rotation.transpose() * translation; 
    // std::cout << "R^t * n:\n" << rotation.transpose() * n << std::endl; 
    Eigen::Matrix3d R1, R2; 
    RotPriorH::Vector3 t1, t2; 
    R1 = str_out.R1; 
    R2 = str_out.R2; 
    t1 = str_out.t1; 
    t2 = str_out.t2;

    // GT homography
    Eigen::Matrix3d Hgt; 
    computeHomography(n, d_plane, R1, R2, t1, t2, Hgt); 


    RotPriorH::Vector3 Trel = t2 - rotation * t1; 
            
    RotPriorH::Vector3 mm = R1 * n; 
    double d11 = d_plane + mm.dot(t1);

    RotPriorH::Vector3 TT = Trel / d11;         
                             
      
      
    std::vector<Matrix2> rot_dlt_s; 
    std::vector<Vector3> trans_dlt_s, n_dlt_s; 
    std::vector<Matrix3> H_dlt_s; 
    computeInitPattern(p0, p1, weights, &H_dlt_s, &rot_dlt_s, &trans_dlt_s, &n_dlt_s);   
      
  
    double best_cost = 100;
    Matrix2 R_opt_iter; 
    Vector3 t_opt_iter; 
    Vector3 n_opt_iter; 
    Matrix3 H_opt_iter;
          
      
    for (size_t ii=0; ii<H_dlt_s.size(); ii++)
    {

    Eigen::Matrix2d R22 = rot_dlt_s[ii]; 
    RotPriorH::Vector3 t22 = trans_dlt_s[ii], n2 = n_dlt_s[ii]; 

    RotPriorH::Matrix34 Rtinit; 
    Rtinit.setZero();

    Rtinit.block<2,2>(0,0) = R22;
    Rtinit.block<3,1>(0,2) = t22;
    Rtinit.block<3,1>(0,3) = n2.normalized();

    /*
    Rtinit.block<2,2>(0,0) << rotation(0,0), rotation(0, 2), rotation(2, 0), rotation(2,2);
    Rtinit.block<3,1>(0,2) = translation * TT.norm();
    Rtinit.block<3,1>(0,3) = R1 * n;
    */
        
    HomographyEstimationOptions options;
    options.verbose = 1;
    options.estimation_verbose=0;
    options.stepsize_tol = 1e-07;
               

    /* Run essential with prior estimation and DLT + pattern*/       
    RotPriorHClass h_prior_est_pattern(p0, p1, weights, options); 
    // run the actual estimation
    HomographyEstimationResult my_result_prior_pattern = h_prior_est_pattern.getResults(Rtinit); 
     
   
    std::cout << "Ground truth rotation:\n" << rotation << std::endl;
    std::cout << "Ground truth translation Tgt:\n" << translation << std::endl;
    std::cout << "Ground truth normal:\n" << R1 * n << std::endl; 

    std::cout << "+++++++++++++++\n\nResults from E with prior estimation PATTERN\n"; 
    // h_prior_est_pattern.printResult(my_result_prior_pattern);
    std::cout << "Rotation:\n" << my_result_prior_pattern.R_opt << std::endl; 
    std::cout << "Traslation:\n" << my_result_prior_pattern.t_opt << std::endl;   
    std::cout << "Normal:\n" << my_result_prior_pattern.n_opt << std::endl;  
    // std::cout << "Rtinit:\n" << Rtinit << std::endl;
    // std::cout << "Errores for mani:\n"; 
    // std::cout << "Error rotation: " << distR(rotation, my_result_prior_pattern.R_opt) << std::endl;
    // std::cout << "Error translation: " << distT(translation, my_result_prior_pattern.t_opt) << std::endl;
      
      
    if (my_result_prior_pattern.f_hat < best_cost)
    {
    std::cout << "Update best solution\n"; 
    best_cost = my_result_prior_pattern.f_hat;                     
    // save bets solution                            
    R_opt_iter     = my_result_prior_pattern.R_opt; 
    t_opt_iter     = my_result_prior_pattern.t_opt; 
    n_opt_iter     = my_result_prior_pattern.n_opt; 
    H_opt_iter     = my_result_prior_pattern.H_opt;   
    }  
    }
    
    // onyl check the solution with lowest cost as we may have more than one init
    
          
    /* Checking optimality */
    std::cout << "########################\nFORMULATION HR\n";
    RotPriorH::Vector12 x_hr; 
    Eigen::Matrix3d H_opt = H_opt_iter; 
    Eigen::Matrix2d R_opt = R_opt_iter; 
    double norm_H = H_opt.norm(); 
    x_hr << vec(H_opt) / norm_H, R_opt(0,0) / norm_H, R_opt(0,1) / norm_H, 1 / norm_H; 
    Vector11 mult_hr_init = Vector11::Zero();

    std::vector<Matrix12> Ahr; 
    RotPriorH::createHRConstraints(Ahr);      

    Eigen::MatrixXd Jhr(12, Ahr.size()); 

    std::cout << "Checking constraints\n"; 
    for (int ii=0; ii < Ahr.size(); ii++)
    {
    // std::cout << "Constraint #"<< ii << " with value: " << x_hr.dot(Ahr[ii] * x_hr) << std::endl;      
    Jhr.col(ii) = Ahr[ii] * x_hr;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_jac(Jhr);       
    std::cout << "Singular values jacobian HR: " << svd_jac.singularValues().transpose() << std::endl;   // rank 6
                                
   
    //compute Ch

    RotPriorH::Matrix9 C = constructDataMatrix(p0,p1,weights);

    RotPriorH::Matrix12 Chr = Matrix12::Zero();
    Chr.block<9,9>(0,0) = C;  
    SymmCert::SymmCertOptions options_symm;
    options_symm.verbose = 1; 
    options_symm.estimation_verbose = 1; 
    options_symm.max_RTR_iterations = 10;

    SymmCert::SymmCertClass<Matrix12, Vector12, Vector11> cert_hr_symm(options_symm);


    SymmCert::SymmCertResult<Matrix12, Vector12, Vector11> res_hr = cert_hr_symm.getResults(x_hr, Chr, Ahr, mult_hr_init, Chr);

           

    // cert_hr_symm.printResult(res_hr);      

    Eigen::SelfAdjointEigenSolver<Matrix12> eig_hr(res_hr.opt_Hessian);
    std::cout << "Eigenvalues Hessian HR:\n" << eig_hr.eigenvalues() << std::endl; 

    const int r_hr = 8;
    Eigen::Matrix<double, 12, r_hr> Yhr; 
    RotPriorH::Matrix12 U = eig_hr.eigenvectors();  
    Eigen::Matrix<double, 12, r_hr> Ured = U.block<12,r_hr>(0,4); 
    Matrix8 Dhr = (eig_hr.eigenvalues().block<r_hr,1>(4,0)).cwiseAbs().cwiseSqrt().asDiagonal();
    Yhr = Ured * Dhr;


    RankYCert::RankYCertOptions options_ranky;
    //  options_ranky.verbose = 1; 
    options_ranky.estimation_verbose = 1; 
    options_ranky.max_RTR_iterations = 10;

     
    RankYCert::RankYCertClass<Matrix12, Eigen::Matrix<double, 12, r_hr>, Vector12, Vector11> cert_hr_ranky(options_ranky);


    RankYCert::RankYCertResult<Matrix12, Eigen::Matrix<double, 12, r_hr>, Vector11> res_hr_ranky = cert_hr_ranky.getResults(x_hr, 
                                                                                                      Chr, Ahr, res_hr.opt_mult, Yhr);



    // Eigen::SelfAdjointEigenSolver<Matrix12> eig_hr_ranky(res_hr_ranky.opt_Hessian);
    // std::cout << "Eigenvalues Hessian HR Rank Y:\n" << eig_hr_ranky.eigenvalues() << std::endl; 
    // cert_hr_ranky.printResult(res_hr_ranky); 


    /** Second formulation: R, t, n **/

    std::cout << "########################\nFORMULATION Rtn\n";
    Vector12 x_rtn; 
    RotPriorH::Vector3 t_opt = t_opt_iter; 
    RotPriorH::Vector3 n_opt = n_opt_iter; 
    Eigen::Matrix3d Q = t_opt * n_opt.transpose();
    x_rtn << vec(Q) / norm_H, R_opt(0,0) / norm_H, R_opt(0,1) / norm_H, 1 / norm_H;
    Vector11 mult_rtn_init = Vector11::Zero();
       
      
    std::vector<Matrix12> Artn; 
    RotPriorH::createRtnConstraints(Artn);     

    /*
    std::cout << "Checking constraints\n"; 
    for (int ii=0; ii < Artn.size(); ii++)
    {
    std::cout << "Constraint #"<< ii << " with value: " << x_rtn.dot(Artn[ii] * x_rtn) << std::endl;

    }      

    */

    Eigen::MatrixXd Jrtn(12, Artn.size()); 

    std::cout << "Checking constraints\n"; 
    for (int ii=0; ii < Artn.size(); ii++)
    {
    // std::cout << "Constraint #"<< ii << " with value: " << x_hr.dot(Ahr[ii] * x_hr) << std::endl;
    Jrtn.col(ii) = Artn[ii] * x_rtn;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_jac2(Jrtn); 

    std::cout << "Singular values jacobian RTN: " << svd_jac2.singularValues().transpose() << std::endl; 
      
      
                               

    //compute Crtn
    Eigen::Matrix3d Qr; 
    RotPriorH::Matrix39 Qtr; 
    createDataCReduced(C, Qr, Qtr);

    RotPriorH::Matrix12 Crtn = Matrix12::Zero();
    Crtn.block<9,9>(0,0) = C; 
    Crtn.block<9,3>(0,9) = 2 * Qtr.transpose(); 
    Crtn.block<3,3>(9,9) = Qr; 
    RotPriorH::Matrix12 Cc = 0.5 * (Crtn + Crtn.transpose()); 
    Cc /= Cc.trace();     

    SymmCert::SymmCertResult<Matrix12, Vector12, Vector11> res_rtn = cert_hr_symm.getResults(x_rtn, Cc, Artn, mult_rtn_init, Cc);


    // cert_hr_symm.printResult(res_hr); 
      

    Eigen::SelfAdjointEigenSolver<Matrix12> eig_rtn(res_rtn.opt_Hessian);
    // std::cout << "Eigenvalues Hessian HR Formulation #2:\n" << eig_rtn.eigenvalues() << std::endl; 

    const int r_rtn = 8;
    Eigen::Matrix<double, 12, r_rtn> Yrtn; 
    RotPriorH::Matrix12 Urtn = eig_rtn.eigenvectors();  
    Eigen::Matrix<double, 12, r_rtn> Ured_rtn = Urtn.block<12,r_rtn>(0,4); 
    Matrix8 Dhr_rtn = (eig_rtn.eigenvalues().block<r_rtn,1>(4,0)).cwiseAbs().cwiseSqrt().asDiagonal();
    Yrtn = Ured_rtn * Dhr_rtn;


    RankYCert::RankYCertResult<Matrix12, Eigen::Matrix<double, 12, r_hr>, Vector11> res_rtn_ranky = cert_hr_ranky.getResults(x_rtn, 
                                                                                                        Cc, Artn, res_rtn.opt_mult, Yrtn);

      
      
      
      
    /** Run redundant certifier **/

    std::cout << "########################\nFORMULATION REDUNDANT\n";
    Vector21 x_red; 

    norm_H = H_opt.norm();        
    x_red << vec(H_opt) / norm_H, R_opt(0,0) / norm_H, R_opt(0,1) / norm_H, 1 / norm_H, vec(Q) / norm_H; 


    std::vector<Matrix21> Ared; 
    RotPriorH::createHRtnConstraints(Ared);  
    const size_t n_const = Ared.size();       
    Vector69 mult_red_init; 
    mult_red_init.setZero(); 

    std::cout << "Number constraints for redundant: " << n_const << std::endl;
     
 
       
    Eigen::MatrixXd Jred(21, Ared.size()); 

    std::cout << "Checking constraints\n"; 
    for (int ii=0; ii < Ared.size(); ii++)
    {
    std::cout << "Constraint #"<< ii << " with value: " << x_red.dot(Ared[ii] * x_red) << std::endl;

    Jred.col(ii) = Ared[ii] * x_red;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_jac3(Jred); 

    std::cout << "Size Jacobian: " << Jred.rows() << " x " << Jred.cols() << std::endl; 
    std::cout << "Singular values jacobian RRED: " << svd_jac3.singularValues().transpose() << std::endl; 

      
    //compute Crtn
    Matrix21 Q2 = Matrix21::Zero(); 
    Q2.block<9,9>(0,0) = C; 
    Q2.block<3,3>(9,9) = Qr; 
    Q2.block<9,9>(12,12)= C; 
    Q2.block<3,9>(9,12) = 2 * Qtr; 
    Matrix21 Qq = 0.5 * (Q2 + Q2.transpose()); 
    Qq /= Qq.trace();     
      

    SymmCert::SymmCertClass<Matrix21, Vector21, Vector69> cert_red_symm(options_symm);

    SymmCert::SymmCertResult<Matrix21, Vector21, Vector69> res_red = cert_red_symm.getResults(x_red, Qq, Ared, mult_red_init, Qq);
      
      
    Eigen::SelfAdjointEigenSolver<Matrix21> eig_red(res_red.opt_Hessian);
    std::cout << "Eigenvalues Hessian HR Formulation #RED:\n" << eig_red.eigenvalues() << std::endl; 

    const int r_red = 16;
    Eigen::Matrix<double, 21, r_red> Yred; 
    Matrix21 Ured_r = eig_red.eigenvectors();  
    Eigen::Matrix<double, 21, r_red> Ured_red = Ured_r.block<21,r_red>(0,5); 
    Matrix16 Dhr_red = (eig_red.eigenvalues().block<r_red,1>(5,0)).cwiseAbs().cwiseSqrt().asDiagonal();
    Yred = Ured_red * Dhr_red;
    // std::cout << "norm eigenvalues:\n" << eig_red.eigenvalues() / eig_red.eigenvalues()(20) << std::endl;

    RankYCert::RankYCertClass<Matrix21, Eigen::Matrix<double, 21, r_red>, Vector21, Vector69> cert_red_ranky(options_ranky);
    RankYCert::RankYCertResult<Matrix21, Eigen::Matrix<double, 21, r_red>, Vector69> res_red_ranky = cert_red_ranky.getResults(x_red, 
                                                                                                         Qq, Ared, res_red.opt_mult, Yred);

    // check independence constraints
    // cert_red_ranky.checkConstIndp(Ared); 
      

    /** Run redundant certifier with simpler cost **/

    std::cout << "########################\nFORMULATION REDUNDANT WITH RANK DEFICIENT\n";
      
    //compute Crtn
    Matrix21 Q2_simpl = Matrix21::Zero(); 
    Q2_simpl.block<9,9>(0,0) = C; 
    Matrix21 Qq_simpl = 0.5 * (Q2_simpl + Q2_simpl.transpose()); 
    Qq_simpl /= Qq_simpl.trace();     

    // std::cout << "Cost rank def: " << x_red.dot(Qq_simpl * x_red) << std::endl;
    SymmCert::SymmCertClass<Matrix21, Vector21, Vector69> cert_red_symm_simpl(options_symm);

    SymmCert::SymmCertResult<Matrix21, Vector21, Vector69> res_red_simpl = cert_red_symm_simpl.getResults(x_red, Qq_simpl, Ared, mult_red_init, Qq_simpl);

            

    Eigen::SelfAdjointEigenSolver<Matrix21> eig_red_simpl(res_red_simpl.opt_Hessian);
    // std::cout << "Eigenvalues Hessian HR Formulation #RED SIMPL:\n" << eig_red_simpl.eigenvalues() << std::endl; 

    const int r_red_simpl = 8;
    Eigen::Matrix<double, 21, r_red_simpl> Yred_simpl; 
    Matrix21 Ured_r_simpl = eig_red_simpl.eigenvectors();  
    Eigen::Matrix<double, 21, r_red_simpl> Ured_red_simpl = Ured_r_simpl.block<21,r_red_simpl>(0,13); 
    Matrix8 Dhr_red_simpl = (eig_red_simpl.eigenvalues().block<r_red_simpl,1>(13,0)).cwiseAbs().cwiseSqrt().asDiagonal();
    Yred_simpl = Ured_red_simpl * Dhr_red_simpl;

    std::cout << "norm eigenvalues:\n" << eig_red_simpl.eigenvalues() / eig_red_simpl.eigenvalues()(20) << std::endl;

    // std::cout << "Dhr red:\n" << Dhr_red_simpl << std::endl;
    options_ranky.estimation_verbose = 1; 
    options_ranky.verbose = 1;
    RankYCert::RankYCertClass<Matrix21, Eigen::Matrix<double, 21, r_red_simpl>, Vector21, Vector69> cert_red_ranky_simpl(options_ranky);
    RankYCert::RankYCertResult<Matrix21, Eigen::Matrix<double, 21, r_red_simpl>, Vector69> res_red_ranky_simpl = cert_red_ranky_simpl.getResults(x_red, 
                                                                                                        Qq_simpl, Ared, res_red_simpl.opt_mult, Yred_simpl);

      

    /** Run redundant certifier with redundant (97) constraints **/
       
    std::cout << "########################\nFORMULATION REDUNDANT 97 CONSTRAINTS\n";    

    std::vector<Matrix21> Ared_l; 
    RotPriorH::createHRtnConstraintsRed(Ared_l);  
    const size_t n_const_l = Ared_l.size();       
    Vector97 mult_red_init_l; 
    mult_red_init_l.setZero(); 

    std::cout << "Number constraints for redundant (should be 97): " << n_const_l << std::endl;
    
    
           

    SymmCert::SymmCertClass<Matrix21, Vector21, Vector97> cert_red_symm_l(options_symm);

    SymmCert::SymmCertResult<Matrix21, Vector21, Vector97> res_red_l = cert_red_symm_l.getResults(x_red, Qq, Ared_l, mult_red_init_l, Qq);


    Eigen::SelfAdjointEigenSolver<Matrix21> eig_red_l(res_red_l.opt_Hessian);
    std::cout << "Eigenvalues Hessian HR Formulation #RED:\n" << eig_red_l.eigenvalues() << std::endl; 

    const int r_red_l = 16;
    Eigen::Matrix<double, 21, r_red_l> Yred_l; 
    Matrix21 Ured_r_l = eig_red_l.eigenvectors();  
    Eigen::Matrix<double, 21, r_red_l> Ured_red_l = Ured_r_l.block<21,r_red_l>(0,5); 
    Matrix16 Dhr_red_l = (eig_red_l.eigenvalues().block<r_red_l,1>(5,0)).cwiseAbs().cwiseSqrt().asDiagonal();
    Yred_l = Ured_red_l * Dhr_red_l;
      
    // std::cout << "norm eigenvalues:\n" << eig_red_l.eigenvalues() / eig_red_l.eigenvalues()(20) << std::endl;

    options_ranky.estimation_verbose = 1; 
    options_ranky.verbose = 1;
    RankYCert::RankYCertClass<Matrix21, Eigen::Matrix<double, 21, r_red_l>, Vector21, Vector97> cert_red_ranky_l(options_ranky);

    // check independence constraints
    // cert_red_ranky_l.checkConstIndp(Ared_l); 

    Eigen::MatrixXd Jred2(21, Ared_l.size()); 

    std::cout << "Checking constraints\n"; 
    for (int ii=0; ii < Ared_l.size(); ii++)
    {
    std::cout << "Constraint #"<< ii << " with value: " << x_red.dot(Ared_l[ii] * x_red) << std::endl;

    Jred2.col(ii) = Ared_l[ii] * x_red;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_jac4(Jred2); 

    std::cout << "Size Jacobian: " << Jred2.rows() << " x " << Jred2.cols() << std::endl; 
    std::cout << "Singular values jacobian RRED 97: " << svd_jac4.singularValues().transpose() << std::endl; 



    RankYCert::RankYCertResult<Matrix21, Eigen::Matrix<double, 21, r_red_l>, Vector97> res_red_ranky_l = cert_red_ranky_l.getResults(x_red, 
                                                                                                        Qq, Ared_l, res_red_l.opt_mult, Yred_l);

    // check independence constraints
    // cert_red_ranky_l.checkConstIndp(Ared_l); 

    /** Run redundant certifier with redundant (97) constraints: low rank C **/

    std::cout << "########################\nFORMULATION REDUNDANT 97 CONSTRAINTS LOW RANK C\n";    

    mult_red_init_l.setZero(); 


    SymmCert::SymmCertClass<Matrix21, Vector21, Vector97> cert_red_symm_l_simpl(options_symm);

    SymmCert::SymmCertResult<Matrix21, Vector21, Vector97> res_red_l_simpl = cert_red_symm_l_simpl.getResults(x_red, Qq_simpl, Ared_l, mult_red_init_l, Qq_simpl);


    Eigen::SelfAdjointEigenSolver<Matrix21> eig_red_l_simpl(res_red_l_simpl.opt_Hessian);
    // std::cout << "Eigenvalues Hessian HR Formulation #RED:\n" << eig_red_l_simpl.eigenvalues() << std::endl; 

    const int r_red_l_simpl = 8;
    Eigen::Matrix<double, 21, r_red_l_simpl> Yred_l_simpl; 
    Matrix21 Ured_r_l_simpl = eig_red_l_simpl.eigenvectors();  
    Eigen::Matrix<double, 21, r_red_l_simpl> Ured_red_l_simpl = Ured_r_l_simpl.block<21,r_red_l_simpl>(0,13); 
    Matrix8 Dhr_red_l_simpl = (eig_red_l_simpl.eigenvalues().block<r_red_l_simpl,1>(13,0)).cwiseAbs().cwiseSqrt().asDiagonal();
    Yred_l_simpl = Ured_red_l_simpl * Dhr_red_l_simpl;

    // std::cout << "norm eigenvalues:\n" << eig_red_l_simpl.eigenvalues() / eig_red_l_simpl.eigenvalues()(20) << std::endl;
    options_ranky.estimation_verbose = 1; 
    options_ranky.verbose = 1; 
    RankYCert::RankYCertClass<Matrix21, Eigen::Matrix<double, 21, r_red_l_simpl>, Vector21, Vector97> cert_red_ranky_l_simpl(options_ranky);
    RankYCert::RankYCertResult<Matrix21, Eigen::Matrix<double, 21, r_red_l_simpl>, Vector97> res_red_ranky_l_simpl = cert_red_ranky_l_simpl.getResults(x_red, 
                                                                                                        Qq_simpl, Ared_l, res_red_l_simpl.opt_mult, Yred_l_simpl);
      
      
     
 
  return 0;

}
