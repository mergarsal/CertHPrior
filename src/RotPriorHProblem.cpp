#include "RotPriorHProblem.h"


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


namespace RotPriorH{
        RotPriorHProblem::RotPriorHProblem(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                           const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                           const Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations) 
        {

                p0_ = bearingVectors0; 
                p1_ = bearingVectors1; 
                w_ = weights_observations; 
                
                data_matrix_C_ = constructDataMatrix(bearingVectors0, bearingVectors1, weights_observations);
                number_points_ = weights_observations.cols();

                /*
                N_.setZero();
                T_.setZero();
                QN_.setZero();
                QT_.setZero();
                Qr_.setZero();
                */


        }; //end of constructor
        RotPriorHProblem::RotPriorHProblem(const Matrix9& data_C)
        {
                data_matrix_C_ = data_C;
                /*
                N_.setZero();
                T_.setZero();
                QN_.setZero();
                QT_.setZero();
                Qr_.setZero();
                */
        }; //end of constructor

         
        Eigen::Matrix<double, 1, Eigen::Dynamic> RotPriorHProblem::computeResiduals(const Matrix34 & Rt) const
        {
            Eigen::Matrix<double, 1, Eigen::Dynamic> residual;
            residual.resize(number_points_);

            Matrix3 H = computeHfromRt(convertRt2ExpandRt(Rt, 1.));
            Matrix3 T1;
            
            for(size_t i = 0; i < number_points_; i++)
            {
                    // std::cout << "Reading point with id " << i << std::endl;
                    Vector3 v1 = p1_.col(i);
                    Vector3 v0 = p0_.col(i);
                    T1.setZero();
                    T1(0,1) = -v1(2);
                    T1(0,2) = v1(1); 
                    T1(1,0) = v1(2); 
                    T1(1,2) = -v1(0); 
                    T1(2,0)=-v1(1);
                    T1(2,1)=v1(0);
                    residual(i) = (T1 * H * v0).norm();
                    // std::cout << "Residual computed " << residual(i) << std::endl;
            }
    
            return residual;
        }        
        
        
        // update points with weights
        void RotPriorHProblem::updateWeights(const Eigen::Matrix<double, 1, Eigen::Dynamic> & new_weights)
        {

            assert(new_weights.cols() == number_points_);
            for(size_t i = 0; i < number_points_; i++)   w_[i] = new_weights(i);


            // update data matrix
            data_matrix_C_ = constructDataMatrix(p0_, p1_, w_);
        }

        
        // Compute pseudo Jacobi preconditioner
        Matrix3 RotPriorHProblem::computePseudoJacobiPrecon(void)
        {
                // If you use the 8 points algorithm, please
                // select the use_precon option

                Eigen::JacobiSVD<Matrix9> svd(data_matrix_C_, Eigen::ComputeFullU | Eigen::ComputeFullV);  //faster

                Vector9 svd_diag = svd.singularValues();
                double scale_c = svd_diag(0) + svd_diag(1) + svd_diag(2);
                return (Matrix3::Identity() / scale_c);
        }

        

      
        
        // apply precon to the full X
        Matrix34 RotPriorHProblem::precondition(const Matrix34 & X, const Matrix34 & Xdot) const
        {
            return tangent_space_projection(X,  Matrix_precon_ * Xdot);
        }
       
        
       

        double RotPriorHProblem::evaluate_objective(const Matrix34 &Y) const {
            // save this for latter (gradient & hessian)
            Matrix3 H = computeHfromRtPrior(Y);            
                        

            return (0.5 * vec(H).transpose() * data_matrix_C_ * vec(H));

        }
       
         double RotPriorHProblem::evaluate_objective(const Matrix34 &Y, ProblemCachedMatrices & problem_matrices) const {

            Vector3 t,n; 
            Matrix3 R;
            Matrix2 rot; 
            n = Y.block<3,1>(0,3); 
            t = Y.block<3,1>(0,2);
            rot=Y.block<2,2>(0,0); 
            // expand matrix
            R.setZero();
            R(0,0)=rot(0,0); 
            R(0,2)=rot(0,1); 
            R(2,0)=rot(1,0);
            R(2,2)=rot(1,1);
            R(1,1) = 1; 
            
            Matrix3 H = computeHfromRtPrior(Y); 
            // std::cout << "H from cost:\n" << H << std::endl; 
            // std::cout << "H from cost normalized:\n" << H / H.norm() << std::endl; 
            // std::cout << "Point Y in cost:\n" << Y << std::endl;
                    
            // problem_matrices.N_ = takeKronProduct(n, Matrix3::Identity());
            Matrix93 temp;
            // std::cout << "n in cost:\n" << n << std::endl;  
            temp.setZero(); 
            temp.block<3,3>(0,0) = n(0) * Matrix3::Identity();
            temp.block<3,3>(3,0) = n(1) * Matrix3::Identity();
            temp.block<3,3>(6,0) = n(2) * Matrix3::Identity();
            problem_matrices.N_ = temp;
            
            // problem_matrices.T_ = takeKronProduct(Matrix3::Identity(), t);
            temp.setZero(); 
            // std::cout << "t in cost:\n" << t << std::endl; 
            temp.block<3,1>(0,0) = t;
            temp.block<3,1>(3,1) = t;
            temp.block<3,1>(6,2) = t;
            problem_matrices.T_ = temp;
            
            problem_matrices.QN_ = data_matrix_C_ * problem_matrices.N_; 
            problem_matrices.QT_ = data_matrix_C_ * problem_matrices.T_; 
            problem_matrices.Qr_ = data_matrix_C_ * vec(R); 
            // std::cout << "Matrix C from cost:\n" << data_matrix_C_ << std::endl;
            // std::cout << "Matrix C from cost normalized:\n" << data_matrix_C_ / data_matrix_C_.trace() << std::endl;
            
            return (0.5 * vec(H).transpose() * data_matrix_C_ * vec(H));
        }



             Matrix34 RotPriorHProblem::Euclidean_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const
            {
                   Vector3 t,n; 
                   Matrix3 R;
                   Matrix2 rot; 
                   n = Y.block<3,1>(0,3); 
                    t = Y.block<3,1>(0,2);
                    rot=Y.block<2,2>(0,0); 
                    // expand matrix
                    R.setZero();
                    R(0,0)=rot(0,0); 
                    R(0,2)=rot(0,1); 
                    R(2,0)=rot(1,0);
                    R(2,2)=rot(1,1);
                    R(1,1) = 1; 
                   
           
                Matrix34 G;
                G.setZero();   
                // std::cout << "Matrix Qr:\n" << problem_matrices.Qr_ << std::endl; 
                // std::cout << "Matrix QT:\n" << problem_matrices.QT_ << std::endl; 
                
                Vector9 mr = problem_matrices.Qr_ + problem_matrices.QT_ * n;  
                Matrix3 Gr = Eigen::Map<Matrix3> (mr.data(), 3, 3);               
                
                G(0,0)=Gr(0,0);
                G(0,1)=Gr(0,2);
                G(1,0)=Gr(2,0);
                G(1,1)=Gr(2,2);
                
                // std::cout << "Gradient for rotation:\n" << Gr << std::endl;
                // translation
                G.block<3,1>(0,2) = problem_matrices.N_.transpose() * problem_matrices.QN_ * t + problem_matrices.N_.transpose() * problem_matrices.Qr_;
                
                // normal
                G.block<3,1>(0,3) = problem_matrices.T_.transpose() * problem_matrices.QT_ * n + problem_matrices.T_.transpose() * problem_matrices.Qr_;
                
                // std::cout << "gradient:\n" << G << std::endl; 
                return G;
            }

             Matrix34 RotPriorHProblem::Riemannian_gradient(const Matrix34 &Y, const Matrix34 &nablaF_Y) const
             {
              return tangent_space_projection(Y, nablaF_Y);
             }
             

            Matrix34 RotPriorHProblem::Riemannian_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const {
            
            Matrix34 Gg = Euclidean_gradient(Y, problem_matrices); 
            
              return tangent_space_projection(Y, Gg);
            }


       /** Given a matrix Y in the domain D of the SE-Sync optimization problem, and
           * a tangent vector dotY in T_D(Y), the tangent space of the domain of the
           * optimization problem at Y, this function computes and returns Hess
           * F(Y)[dotY], the action of the Riemannian Hessian on dotX */                                               

     Matrix34 RotPriorHProblem::Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                                   const ProblemCachedMatrices & problem_matrices,
                                                   const Matrix34 &dotY) const
                                                   {
                                                   // Euclidean Hessian-vector product
                                                    Matrix34 HessRiemannian;
                                                    HessRiemannian.setZero();

                                                 
                                                    
                                                   Vector3 t,n; 
                                                   Matrix3 R;
                                                   Matrix2 rot; 
                                                   n = Y.block<3,1>(0,3); 
                                                   t = Y.block<3,1>(0,2);
                                                   rot=Y.block<2,2>(0,0); 
                                                   // expand matrix
                                                   R.setZero();
                                                   R(0,0)=rot(0,0); 
                                                   R(0,2)=rot(0,1); 
                                                   R(2,0)=rot(1,0);
                                                   R(2,2)=rot(1,1);
                                                   R(1,1) = 1; 
                                                    
                                                                                                        
                                                    Matrix35 dotYext = convertRt2ExpandRt(dotY, 0.0);
                                                    
                                                    Vector3 Vn = dotYext.block<3, 1>(0, 4);
                                                    Vector3 Vt = dotYext.block<3, 1>(0, 3);
                                                    Matrix3 VR = dotYext.block<3, 3>(0, 0);
                                                    // std::cout << "dotYExt:\n" << dotYext << std::endl; 
                                                    // std::cout << "dotRot:\n" << VR << std::endl;
                                                    
                                                    
                                                    // std::cout << "MIxed term for hessian:\n" << mixed_terms_der << std::endl;
                                                    // Compute Euclidean Hessian
                                                    Matrix9By15 coeff_R;
                                                    coeff_R.setZero(); 
                                                    coeff_R.block<9, 9>(0, 0) = data_matrix_C_;
                                                    coeff_R.block<9, 3>(0, 9) = problem_matrices.QN_;
                                                    coeff_R.block<9, 3>(0, 12) = problem_matrices.QT_;


                                                    Vector15 Vx;
                                                    Vx.setZero();
                                                    Vx.block<9, 1>(0, 0) = vec(VR);
                                                    Vx.block<3, 1>(9, 0) = Vt;
                                                    Vx.block<3, 1>(12, 0) = Vn;


                                                    Vector9 mr = coeff_R * Vx;
                                                    Matrix3 HessR = Eigen::Map<Matrix3> (mr.data(), 3, 3);
                                                    
                                                    
                                                    Eigen::Matrix<double, 1, 9> row_qr;
                                                    row_qr = problem_matrices.Qr_.transpose();
                                                    
                                                    // Matrix3 RQ = Eigen::Map<Matrix3> (row_qr.data(), 3, 3);    
                                                    Matrix3 Cs = Eigen::Map<Matrix3> (row_qr.data(), 3, 3); 
                                                    /*
                                                    Vector3 temp1;
                                                    temp1 << 1, 0, 0;
                                                    Cs.row(0) = temp1.transpose() * (RQ);   
                                                    temp1 << 0, 1, 0;  
                                                    Cs.row(1) = temp1.transpose() * (RQ); 
                                                    temp1 << 0, 0, 1;
                                                    Cs.row(2) = temp1.transpose() * (RQ); 
                                                    */
                                                              
                                                    Matrix3By15 coeff_T;
                                                    coeff_T.setZero(); 
                                                    coeff_T.block<3, 9>(0, 0) = problem_matrices.QN_.transpose();
                                                    coeff_T.block<3, 3>(0, 9) = problem_matrices.N_.transpose() * problem_matrices.QN_;
                                                    coeff_T.block<3, 3>(0, 12) = 2 * problem_matrices.N_.transpose() * problem_matrices.QT_ + Cs;

                                                    Vector3 HessT = coeff_T * Vx;
                                                    
                                                    
                                                    Matrix3By15 coeff_N;
                                                    coeff_N.setZero(); 
                                                    coeff_N.block<3, 9>(0, 0) = problem_matrices.QT_.transpose();
                                                    coeff_N.block<3, 3>(0, 9) = coeff_T.block<3, 3>(0, 12).transpose();
                                                    coeff_N.block<3, 3>(0, 12) = problem_matrices.T_.transpose() * problem_matrices.QT_;

                                                    Vector3 HessN = coeff_N * Vx;
                                                    

                                                    // recover NablaF(Y)
                                                    
                                                    Matrix35 nablaExt = convertRt2ExpandRt(problem_matrices.NablaF_Y, 0.0);
                                                                                                                                                            
                                                    Vector3 nabla_Yn = (nablaExt).block<3,1>(0,4);
                                                    Vector3 nabla_Yt = (nablaExt).block<3,1>(0,3);
                                                    Matrix3 nabla_YR = (nablaExt).block<3,3>(0,0);
                                                    
                                                    // std::cout << "Saved nabla (should be 3x5):\n" << nablaExt << std::endl;

                                                    // clean
                                                    HessRiemannian.setZero();

                                                    // Riemannain Hessian for t (sphere)
                                                    HessRiemannian.block<3,1>(0, 3) = domain_.ProjSphere(n, HessN - (n.dot(nabla_Yn)) * Vn);
                                                    HessRiemannian.block<3,1>(0, 2) = HessT;

                                                    // Riemannain Hessian for R (rotation)
                                                    // Riemannain Hessian for R (rotation)
                                                    // ROTATION 
                                                    /**
                                                    Matrix3 SR = R * VR;  // tangent to ambient fcn
                                                    
                                                    Matrix3 HessRR_mult = ( R.transpose() * HessR - SymProduct(SR, R, nabla_YR) );
                                                    Matrix3 HessRR = 0.5 * (HessRR_mult - HessRR_mult.transpose());
                                                    
                                                    Matrix2 HessRred = Matrix2::Zero(); 
                                                    HessRred(0,0)=HessRR(0,0); 
                                                    HessRred(0,1)=HessRR(0,2);
                                                    HessRred(1,0)=HessRR(2,0);
                                                    HessRred(1,1)=HessRR(2,2);
                                                    
                                                    // HessRiemannian.block<2, 2>(0, 0) = domain_.ProjRotation(rot, HessRred);
                                                    HessRiemannian.block<2, 2>(0, 0) = HessRred;
                                                    // Riemannian hessian
                                                    // std::cout << "Hessian for R:\n" << HessRR << std::endl;
                                                    // std::cout << "Hessian:\n" << HessRiemannian << std::endl;
                                                    return HessRiemannian;
                                                    **/
                                                    
                                                     
                                                    // Riemannain Hessian for R (rotation)
                                                    Matrix3 HessRR = HessR - SymProduct(VR, R, nabla_YR);
                                                    
                                                    Matrix2 HessRred = Matrix2::Zero(); 
                                                    HessRred(0,0)=HessRR(0,0); 
                                                    HessRred(0,1)=HessRR(0,2);
                                                    HessRred(1,0)=HessRR(2,0);
                                                    HessRred(1,1)=HessRR(2,2);
                                                    HessRiemannian.block<2, 2>(0, 0) = domain_.ProjRotation(rot, HessRred);
                                                    return HessRiemannian;
                           }


     

          Matrix34 RotPriorHProblem::tangent_space_projection(const Matrix34 &Y,
                                                            const Matrix34 &dotY) const 
                                                            { return domain_.Proj(Y, dotY); }


            Matrix34 RotPriorHProblem::retract(const Matrix34 &Y, const Matrix34 &dotY) const
            {
                return domain_.retract(Y, dotY);

            }

      

} // end of Essential namespace
