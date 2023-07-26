#include "RotPriorHUtils.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


#include <Eigen/Eigenvalues> 


#define SQRT2 1.41421356237
#define R180PI 57.2957795131

namespace RotPriorH
{

        
        Matrix3 computeHfromRtPrior(const Matrix2 & rot, const Vector3 trans, const Vector3 n)
        {
                Matrix34 rt; 
                rt.block<2, 2>(0, 0) = rot; 
                rt.block<3, 1>(0, 2) = trans;
                rt.block<3, 1>(0, 3) = n;
                return (computeHfromRtPrior(rt));
        }

        Matrix3 computeHfromRtPrior(const Matrix34 & rt)
        {
                return (computeHfromRt(convertRt2ExpandRt(rt, 1.0)));
        }
        
        
        Matrix3 computeHfromRt(const Matrix3 & R, const Vector3 & t, const Vector3 & n)
        {
            Matrix3 H = R + t*n.transpose();
            
            return H; 
        }
        

        Matrix3 computeHfromRt(const Matrix35& Rt)
        {
            Vector3 t = Rt.block<3, 1>(0, 3);
            Matrix3 R = Rt.block<3, 3>(0, 0);
            Vector3 n = Rt.block<3, 1>(0, 4); 
            return (computeHfromRt(R, t, n));
        }


        Matrix35 convertRt2ExpandRt(const Matrix34 & Rt, double pad_int)
        {
                Matrix35 Rt_exp; 
                   
                
                Rt_exp.setZero(); 
                Rt_exp(0,0)=Rt(0,0); 
                Rt_exp(0,2)=Rt(0,1);
                Rt_exp(2,0)=Rt(1,0);
                Rt_exp(2,2)=Rt(1,1);                 
                Rt_exp(1, 1) = pad_int;
                Rt_exp.block<3, 1>(0, 3) = Rt.block<3, 1>(0, 2);
                Rt_exp.block<3, 1>(0, 4) = Rt.block<3, 1>(0, 3);
                
                return Rt_exp;
        }

       
        Matrix2 projectRtoPrior(const Matrix3 & R)
        {
                Matrix3 Rt = R; 

                if (Rt.determinant() < 0)
                    Rt.col(2) *=-1;
                Eigen::Quaternion<double> eq = Eigen::Quaternion<double>(Rt); 
                
                      
                Vector4 init_rot = Vector4(); 
                init_rot << eq.w(), -eq.x(), -eq.y(), -eq.z();
                

                // let be a Z-rotation
                init_rot(1) = 0; 
                init_rot(3) = 0; 
                init_rot.normalize();  
                
                // NOTE: use conjugate 
                Eigen::Quaternion<double> rot_quat = Eigen::Quaternion<double>(init_rot[0], -init_rot[1], -init_rot[2], -init_rot[3]); 
                Matrix3 rot = rot_quat.toRotationMatrix(); 

                
                Matrix2 rot_r = Matrix2::Zero(); 
                rot_r(0,0)=rot(0,0);
                rot_r(0,1)=rot(0,2); 
                rot_r(1,0)=rot(2,0);
                rot_r(1,1)=rot(2,2);
                               
                
                // convert back to rotation matrix                          
                return (rot_r); 
        }
        

        Matrix2 obtainRot(const Matrix2 & R)
            {
                Matrix2 R_p;
                Eigen::JacobiSVD<Matrix2> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

                R_p =  svd.matrixU() * svd.matrixV().transpose();
                
                if (R_p.determinant() < 0)                                   
                    R_p.col(2) *= -1;
                
            
                return R_p;
            
            }

          Matrix3 SymProduct(const Matrix3 & Ydot, const Matrix3 & Y, const Matrix3 & nabla_Y) //const
         {
                Matrix3 P = Y.transpose() * nabla_Y;
                
                return (Ydot * 0.5 * (P + P.transpose()));
         }
         
        
         void createDataCReduced(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 Matrix9& C, Matrix3& Qr, Matrix39 & Qtr)
         {
         
               C = constructDataMatrix(bearingVectors0, bearingVectors1, weights_observations);
                                
               // auxiliar matrix        
               
               Qr.setZero(); 
               Qtr.setZero(); 
               
               Qr(0,0) = C(0,0) + C(0,8) + C(8,0) + C(8,8); 
               Qr(0,1) = C(0,6) - C(0,2) + - C(2,8) + C(6, 8); 
               Qr(0,2) = C(0,4) + C(4,8); 
               
               Qr(1,0) = C(6,0) - C(2,0) - C(8,2) + C(8,6); 
               Qr(1,1) = C(2,2) - C(2,6) - C(6,2) + C(6,6); 
               Qr(1,2) = C(4,6) - C(4,2);
               
               Qr(2,0) = C(8,4) + C(4,0); 
               Qr(2,1) = C(6,4) - C(2,4); 
               Qr(2,2) = C(4,4); 
                              
               
               Qtr.row(0) = C.row(0) + C.row(8); 
               Qtr.row(1) = C.row(6) - C.row(2); 
               Qtr.row(2) = C.row(4); 
                              
               return;   
         }
         
         
        Vector9 vec(Matrix3 & M)
        {    return Eigen::Map<Vector9> (M.data(), 9, 1);       }

       
         
       
        
        Eigen::MatrixXd skew(const Eigen::MatrixXd & M)
        {
            return (0.5 * (M - M.transpose()));
        }
        
        Eigen::Matrix3d skew(const Eigen::Vector3d & m)
        {
            Eigen::Matrix3d s = Eigen::Matrix3d::Zero(); 
            s(0,1) = -m(2); 
            s(0,2) = m(1); 
            s(1,0) = m(2); 
            s(1,2) = -m(0); 
            s(2,0) = -m(1); 
            s(2,1) = m(0);
            
            return s; 
        }        
        
         Matrix9 constructDataMatrix(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                     const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                     const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations)
        {

                Matrix9 C;
                Matrix9 temp;
                Vector3 vtemp;
                // clean output matrix
                C.setZero();
                
              
                for (int i = 0; i < bearingVectors0.cols(); i++)
                {
                        // clean data
                        temp.setZero();
                        const Vector3 v1 = bearingVectors1.col(i);
                        const Vector3 v0 = bearingVectors0.col(i);
                        const double weight = weights_observations[i];

                        vtemp.setZero(); 
                        
                        // B1 * v1
                        vtemp(1) = v1(2);
                        vtemp(2) = -v1(1);                     
                        for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v0[j] * vtemp;

                        C += weight * temp * temp.transpose();
                        
                        /* Second */
                        vtemp.setZero(); 
                        // B2 * v1
                        vtemp(0) = -v1(2);
                        vtemp(2) = v1(0);  
                                             
                        temp.setZero();                
                        for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v0[j] * vtemp;

                        C += weight * temp * temp.transpose();
                        
                        /* Third*/
                        vtemp.setZero(); 
                        // B3 * v1
                        vtemp(0) = -v1(1);
                        vtemp(1) = v1(0); 
                        
                                             
                        temp.setZero();                
                        for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v0[j] * vtemp;
                        

                }
                
                Matrix9 Q = 0.5 * (C + C.transpose()); 

             
                return (Q);
        }
        
          void computeHomography(Vector3 & n, double d_plane, 
                               Matrix3 & R1, Matrix3 & R2, 
                               Vector3 & t1, Vector3 & t2, 
                               Matrix3 & H)
                               
                               
        {
                Matrix3 Rrel = R2 * R1.transpose(); 
                Vector3 Trel = t2 - Rrel * t1; 
                
                Vector3 m = R1 * n; 
                double d1 = d_plane + m.dot(t1);
                
                H.setZero(); 
                H = Rrel + Trel * m.transpose() / d1; 
        
        }
        
         void createDataCReduced(const Matrix9& C, Matrix3& Qr, Matrix39 & Qtr)
         {
                
               // auxiliar matrix        
               
               Qr.setZero(); 
               Qtr.setZero(); 
               
               Qr(0,0) = C(0,0) + C(0,8) + C(8,0) + C(8,8); 
               Qr(0,1) = C(0,6) - C(0,2) + - C(2,8) + C(6, 8); 
               Qr(0,2) = C(0,4) + C(4,8); 
               
               Qr(1,0) = C(6,0) - C(2,0) - C(8,2) + C(8,6); 
               Qr(1,1) = C(2,2) - C(2,6) - C(6,2) + C(6,6); 
               Qr(1,2) = C(4,6) - C(4,2);
               
               Qr(2,0) = C(8,4) + C(4,0); 
               Qr(2,1) = C(6,4) - C(2,4); 
               Qr(2,2) = C(4,4); 
                              
               Qtr.row(0) = C.row(0) + C.row(8); 
               Qtr.row(1) = C.row(6) - C.row(2); 
               Qtr.row(2) = C.row(4); 
               

               
               return;      
         
         }
         
                 /* From colmap:
                 src/base/homography_matrix
                 */
                double ComputeOppositeOfMinor(const Eigen::Matrix3d& matrix, const size_t row,
                                              const size_t col) {
                  const size_t col1 = col == 0 ? 1 : 0;
                  const size_t col2 = col == 2 ? 1 : 2;
                  const size_t row1 = row == 0 ? 1 : 0;
                  const size_t row2 = row == 2 ? 1 : 2;
                  return (matrix(row1, col2) * matrix(row2, col1) -
                          matrix(row1, col1) * matrix(row2, col2));
                }
                /* From colmap:
                 src/base/homography_matrix
                 */
               Matrix3 ComputeHomographyRotation(const Matrix3& H_normalized,
                                          const Vector3& tstar,
                                          const Vector3& n,
                                          const double v) {
                          return H_normalized *
                                 (Matrix3::Identity() - (2.0 / v) * tstar * n.transpose());
                        }

                 /* From colmap:
                 src/util/math
                 */
                template <typename T>
                int SignOfNumber(const T val) {
                  return (T(0) < val) - (val < T(0));
                }
                           
            
         void DecomposeHomographyMatrix3D(const Matrix3& H,
                                        std::vector<Matrix3>* R,
                                        std::vector<Vector3>* t,
                                        std::vector<Vector3>* n) 
         {
                                        
                  Matrix3 H_normalized = H;
                  // Remove scale from normalized homography.
                  Eigen::JacobiSVD<Eigen::Matrix3d> hmatrix_norm_svd(H_normalized);
                  H_normalized.array() /= hmatrix_norm_svd.singularValues()[1];
                  

                  const Matrix3 S =
                      H_normalized.transpose() * H_normalized;

                  // Check if H is rotation matrix.
                  const double kMinInfinityNorm = 1e-2; 
                  if ((S - Eigen::Matrix3d::Identity()).lpNorm<Eigen::Infinity>() < kMinInfinityNorm) {
            
                    double sR = SignOfNumber(H_normalized(1,1));    
                    *R = {sR * H_normalized};
                    *t = {Vector3::Zero()};
                    Vector3 ss = Vector3::Zero(); 
                    ss(0) = 1; 
                    *n = {ss};
                    return;
                  }         
                  
                    
                 
                 Eigen::JacobiSVD<Eigen::Matrix3d> svd_hh(S, Eigen::ComputeFullV);
                 
                 Matrix3 V = svd_hh.matrixV();
                 
                 if (V.determinant() < 0) V *= -1;
                                
                 
                 // svd
                 double s3 = svd_hh.singularValues()(2); 
                 double s1 = svd_hh.singularValues()(0);
                 
                 // svd
                 Eigen::Vector3d v1 = V.col(0), v2 = V.col(1), v3 = V.col(2);
                 
             
                 Eigen::Vector3d u1 = 1 / std::sqrt(s1-s3) * (std::sqrt(1-s3)*v1 + std::sqrt(s1-1)*v3);
                 Eigen::Vector3d u2 = 1 / std::sqrt(s1-s3) * (std::sqrt(1-s3)*v1 - std::sqrt(s1-1)*v3);
                                   
               
                 
                 Eigen::Vector3d Hv2, Hu1, Hu2; 
                 Hv2 = H_normalized*v2; 
                 Hu1 = H_normalized*u1; 
                 Hu2 = H_normalized*u2;
                 Eigen::Matrix3d skewHv2 = skew(Hv2), skewv2 = skew(v2); 
                 
                 Eigen::Matrix3d U1; 
                 U1 << v2, u1, skewv2*u1;         
                        
                 Eigen::Matrix3d W1; 
                 W1 << Hv2, Hu1, skewHv2*Hu1;
                 
                 Eigen::Matrix3d U2; 
                 U2 << v2, u2, skewv2*u2; 
                 Eigen::Matrix3d W2;
                 W2 << Hv2, Hu2, skewHv2*Hu2;
                 
                 Eigen::Matrix3d R1 = W1 * U1.transpose(); 
                 Eigen::Matrix3d R2 = W2 * U2.transpose();
                 
                 // make sure R1(1,1) > 0 (same for R2)
                 double sR1 = SignOfNumber(R1(1,1));                 
                 double sR2 = SignOfNumber(R2(1,1));               
                 
                 
                 Eigen::Vector3d n_t1 = skewv2 * u1;  
                 // if (n_t1(2) < 0)        n_t1 *= -1;
                 Eigen::Vector3d t_t1 = sR1 * (H_normalized - R1) * n_t1;
                 
                 
                 Eigen::Vector3d n_t2 = skewv2 * u2; 
                 // if (n_t2(2)< 0)         n_t2 *= -1;
                 Eigen::Vector3d t_t2 = sR2 *(H_normalized - R2) * n_t2;
                 
                  
                
                  
                 // BORRAR LO DE ARRIBA
                 // double nh = std::sqrt((H.transpose()*H).trace()); 
                 
                 /*
                 Matrix3 Q1 =  H_normalized - (sR1 * R1) ; 
                 
                 Eigen::JacobiSVD<Eigen::Matrix3d> svd_Q1(Q1, Eigen::ComputeFullV | Eigen::ComputeFullU);
                 
                 
                 
                 t_t1 = svd_Q1.singularValues()(0) * nh * svd_Q1.matrixU().col(0); 
                 n_t1 = svd_Q1.matrixV().col(0);
                 
                 Matrix3 Q2 = H_normalized - (sR2 * R2).transpose() ; 
                 Eigen::JacobiSVD<Eigen::Matrix3d> svd_Q2(Q2, Eigen::ComputeFullV | Eigen::ComputeFullU);
                 
                
                 t_t2 = svd_Q2.singularValues()(0) * nh * svd_Q2.matrixU().col(0); 
                 n_t2 = svd_Q2.matrixV().col(0);
                 
                 std::cout << "n_t2 after: " << n_t2.transpose() << std::endl; 
                 
                 */
                
                
                 *R = {sR1 * R1, sR2 * R2};
                 *t = {t_t1, t_t2};
                 *n = {n_t1, n_t2};                               
                                        
         }  
              
                         
         void DecomposeHomographyMatrix3DPM(const Matrix3& H,
                                        std::vector<Matrix3>* R,
                                        std::vector<Vector3>* t,
                                        std::vector<Vector3>* n) 
                                        
         {
                  std::vector<Eigen::Matrix3d> Rs, Rt;
                  std::vector<Eigen::Vector3d> ts, tt;
                  std::vector<Eigen::Vector3d> ns, nt;

                  /* Call for + H */
                  DecomposeHomographyMatrix3D(H, &Rt, &tt, &nt);         
                  Rs = Rt; 
                  ts = tt; 
                  ns = nt;
                  
                  /* Call for -H */
                  DecomposeHomographyMatrix3D(-H, &Rt, &tt, &nt);   
                  int n_rt = Rt.size(); 
                  for (int i=0;i<n_rt;i++)
                  {
                    Rs.push_back(Rt[i]);
                    ts.push_back(tt[i]); 
                    ns.push_back(nt[i]);
                  }      
                  
                  *R = Rs; 
                  *t = ts; 
                  *n = ns;
         
         }    
              
              
         /* From colmap:
         src/base/homography_matrix
         */
         void DecomposeHomographyMatrix(const Matrix3& H,
                                        std::vector<Matrix3>* R,
                                        std::vector<Vector3>* t,
                                        std::vector<Vector3>* n) {
                  
                  Matrix3 H_normalized = H;
                  // Remove scale from normalized homography.
                  Eigen::JacobiSVD<Eigen::Matrix3d> hmatrix_norm_svd(H_normalized);
                  H_normalized.array() /= hmatrix_norm_svd.singularValues()[1];

                  const Matrix3 S =
                      H_normalized.transpose() * H_normalized - Eigen::Matrix3d::Identity();

                  // Check if H is rotation matrix.
                  const double kMinInfinityNorm = 1e-2;
                  if (S.lpNorm<Eigen::Infinity>() < kMinInfinityNorm) {
            
                    *R = {H_normalized};
                    *t = {Vector3::Zero()};
                    Vector3 ss = Vector3::Zero(); 
                    ss(0) = 1; 
                    *n = {ss};
                    return;
                  }

                  const double M00 = ComputeOppositeOfMinor(S, 0, 0);
                  const double M11 = ComputeOppositeOfMinor(S, 1, 1);
                  const double M22 = ComputeOppositeOfMinor(S, 2, 2);

                  const double rtM00 = std::sqrt(M00);
                  const double rtM11 = std::sqrt(M11);
                  const double rtM22 = std::sqrt(M22);

                  const double M01 = ComputeOppositeOfMinor(S, 0, 1);
                  const double M12 = ComputeOppositeOfMinor(S, 1, 2);
                  const double M02 = ComputeOppositeOfMinor(S, 0, 2);

                  const int e12 = SignOfNumber(M12);
                  const int e02 = SignOfNumber(M02);
                  const int e01 = SignOfNumber(M01);

                  const double nS00 = std::abs(S(0, 0));
                  const double nS11 = std::abs(S(1, 1));
                  const double nS22 = std::abs(S(2, 2));

                  const std::array<double, 3> nS{{nS00, nS11, nS22}};
                  const size_t idx =
                      std::distance(nS.begin(), std::max_element(nS.begin(), nS.end()));

                  Vector3 np1;
                  Vector3 np2;
                  if (idx == 0) {
                    np1[0] = S(0, 0);
                    np2[0] = S(0, 0);
                    np1[1] = S(0, 1) + rtM22;
                    np2[1] = S(0, 1) - rtM22;
                    np1[2] = S(0, 2) + e12 * rtM11;
                    np2[2] = S(0, 2) - e12 * rtM11;
                  } else if (idx == 1) {
                    np1[0] = S(0, 1) + rtM22;
                    np2[0] = S(0, 1) - rtM22;
                    np1[1] = S(1, 1);
                    np2[1] = S(1, 1);
                    np1[2] = S(1, 2) - e02 * rtM00;
                    np2[2] = S(1, 2) + e02 * rtM00;
                  } else if (idx == 2) {
                    np1[0] = S(0, 2) + e01 * rtM11;
                    np2[0] = S(0, 2) - e01 * rtM11;
                    np1[1] = S(1, 2) + rtM00;
                    np2[1] = S(1, 2) - rtM00;
                    np1[2] = S(2, 2);
                    np2[2] = S(2, 2);
                  }

                  const double traceS = S.trace();
                  const double v = 2.0 * std::sqrt(1.0 + traceS - M00 - M11 - M22);

                  const double ESii = SignOfNumber(S(idx, idx));
                  const double r_2 = 2 + traceS + v;
                  const double nt_2 = 2 + traceS - v;

                  const double r = std::sqrt(r_2);
                  const double n_t = std::sqrt(nt_2);

                  const Vector3 n1 = np1.normalized();
                  const Vector3 n2 = np2.normalized();

                  const double half_nt = 0.5 * n_t;
                  const double esii_t_r = ESii * r;

                  const Vector3 t1_star = half_nt * (esii_t_r * n2 - n_t * n1);
                  const Vector3 t2_star = half_nt * (esii_t_r * n1 - n_t * n2);

                  const Matrix3 R1 =
                      ComputeHomographyRotation(H_normalized, t1_star, n1, v);
                  const Vector3 t1 = R1 * t1_star;

                  const Matrix3 R2 =
                      ComputeHomographyRotation(H_normalized, t2_star, n2, v);
                  const Vector3 t2 = R2 * t2_star;

                  *R = {R1, R1, R2, R2};
                  *t = {t1, -t1, t2, -t2};
                  *n = {-n1, n1, -n2, n2};
                }
                
                /* 
                Check whether the rotation is a Y-rotation
                Select the closest to a Y-rotation
                 */
                void PoseFromHomographyMatrixYRot(const Matrix3& H,
                                              std::vector<Eigen::Matrix3d>* R, 
                                              std::vector<Eigen::Vector3d>* t,
                                              std::vector<Eigen::Vector3d>* n) {
                                              
                  std::vector<Eigen::Matrix3d> R_cmbs;
                  std::vector<Eigen::Vector3d> t_cmbs;
                  std::vector<Eigen::Vector3d> n_cmbs;

                  DecomposeHomographyMatrix3DPM(H, &R_cmbs, &t_cmbs, &n_cmbs);    
                                
                 
                  *R = R_cmbs;
                  *t = t_cmbs;
                  *n = n_cmbs;
                                           
                  
               }
                              
               
                 /* From colmap:
                 src/base/homography_matrix
                 */
                void PoseFromHomographyMatrix(const Matrix3& H,
                                              Matrix3* R, 
                                              Vector3* t,
                                              Vector3* n) {
                  std::vector<Eigen::Matrix3d> R_cmbs;
                  std::vector<Eigen::Vector3d> t_cmbs;
                  std::vector<Eigen::Vector3d> n_cmbs;
                  DecomposeHomographyMatrix(H, &R_cmbs, &t_cmbs, &n_cmbs);

                  // check R1 and R2
                  Matrix3 Rtemp = R_cmbs[0]; 
                  // get entries for Y rotation
                  double c = Rtemp(0,0), s = Rtemp(0,2); 
                  double err_y = std::fabs(c*c + s*s - 1); 
                  
                  Rtemp = R_cmbs[2]; 
                  c = Rtemp(0,0);
                  s = Rtemp(0,2); 
                  // NOTE: Maybe we should add also the entry (1,1) = 1.0
                  if (err_y < std::fabs(c*c + s*s - 1))
                        {
                                *R = R_cmbs[0];
                                *t = t_cmbs[0];
                                *n = n_cmbs[0];
                        }
                  else
                  {
                                *R = R_cmbs[2];
                                *t = t_cmbs[2];
                                *n = n_cmbs[2];  
                  }                
               }
               

         /** Use `An invitation to 3D vision` algorithm **/
         /* initialization based on DLT + pattern */     
         void computeInitPattern(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 std::vector<Matrix3>* H, std::vector<Matrix2> * rot, std::vector<Vector3>* trans, std::vector<Vector3> * n)
         {
         
               Matrix9 Q; 
               Matrix3 Qr; 
               Matrix39 Qtr;
               createDataCReduced(bearingVectors0, bearingVectors1, weights_observations, Q, Qr, Qtr);
               
               
               
             Eigen::SelfAdjointEigenSolver<Matrix9> eig_qh(Q);
             Vector9 q_eig = eig_qh.eigenvalues(); 
             // std::cout << "Eigenvalues: " << q_eig.transpose() << std::endl; 
            
            
             Vector9 qq = eig_qh.eigenvectors().col(0); 
               
             Matrix3 Hh = Eigen::Map<Matrix3> (qq.data(), 3, 3);  
             
             
             std::vector<Matrix3> Rv = {}; 
             std::vector<Vector3> transv = {}, nv = {};
             
             
             PoseFromHomographyMatrixYRot(Hh, &Rv, &transv, &nv);    
             
             /* select here best */        
             size_t best_id = 0; 
             Matrix3 Hi; 
             Eigen::Matrix<double, 9, 1> hhi;
             
             double cost_best = 1000;
             
             std::vector<Eigen::Matrix2d> rot_s = {};
             std::vector<Eigen::Matrix3d> hom_s = {};
             std::vector<Eigen::Vector3d> trans_s = {};
             std::vector<Eigen::Vector3d> normal_s = {};



             for (size_t ii = 0; ii < Rv.size(); ii++)
             {                
                Eigen::Matrix3d Rog = Rv[ii];  // not nece Y rot
                Vector3 tii = transv[ii];
                Vector3 nii = nv[ii];       
                
                // Rog is Y-rot
                Hi = Rog + tii*nii.transpose();  
                
                hhi = Eigen::Map<Eigen::Matrix<double, 9, 1>> (Hi.data(), 9, 1);  
             
                double cost_i = hhi.transpose() * (Q * hhi); 

                // std::cout << "Cost for solution i: " << cost_i << std::endl; 
                // std::cout << "Rotation: " << Rog << std::endl;                                           
                Matrix2 rot2 = projectRtoPrior(Rog);
                
                if ((rot2(0,0) < -0.5) || (rot2(1,1) < -0.5))
                {
                    rot2 *= -1;
                    tii *= -1;
                }
                             
                if (cost_i < cost_best * 1.05)
                {
                    cost_best = cost_i;  // let some margin hee    
                            
                    // Update bets
                    rot_s.push_back(rot2); 
                    trans_s.push_back(tii); 
                    normal_s.push_back(nii);
                    // Save results
                     Matrix3 Rr = Matrix3::Identity();                   
                    
                     Rr(0,0)=rot2(0,0); 
                     Rr(0,2)=rot2(0,1);
                     Rr(2,0)=rot2(1,0);
                     Rr(2,2)=rot2(1,1);    
                     Matrix3 Hjj = Rr + tii * nii.transpose();  
                     hom_s.push_back(Hjj);             
                }   
               
                
                
             }
             
             *H = hom_s; 
             *rot = rot_s;   
             *trans = trans_s; 
             *n = normal_s;
              
         }
         



         /** Use colmap implementation **/
         /* initialization based on DLT + pattern */     
         void computeInitPattern2(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 Matrix3& H, Matrix2 & rot, Vector3& trans, Vector3 & n)
         {
         
               Matrix9 Q; 
               Matrix3 Qr; 
               Matrix39 Qtr;
               createDataCReduced(bearingVectors0, bearingVectors1, weights_observations, Q, Qr, Qtr);
              
               
             Eigen::SelfAdjointEigenSolver<Matrix9> eig_qh(Q);
             Vector9 q_eig = eig_qh.eigenvalues(); 

            
             int id_min = 0; 
             double val_min = q_eig(0); 
             for (int i=1;i<9;i++)
             {
                if (q_eig(i) < val_min)
                        {val_min = q_eig(i); id_min = i;}
             }
             Vector9 qq = eig_qh.eigenvectors().col(id_min); 
               
             Matrix3 Hh = Eigen::Map<Matrix3> (qq.data(), 3, 3);  
             
             Matrix3 R3, Rproj; 
             PoseFromHomographyMatrix(Hh, &R3, &trans, &n);
             // try + H and - H
                                      
             const double ESii = SignOfNumber(R3(1, 1));
                
             
             // This is the decomposition for -Hh
             R3 *= ESii;
             trans *= ESii;
             
             // obtain feasible solution for the primal problem
             rot = projectRtoPrior(R3);
             
             Matrix3 Rr = Matrix3::Identity(); 
             Rr(0,0) = rot(0,0); 
             Rr(0,2) = rot(0,1); 
             Rr(2,0) = rot(1,0); 
             Rr(2,2) = rot(1,1); 
            
             // Check options
             Matrix3 Htemp = Rr + trans * n.transpose(); 
                         
             H = Htemp;                        
           
         }
         
       
         Matrix3 projectToHomography(const Matrix3 & H_hat, const Matrix2 & rot, Vector3 & trans, Vector3 & n)
         {
         
            Matrix3 R = Matrix3::Identity(); 
            R(0,0)=rot(0,0);
            R(0,2)=rot(0,1);
            R(2,0)=rot(1,0);
            R(2,2)=rot(1,1);
            
            Matrix3 RH; 
            RH = R.transpose() * H_hat; 
            
            Eigen::EigenSolver<Matrix3> eig_rh(RH);

            Vector3 s_eig; 
            s_eig = eig_rh.eigenvalues().real(); 

            
            double m_eig = s_eig(0), M_eig = s_eig(0);
            int id_m = 0, id_M = 0;
            if (s_eig(1) < m_eig)
            {m_eig = s_eig(1);id_m = 1;}
                
            if (s_eig(1) > M_eig)
            {M_eig = s_eig(1);id_M = 1;}
                 
            if (s_eig(2) < m_eig)
            {m_eig = s_eig(2);id_m = 2;}
                
            if (s_eig(2) > M_eig)
            {M_eig = s_eig(2);id_M = 2;}
            
            double m_t = 0;
            for (int i=0;i<3;i++)
            {
                if ((id_m == i) || (id_M == i))
                        continue; 
                // else
                m_t = s_eig(i);
                break;
            }

            
            Matrix3 RHI = RH - m_t * Matrix3::Identity(); 
            
                        
            const Eigen::JacobiSVD<Matrix3> svd(RHI, Eigen::ComputeFullU | Eigen::ComputeFullV);
            
            Vector3 s_svd = svd.singularValues(); 

            double max_svd = s_svd(0); 
            int id_max = 0; 
            if (s_svd(1) > max_svd)
            {
                max_svd = s_svd(1);         id_max = 1;        
            }
            
            if (s_svd(2) > max_svd)
            {
                max_svd = s_svd(2);         id_max = 2;
            }               
            // check determinant
            Matrix3 Ur = svd.matrixU();
            Matrix3 V = svd.matrixV();
            trans.setZero(); 
            n.setZero(); 
            trans = max_svd / m_t * R * Ur.col(id_max); //scaled translation
            n = V.col(id_max);             

            Matrix3 H = Matrix3::Zero(); 
            H = R + trans * n.transpose(); 
            
            return (H);
         }

        
        
}  // end of namespace RotPrior




