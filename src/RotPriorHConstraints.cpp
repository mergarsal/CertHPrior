#include "RotPriorHConstraints.h"


namespace RotPriorH{

 void createHRConstraints(std::vector<Matrix12>& As)
    {
          // auxiliary variables
          size_t h1=0,h4=h1+1,h7=h4+1,h2=h7+1,h5=h2+1,h8=h5+1,h3=h8+1,h6=h3+1,h9=h6+1; 
          size_t c=h9+1,s=c+1,k=s+1;

          Matrix12 Ai; 
          Ai.setZero(); 
          Ai(c,c)=1;
          Ai(s,s)=1;
          Ai(k,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(c,k)=1;
          Ai(h5,c)=-1;
          Ai(h1,h5)=1;
          Ai(h2,h4)=-1;
          Ai(h1,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h2,h6)=1;
          Ai(h3,h5)=-1;
          Ai(h3,k)=1;
          Ai(h5,s)=1;
          Ai(s,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h6)=1;
          Ai(c,h6)=-1;
          Ai(h3,h4)=-1;
          Ai(h4,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h4,h8)=1;
          Ai(h5,h7)=-1;
          Ai(k,h7)=1;
          Ai(s,h5)=-1;
          Ai(s,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(c,k)=1;
          Ai(h5,c)=-1;
          Ai(h5,h9)=1;
          Ai(h6,h8)=-1;
          Ai(h9,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h4,h9)=1;
          Ai(h4,c)=-1;
          Ai(h7,h6)=-1;
          Ai(s,h6)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h8)=1;
          Ai(c,h8)=-1;
          Ai(h2,h7)=-1;
          Ai(h2,s)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h2,h9)=1;
          Ai(h2,c)=-1;
          Ai(h3,h8)=-1;
          Ai(h8,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h9)=1;
          Ai(c,h9)=-1;
          Ai(h1,c)=-1;
          Ai(h3,h7)=-1;
          Ai(h3,s)=-1;
          Ai(s,h7)=1;
          Ai(c,c)=1;
          Ai(s,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          
          // norm on H
          Ai.setZero();
          Ai(h1,h1)=1;
          Ai(h2,h2)=1;
          Ai(h3,h3)=1;
          Ai(h4,h4)=1;
          Ai(h5,h5)=1;
          Ai(h6,h6)=1;
          Ai(h7,h7)=1;
          Ai(h8,h8)=1;
          Ai(h9,h9)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));  
          
          
          // value on k
          /*
          Ai.setZero();
          Ai(k,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose())); 
          */ 
    } 
    
    void createRtnConstraints(std::vector<Matrix12>& As)
    {
          // auxiliary variables
          size_t q1=0,q4=q1+1,q7=q4+1,q2=q7+1,q5=q2+1,q8=q5+1,q3=q8+1,q6=q3+1,q9=q6+1;  
          size_t c=q9+1,s=c+1,k=s+1;
          
          Matrix12 Ai; 
          Ai.setZero();
          Ai(c,c)=1;
          Ai(s,s)=1;
          Ai(k,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q5)=1;
          Ai(q2,q4)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q2,q6)=1;
          Ai(q3,q5)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q6)=1;
          Ai(q3,q4)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q4,q8)=1;
          Ai(q5,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q5,q9)=1;
          Ai(q6,q8)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q4,q9)=1;
          Ai(q6,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q8)=1;
          Ai(q2,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q2,q9)=1;
          Ai(q3,q8)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q9)=1;
          Ai(q3,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          // for k = 1
          /*
          Ai.setZero();
          Ai(k,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));          
          */
          // for norm(R+Q) = 1
          Ai.setZero();
          Ai(c,c)=2;Ai(c,q1)=2;Ai(c,q9)=2;Ai(k,k)=1;Ai(k,q5)=2;
          Ai(q1,q1)=1;Ai(q2,q2)=1;Ai(q3,q3)=1;Ai(q3,s)=2;Ai(q4,q4)=1;Ai(q5,q5)=1;
          Ai(q6,q6)=1;Ai(q7,q7)=1;Ai(q7,s)=-2;Ai(q8,q8)=1;Ai(q9,q9)=1;Ai(s,s)=2;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          
    }


 void createHRtnConstraints(std::vector<Matrix21>& As)
    {
          // auxiliary variables
          size_t h1=0,h4=h1+1,h7=h4+1,h2=h7+1,h5=h2+1,h8=h5+1,h3=h8+1,h6=h3+1,h9=h6+1; 
          size_t c=h9+1,s=c+1,k=s+1;
          size_t q1=k+1,q4=q1+1,q7=q4+1,q2=q7+1,q5=q2+1,q8=q5+1,q3=q8+1,q6=q3+1,q9=q6+1; 

          /** Constraints for HR **/
          Matrix21 Ai; 
          Ai.setZero(); 
          Ai(c,c)=1;
          Ai(s,s)=1;
          Ai(k,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(c,k)=1;
          Ai(h5,c)=-1;
          Ai(h1,h5)=1;
          Ai(h2,h4)=-1;
          Ai(h1,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h2,h6)=1;
          Ai(h3,h5)=-1;
          Ai(h3,k)=1;
          Ai(h5,s)=1;
          Ai(s,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h6)=1;
          Ai(c,h6)=-1;
          Ai(h3,h4)=-1;
          Ai(h4,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h4,h8)=1;
          Ai(h5,h7)=-1;
          Ai(k,h7)=1;
          Ai(s,h5)=-1;
          Ai(s,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(c,k)=1;
          Ai(h5,c)=-1;
          Ai(h5,h9)=1;
          Ai(h6,h8)=-1;
          Ai(h9,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h4,h9)=1;
          Ai(h4,c)=-1;
          Ai(h7,h6)=-1;
          Ai(s,h6)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h8)=1;
          Ai(c,h8)=-1;
          Ai(h2,h7)=-1;
          Ai(h2,s)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h2,h9)=1;
          Ai(h2,c)=-1;
          Ai(h3,h8)=-1;
          Ai(h8,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h9)=1;
          Ai(c,h9)=-1;
          Ai(h1,c)=-1;
          Ai(h3,h7)=-1;
          Ai(h3,s)=-1;
          Ai(s,h7)=1;
          Ai(c,c)=1;
          Ai(s,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h1)=1;
          Ai(h2,h2)=1;
          Ai(h3,h3)=1;
          Ai(h4,h4)=1;
          Ai(h5,h5)=1;
          Ai(h6,h6)=1;
          Ai(h7,h7)=1;
          Ai(h8,h8)=1;
          Ai(h9,h9)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));  
          
          /** Constraints for Rtn **/
          Ai.setZero();
          Ai(q1,q5)=1;
          Ai(q2,q4)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q2,q6)=1;
          Ai(q3,q5)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q6)=1;
          Ai(q3,q4)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q4,q8)=1;
          Ai(q5,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q5,q9)=1;
          Ai(q6,q8)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q4,q9)=1;
          Ai(q6,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q8)=1;
          Ai(q2,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q2,q9)=1;
          Ai(q3,q8)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q9)=1;
          Ai(q3,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          // for k = 1
          /*
          Ai.setZero();
          Ai(k,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));          
          */
          // for norm(R+Q) = 1
          Ai.setZero();
          Ai(c,c)=2;Ai(c,q1)=2;Ai(c,q9)=2;Ai(k,k)=1;Ai(k,q5)=2;
          Ai(q1,q1)=1;Ai(q2,q2)=1;Ai(q3,q3)=1;Ai(q3,s)=2;Ai(q4,q4)=1;Ai(q5,q5)=1;
          Ai(q6,q6)=1;Ai(q7,q7)=1;Ai(q7,s)=-2;Ai(q8,q8)=1;Ai(q9,q9)=1;Ai(s,s)=2;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          /** Constraints for HH **/
          Ai.setZero();
          Ai(h2,h4)=1;
          Ai(c,h1)=-1;
          Ai(h3,h7)=1;
          Ai(h1,q1)=-1;
          Ai(h4,q2)=-1;
          Ai(h7,q3)=-1;
          Ai(h7,s)=-1;
          Ai(h1,h1)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h2)=1,Ai(c,h2)=-1;Ai(h2,h5)=1;
          Ai(h3,h8)=1;Ai(h2,q1)=-1;Ai(h5,q2)=-1;
          Ai(h8,q3)=-1;Ai(h8,s)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h3)=1;Ai(c,h3)=-1;Ai(h2,h6)=1;Ai(h3,h9)=1;Ai(h3,q1)=-1;
          Ai(h6,q2)=-1;Ai(h9,q3)=-1;Ai(h9,s)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
        Ai.setZero();
        Ai.setZero(); Ai(h1,h4)=1;Ai(h4,h5)=1;Ai(h6,h7)=1;
        Ai(h4,k)=-1;Ai(h1,q4)=-1;Ai(h4,q5)=-1;Ai(h7,q6)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
      
        Ai.setZero();
        Ai(h5,h5)=1;Ai(h2,h4)=1;Ai(h6,h8)=1;
        Ai(h5,k)=-1;Ai(h2,q4)=-1;Ai(h5,q5)=-1;Ai(h8,q6)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
      
        Ai.setZero();
        Ai(h3,h4)=1;Ai(h5,h6)=1;Ai(h6,h9)=1;
        Ai(h6,k)=-1;Ai(h3,q4)=-1;Ai(h6,q5)=-1;Ai(h9,q6)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

      
        Ai.setZero();
        Ai(h1,h7)=1;Ai(h4,h8)=1;Ai(h9,h7)=1;Ai(h7,q9)=-1;
        Ai(c,h7)=-1;Ai(h1,q7)=-1;Ai(h4,q8)=-1;Ai(h1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(h2,h7)=1;Ai(h5,h8)=1;Ai(h8,h9)=1;
        Ai(h2,q7)=-1;Ai(h5,q8)=-1;Ai(h8,q9)=-1;Ai(h2,s)=1;Ai(c,h8)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(h9,h9)=1;Ai(h3,h7)=1;Ai(h6,h8)=1;
        Ai(h3,q7)=-1;Ai(h6,q8)=-1;Ai(h9,q9)=-1;Ai(h3,s)=1;Ai(c,h9)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        
        /** Constraints for HQ **/
        /* HQt */
        Ai.setZero();
        Ai(h1, q1) = 1; Ai(h2, q2) = 1; Ai(h3, q3) = 1; 
        Ai(q1, q1) = -1; Ai(c, q1) = -1; Ai(q2, q2) = -1; Ai(s, q3) = -1; Ai(q3, q3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));

 
        Ai.setZero();
        Ai(h1, q4) = 1; Ai(h2, q5) = 1; Ai(h3, q6) = 1; 
        Ai(c, q4) = -1; Ai(q1, q4) = -1; Ai(q2, q5) = -1; Ai(q3, q6) = -1; Ai(q6, s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));

       
        Ai.setZero();
        Ai(h1, q7) = 1; Ai(h2, q8) = 1; Ai(h3, q9) = 1; 
        Ai(c, q7) = -1; Ai(q1, q7) = -1; Ai(q2, q8) = -1; Ai(q3, q9) = -1; Ai(q9, s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));

      
        Ai.setZero();
        Ai(h4, q1) = 1; Ai(h5, q2) = 1; Ai(h6, q3) = 1; 
        Ai(k, q2) = -1; Ai(q1, q4) = -1; Ai(q2, q5) = -1; Ai(q3, q6) = -1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

       
        Ai.setZero();
        Ai(h4, q4) = 1; Ai(h5, q5) = 1; Ai(h6, q6) = 1; 
        Ai(q4, q4) = -1; Ai(q5, q5) = -1; Ai(k, q5) = -1; Ai(q6, q6)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));

      
        Ai.setZero();
        Ai(h4, q7) = 1; Ai(h5, q8) = 1; Ai(h6, q9) = 1; 
        Ai(k, q8) = -1; Ai(q4, q7) = -1; Ai(q5, q8) = -1; Ai(q6, q9) = -1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        
        Ai.setZero();
        Ai(h7, q1) = 1; Ai(h8, q2) = 1; Ai(h9, q3) = 1; 
        Ai(c, q3) = -1; Ai(q1, q7) = -1; Ai(q2, q8) = -1; Ai(q3, q9) = -1;  Ai(q1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

      
        Ai.setZero();
        Ai(h7, q4) = 1; Ai(h8, q5) = 1; Ai(h9, q6) = 1; 
        Ai(c, q6) = -1; Ai(q4, q7) = -1; Ai(q5, q8) = -1; Ai(q6, q9)=-1; Ai(q4,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

       
        Ai.setZero();
        Ai(h7, q7) = 1; Ai(h8, q8) = 1; Ai(h9, q9) = 1; 
        Ai(q7, q7) = -1; Ai(q8, q8) = -1; Ai(q9, q9) = -1; Ai(c, q9) = -1; Ai(q7,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
                  
                  
        /* QH */
        Ai.setZero();
        Ai(h1, q1) = 1; Ai(q2, h4) = 1; Ai(q3, h7) = 1; 
        Ai(q1, q1) = -1; Ai(c, q1) = -1; Ai(q2, q4) = -1; Ai(q3, q7) = -1; Ai(q3,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
    
        Ai.setZero();
        Ai(q1, h2) = 1; Ai(q2, h5) = 1; Ai(q3, h8) = 1; 
        Ai(k, q2) = -1; Ai(q1, q2) = -1; Ai(q2, q5) = -1; Ai(q3, q8) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
     
        Ai.setZero();
        Ai(q1, h3) = 1; Ai(q2, h6) = 1; Ai(q3, h9) = 1; 
        Ai(c, q3) = -1; Ai(q1, q3) = -1; Ai(q2, q6) = -1; Ai(q3, q9) = -1; Ai(q1,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
 
        Ai.setZero();
        Ai(q4, h1) = 1; Ai(q5, h4) = 1; Ai(q6, h7) = 1; 
        Ai(c, q4) = -1; Ai(q1, q4) = -1; Ai(q4, q5) = -1; Ai(q6, q7) = -1; Ai(q6,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(q4, h2) = 1; Ai(q5, h5) = 1; Ai(q6, h8) = 1; 
        Ai(q5, q5) = -1; Ai(k, q5) = -1; Ai(q2, q4) = -1; Ai(q6, q8) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
   
        Ai.setZero();
        Ai(q4, h3) = 1; Ai(q5, h6) = 1; Ai(q6, h9) = 1; 
        Ai(c, q6) = -1; Ai(q3, q4) = -1; Ai(q5, q6) = -1; Ai(q6, q9) = -1; Ai(q4,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
       
        Ai.setZero();
        Ai(q7, h1) = 1; Ai(q8, h4) = 1; Ai(q9, h7) = 1; 
        Ai(c, q7) = -1; Ai(q1, q7) = -1; Ai(q4, q8) = -1; Ai(q7, q9) = -1; Ai(s, q9)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
    
        Ai.setZero();
        Ai(q7, h2) = 1; Ai(q8, h5) = 1; Ai(q9, h8) = 1; 
        Ai(k, q8) = -1; Ai(q2, q7) = -1; Ai(q5, q8) = -1; Ai(q8, q9) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(q7, h3) = 1; Ai(q8, h6) = 1; Ai(q9, h9) = 1; 
        Ai(q9, q9) = -1; Ai(c, q9) = -1; Ai(q3, q7) = -1; Ai(q6, q8) = -1; Ai(q7,s)=-1;          
        As.push_back(0.5 * (Ai + Ai.transpose()));                 
                  
                  
        /** Constraints for HR **/
        /* HRt */  
        // Ai.setZero();        
        // Ai(c, h1) = 1; Ai(h3, s)=1; 
        // Ai(q1,c)=-1;Ai(c,c)=-1;Ai(q3,s)=-1;Ai(s,s)=-1;
        // As.push_back(0.5 * (Ai + Ai.transpose())); // 48
        
       // Ai.setZero(); Ai(k, h2) = 1; 
       // Ai(q2,k)=-1; 
       // As.push_back(0.5 * (Ai + Ai.transpose()));  // 49
      
       // Ai.setZero(); Ai(c, h3)=1;Ai(c,q3)=-1;Ai(h1,s)=-1;Ai(q1,s)=1;
       // As.push_back(0.5 * (Ai + Ai.transpose())); // 50

       
        Ai.setZero(); Ai(c, h4) = 1; Ai(h6, s)=1; 
        Ai(q4,c)=-1;Ai(q6,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
       
        Ai.setZero(); Ai(k, h5) = -1;
        Ai(q5,k)=1;Ai(k,k)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
      
         Ai.setZero(); Ai(c, h6) = 1;Ai(c,q6)=-1;Ai(h4,s)=-1;Ai(q4,s)=1;
         As.push_back(0.5 * (Ai + Ai.transpose()));  

       
        // Ai.setZero(); Ai(c, h7) = 1; Ai(h9, s)=1; Ai(c,q7)=-1;Ai(q9,s)=-1;
        // As.push_back(0.5 * (Ai + Ai.transpose())); // 54
       
         Ai.setZero(); Ai(k, h8) = 1;Ai(q8,k)=-1;
         As.push_back(0.5 * (Ai + Ai.transpose())); 
       
        // Ai.setZero(); Ai(c, h9) = 1;Ai(q9,c)=-1;Ai(h7,s)=-1;Ai(q7,s)=1;
        // Ai(c,c)=-1;Ai(s,s)=-1;         
        // As.push_back(0.5 * (Ai + Ai.transpose()));  // 56                  
                  
        /* RtH */
        Ai.setZero(); 
        Ai(c,h1)=1;Ai(h7,s)=-1;Ai(q7,s)=1;Ai(q1,c)=-1;
        Ai(c,c)=-1;Ai(s,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
      
        Ai.setZero(); Ai(c,h2)=1;Ai(h8,s)=-1; 
        Ai(q2,c)=-1;Ai(q8,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
       
        Ai.setZero(); Ai(c,h3)=1;Ai(h9,s)=-1;
        Ai(q3,c)=-1;Ai(q9,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 

    
        Ai.setZero(); Ai(k,h4)=1;Ai(q4,k)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
       
       // Ai.setZero(); Ai(k,h5)=-1;Ai(k,k)=1;Ai(q5,k)=1;
       // As.push_back(0.5 * (Ai + Ai.transpose())); // 61
       
        Ai.setZero(); Ai(k,h6)=1;Ai(q6,k)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 

       
        Ai.setZero(); Ai(c,h7)=1;Ai(q7,c)=-1;Ai(h1,s)=1;Ai(q1,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
      
        Ai.setZero(); Ai(c,h8)=1;Ai(q8,c)=-1;Ai(h2,s)=1;Ai(q2,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
       
        // Ai.setZero(); Ai(c,h9)=1;Ai(q9,c)=-1;Ai(q3,s)=-1;Ai(h3,s)=1;
        // Ai(c,c)=-1;Ai(s,s)=-1;        
        // As.push_back(0.5 * (Ai + Ai.transpose()));   // 65 
                  
                  
        /* HR */
         Ai.setZero(); 
         Ai(c,h1)=1;Ai(h3,s)=-1;
         Ai(c,c)=-1;Ai(s,s)=1;Ai(q1,c)=-1;Ai(q3,s)=1;
         As.push_back(0.5 * (Ai + Ai.transpose())); 
        
         Ai.setZero();  Ai(k,h2)=1;Ai(k,q2)=-1;
         As.push_back(0.5 * (Ai + Ai.transpose()));  
       
       // Ai.setZero();  Ai(c,h3)=1;Ai(c,q3)=-1;Ai(c,s)=-2;Ai(h1,s)=1;Ai(q1,s)=-1;
       //  As.push_back(0.5 * (Ai + Ai.transpose()));  // 68

       
         Ai.setZero();  Ai(c,h4)=1;Ai(h6,s)=-1;
         Ai(c,q4)=-1;Ai(q6,s)=1;
         As.push_back(0.5 * (Ai + Ai.transpose())); 
     
        // Ai.setZero();  Ai(k,h5)=-1;Ai(q5,k)=1;Ai(k,k)=1;
        // As.push_back(0.5 * (Ai + Ai.transpose()));  // 70
       
         Ai.setZero();  Ai(c,h6)=1;Ai(c,q6)=-1;Ai(h4,s)=1;Ai(q4,s)=-1;
         As.push_back(0.5 * (Ai + Ai.transpose()));  

       
        //  Ai.setZero();  Ai(c,h7)=1;Ai(h9,s)=-1;Ai(q7,c)=-1;Ai(q9,s)=1;Ai(c,s)=2;
        // As.push_back(0.5 * (Ai + Ai.transpose())); // 72
       
        // Ai.setZero();  Ai(k,h8)=1;Ai(k,q8)=-1;
        // As.push_back(0.5 * (Ai + Ai.transpose()));  // 73
    
        // Ai.setZero();  Ai(c,h9)=1;Ai(c,q9)=-1;Ai(h7,s)=1;Ai(q7,s)=-1;
        // Ai(c,c)=-1;Ai(s,s)=1;    
        // As.push_back(0.5 * (Ai + Ai.transpose()));       // 74    
                  
        /* RH */
        Ai.setZero(); 
        Ai(c,h1)=1;Ai(h7,s)=1;
        Ai(c,c)=-1;Ai(s,s)=1;Ai(q1,c)=-1;Ai(q7,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        Ai.setZero(); Ai(c,h2)=1;Ai(h8,s)=1;Ai(q2,c)=-1;Ai(q8,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
         Ai.setZero(); Ai(c,h3)=1;Ai(h9,s)=1;
         Ai(q3,c)=-1;Ai(q9,s)=-1;Ai(c,s)=-2;
         As.push_back(0.5 * (Ai + Ai.transpose())); 

        
        // Ai.setZero(); Ai(k,h4)=1;Ai(k,q4)=-1;
        // As.push_back(0.5 * (Ai + Ai.transpose()));   // 78
       
        // Ai.setZero(); Ai(k,h5)=-1;Ai(k,k)=1;Ai(q5,k)=1;
        // As.push_back(0.5 * (Ai + Ai.transpose())); // 79
       
        // Ai.setZero(); Ai(k,h6)=1;Ai(q6,k)=-1;
        // As.push_back(0.5 * (Ai + Ai.transpose())); // 80
        

        
        Ai.setZero(); Ai(c,h7)=1;Ai(q7,c)=-1;Ai(c,s)=2;Ai(h1,s)=-1;Ai(q1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        Ai.setZero(); Ai(c,h8)=1;Ai(q8,c)=-1;Ai(h2,s)=-1;Ai(q2,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        Ai.setZero(); Ai(c,h9)=1;Ai(c,q9)=-1;Ai(h3,s)=-1;Ai(q3,s)=1;
        Ai(c,c)=-1;Ai(s,s)=1;          
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 83 (index from 0)           
                       
                  
    } 
    
    
     void createHRtnConstraintsRed(std::vector<Matrix21>& As)
    {
          // auxiliary variables
          size_t h1=0,h4=h1+1,h7=h4+1,h2=h7+1,h5=h2+1,h8=h5+1,h3=h8+1,h6=h3+1,h9=h6+1; 
          size_t c=h9+1,s=c+1,k=s+1;
          size_t q1=k+1,q4=q1+1,q7=q4+1,q2=q7+1,q5=q2+1,q8=q5+1,q3=q8+1,q6=q3+1,q9=q6+1; 

          /** Constraints for HR **/
          Matrix21 Ai; 
          Ai.setZero(); 
          Ai(c,c)=1;
          Ai(s,s)=1;
          Ai(k,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(c,k)=1;
          Ai(h5,c)=-1;
          Ai(h1,h5)=1;
          Ai(h2,h4)=-1;
          Ai(h1,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h2,h6)=1;
          Ai(h3,h5)=-1;
          Ai(h3,k)=1;
          Ai(h5,s)=1;
          Ai(s,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h6)=1;
          Ai(c,h6)=-1;
          Ai(h3,h4)=-1;
          Ai(h4,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h4,h8)=1;
          Ai(h5,h7)=-1;
          Ai(k,h7)=1;
          Ai(s,h5)=-1;
          Ai(s,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(c,k)=1;
          Ai(h5,c)=-1;
          Ai(h5,h9)=1;
          Ai(h6,h8)=-1;
          Ai(h9,k)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));  // 5
          
          Ai.setZero();
          Ai(h4,h9)=1;
          Ai(h4,c)=-1;
          Ai(h7,h6)=-1;
          Ai(s,h6)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h8)=1;
          Ai(c,h8)=-1;
          Ai(h2,h7)=-1;
          Ai(h2,s)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h2,h9)=1;
          Ai(h2,c)=-1;
          Ai(h3,h8)=-1;
          Ai(h8,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
           
          Ai.setZero();
          Ai(h1,h9)=1;
          Ai(c,h9)=-1;
          Ai(h1,c)=-1;
          Ai(h3,h7)=-1;
          Ai(h3,s)=-1;
          Ai(s,h7)=1;
          Ai(c,c)=1;
          Ai(s,s)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(h1,h1)=1;
          Ai(h2,h2)=1;
          Ai(h3,h3)=1;
          Ai(h4,h4)=1;
          Ai(h5,h5)=1;
          Ai(h6,h6)=1;
          Ai(h7,h7)=1;
          Ai(h8,h8)=1;
          Ai(h9,h9)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));   // 10
          
          /** Constraints for Rtn **/
          Ai.setZero();
          Ai(q1,q5)=1;
          Ai(q2,q4)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q2,q6)=1;
          Ai(q3,q5)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q6)=1;
          Ai(q3,q4)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q4,q8)=1;
          Ai(q5,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q5,q9)=1;
          Ai(q6,q8)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));  //15
          
          Ai.setZero();
          Ai(q4,q9)=1;
          Ai(q6,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q8)=1;
          Ai(q2,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q2,q9)=1;
          Ai(q3,q8)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          Ai.setZero();
          Ai(q1,q9)=1;
          Ai(q3,q7)=-1;
          As.push_back(0.5 * (Ai + Ai.transpose()));
          
          // for k = 1
          /*
          Ai.setZero();
          Ai(k,k)=1;
          As.push_back(0.5 * (Ai + Ai.transpose()));          
          */
          // for norm(R+Q) = 1
          Ai.setZero();
          Ai(c,c)=2;Ai(c,q1)=2;Ai(c,q9)=2;Ai(k,k)=1;Ai(k,q5)=2;
          Ai(q1,q1)=1;Ai(q2,q2)=1;Ai(q3,q3)=1;Ai(q3,s)=2;Ai(q4,q4)=1;Ai(q5,q5)=1;
          Ai(q6,q6)=1;Ai(q7,q7)=1;Ai(q7,s)=-2;Ai(q8,q8)=1;Ai(q9,q9)=1;Ai(s,s)=2;
          As.push_back(0.5 * (Ai + Ai.transpose()));  // 20
          
        
        
        /** Constraints for HQ **/
        /* HQt */
        /*
        Ai.setZero();
        Ai(h1, q1) = 1; Ai(h2, q2) = 1; Ai(h3, q3) = 1; 
        Ai(q1, q1) = -1; Ai(c, q1) = -1; Ai(q2, q2) = -1; Ai(s, q3) = -1; Ai(q3, q3)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 21
        */
 
        Ai.setZero();
        Ai(h1, q4) = 1; Ai(h2, q5) = 1; Ai(h3, q6) = 1; 
        Ai(c, q4) = -1; Ai(q1, q4) = -1; Ai(q2, q5) = -1; Ai(q3, q6) = -1; Ai(q6, s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 22

        /*
        Ai.setZero();
        Ai(h1, q7) = 1; Ai(h2, q8) = 1; Ai(h3, q9) = 1; 
        Ai(c, q7) = -1; Ai(q1, q7) = -1; Ai(q2, q8) = -1; Ai(q3, q9) = -1; Ai(q9, s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   //23
        */
      
        Ai.setZero();
        Ai(h4, q1) = 1; Ai(h5, q2) = 1; Ai(h6, q3) = 1; 
        Ai(k, q2) = -1; Ai(q1, q4) = -1; Ai(q2, q5) = -1; Ai(q3, q6) = -1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        /*
        Ai.setZero();
        Ai(h4, q4) = 1; Ai(h5, q5) = 1; Ai(h6, q6) = 1; 
        Ai(q4, q4) = -1; Ai(q5, q5) = -1; Ai(k, q5) = -1; Ai(q6, q6)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   //25
        */
        
        /*
        Ai.setZero();
        Ai(h4, q7) = 1; Ai(h5, q8) = 1; Ai(h6, q9) = 1; 
        Ai(k, q8) = -1; Ai(q4, q7) = -1; Ai(q5, q8) = -1; Ai(q6, q9) = -1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //26
        */
        
        /*
        Ai.setZero();
        Ai(h7, q1) = 1; Ai(h8, q2) = 1; Ai(h9, q3) = 1; 
        Ai(c, q3) = -1; Ai(q1, q7) = -1; Ai(q2, q8) = -1; Ai(q3, q9) = -1;  Ai(q1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //27
        */
        
      
        Ai.setZero();
        Ai(h7, q4) = 1; Ai(h8, q5) = 1; Ai(h9, q6) = 1; 
        Ai(c, q6) = -1; Ai(q4, q7) = -1; Ai(q5, q8) = -1; Ai(q6, q9)=-1; Ai(q4,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        /*
        Ai.setZero();
        Ai(h7, q7) = 1; Ai(h8, q8) = 1; Ai(h9, q9) = 1; 
        Ai(q7, q7) = -1; Ai(q8, q8) = -1; Ai(q9, q9) = -1; Ai(c, q9) = -1; Ai(q7,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //29
        */          
        
        
        
        
        /* HQ - R*Q - Q*Q */
        Ai.setZero();
        Ai(h1, q1) = 1; Ai(c, q1) = -1; Ai(q4, h2) = 1; Ai(q7, h3) = 1; 
        Ai(q1, q1) = -1; Ai(q2, q4) = -1; Ai(q3, q7) = -1; Ai(q7,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //30
        
        Ai.setZero();
        Ai(h4, q1) = 1; Ai(q4, h5) = 1; Ai(q7, h6) = 1; 
        Ai(k, q4) = -1; Ai(q1, q4) = -1; Ai(q5, q4) = -1; Ai(q7,q6)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        // h7*q1 - c*q7 + h8*q4 + h9*q7 - q1*q7 - q4*q8 - q7*q9 + q1*s
        Ai.setZero();
        Ai(h7, q1) = 1; Ai(c, q7)=-1; Ai(h8,q4)=1; 
        Ai(h9,q7) = 1; Ai(q1, q7)=-1; Ai(q4, q8)=-1; 
        Ai(q7,q9) = -1; Ai(q1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
        // h1*q2 - c*q2 + h2*q5 + h3*q8 - q1*q2 - q2*q5 - q3*q8 - q8*s
        Ai.setZero();
        Ai(h1, q2) = 1; Ai(c, q2) = -1; Ai(q5, h2) = 1; Ai(h3,q8) = 1;
        Ai(q1, q2) = -1; Ai(q2, q5) = -1; Ai(q3, q8) = -1; Ai(q8,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
        // h4*q2 + h5*q5 + h6*q8 - k*q5 - q2*q4 - q6*q8 - q5^2
        Ai.setZero();
        Ai(h4, q2) = 1; Ai(h5, q5) = 1; Ai(q8, h6) = 1; Ai(k,q5) = -1;
        Ai(q2, q4) = -1; Ai(q6, q8) = -1; Ai(q5, q5) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
        // h7*q2 - c*q8 + h8*q5 + h9*q8 - q2*q7 - q5*q8 - q8*q9 + q2*s
        Ai.setZero();
        Ai(h7, q2) = 1; Ai(h8, q5) = 1; Ai(q8, c) = -1; Ai(h9,q8) = 1;
        Ai(q2, q7) = -1; Ai(q5, q8) = -1; Ai(q8, q9) = -1; Ai(q2,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //35
        
        // h1*q3 - c*q3 + h2*q6 + h3*q9 - q1*q3 - q2*q6 - q3*q9 - q9*s
        Ai.setZero();
        Ai(h1, q3) = 1; Ai(c, q3) = -1; Ai(q6, h2) = 1; Ai(h3,q9) = 1;
        Ai(q1, q3) = -1; Ai(q2, q6) = -1; Ai(q3, q9) = -1; Ai(q9,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h4*q3 + h5*q6 + h6*q9 - k*q6 - q3*q4 - q5*q6 - q6*q9
        Ai.setZero();
        Ai(h4, q3) = 1; Ai(h5, q6) = 1; Ai(q9, h6) = 1; Ai(k,q6) = -1;
        Ai(q4, q3) = -1; Ai(q5, q6) = -1; Ai(q6, q9) = -1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   
        
        // h7*q3 - c*q9 + h8*q6 + h9*q9 - q3*q7 - q6*q8 + q3*s - q9^2
        Ai.setZero();
        Ai(h7, q3) = 1; Ai(c, q9) = -1; Ai(q6, h8) = 1; Ai(h9,q9) = 1;
        Ai(q7, q3) = -1; Ai(q8, q6) = -1; Ai(q9, q9) = -1; Ai(q3,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        
        /* QtH */
        // h1*q1 - c*q1 + h4*q4 + h7*q7 + q7*s - q1^2 - q4^2 - q7^2
        /*
        Ai.setZero();
        Ai(h1, q1) = 1; Ai(c, q1) = -1;  Ai(q4, h4) = 1; Ai(q7, h7) = 1; 
        Ai(q1, q1) = -1; Ai(q4, q4) = -1; Ai(q7, q7) = -1; Ai(q7,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //39
        */
              
        // h1*q2 - c*q2 + h4*q5 + h7*q8 - q1*q2 - q4*q5 - q7*q8 + q8*s
        /*
        Ai.setZero();
        Ai(h1, q2) = 1; Ai(c, q2) = -1;  Ai(q5, h4) = 1; Ai(q8, h7) = 1; 
        Ai(q1, q2) = -1; Ai(q4, q5) = -1; Ai(q7, q8) = -1; Ai(q8,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //40
        */
        
        // h1*q3 - c*q3 + h4*q6 + h7*q9 - q1*q3 - q4*q6 - q7*q9 + q9*s
        Ai.setZero();
        Ai(h1, q3) = 1; Ai(c, q3) = -1;  Ai(q6, h4) = 1; Ai(q9, h7) = 1; 
        Ai(q1, q3) = -1; Ai(q4, q6) = -1; Ai(q7, q9) = -1; Ai(q9,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        // h2*q1 + h5*q4 + h8*q7 - k*q4 - q1*q2 - q4*q5 - q7*q8
       /*
        Ai.setZero();
        Ai(h2, q1) = 1; Ai(h5, q4) = 1;  Ai(q7, h8) = 1; Ai(k, q4) = -1; 
        Ai(q1, q2) = -1; Ai(q4, q5) = -1; Ai(q7, q8) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 42
        */
        
        // h2*q2 + h5*q5 + h8*q8 - k*q5 - q2^2 - q5^2 - q8^2
        /*
        Ai.setZero();
        Ai(h2, q2) = 1; Ai(h5, q5) = 1;  Ai(q8, h8) = 1; Ai(k, q5) = -1; 
        Ai(q2, q2) = -1; Ai(q5, q5) = -1; Ai(q8, q8) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  //43
        */
        
        // h2*q3 + h5*q6 + h8*q9 - k*q6 - q2*q3 - q5*q6 - q8*q9
        /*
        Ai.setZero();
        Ai(h2, q3) = 1; Ai(h5, q6) = 1;  Ai(q9, h8) = 1; Ai(k, q6) = -1; 
        Ai(q2, q3) = -1; Ai(q5, q6) = -1; Ai(q8, q9) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   //44
        */
        
        
        // h3*q1 - c*q7 + h6*q4 + h9*q7 - q1*q3 - q4*q6 - q7*q9 - q1*s
        Ai.setZero();
        Ai(q1, h3) = 1; Ai(c, q7) = -1; Ai(h6, q4) = 1;  Ai(q7, h9) = 1; 
        Ai(q1, q3) = -1; Ai(q4, q6) = -1; Ai(q7, q9) = -1; Ai(q1, s) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  //45
        
        // h3*q2 - c*q8 + h6*q5 + h9*q8 - q2*q3 - q5*q6 - q8*q9 - q2*s
        /*
        Ai.setZero();
        Ai(h3, q2) = 1; Ai(c, q8) = -1; Ai(h6, q5) = 1;  Ai(q8, h9) = 1; 
        Ai(q2, q3) = -1; Ai(q5, q6) = -1; Ai(q8, q9) = -1; Ai(q2, s) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  //46
        */
            
        // h3*q3 - c*q9 + h6*q6 + h9*q9 - q3*s - q3^2 - q6^2 - q9^2
        /*
        Ai.setZero();
        Ai(h3, q3) = 1; Ai(c, q9) = -1; Ai(h6, q6) = 1;  Ai(q9, h9) = 1; 
        Ai(q3, q3) = -1; Ai(q6, q6) = -1; Ai(q9, q9) = -1; Ai(q3, s) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));    //47
        */
              
              
        /* QH */
        Ai.setZero();
        Ai(h1, q1) = 1; Ai(q2, h4) = 1; Ai(q3, h7) = 1; 
        Ai(q1, q1) = -1; Ai(c, q1) = -1; Ai(q2, q4) = -1; Ai(q3, q7) = -1; Ai(q3,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
    
        Ai.setZero();
        Ai(q1, h2) = 1; Ai(q2, h5) = 1; Ai(q3, h8) = 1; 
        Ai(k, q2) = -1; Ai(q1, q2) = -1; Ai(q2, q5) = -1; Ai(q3, q8) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));
     
        Ai.setZero();
        Ai(q1, h3) = 1; Ai(q2, h6) = 1; Ai(q3, h9) = 1; 
        Ai(c, q3) = -1; Ai(q1, q3) = -1; Ai(q2, q6) = -1; Ai(q3, q9) = -1; Ai(q1,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //50
 
        Ai.setZero();
        Ai(q4, h1) = 1; Ai(q5, h4) = 1; Ai(q6, h7) = 1; 
        Ai(c, q4) = -1; Ai(q1, q4) = -1; Ai(q4, q5) = -1; Ai(q6, q7) = -1; Ai(q6,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        /*
        Ai.setZero();
        Ai(q4, h2) = 1; Ai(q5, h5) = 1; Ai(q6, h8) = 1; 
        Ai(q5, q5) = -1; Ai(k, q5) = -1; Ai(q2, q4) = -1; Ai(q6, q8) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   // 52
        */
        
        Ai.setZero();
        Ai(q4, h3) = 1; Ai(q5, h6) = 1; Ai(q6, h9) = 1; 
        Ai(c, q6) = -1; Ai(q3, q4) = -1; Ai(q5, q6) = -1; Ai(q6, q9) = -1; Ai(q4,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
       
        Ai.setZero();
        Ai(q7, h1) = 1; Ai(q8, h4) = 1; Ai(q9, h7) = 1; 
        Ai(c, q7) = -1; Ai(q1, q7) = -1; Ai(q4, q8) = -1; Ai(q7, q9) = -1; Ai(s, q9)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
    
        Ai.setZero();
        Ai(q7, h2) = 1; Ai(q8, h5) = 1; Ai(q9, h8) = 1; 
        Ai(k, q8) = -1; Ai(q2, q7) = -1; Ai(q5, q8) = -1; Ai(q8, q9) = -1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));  //55
        
        Ai.setZero();
        Ai(q7, h3) = 1; Ai(q8, h6) = 1; Ai(q9, h9) = 1; 
        Ai(q9, q9) = -1; Ai(c, q9) = -1; Ai(q3, q7) = -1; Ai(q6, q8) = -1; Ai(q7,s)=-1;          
        As.push_back(0.5 * (Ai + Ai.transpose()));                 
                  
                  
        /** Constraints for HR **/
        /* HRt */  
        /*
        Ai.setZero();        
        Ai(c, h1) = 1; Ai(h3, s)=1; 
        Ai(q1,c)=-1;Ai(c,c)=-1;Ai(q3,s)=-1;Ai(s,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //57
        */
        
        /*
        Ai.setZero(); Ai(k, h2) = 1; 
        Ai(q2,k)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   //58
       */
       
        /*
        Ai.setZero(); Ai(c, h3)=1;Ai(c,q3)=-1;Ai(h1,s)=-1;Ai(q1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //59
        */
       
        Ai.setZero(); Ai(c, h4) = 1; Ai(h6, s)=1; 
        Ai(q4,c)=-1;Ai(q6,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));    //60
       
       /*
        Ai.setZero(); Ai(k, h5) = -1;
        Ai(q5,k)=1;Ai(k,k)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //61
      */
      
        Ai.setZero(); Ai(c, h6) = 1;Ai(c,q6)=-1;Ai(h4,s)=-1;Ai(q4,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  

       /*
        Ai.setZero(); Ai(c, h7) = 1; Ai(h9, s)=1; Ai(c,q7)=-1;Ai(q9,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //63
       */
       
       /*
        Ai.setZero(); Ai(k, h8) = 1;Ai(q8,k)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  //64
       */
       /*
        Ai.setZero(); Ai(c, h9) = 1;Ai(q9,c)=-1;Ai(h7,s)=-1;Ai(q7,s)=1;
        Ai(c,c)=-1;Ai(s,s)=-1;         
        As.push_back(0.5 * (Ai + Ai.transpose()));     //65               
         */
                  
        /* RtH */
        /*
        Ai.setZero(); 
        Ai(c,h1)=1;Ai(h7,s)=-1;Ai(q7,s)=1;Ai(q1,c)=-1;
        Ai(c,c)=-1;Ai(s,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //66
      */
      
        Ai.setZero(); Ai(c,h2)=1;Ai(h8,s)=-1; 
        Ai(q2,c)=-1;Ai(q8,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
       
        Ai.setZero(); Ai(c,h3)=1;Ai(h9,s)=-1;
        Ai(q3,c)=-1;Ai(q9,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 

       /*
        Ai.setZero(); Ai(k,h4)=1;Ai(q4,k)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //69
       */
       
        Ai.setZero(); Ai(k,h5)=-1;Ai(k,k)=1;Ai(q5,k)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //70
       
        /*
        Ai.setZero(); Ai(k,h6)=1;Ai(q6,k)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //71
        */        
        
        /*       
        Ai.setZero(); Ai(c,h7)=1;Ai(q7,c)=-1;Ai(h1,s)=1;Ai(q1,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //72
        */
      
        Ai.setZero(); Ai(c,h8)=1;Ai(q8,c)=-1;Ai(h2,s)=1;Ai(q2,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
       
       /*
        Ai.setZero(); Ai(c,h9)=1;Ai(q9,c)=-1;Ai(q3,s)=-1;Ai(h3,s)=1;
        Ai(c,c)=-1;Ai(s,s)=-1;        
        As.push_back(0.5 * (Ai + Ai.transpose()));    //74
                  */
                  
        /* HR */
        Ai.setZero(); 
        Ai(c,h1)=1;Ai(h3,s)=-1;
        Ai(c,c)=-1;Ai(s,s)=1;Ai(q1,c)=-1;Ai(q3,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
       
       /*
        Ai.setZero();  Ai(k,h2)=1;Ai(k,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //76
       */
       
       /*
        Ai.setZero();  Ai(c,h3)=1;Ai(c,q3)=-1;Ai(c,s)=-2;Ai(h1,s)=1;Ai(q1,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));    // 77
        */
       
        Ai.setZero();  Ai(c,h4)=1;Ai(h6,s)=-1;
        Ai(c,q4)=-1;Ai(q6,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
     
        /*
        Ai.setZero();  Ai(k,h5)=-1;Ai(q5,k)=1;Ai(k,k)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   // 79
        */
       
        Ai.setZero();  Ai(c,h6)=1;Ai(c,q6)=-1;Ai(h4,s)=1;Ai(q4,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); //80
       
        /*
        Ai.setZero();  Ai(c,h7)=1;Ai(h9,s)=-1;Ai(q7,c)=-1;Ai(q9,s)=1;Ai(c,s)=2;
        As.push_back(0.5 * (Ai + Ai.transpose()));   // 81
        */
       
        
        Ai.setZero();  Ai(k,h8)=1;Ai(k,q8)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); // 82
        
      
        
        Ai.setZero();  Ai(c,h9)=1;Ai(c,q9)=-1;Ai(h7,s)=1;Ai(q7,s)=-1;
        Ai(c,c)=-1;Ai(s,s)=1;    
        As.push_back(0.5 * (Ai + Ai.transpose()));      // 83
        
                  
        /* RH */
        
        /*
        Ai.setZero(); 
        Ai(c,h1)=1;Ai(h7,s)=1;
        Ai(c,c)=-1;Ai(s,s)=1;Ai(q1,c)=-1;Ai(q7,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   // 84
        */
        
        
        Ai.setZero(); Ai(c,h2)=1;Ai(h8,s)=1;Ai(q2,c)=-1;Ai(q8,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); //85
        
        Ai.setZero(); Ai(c,h3)=1;Ai(h9,s)=1;
        Ai(q3,c)=-1;Ai(q9,s)=-1;Ai(c,s)=-2;
        As.push_back(0.5 * (Ai + Ai.transpose())); 

        
        Ai.setZero(); Ai(k,h4)=1;Ai(k,q4)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
       
        
        /*
        Ai.setZero(); Ai(k,h5)=-1;Ai(k,k)=1;Ai(q5,k)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //88
        */
        
        
        Ai.setZero(); Ai(k,h6)=1;Ai(q6,k)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 89
        
        
        Ai.setZero(); Ai(c,h7)=1;Ai(q7,c)=-1;Ai(c,s)=2;Ai(h1,s)=-1;Ai(q1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); //90
        
        Ai.setZero(); Ai(c,h8)=1;Ai(q8,c)=-1;Ai(h2,s)=-1;Ai(q2,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
     
        Ai.setZero(); Ai(c,h9)=1;Ai(c,q9)=-1;Ai(h3,s)=-1;Ai(q3,s)=1;
        Ai(c,c)=-1;Ai(s,s)=1;          
        As.push_back(0.5 * (Ai + Ai.transpose())); 
                
          
        /** Constraints with HH **/
        /* HHt */
        // h1^2 - h1*q1 - h2*q2 - h3*q3 - h3*s - c*h1 + h2^2 + h3^2
        Ai.setZero(); 
        Ai(h1,h1)=1; Ai(h1,q1)=-1; Ai(h2,q2)=-1; Ai(h3,q3)=-1;
        Ai(h3,s)=-1; Ai(c,h1)=-1; Ai(h2,h2)=1; Ai(h3,h3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   
        
        // h1*h4 + h2*h5 + h3*h6 - h2*k - h1*q4 - h2*q5 - h3*q6
        Ai.setZero(); 
        Ai(h1,h4)=1; Ai(h2,h5)=1; Ai(h3,h6)=1; Ai(h2,k)=-1;
        Ai(h1,q4)=-1; Ai(h2,q5)=-1; Ai(h3,q6)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   
        
        // h1*h7 - c*h3 + h2*h8 + h3*h9 - h1*q7 - h2*q8 - h3*q9 + h1*s
        Ai.setZero(); 
        Ai(h1,h7)=1; Ai(c,h3)=-1; Ai(h2,h8)=1; Ai(h3,h9)=1;
        Ai(h1,q7)=-1; Ai(h2,q8)=-1; Ai(h3,q9)=-1; Ai(h1,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));       //95      
         
        // h1*h4 - c*h4 + h2*h5 + h3*h6 - h4*q1 - h5*q2 - h6*q3 - h6*s
        Ai.setZero(); 
        Ai(h1,h4)=1; Ai(c,h4)=-1; Ai(h2,h5)=1; Ai(h3,q6)=1;
        Ai(h4,q1)=-1; Ai(h5,q2)=-1; Ai(h6,q3)=-1; Ai(h6,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));              
                  
        // h4^2 - h4*q4 - h5*q5 - h6*q6 - h5*k + h5^2 + h6^2
        Ai.setZero(); 
        Ai(h4,h4)=1; Ai(q4,h4)=-1; Ai(q5,h5)=-1; Ai(h6,q6)=-1;
        Ai(h5,k)=-1; Ai(h5,h5)=1; Ai(h6,h6)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));          
              
        // h4*h7 - c*h6 + h5*h8 + h6*h9 - h4*q7 - h5*q8 - h6*q9 + h4*s
        Ai.setZero(); 
        Ai(h4,h7)=1; Ai(c,h6)=-1; Ai(q8,h5)=1; Ai(h6,h9)=1;
        Ai(h4,q7)=-1; Ai(h5,q8)=-1; Ai(h6,q9)=-1; Ai(h4,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));            
         
        // h1*h7 - c*h7 + h2*h8 + h3*h9 - h7*q1 - h8*q2 - h9*q3 - h9*s
        Ai.setZero(); 
        Ai(h1,h7)=1; Ai(c,h7)=-1; Ai(h8,h2)=1; Ai(h3,h9)=1;
        Ai(h7,q1)=-1; Ai(h8,q2)=-1; Ai(h9,q3)=-1; Ai(h9,s)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));           
          
        // h4*h7 + h5*h8 + h6*h9 - h8*k - h7*q4 - h8*q5 - h9*q6
        Ai.setZero(); 
        Ai(h4,h7)=1; Ai(h5,h8)=1; Ai(h9,h6)=1; Ai(h8,k)=-1;
        Ai(h7,q4)=-1; Ai(h8,q5)=-1; Ai(h9,q6)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //100    
           
        // h7*s - h7*q7 - h8*q8 - h9*q9 - c*h9 + h7^2 + h8^2 + h9^2
 
        Ai.setZero(); 
        Ai(s, h7)=1; Ai(q7,h7)=-1; Ai(h8,q8)=-1; Ai(q9,h9)=-1;
        Ai(h9,c)=-1; Ai(h7,h7)=1; Ai(h8,h8)=1; Ai(h9,h9)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));    //101
 
        
                
        /* HtH */
        // h7*s - h1*q1 - h4*q4 - h7*q7 - c*h1 + h1^2 + h4^2 + h7^2
        /*
        Ai.setZero(); 
        Ai(h7,s)=1; Ai(h1,q1)=-1; Ai(h4,q4)=-1; Ai(h7,q7)=-1;
        Ai(c,h1)=-1; Ai(h1,h1)=1; Ai(h4,h4)=1; Ai(h7,h7)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); //102
        */
        
        // h1*h2 - c*h2 + h4*h5 + h7*h8 - h2*q1 - h5*q4 - h8*q7 + h8*s
        Ai.setZero(); 
        Ai(h1,h2)=1; Ai(c,h2)=-1; Ai(h4,h5)=1; Ai(h7,h8)=1;
        Ai(q1,h2)=-1; Ai(h5,q4)=-1; Ai(h8,q7)=-1; Ai(h8,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h1*h3 - c*h3 + h4*h6 + h7*h9 - h3*q1 - h6*q4 - h9*q7 + h9*s
        /*
        Ai.setZero(); 
        Ai(h1,h3)=1; Ai(c,h3)=-1; Ai(h4,h6)=1; Ai(h7,h9)=1;
        Ai(q1,h3)=-1; Ai(h6,q4)=-1; Ai(h9,q7)=-1; Ai(h9,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   //104
        */
        
        // h1*h2 + h4*h5 + h7*h8 - h4*k - h1*q2 - h4*q5 - h7*q8
        Ai.setZero(); 
        Ai(h1,h2)=1; Ai(h4,h5)=1; Ai(h7,h8)=1; Ai(h4,k)=-1;
        Ai(q2,h1)=-1; Ai(h4,q5)=-1; Ai(h7,q8)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); //105
        
        // h2^2 - h2*q2 - h5*q5 - h8*q8 - h5*k + h5^2 + h8^2
        Ai.setZero(); 
        Ai(h2,h2)=1; Ai(h2,q2)=-1; Ai(h5,q5)=-1; Ai(h8,q8)=-1;
        Ai(k,h5)=-1; Ai(h5,h5)=1; Ai(h8,h8)=1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h2*h3 + h5*h6 + h8*h9 - h6*k - h3*q2 - h6*q5 - h9*q8
        Ai.setZero(); 
        Ai(h2,h3)=1; Ai(h5,h6)=1; Ai(h8,h9)=1; Ai(h6,k)=-1;
        Ai(q2,h3)=-1; Ai(h6,q5)=-1; Ai(h9,q8)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h1*h3 - c*h7 + h4*h6 + h7*h9 - h1*q3 - h4*q6 - h7*q9 - h1*s
        Ai.setZero(); 
        Ai(h1,h3)=1; Ai(c,h7)=-1; Ai(h4,h6)=1; Ai(h1,s)=-1;
        Ai(h7,h9)=1; Ai(h1,q3)=-1; Ai(h4,q6)=-1; Ai(h7,q9)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h2*h3 - c*h8 + h5*h6 + h8*h9 - h2*q3 - h5*q6 - h8*q9 - h2*s
        Ai.setZero(); 
        Ai(h2,h3)=1; Ai(c,h8)=-1; Ai(h5,h6)=1; Ai(h8,h9)=1;
        Ai(h2,q3)=-1; Ai(h5,q6)=-1; Ai(h8,q9)=-1; Ai(h2,s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose()));   
        
        // h3^2 - h3*q3 - h6*q6 - h9*q9 - h3*s - c*h9 + h6^2 + h9^2
        Ai.setZero(); 
        Ai(h3,h3)=1; Ai(h3,q3)=-1; Ai(h6,q6)=-1; Ai(h9,q9)=-1;
        Ai(h3,s)=-1; Ai(h9,c)=-1; Ai(h6,h6)=1; Ai(h9,h9)=1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); //110
        
        
        /* HH */
        // h2*h4 - c*h1 + h3*h7 - h1*q1 - h4*q2 - h7*q3 - h7*s + h1^2
        Ai.setZero(); 
        Ai(h2,h4)=1; Ai(c,h1)=-1; Ai(h3,h7)=1; Ai(h1,q1)=-1;
        Ai(h4,q2)=-1; Ai(h7,q3)=-1; Ai(h7,s)=-1; Ai(h1,h1)=1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h1*h2 - c*h2 + h2*h5 + h3*h8 - h2*q1 - h5*q2 - h8*q3 - h8*s
        Ai.setZero(); 
        Ai(h1,h2)=1; Ai(c,h2)=-1; Ai(h2,h5)=1; Ai(h3,h8)=1;
        Ai(h2,q1)=-1; Ai(h5,q2)=-1; Ai(h8,q3)=-1; Ai(h8,s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h1*h3 - c*h3 + h2*h6 + h3*h9 - h3*q1 - h6*q2 - h9*q3 - h9*s
        Ai.setZero(); 
        Ai(h1,h3)=1; Ai(c,h3)=-1; Ai(h2,h6)=1; Ai(h3,h9)=1;
        Ai(h3,q1)=-1; Ai(h6,q2)=-1; Ai(h9,q3)=-1; Ai(h9,s)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h1*h4 + h4*h5 + h6*h7 - h4*k - h1*q4 - h4*q5 - h7*q6
        Ai.setZero(); 
        Ai(h1,h4)=1; Ai(h4,h5)=1; Ai(h7,h6)=1; Ai(h4,k)=-1;
        Ai(h1,q4)=-1; Ai(h4,q5)=-1; Ai(h7,q6)=-1; 
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h2*h4 + h6*h8 - h5*k - h2*q4 - h5*q5 - h8*q6 + h5^2
        Ai.setZero(); 
        Ai(h2,h4)=1; Ai(h6,h8)=1; Ai(h5,k)=-1;
        Ai(h2,q4)=-1; Ai(h5,q5)=-1; Ai(h8,q6)=-1; Ai(h5,h5)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); //115
        
        // h3*h4 + h5*h6 + h6*h9 - h6*k - h3*q4 - h6*q5 - h9*q6
        Ai.setZero(); 
        Ai(h3,h4)=1; Ai(h6,h5)=1; Ai(h6,h9)=1;
        Ai(h6,k)=-1; Ai(h3,q4)=-1; Ai(h6,q5)=-1; Ai(h9,q6)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
        // h1*h7 - c*h7 + h4*h8 + h7*h9 - h1*q7 - h4*q8 - h7*q9 + h1*s
        Ai.setZero(); 
        Ai(h1,h7)=1; Ai(c,h7)=-1; Ai(h4,h8)=1; Ai(h1,s)=1;
        Ai(h7,h9)=1; Ai(h1,q7)=-1; Ai(h4,q8)=-1; Ai(h7,q9)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h2*h7 - c*h8 + h5*h8 + h8*h9 - h2*q7 - h5*q8 - h8*q9 + h2*s
        Ai.setZero(); 
        Ai(h2,h7)=1; Ai(c,h8)=-1; Ai(h5,h8)=1; Ai(h8,h9)=1;
        Ai(q7,h2)=-1; Ai(h5,q8)=-1; Ai(h8,q9)=-1; Ai(h2,s)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h3*h7 - c*h9 + h6*h8 - h3*q7 - h6*q8 - h9*q9 + h3*s + h9^2
        Ai.setZero(); 
        Ai(h3,h7)=1; Ai(c,h9)=-1; Ai(h6,h8)=1; Ai(h3,q7)=-1;
        Ai(h6,q8)=-1; Ai(h9,q9)=-1; Ai(h3,s)=1; Ai(h9,h9)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
        
        
        /** Constraints for rotation **/
        // K = H - Q
        
        /* trace(K*K') - 3 * k*k */
        //  h1^2 - 2*h1*q1 + h2^2 - 2*h2*q2 + h3^2 - 2*h3*q3 
        // + h4^2 - 2*h4*q4 + h5^2 - 2*h5*q5 + h6^2 
        // - 2*h6*q6 + h7^2 - 2*h7*q7 + h8^2 - 2*h8*q8 + h9^2 
        // - 2*h9*q9 - 3*k^2 + q1^2 + q2^2 + q3^2 + q4^2 
        // + q5^2 + q6^2 + q7^2 + q8^2 + q9^2
        Ai.setZero(); 
        Ai(h1,h1)=1; Ai(h1,q1)=-2; Ai(h2,h2)=1; Ai(h2,q2)=-2;
        Ai(h3,h3)=1; Ai(h3,q3)=-2; Ai(h4,h4)=1; Ai(h4,q4)=-2; 
        Ai(h5,h5)=1; Ai(h5,q5)=-2; Ai(h6,h6)=1; Ai(h6,q6)=-2; 
        Ai(h7,h7)=1; Ai(h7,q7)=-2; Ai(h8,h8)=1; Ai(h8,q8)=-2;  
        Ai(h9,h9)=1; Ai(h9,q9)=-2; Ai(k,k)=-3; Ai(q1,q1)=1; 
        Ai(q2,q2)=1; Ai(q3,q3)=1; Ai(q4,q4)=1; Ai(q5,q5)=1; 
        Ai(q6,q6)=1; Ai(q7,q7)=1; Ai(q8,q8)=1; Ai(q9,q9)=1;         
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 120
        
        /* K*K' - k*k * eye(3) */
        // h1^2 - 2*h1*q1 + h2^2 - 2*h2*q2 + h3^2 - 2*h3*q3 - k^2 + q1^2 + q2^2 + q3^2
        Ai.setZero(); 
        Ai(h1,h1)=1; Ai(h1,q1)=-2; Ai(h2,h2)=1; Ai(h2,q2)=-2;
        Ai(h3,h3)=1; Ai(h3,q3)=-2; Ai(k,k)=-1; Ai(q1,q1)=1;
        Ai(q2,q2)=1; Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));  
         
        // h1*h4 + h2*h5 + h3*h6 - h1*q4 - h4*q1 - h2*q5 
        // - h5*q2 - h3*q6 - h6*q3 + q1*q4 + q2*q5 + q3*q6
        Ai.setZero(); 
        Ai(h1,h4)=1; Ai(h2,h5)=1; Ai(h3,h6)=1; Ai(h1,q4)=-1;
        Ai(h4,q1)=-1; Ai(h2,q5)=-1; Ai(h5,q2)=-1; Ai(h3,q6)=-1;
        Ai(h6,q3)=-1; Ai(q1,q4)=1; Ai(q2,q5)=1; Ai(q3,q6)=1;
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h1*h7 + h2*h8 + h3*h9 - h1*q7 - h7*q1 - h2*q8 
        // - h8*q2 - h3*q9 - h9*q3 + q1*q7 + q2*q8 + q3*q9
        Ai.setZero(); 
        Ai(h1,h7)=1; Ai(h2,h8)=1; Ai(h3,h9)=1; Ai(h1,q7)=-1;
        Ai(h7,q1)=-1; Ai(h2,q8)=-1; Ai(h8,q2)=-1; Ai(h3,q9)=-1;
        Ai(h9,q3)=-1; Ai(q1,q7)=1; Ai(q2,q8)=1; Ai(q3,q9)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));   
        
        // h4^2 - 2*h4*q4 + h5^2 - 2*h5*q5 + h6^2 - 2*h6*q6 - k^2 + q4^2 + q5^2 + q6^2
        Ai.setZero(); 
        Ai(h4,h4)=1; Ai(q4,h4)=-2; Ai(h5,h5)=1; Ai(h5,q5)=-2;
        Ai(h6,h6)=1; Ai(h6,q6)=-2; Ai(k,k)=-1; Ai(q4,q4)=1;
        Ai(q5,q5)=1; Ai(q6,q6)=1;        
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        
        // h4*h7 + h5*h8 + h6*h9 - h4*q7 - h7*q4 - h5*q8 
        // - h8*q5 - h6*q9 - h9*q6 + q4*q7 + q5*q8 + q6*q9        
        Ai.setZero(); 
        Ai(h7,h4)=1; Ai(h5,h8)=1; Ai(h6,h9)=1; Ai(h4,q7)=-1;
        Ai(h7,q4)=-1; Ai(h5,q8)=-1; Ai(h8,q5)=-1; Ai(h6,q9)=-1;
        Ai(h9,q6)=-1; Ai(q4,q7)=1; Ai(q5,q8)=1; Ai(q6,q9) = 1;        
        As.push_back(0.5 * (Ai + Ai.transpose()));   //125
        
        
        // h7^2 - 2*h7*q7 + h8^2 - 2*h8*q8 + h9^2 - 2*h9*q9 - k^2 + q7^2 + q8^2 + q9^2
        /*
        Ai.setZero(); 
        Ai(h7,h7)=1; Ai(h7,q7)=-2; Ai(h8,h8)=1; Ai(h8,q8)=-2;
        Ai(h9,h9)=1; Ai(h9,q9)=-2; Ai(k,k)=-1; Ai(q7,q7)=1;
        Ai(q8,q8)=1; Ai(q9,q9)=1;        
        As.push_back(0.5 * (Ai + Ai.transpose())); // 126
        */
        
        /** KtK **/
        /* K'*K - k*k*eye(3) */
        // h1^2 - 2*h1*q1 + h4^2 - 2*h4*q4 + h7^2 - 2*h7*q7 - k^2 + q1^2 + q4^2 + q7^2

        Ai.setZero(); 
        Ai(h1,h1)=1; Ai(h1,q1)=-2; Ai(h4,h4)=1; Ai(h4,q4)=-2;
        Ai(h7,h7)=1; Ai(h7,q7)=-2; Ai(k,k)=-1; Ai(q1,q1)=1;
        Ai(q4,q4)=1; Ai(q7,q7)=1;        
        As.push_back(0.5 * (Ai + Ai.transpose())); //127
    
        
        // h1*h2 + h4*h5 + h7*h8 - h1*q2 - h2*q1 - h4*q5 - h5*q4 - h7*q8 - h8*q7 + q1*q2 + q4*q5 + q7*q8
        Ai.setZero(); 
        Ai(h1,h2)=1; Ai(h4,h5)=1; Ai(h7,h8)=1; Ai(h1,q2)=-1;
        Ai(h2,q1)=-1; Ai(h4,q5)=-1; Ai(h5,q4)=-1; Ai(h7,q8)=-1; 
        Ai(h8,q7)=-1; Ai(q1,q2)=1;Ai(q4,q5)=1;Ai(q7,q8)=1;       
        As.push_back(0.5 * (Ai + Ai.transpose()));   
        
        // h1*h3 + h4*h6 + h7*h9 - h1*q3 - h3*q1 - h4*q6 - h6*q4 - h7*q9 - h9*q7 + q1*q3 + q4*q6 + q7*q9
        Ai.setZero(); 
        Ai(h1,h3)=1; Ai(h4,h6)=1; Ai(h7,h9)=1; Ai(h1,q3)=-1;
        Ai(h3,q1)=-1; Ai(h4,q6)=-1; Ai(h6,q4)=-1; Ai(h7,q9)=-1; 
        Ai(h9,q7)=-1; Ai(q1,q3)=1;Ai(q4,q6)=1;Ai(q7,q9)=1;       
        As.push_back(0.5 * (Ai + Ai.transpose())); 
        
        // h2^2 - 2*h2*q2 + h5^2 - 2*h5*q5 + h8^2 - 2*h8*q8 - k^2 + q2^2 + q5^2 + q8^2
        /*
        Ai.setZero(); 
        Ai(h2,h2)=1; Ai(h2,q2)=-2; Ai(h5,h5)=1; Ai(h5,q5)=-2;
        Ai(h8,h8)=1; Ai(h8,q8)=-2; Ai(k,k)=-1; Ai(q2,q2)=1;
        Ai(q5,q5)=1; Ai(q8,q8)=1;        
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 130
        */
        
        // h2*h3 + h5*h6 + h8*h9 - h2*q3 - h3*q2 - h5*q6 - h6*q5 - h8*q9 - h9*q8 + q2*q3 + q5*q6 + q8*q9
        Ai.setZero(); 
        Ai(h2,h3)=1; Ai(h5,h6)=1; Ai(h8,h9)=1; Ai(h2,q3)=-1;
        Ai(h3,q2)=-1; Ai(h5,q6)=-1; Ai(h6,q5)=-1; Ai(h8,q9)=-1; 
        Ai(h9,q8)=-1; Ai(q2,q3)=1;Ai(q5,q6)=1;Ai(q8,q9)=1;       
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        // h3^2 - 2*h3*q3 + h6^2 - 2*h6*q6 + h9^2 - 2*h9*q9 - k^2 + q3^2 + q6^2 + q9^2
        Ai.setZero(); 
        Ai(h3,h3)=1; Ai(h3,q3)=-2; Ai(h6,h6)=1; Ai(h6,q6)=-2;
        Ai(h9,h9)=1; Ai(h9,q9)=-2; Ai(k,k)=-1; Ai(q3,q3)=1;
        Ai(q6,q6)=1; Ai(q9,q9)=1;        
        As.push_back(0.5 * (Ai + Ai.transpose()));  // 132 (index from 0)
        
        /* Redundant 
         21
         23
         25
         26
         27
         29
         39
         40
         42
         43
         44
         46
         47         
         52
         57
         58 
         59
         61
         63
         64
         65
         66
         69              
         71
         72
         74
         76
         77
         79
         81
         84
         88         
         102
         104
         126
         130
        */
        
        /*
         21
         23
         24
         25
         26
         27
         29         
         39
         40
         44
         46
         47
         56
         57
         59
         61
         63
         64
         65
         66
         69
         70
         71
         72
         74
         76
         77
         79
         81
         84
         87
         88         
         102
         104
         126
         127              
        */
                         
    } 
    
}  // end of namespace Rot Prior H
