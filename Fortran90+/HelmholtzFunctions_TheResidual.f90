Module HelmholtzFunctions_TheResidual
   
   Use Toolbox
   Use HFEConstant, Only: c1,d1,t,n,alpha,beta,betaInv,gamma,epsilon,a1,b1,c2,d2,a2,b2

   Implicit None


   Private GetTheta,GetDel


Contains


   Elemental Function Driver_HelmholtzResidual_d(delta,tau) Result(Helm_d)

      Real(Kind=ExDouble), Intent(in)  :: delta,tau
      Real(Kind=ExDouble)              :: Helm_d

      Integer :: k,m,p
      Real(Kind=ExDouble) :: Part,deltaMod,Theta,Del,Psi,Psi_d,Delbi_d

      !Initialize local variables
      Helm_d = zero;
      deltaMod = (delta-one)**itwo + SmallNumber


      Do k = 1,7
         Part = n(k) * d1(k) * delta**(d1(k)-one) * tau**(t(k))
         Helm_d = Helm_d + Part
      End Do


      Do k =  8,51
         Part = d1(k) - c1(k) * delta**(c1(k))
         Part = delta**(d1(k)-one) * tau**(t(k)) * Part
         Part = n(k) * Exp(-delta**(c1(k))) * Part
         Helm_d = Helm_d + Part
      End Do



      Do k =  52, 54

         m = k - 51

         Part = -alpha(m)*(delta-Epsilon(m))**itwo - beta(m)*(tau-gamma(m))**itwo
         Part = n(k) * delta**(d1(k)) * tau**(t(k)) * Exp(Part)
         Part = Part * ( d1(k)/delta - two*alpha(m)*(delta-epsilon(m)))
         Helm_d = Helm_d + Part

      End Do

      Do k = 55,56
         m = k-51
         p = k-54


         Theta   = GetTheta(deltaMod, tau  , a2(p) , betaInv(m))
         Del     = GetDel  (deltaMod, Theta, b1(p) , a2(p))
         Psi     = GetPsi  (deltaMod,tau,c2(p),d2(p))
         Psi_d   = GetPsi_d(delta,Psi,c2(p))
         Delbi_d = GetDelbi_d(delta,deltaMod,Del,Theta,a1(p),b1(p),a2(p),b2(p),betaInv(m))

         Part = Del**(b2(p)) * (Psi + delta * Psi_d)
         Part = n(k)*(Part + Delbi_d * delta * Psi)
         Helm_d = Helm_d + Part
      End Do

   End Function Driver_HelmholtzResidual_d







   Elemental Function HelmholtzResidual_dWithKahanSum(delta,tau) Result(Helm_d)

      Real(Kind=ExDouble), Intent(in)  :: delta,tau
      Real(Kind=ExDouble)              :: Helm_d

      Integer :: k,m,p
      Real(Kind=ExDouble) :: Part,deltaMod,Theta,Del,Psi,Psi_d,Delbi_d
      Real(Kind=ExDouble) :: SumErrorIn,SumErrorOut,Helm_dWork

      !Initialize local variables
      Helm_d      = zero
      Helm_dWork  = zero
      SumErrorIn  = zero
      SumErrorOut = zero
      deltaMod    = (delta-one)**itwo + SmallNumber

      Do k = 1,7
         Part = n(k) * d1(k) * delta**(d1(k)-one) * tau**(t(k))

         Call KahanSum(Helm_d,Part,SumErrorIn,Helm_dWork,SumErrorOut)
         SumErrorIn = SumErrorOut
         Helm_d     = Helm_dWork
      End Do
      
      Do k =  8,51
         Part = d1(k) - c1(k) * delta**(c1(k))
         Part = delta**(d1(k)-one) * tau**(t(k)) * Part
         Part = n(k) * Exp(-delta**(c1(k))) * Part

         Call KahanSum(Helm_d,Part,SumErrorIn,Helm_dWork,SumErrorOut)
         SumErrorIn = SumErrorOut
         Helm_d     = Helm_dWork
      End Do


      
      Do k =  52, 54
         m = k - 51

         Part = -alpha(m)*(delta-Epsilon(m))**itwo - beta(m)*(tau-gamma(m))**itwo
         Part = n(k) * delta**(d1(k)) * tau**(t(k)) * Exp(Part)
         Part = Part * ( d1(k)/delta - two*alpha(m)*(delta-epsilon(m)))

         Call KahanSum(Helm_d,Part,SumErrorIn,Helm_dWork,SumErrorOut)
         SumErrorIn = SumErrorOut
         Helm_d     = Helm_dWork
      End Do

      Do k = 55,56
         m = k-51
         p = k-54


         Theta   = GetTheta(deltaMod, tau  , a2(p) , 1/beta(m))
         Del     = GetDel  (deltaMod, Theta, b1(p) , a2(p))
         Psi     = GetPsi  (deltaMod,tau,c2(p),d2(p))
         Psi_d   = GetPsi_d(delta,Psi,c2(p))
         Delbi_d = GetDelbi_d(delta,deltaMod,Del,Theta,a1(p),b1(p),a2(p),b2(p),1/beta(m))

         Part = Del**(b2(p)) * (Psi + delta * Psi_d)
         Part = n(k)*(Part + Delbi_d * delta * Psi)

         Call KahanSum(Helm_d,Part,SumErrorIn,Helm_dWork,SumErrorOut)
         SumErrorIn = SumErrorOut
         Helm_d     = Helm_dWork 
     End Do

  End Function HelmholtzResidual_dWithKahanSum






   ! =================================================================================== !
   !                               Private Helper Functions                              !
   ! =================================================================================== !

   Elemental Function GetTheta(deltaMod,tau,a1elem,betaInv) Result(Theta)
      Real(Kind=ExDouble),Intent(in) :: deltaMod,tau,a1elem,betaInv
      Real(Kind=ExDouble)            :: Theta

      Theta = (one-tau) + a1elem * deltaMod**(half*betainv)

   End Function GetTheta



   Elemental Function GetDel(deltaMod,Theta,b1elem,a2elem) Result(Del)
      Real(Kind=ExDouble),Intent(in) :: deltaMod,Theta,b1elem,a2elem
      Real(Kind=ExDouble)            :: Del

      Del = Theta**2 + b1elem * deltaMod**a2elem

   End Function GetDel



   Elemental Function GetPsi(deltaMod,tau,c2elem,d2elem) Result(Psi)
      Real(Kind=ExDouble),Intent(in) :: deltaMod,tau,c2elem,d2elem
      Real(Kind=ExDouble)            :: Psi

      Psi = Exp(-c2elem*deltaMod - d2elem*(tau-one)**itwo)

   End Function GetPsi



   Elemental Function GetPsi_d(delta,Psi,c2elem) Result(Psi_d)
      Real(Kind=ExDouble),Intent(in) :: delta,Psi,c2elem
      Real(Kind=ExDouble)            :: Psi_d

      Psi_d = -two*c2elem *(delta-one) * Psi

   End Function GetPsi_d



   Elemental Function GetDelbi_d(delta,deltaMod,Del,Theta,a1elem,b1elem, & 
                                   a2elem,b2elem,betaInv) Result(Delbi_d)
      
      Real(Kind=ExDouble), Intent(in) :: delta,deltaMod,Del,Theta,a1elem,b1elem
      Real(Kind=ExDouble), Intent(in) :: a2elem,b2elem,betaInv
      Real(Kind=ExDouble)             :: Delbi_d

      Real(Kind=ExDouble) :: Del_d

      Del_d   = GetDel_d(delta,deltaMod,Theta,a1elem,b1elem,a2elem,betaInv)
      Delbi_d = b2elem * Del**(b2elem-one) * Del_d

   End Function GetDelbi_d



   Elemental Function GetDel_d(delta,deltaMod,Theta,a1elem,b1elem,a2elem, & 
                               betaInv) Result(Del_d)

      Real(Kind=ExDouble), Intent(in) :: delta,deltaMod,Theta,a1elem,b1elem
      Real(Kind=ExDouble), Intent(in) :: a2elem,betaInv
      Real(Kind=ExDouble)             :: Del_d

      Real(Kind=ExDouble) :: Part

      Part  = two*betaInv*a1elem*Theta*deltaMod**(half*betaInv-one)
      Part  = Part + two*b1elem*a2elem*deltaMod**(a2elem-one)
      Del_d = (delta-one)*Part
      
   End Function GetDel_d

End Module HelmholtzFunctions_TheResidual
