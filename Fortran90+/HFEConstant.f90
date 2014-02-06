Module HFEConstant

   ! =================================================================================== ! 
   !                                        Module Loads                                 !
   ! =================================================================================== !
!   Use Iso_Fortran_Env, Only: Real128
   Use Toolbox, Only: ExDouble


   Implicit None

!   Integer,Parameter :: Double = 8, ExDouble = 10


   ! =================================================================================== !
   !                     Helmholtz Free Energy Ideal Gas Constants                       !
   ! =================================================================================== !

   !   These commented coefficients are the exact values from the IAPWS-95 
   !   reference.  The coefficients actually used have been adjusted to 
   !   yield exact 0s for the saturated liquid internal energy and entropy at
   !   the triple point temperature.
   !
   !     n_o_Exact(1) = -8.3204464837497;
   !     n_o_Exact(2) =  6.6832105275932;
   !
   
   ! Ideal Gas Coefficient Array 1
   Real(Kind=ExDouble), Dimension(8), Parameter :: &
         n_o = (/ -8.320446483749712_ExDouble , & 
                   6.683210527593234_ExDouble , &
                   3.00632_ExDouble           , &
                   0.012436_ExDouble          , &
                   0.97315_ExDouble           , &
                   1.27950_ExDouble           , &
                   0.96956_ExDouble           , &
                   0.24873_ExDouble            /) 


   ! Ideal Gas Coefficient Array 2
   Real(Kind=ExDouble), Dimension(8),Parameter :: &
         gamma_o = (/  0.0_ExDouble        , & 
                       0.0_ExDouble        , & 
                       0.0_ExDouble        , & 
                       1.28728967_ExDouble , & 
                       3.53734222_ExDouble , & 
                       7.74073708_ExDouble , &
                       9.24437796_ExDouble , & 
                      27.50751050_ExDouble  /) 




   ! =================================================================================== !
   !                     Helmholtz Free Energy Residual Constants                        !
   ! =================================================================================== !

   ! --------------------------------------------------------------- !
   !                             Group 1                             !
   ! --------------------------------------------------------------- !

   ! c1:Group1:Array1
    Real(Kind=ExDouble),Dimension(51), Parameter :: &
       c1 = (/ 0.0_ExDouble,&
               0.0_ExDouble,&
               0.0_ExDouble,&
               0.0_ExDouble,&
               0.0_ExDouble,&
               0.0_ExDouble,&
               0.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               1.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               2.0_ExDouble,&
               3.0_ExDouble,&
               3.0_ExDouble,&
               3.0_ExDouble,&
               3.0_ExDouble,&
               4.0_ExDouble,&
               6.0_ExDouble,&
               6.0_ExDouble,&
               6.0_ExDouble,&
               6.0_ExDouble  /)

 
    ! d:Group1:Array2
    Real(Kind=ExDouble),Dimension(54),Parameter ::   &
       d1 = (/  1.0_ExDouble,&
                1.0_ExDouble,&
                1.0_ExDouble,&
                2.0_ExDouble,&
                2.0_ExDouble,&
                3.0_ExDouble,&
                4.0_ExDouble,&
                1.0_ExDouble,&
                1.0_ExDouble,&
                1.0_ExDouble,&
                2.0_ExDouble,&
                2.0_ExDouble,&
                3.0_ExDouble,&
                4.0_ExDouble,&
                4.0_ExDouble,&
                5.0_ExDouble,&
                7.0_ExDouble,&
                9.0_ExDouble,&
               10.0_ExDouble,&
               11.0_ExDouble,&
               13.0_ExDouble,&
               15.0_ExDouble,&
                1.0_ExDouble,&
                2.0_ExDouble,&
                2.0_ExDouble,&
                2.0_ExDouble,&
                3.0_ExDouble,&
                4.0_ExDouble,&
                4.0_ExDouble,&
                4.0_ExDouble,&
                5.0_ExDouble,&
                6.0_ExDouble,&
                6.0_ExDouble,&
                7.0_ExDouble,&
                9.0_ExDouble,&
                9.0_ExDouble,&
                9.0_ExDouble,&
                9.0_ExDouble,&
                9.0_ExDouble,&
               10.0_ExDouble,&
               10.0_ExDouble,&
               12.0_ExDouble,&
                3.0_ExDouble,&
                4.0_ExDouble,&
                4.0_ExDouble,&
                5.0_ExDouble,&
               14.0_ExDouble,&
                3.0_ExDouble,&
                6.0_ExDouble,&
                6.0_ExDouble,&
                6.0_ExDouble,&
                3.0_ExDouble,&
                3.0_ExDouble,&
                3.0_ExDouble /)

  
    ! t:Group1:Array3
    Real(Kind=ExDouble),Dimension(54), Parameter ::   &
        t = (/ -0.5_ExDouble,&
              0.875_ExDouble,&
                1.0_ExDouble,&
                0.5_ExDouble,&
               0.75_ExDouble,&
              0.375_ExDouble,&
                1.0_ExDouble,&
                4.0_ExDouble,&
                6.0_ExDouble,&
               12.0_ExDouble,&
                1.0_ExDouble,&
                5.0_ExDouble,&
                4.0_ExDouble,&
                2.0_ExDouble,&
               13.0_ExDouble,&
                9.0_ExDouble,&
                3.0_ExDouble,&
                4.0_ExDouble,&
               11.0_ExDouble,&
                4.0_ExDouble,&
               13.0_ExDouble,&
                1.0_ExDouble,&
                7.0_ExDouble,&
                1.0_ExDouble,&
                9.0_ExDouble,&
               10.0_ExDouble,&
               10.0_ExDouble,&
                3.0_ExDouble,&
                7.0_ExDouble,&
               10.0_ExDouble,&
               10.0_ExDouble,&
                6.0_ExDouble,&
               10.0_ExDouble,&
               10.0_ExDouble,&
                1.0_ExDouble,&
                2.0_ExDouble,&
                3.0_ExDouble,&
                4.0_ExDouble,&
                8.0_ExDouble,&
                6.0_ExDouble,&
                9.0_ExDouble,&
                8.0_ExDouble,&
               16.0_ExDouble,&
               22.0_ExDouble,&
               23.0_ExDouble,&
               23.0_ExDouble,&
               10.0_ExDouble,&
               50.0_ExDouble,&
               44.0_ExDouble,&
               46.0_ExDouble,&
               50.0_ExDouble,&
                0.0_ExDouble,&
                1.0_ExDouble,&
                4.0_ExDouble /)

    ! n:Group1:Array4 
    Real(Kind=ExDouble),Dimension(56), Parameter ::   &
        n = (/  0.12533547935523E-1_ExDouble,&
                0.78957634722828E+1_ExDouble,&
               -0.87803203303561E+1_ExDouble,&
                   0.31802509345418_ExDouble,&
                  -0.26145533859358_ExDouble,&
               -0.78199751687981E-2_ExDouble,&
                0.88089493102134E-2_ExDouble,&
                  -0.66856572307965_ExDouble,&
                   0.20433810950965_ExDouble,&
               -0.66212605039687E-4_ExDouble,&
                  -0.19232721156002_ExDouble,&
                  -0.25709043003438_ExDouble,&
                   0.16074868486251_ExDouble,&
               -0.40092828925807E-1_ExDouble,&
                0.39343422603254E-6_ExDouble,&
               -0.75941377088144E-5_ExDouble,&
                0.56250979351888E-3_ExDouble,&
               -0.15608652257135E-4_ExDouble,&
                0.11537996422951E-8_ExDouble,&
                0.36582165144204E-6_ExDouble,&
              -0.13251180074668E-11_ExDouble,&
               -0.62639586912454E-9_ExDouble,&
                  -0.10793600908932_ExDouble,&
                0.17611491008752E-1_ExDouble,&
                   0.22132295167546_ExDouble,&
                  -0.40247669763528_ExDouble,&
                   0.58083399985759_ExDouble,&
                0.49969146990806E-2_ExDouble,&
               -0.31358700712549E-1_ExDouble,&
                  -0.74315929710341_ExDouble,&
                   0.47807329915480_ExDouble,&
                0.20527940895948E-1_ExDouble,&
                  -0.13636435110343_ExDouble,&
                0.14180634400617E-1_ExDouble,&
                0.83326504880713E-2_ExDouble,&
               -0.29052336009585E-1_ExDouble,&
                0.38615085574206E-1_ExDouble,&
               -0.20393486513704E-1_ExDouble,&
               -0.16554050063734E-2_ExDouble,&
                0.19955571979541E-2_ExDouble,&
                0.15870308324157E-3_ExDouble,&
               -0.16388568342530E-4_ExDouble,&
                0.43613615723811E-1_ExDouble,&
                0.34994005463765E-1_ExDouble,&
               -0.76788197844621E-1_ExDouble,&
                0.22446277332006E-1_ExDouble,&
               -0.62689710414685E-4_ExDouble,&
               -0.55711118565645E-9_ExDouble,&
                  -0.19905718354408_ExDouble,&
                   0.31777497330738_ExDouble,&
                  -0.11841182425981_ExDouble,&
                -0.31306260323435E2_ExDouble,&
                 0.31546140237781E2_ExDouble,&
                -0.25213154341695E4_ExDouble,&
                  -0.14874640856724_ExDouble,&
                   0.31806110878444_ExDouble /)



   ! --------------------------------------------------------------- !
   !                             Group 2                             !
   ! --------------------------------------------------------------- !
 
    ! alpha:Group2:Array1
    Real(Kind=ExDouble), Dimension(3), Parameter :: &
          alpha = 20.0_ExDouble


    ! beta:Group2:Array2
    Real(Kind=ExDouble), Dimension(5), Parameter :: &
          beta = (/150.0_ExDouble, 150.0_ExDouble, 250.0_ExDouble, 0.3_ExDouble, 0.3_ExDouble  /)

    ! beta:Group2:Array2
    Real(Kind=ExDouble), Dimension(5), Parameter :: &
          betaInv = beta**(-1)


    ! gamma:Group2:Array3  
    Real(Kind=ExDouble), Dimension(3), Parameter :: &
          gamma = (/ 1.21_ExDouble, 1.21_ExDouble, 1.25_ExDouble /)


    ! epsilon:Group2:Array4
    Real(Kind=ExDouble), Dimension(3), Parameter :: &
          epsilon = 1.0_ExDouble



   ! --------------------------------------------------------------- !
   !                             Group 3                             !
   ! --------------------------------------------------------------- !
 
    ! a1:Group3:Array1
    Real(Kind=ExDouble), Dimension(2), Parameter :: &
          a1 = 0.32_ExDouble


    ! b1:Group3:Array2
    Real(Kind=ExDouble), Dimension(2), Parameter :: &
          b1 = 0.2_ExDouble


    ! c2:Group3:Array3
    Real(Kind=ExDouble), Dimension(2), Parameter :: &
          c2 = (/ 28.0_ExDouble, 32.0_ExDouble /)


    ! d1:Group3:Array4
    Real(Kind=ExDouble), Dimension(2), Parameter :: &
          d2 = (/ 700.0_ExDouble, 800.0_ExDouble /)


    ! a2:Group3:Array5
    Real(Kind=ExDouble), Dimension(2), Parameter :: &
          a2 = 3.5_ExDouble

    ! b2:Group3:Array6
    Real(Kind=ExDouble), Dimension(2), Parameter :: &
          b2 = (/ 0.85_ExDouble , 0.95_ExDouble /)


End Module HFEConstant
