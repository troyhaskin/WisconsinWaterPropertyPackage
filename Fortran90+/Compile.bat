gfortran -c Toolbox.f90
gfortran -c HFEConstant.f90 Toolbox.f90
gfortran -c HelmholtzFunctions_TheResidual.f90 HFEConstant.f90 Toolbox.f90
gfortran -o HelmholtzResidual_d.exe HelmholzResidual_d.f90 HelmholtzFunctions_TheResidual.f90 HFEConstant.f90 Toolbox.f90