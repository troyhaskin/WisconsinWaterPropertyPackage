Module Toolbox

   Implicit None

   ! These will be unneeded when the ISO_Fortran_Env module is fully supported
   ! Real Kind parameters
   Integer, Parameter :: Single   = Selected_Real_Kind(p=6)
   Integer, Parameter :: Double   = Selected_Real_Kind(p=15)
   Integer, Parameter :: ExDouble = Selected_Real_Kind(p=18)

   ! Integer Kind Parameters
   Integer, Parameter :: ShortInt = Selected_Int_Kind(8)
   Integer, Parameter :: LongInt  = Selected_Int_Kind(16)


   ! Shortcuts to high-precision constants
   Real(Kind=ExDouble), Parameter :: zero = 0.0_ExDouble
   Real(Kind=ExDouble), Parameter :: one  = 1.0_ExDouble
   Real(Kind=ExDouble), Parameter :: two  = 2.0_ExDouble
   Real(Kind=ExDouble), Parameter :: half = 0.5_ExDouble

   Integer(Kind=ShortInt), Parameter :: izero = 0
   Integer(Kind=ShortInt), Parameter :: ione  = 1
   Integer(Kind=ShortInt), Parameter :: itwo  = 2


   ! Just an approximation of absolute  machine zero
   Real(Kind=ExDouble), Parameter :: SmallNumber = 1E-16_ExDouble



   ! Default File names
   Character(Len=100),Parameter :: DefaultSizeFileName = 'Nelem.txt'
   Character(Len=100),Parameter :: DefaultDataFileName = 'DeltaTau.txt'
   Character(Len=100),Parameter :: DefaultOutputFileName = 'Output.txt'


Contains


   Elemental Subroutine KahanSum(SumN,Part,ErrorIn  ,  SumNp1,ErrorOut)
      Real(Kind=ExDouble), Intent(in)  :: SumN,Part, ErrorIn
      Real(Kind=ExDouble), Intent(out) :: SumNp1   , ErrorOut

      Real(Kind=ExDouble) :: Shift


      Shift    = Part - ErrorIn
      SumNp1   = SumN + Shift
      ErrorOut = (SumNp1 - SumN) - Shift

   End Subroutine KahanSum


   Subroutine ReadElementCount(Nelem)
      Integer,Intent(out) :: Nelem

      Character(Len=100) :: FileName = DefaultSizeFileName
      Integer :: UnitNumber = 50
      

      Open(Unit=UnitNumber,File=FileName,Status='Old',Action='Read');
      Read(UnitNumber,*) Nelem
      Close(UnitNumber)

      Return

   End Subroutine ReadElementCount

   Subroutine ReadDeltaAndTau(delta,tau)
      Real(Kind=ExDouble),Dimension(:),Intent(out) :: delta,tau

      Character(Len=100) :: FileName = DefaultDataFileName
      Integer :: UnitNumber = 50
      

      Open(Unit=UnitNumber,File=FileName,Status='Old',Action='Read');
      Read(UnitNumber,*) delta
      Read(UnitNumber,*) tau
      Close(UnitNumber)

      Return

   End Subroutine ReadDeltaAndTau


   Subroutine WriteOutput(Output)
      Real(Kind=ExDouble),Dimension(:),Intent(in)::Output

      Character(Len=100) :: FileName = DefaultOutputFileName
      Integer :: UnitNumber = 50


      Open(Unit=UnitNumber,File=FileName,Status='Replace',Action='Write');
      Write(UnitNumber,*) Output
      Close(UnitNumber)

      Return


   End Subroutine WriteOutput


End Module Toolbox
