module parameters
  IMPLICIT NONE

  !-- common parameters and variables ------------------------------
  ! THIS IS ALMOST PROJECT-INDEPENDENT 
  double precision, parameter :: tm32   = 1.d0/(2.d0**32.d0)
  double precision, parameter :: eps    = 1.d-14            ! very small number
  double precision, parameter :: tol    = 0.15d0            ! tolerance for Cor
  logical                     :: prt                        ! flag for write2file
  integer,          parameter :: Mxint  = 2147483647        ! maximum integer
  integer,          parameter :: Mnint  =-2147483647        ! minimum integer
  double precision :: pi=3.1415926

  !-- Parameters -------------------------------------------------
  integer, parameter :: D=2  
  integer, parameter :: UP=1
  integer, parameter :: DOWN=0
  integer, parameter :: MxL=512     !Max size of the system
  integer          :: PID      ! the ID of this job
  integer :: L     !the actual size of the system
  integer          :: Order
  double precision :: Beta
  double precision :: J
  double precision :: TotalStep  !total steps of this MC simulation
  double precision :: Step    ! a counter to keep track of the current step number
  double precision, parameter :: UpdateNum=6    ! number of updates
  integer                      :: Seed                   ! random-number seed

  !-- Diagram Permutation Table ----------------------------------
  integer, parameter :: PHYSICAL=1
  integer, parameter :: NORMALIZATION=2
  integer, parameter :: MaxOrder=8 ! Max diagram order
  integer, parameter :: MaxDiagNum=10000 ! Max diagram number 

  ! all diagram-related arraies have two copies: 1 is for physical diagrams, 2 is normalization diagrams
  integer, dimension(2*MaxOrder, MaxDiagNum, 2) :: Permutation 
  double precision,dimension(MaxDiagNum, 2) ::SymFactor  ! Symmetry Factor (includes diagram sign)
  integer, dimension(2*MaxOrder, MaxOrder+1, MaxDiagNum, 2) :: LoopBases ! Bases for loops
  integer, dimension(MaxOrder+1, 2) :: LoopValue ! values to attach to each loop basis
  complex, dimension(MaxDiagNum, 2) :: Weight  !weight of each diagram

  !-- Configuration -----------------------------------------------
  double precision, dimension(2*MaxOrder) :: TimeTable  !time variable for each vertex, all tau are between [0, beta)
  integer, dimension(D, 2*MaxOrder) :: CoordinateTable  !space variable for each vertex, all coordinates are between [1, L]
  integer, dimension(2, 2*MaxOrder) :: SpinTable  !spin variable for each vertex function, 1: spin for in-leg, 2:spin for out-leg
  integer                           :: Sector  !a flag, 1: physical diagram, 2:  normalization diagram
  double precision                  :: AbsWeight !the absolute value of the weight of all diagrams in current sector
  complex                           :: Phase     !the phase factor of the weight of all diagrams in current sector
  !the weight of all diagrams in current sector=AbsWeight*Phase

  !-- Measurement  ------------------------------------------------
  complex :: NormWeight    !the accumulated weight of the normalization diagram
  integer, parameter          :: TauBinNum=128     !number of tau bins
  complex, dimension(MxL, MxL, TauBinNum) :: Polarization  !the accumulated weight of the Spin-zz polarization
end module

INCLUDE "rng.f90"

program main
  use mt19937
  use parameters
  implicit none
  double precision :: x
  integer :: PrintCounter, SaveCounter

  print *, 'L, Beta, J, Order, TotalStep(*1e6), Seed, PID'
  read(*,*)  L, Beta, J, Order, TotalStep, Seed, PID
  Step=0.0
  PrintCounter=0
  SaveCounter=0
  call sgrnd(Seed) 
  !initialize rng seed
  
  call Test() !call test first to make sure all subroutines work as expected

  TotalStep=TotalStep*1.0e6
  print *, "Start simulation..."
  do while (Step<TotalStep)
    Step=Step+1.0
    x=grnd()
    if (x<1.0/UpdateNum) then
      call IncreaseOrder()
    else if (x<2.0/UpdateNum) then
      call DecreaseOrder()
    else if (x<3.0/UpdateNum) then
      call ChangeTau()
    else if (x<4.0/UpdateNum) then
      call ChangeCoordinate()
    else if (x<5.0/UpdateNum) then
      call ChangeSpin()
    endif
    call Measure()

    PrintCounter=PrintCounter+1
    if (PrintCounter==1e6)  then
      print *, Step, "million steps"
      PrintCounter=0
    endif

    SaveCounter=SaveCounter+1
    if (SaveCounter==1e8) then
      !write estimator to disk for every 1e6 steps
      call SaveToDisk()
      SaveCounter=0
    endif
  end do
  print *, "End simulation."

  CONTAINS

  subroutine Test()
    implicit none
    integer :: i
    ! Put all tests here
    if (cabs(Green(0.d0, -1.d0)+Green(0.d0, Beta-1.d0))>1e-6) then
      print *, "Green's function is not anti-periodic"
      stop
    endif
  end subroutine

  subroutine ReadDiagram()
    implicit none
    !Read diagrams from file to memory
    !Initialize physical diagrams
    Permutation(1:2*Order,1,PHYSICAL)=(/3,2,0,1/) !an example of a 2-Order diagram
    SymFactor(1,PHYSICAL)=-1
    LoopBases(1:2*Order, 1, 1, PHYSICAL)=(/1,0,1,0/)
    LoopBases(1:2*Order, 2, 1, PHYSICAL)=(/0,1,0,1/)
    LoopBases(1:2*Order, 3, 1, PHYSICAL)=(/1,0,0,1/)
    !Initialize normalization diagrams
    Permutation(1:2*(Order-1),1,NORMALIZATION)=(/1,0/) !an example of a 1-Order diagram
    SymFactor(1,NORMALIZATION)=1
    LoopBases(1:2*(Order-1), 1, 1, NORMALIZATION)=(/1,1/)
    LoopBases(1:2*(Order-1), 2, 1, NORMALIZATION)=(/1,0/)
  end subroutine

  complex function Green(t_in, t_out)
    !calculate Green's function
    implicit none
    double precision :: t_in, t_out, t
    integer :: sig
    t=t_out-t_in
    sig=1.0
    if (t<0) then
      t=t+Beta
      sig=-sig
    endif
    Green=sig*exp(cmplx(0.d0, -2.d0*pi/Beta*t))/cmplx(1.d0,1.d0)
    return
  end function Green

  complex function Interaction(IsBare, r_in, r_out, t_in, t_out, spin_in, spin_out)
    implicit none
    !calculate the weight of an Interaction line
    !Notice there are two types of interaction line:
    !IsBare==True: bare interaction, in this case, we expect t_in=t_out
    !IsBare==False: dressed interaction, in this case, we expect t_in!=t_out in general
    logical :: IsBare
    double precision :: t_in, t_out
    integer, dimension(D) :: r_in, r_out
    integer, dimension(2) :: spin_in, spin_out !there are spin indexes on both sides of the interaction line
    if (IsBare .eqv. .true.) then
      !calculate the weight of bare interaction
      Interaction=0.0
    else
      !calculate the weight of a dressed interaction
      Interaction=0.0
    endif
    return
  end function Interaction

  complex function CalcWeight(sector)
    !calculate the weight for ALL diagrams in a given sector
    implicit none
    integer :: sector

    CalcWeight=0.0
    return
  end function CalcWeight

  subroutine Measure()
    implicit none
    integer, dimension(D) :: dR
    integer :: dT
    if (Sector==NORMALIZATION) then
      NormWeight=NormWeight+Phase
    else
      call DeltaCoordinate(CoordinateTable(:,1), CoordinateTable(:,2), dR)
      dT=DeltaTime(TimeTable(1), TimeTable(2))
      if ((SpinTable(1, 1)==UP .and. SpinTable(2, 1)==UP) .and. (SpinTable(1,2)==UP .and. SpinTable(2,2)==UP)) then
        !we only measure Pi(up, up, up, up) here. Add other spin indexes if necessary
        Polarization(dR(1), dR(2), dT)=Polarization(dR(1), dR(2), dT)+Phase
      endif
    endif
    return
  end subroutine

  subroutine SaveToDisk()
    implicit none
    !Save NormWeight and Polarization to disk
    return
  end subroutine

  subroutine IncreaseOrder()
  !increase diagram order by one/change normalization diagram to physical diagram
    implicit none
    if (Sector==PHYSICAL) return
    !if the current diagrams are already in physical sector, then return
    return
  end subroutine

  subroutine DecreaseOrder()
  !decrease diagram order by one/change physical diagram to normalization diagram
    implicit none
    if (Sector==NORMALIZATION) return
    !if the current diagrams are already in normalization sector, then return
    return
  end subroutine

  subroutine ChangeTau()
  !randomly choose a vertex, change the time variable
    implicit none
    return
  end subroutine

  subroutine ChangeCoordinate()
  !randomly choose a vertex, change the space variable
    implicit none
    return
  end subroutine

  subroutine ChangeSpin()
  !randomly choose a loop basis, change the spin variable
    implicit none
    return
  end subroutine

  subroutine DeltaCoordinate(coord_in, coord_out, delta_coord)
    implicit none
    !calculate coord_out-coord_in
    integer :: i, dR
    integer, dimension(D) :: coord_in, coord_out, delta_coord
    do i=1,2
      dR=coord_out(i)-coord_in(i)
      if (dR<0) dR=dR+L
      delta_coord(i)=dR
    enddo
  end subroutine

  integer function DeltaTime(t_in, t_out)
    implicit none
    !calculate t_out-t_in
    double precision :: t_in, t_out, dt
    integer :: dt_bin
    dt=t_out-t_in
    if (dt<0) dt=dt+Beta
    DeltaTime=int(dt/Beta*TauBinNum)
    return
  end function

end program main
