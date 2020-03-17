module material_module
!
!   Purpose:
!
!
!   Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!    30/03/15  B. Y. Chen            Original code
!    11/04/15  B. Y. Chen            Added lamina_ig_point object and its assoc.
!                                    procedures
!	 19/11/18  G.TostiBalducci 		 Modified to material_ELI_module for isotropics
!    30/11/18  G.TostiBalducci		 Removed G as a material type field
!	 01/12/18  G.TostiBalducci		 Added 'extract_ELI' subroutine to be used in
!									 TUBA3_toolkit_module
!    05/01/19  G.TostiBalducci		 Added type ply_material
!									 Added set_ply, extract_ply
!	 16/01/19  G.TostiBalducci	     Added cohesive_material type and set_coh
!    10/03/20  X. P. Ai              Removed type ply_material
!                                    Removed type cohesive_material
!                                    Removed empty_coh 
!                                    Removed set_ply
!                                    Removed set_coh
!                                    Removed extract_ply
!                                    Removed extract_coh
!
!	Debugging Log:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!	 02/12/18  G.TostiBalducci	     removed compile time errors
!	 02/12/18  G.TostiBalducci	     ELI_type, set, extract work as expected
!
use parameter_module, only : DP, ZERO, ONE, TWO

implicit none
private


!**** actual objects of the module ****
! with encapsulation and allowed procedures

! elastic material, standard notations
type, public :: ELI_material
  private
  real(DP) :: E  = ZERO
  real(DP) :: nu = ZERO
end type ELI_material



!~type, public :: ply_material
! 1,2: material axes
! ply: transverse isotropy + plane stress
!~  private
!~ real(DP) :: E1   = ZERO
!~  real(DP) :: E2   = ZERO
!~  real(DP) :: nu12 = ZERO
!~  real(DP) :: G12  = ZERO
!~end type ply_material



!~type, public :: cohesive_material
! material constants for cohe el
!~  private
!~  real(DP) :: E3   = ZERO
!~ real(DP) :: GcI   = ZERO
!~  real(DP) :: GcII = ZERO
!~  real(DP) :: tau0I  = ZERO
!~  real(DP) :: tau0II  = ZERO
!~  real(DP) :: BKcoeff  = ZERO
!~end type cohesive_material
! associated procedures: empty, set, display, ddsdde
! all other operations on the object are performed within this module

! ELI material integration point object
! stores everything needed for the integration of ELI material in elements
! G: The sdvs are allocated at the UEL level, in the declaration part
! G: Don't need to define a type specifically for the sdv in the ELI case:
! G: the variables are all real and are stresses and strains
!~type, public :: ELI_ig_point
!~  private
!~  real(DP), allocatable :: x(:)
!~  real(DP), allocatable :: u(:)
!~  real(DP) :: sdv_previous(2*NST) !include strains and stresses
!~  real(DP) :: sdv_current(2*NST) !include strains and stresses
!~end type ELI_ig_point
! associated procedures: empty, display, update, extract
! all other operations on the object are performed within this module


interface empty
 	module procedure empty_ELI
!~	module procedure empty_coh
!~	module procedure empty_ELI_ig_point
end interface empty

interface set
	module procedure set_ELI
!~	module procedure set_ply
!~	module procedure set_coh
end interface set

interface display
 	module procedure display_ELI
!~	module procedure display_ELI_ig_point
end interface display

!~interface ddsdde
!~    module procedure ddsdde_ELI
!~end interface ddsdde

!~interface update
!~    module procedure update_ELI_ig_point
!~end interface update

interface extract
	module procedure extract_ELI
!~	module procedure extract_ply
!~	module procedure extract_coh
end interface extract



public :: set, empty, display, extract !update, ddsdde







contains







pure subroutine empty_ELI (this)
  ! Purpose:
  ! to reset this ELI object's components into their default values (ZERO)

    type(ELI_material), intent(inout) :: this

    this%E         = ZERO
    this%nu        = ZERO

end subroutine empty_ELI




!~pure subroutine empty_coh(this)
!~
!~	type(cohesive_material),	intent(inout)	::	this
!~	
!~	this%E3=ZERO
!~	this%GcI=ZERO
!~	this%GcII=ZERO
!~	this%tau0I=ZERO
!~	this%tau0II=ZERO
!~	this%BKcoeff=ZERO
!~	
!~end subroutine empty_coh




pure subroutine set_ELI (this, E, nu)
  ! Purpose:
  ! to set this ELI object's components during preprocessing before analysis

    type(ELI_material),       	 intent(inout) :: this
    real(DP),        			 intent(in)    :: E
    real(DP),                    intent(in)    :: nu

    this%E        = E
    this%nu       = nu

end subroutine set_ELI



!~pure subroutine set_ply (this, E1, E2, nu12, G12)
  ! Purpose:
  ! to set this ply object's components during preprocessing before analysis

!~    type(ply_material),       		 intent(inout) :: this
!~    real(DP),        			 	 intent(in)    :: E1
!~    real(DP),        			 	 intent(in)    :: E2
!~    real(DP),                    	 intent(in)    :: nu12
!~    real(DP),                    	 intent(in)    :: G12
!~
!~    this%E1        = E1
!~    this%E2        = E2
!~   this%nu12      = nu12
!~    this%G12       = G12

!~end subroutine set_ply



!~pure subroutine set_coh (this, E3, GcI, GcII, tau0I, tau0II, BKcoeff)
!~  ! Purpose:
!~  ! to set this ply object's components during preprocessing before analysis
!~
!~    type(cohesive_material),       	 intent(inout) :: this
!~    real(DP),        			 intent(in)    :: E3
!~    real(DP), optional,       	 intent(in)    :: GcI
!~    real(DP), optional,      	 intent(in)    :: GcII
!~    real(DP), optional,       	 intent(in)    :: tau0I
!~    real(DP), optional,       	 intent(in)    :: tau0II
!~    real(DP), optional,       	 intent(in)    :: BKcoeff
!~
!~    this%E3        = E3
!~    if (present(GcI)) this%GcI=GcI
!~    if (present(GcII)) this%GcII=GcII
!~    if (present(tau0I)) this%tau0I=tau0I
!~    if (present(tau0II)) this%tau0II=tau0II
!~    if (present(BKcoeff)) this%BKcoeff=BKcoeff
!~
!~end subroutine set_coh



pure subroutine extract_ELI(E, nu, material)
!	Purpose: 
!	takes the fields of ELI_material type object and stores them in variables
!	to be used by modules different than this one

	real(DP), 			optional,   intent(out) :: E
	real(DP), 			optional,   intent(out) :: nu
	type(ELI_material), optional,  	intent(in)  :: material
	
	if(present(E)) then
		E=material%E
	end if
	
	if(present(nu)) then
		nu=material%nu
	end if
	
end subroutine extract_ELI



!~pure subroutine extract_ply(E1, E2, nu12, G12, material)
!~!	Purpose: 
!~!	takes the fields of ply_material type object and stores them in variables
!~!	to be used by modules different than this one
!~
!~	real(DP), 			optional,   intent(out) :: E1
!~	real(DP), 			optional,   intent(out) :: E2
!~	real(DP), 			optional,   intent(out) :: nu12
!~	real(DP), 			optional,   intent(out) :: G12
!~	type(ply_material), optional,  	intent(in)  :: material
!~	
!~	if(present(E1)) then
!~		E1=material%E1
!~	end if
!~	
!~	if(present(E2)) then
!~		E2=material%E2
!~	end if
!~	
!~	if(present(nu12)) then
!~		nu12=material%nu12
!~	end if
!~	
!~	if(present(G12)) then
!~		G12=material%G12
!~	end if
!~	
!~end subroutine extract_ply




!~pure subroutine extract_coh(E3, GcI, GcII, tau0I, tau0II, BKcoeff, material)
!~!	Purpose: 
!~!	takes the fields of coh_material type object and stores them in variables
!~!	to be used by modules different than this one
!~
!~	real(DP), 			optional,   intent(out) :: E3
!~	real(DP), 			optional,   intent(out) :: GcI
!~	real(DP), 			optional,   intent(out) :: GcII
!~	real(DP), 			optional,   intent(out) :: tau0I
!~	real(DP), 			optional,   intent(out) :: tau0II
!~	real(DP), 			optional,   intent(out) :: BKcoeff
!~	type(cohesive_material), 		intent(in)  :: material
!~	
!~	if(present(E3)) then
!~		E3=material%E3
!~	end if
!~	
!~	if(present(GcI)) then
!~		GcI=material%GcI
!~	end if
!~	
!~	if(present(GcII)) then
!~		GcII=material%GcII
!~	end if
!~	
!~	if(present(tau0I)) then
!~		tau0I=material%tau0I
!~	end if
!~	
!~	if(present(tau0II)) then
!~		tau0II=material%tau0II
!~	end if
!~	
!~	if(present(BKcoeff)) then
!~		BKcoeff=material%BKcoeff
!~	end if
!~	
!~end subroutine extract_coh



subroutine display_ELI (this)
!	G: cannot define the subroutine as 'pure', because pure SRs 
!	G: can't do I/O (no write allowed)

  ! Purpose:
  ! to display this ELI object's components on cmd window
  ! this is useful for debugging

    type(ELI_material), intent(in) :: this

    ! local variable to set the output format
    character(len=20) :: display_fmt

    ! initialize local variable
    display_fmt = ''

    ! set display format for string and real
    ! A for string, ES for real (scientific notation)
    ! 10 is width, 3 is no. of digits aft decimal point
    ! note that for scientific real, ESw.d, w>=d+7
    display_fmt = '(1X, A, ES10.3)'

    write(*,'(1X, A)') ''
    write(*,'(1X, A)') 'Display the properties of the inquired ELI object :'
    write(*,display_fmt) 'lamina E   is: ', this%E
    write(*,display_fmt) 'lamina nu  is: ', this%nu

end subroutine display_ELI
  
  
  
!~pure subroutine ddsdde_ELI(this_mat, dee, stress, strain)
!~! Purpose
!~!
!~!
!~  !dummy argument list
!~  !
!~  !
!~  !
!~  type(ELI_material),	intent(in)	 :: this_mat
!~  real(DP),				intent(inout):: dee(:,:)
!~  real(DP),				intent(inout):: stress(:)
!~  real(DP),				intent(in)   :: strain(:)
!~  
!~  !**** check validity of non-pass dummy arguments with intent(in/inout) ****
!~  ! they include : dee, stress and strain
!~
!~  ! dee and stress input values are not used; they can be intent(out).
!~  ! they are defined as intent(inout) to avoid any potential memory leak.
!~  ! so no need to check their input values
!~
!~  ! strain components can take any real value, nothing to check
!~  
!~  !**** MAIN CALCULATIONS ****
!~
!~  ! calculate dee using original material properties only
!~  call deemat_2d (this_mat, dee)
!~
!~  ! calculate stress
!~  stress = matmul(dee, strain)
!~
!~end subroutine ddsdde_ELI



!~pure subroutine deemat_2d (this_mat, dee)
!~	!Purpose:
!~
!~	type(ELI_material), intent(in)	    :: this_mat
!~	real(DP)		  , intent(inout)	:: dee(:,:)
!~	
!~	!local var list
!~	real(DP) :: E, nu
!~	real(DP) :: comm_fac
!~	
!~	!G: ddsdde_ELI does not know the interface of deemat_2d
!~	! for private procedures used in the main procedure (ddsdde_ELI),
!~    ! the values of the input arguments are not checked.
!~    ! they are assumed to be of valid values.
!~	
!~	
!~	!initialize locals
!~	E = ZERO
!~	nu = ZERO
!~	
!~	comm_fac = ZERO
!~	
!~	!extract ELI material properties
!~	E = this_mat%E
!~	nu = this_mat%nu
!~	
!~	!calculate D matrix terms
!~	!zero all terms before
!~	dee = ZERO
!~	
!~	comm_fac = (E*(ONE-nu))/((ONE+nu)*(ONE-TWO*nu))
!~	
!~	dee(1,1) = comm_fac
!~	dee(1,2) = comm_fac*(nu/(ONE-nu))
!~	dee(2,1) = comm_fac*(nu/(ONE-nu))
!~	dee(2,2) = comm_fac
!~	dee(3,3) = comm_fac*((ONE-TWO*nu)/(TWO*(ONE-nu)))
!~	
!~end subroutine deemat_2d
  


!~pure subroutine empty_ELI_ig_point (ig_point)
!~
!~    type(ELI_ig_point), intent(inout) :: ig_point
!~
!~    ! local varibale; derived type is initialized upon declaration
!~    type(ELI_ig_point) :: ig_point_lcl
!~
!~    ! reset the dummy arg. to initial values
!~    ig_point = ig_point_lcl
!~
!~end subroutine empty_ELI_ig_point
!~
!~
!~
!~pure subroutine update_ELI_ig_point(ig_point, x, u, sdv_previous, sdv_current)
!~	! Purpose:
!~    ! to update ELI ig point components
!~    ! this is an outbound procedure, so its inputs should be checked for validity
!~    ! x, u, traction, separation can take any value
!~    ! ELI sdv objects are only modified within this module
!~    ! so, checking is spared
!~	!G: for ELI mat. ig_point --> the only sdvs are stresses/strains
!~		
!~	type(ELI_ig_point),intent(inout) :: ig_point
!~	real(DP), optional, intent(in)	 :: x(NDIM)
!~	real(DP), optional, intent(in)	 :: u(NDIM)
!~	real(DP), optional, intent(in)	 :: sdv_previous(2*NST)
!~	real(DP), optional, intent(in)	 :: sdv_current(2*NST)
!~	
!~	if (present(x)) then
!~		if(.not. allocated(ig_point%x)) allocate(ig_point%x(NDIM))
!~		ig_point%x = x
!~	end if
!~	
!~	if (present(u)) then
!~		if(.not. allocated(ig_point%u)) allocate(ig_point%u(NDIM))
!~		ig_point%u = u
!~	end if
!~	
!~	if (present(sdv_previous)) ig_point%sdv_previous = sdv_previous
!~	if (present(sdv_current))  ig_point%sdv_current  = sdv_current
!~	
!~end subroutine update_ELI_ig_point
!~
!~
!~
!~pure subroutine extract_ELI_ig_point(x, u, sdv_previous, sdv_current, ig_point)
!~	!
!~	real(DP), optional,     intent(out) :: x(NDIM)
!~	real(DP), optional,     intent(out) :: u(NDIM)
!~	real(DP), optional,     intent(out) :: sdv_previous(2*NST)
!~	real(DP), optional,     intent(out) :: sdv_current(2*NST)
!~	type(ELI_ig_point), intent(in)  :: ig_point
!~	
!~	if(present(x)) then
!~		if (allocated(ig_point%x)) x=ig_point%x
!~	end if
!~	
!~	if(present(u)) then
!~		if (allocated(ig_point%u)) u=ig_point%u
!~	end if
!~	
!~	if(present(sdv_previous)) sdv_previous=ig_point%sdv_previous
!~	if(present(sdv_current)) sdv_current=ig_point%sdv_current
!~	
!~end subroutine extract_ELI_ig_point
!~
!~
!~
!~subroutine display_ELI_ig_point(this)
!~	!G: cannot define the subroutine as 'pure', because pure SRs 
!~	!G: can't do I/O (no write allowed)
!~	
!~	type(ELI_ig_point), intent(in)  :: this
!~	!declare local variables
!~	character(len=20) :: display_fmt
!~	integer			  :: i
!~	
!~	!initialize local variables
!~	display_fmt = ''
!~	i=0
!~
!~    ! set display format for real
!~    ! ES for real (scientific notation)
!~    ! 10 is width, 3 is no. of digits aft decimal point
!~    ! note that for scientific real, ESw.d, w>=d+7
!~    display_fmt = '(ES10.3)'
!~
!~    write(*,'(1X, A)') ''
!~    write(*,'(1X, A)') 'Display the components of the inquired ELI_ig_point &
!~                       &object :'
!~    write(*,'(1X, A)') ''
!~	
!~    if (allocated(this%x)) then
!~      write(*,'(1X, A)') '- x of this ELI_ig_point is: '
!~      do i = 1, NDIM
!~        write(*,display_fmt,advance="no") this%x(i)
!~      end do
!~      write(*,'(1X, A)') ''
!~      write(*,'(1X, A)') ''
!~    end if
!~	
!~	if (allocated(this%u)) then
!~      write(*,'(1X, A)') '- u of this ELI_ig_point is: '
!~      do i = 1, NDIM
!~        write(*,display_fmt,advance="no") this%u(i)
!~      end do
!~      write(*,'(1X, A)') ''
!~      write(*,'(1X, A)') ''
!~    end if
!~	
!~    write(*,'(1X, A)') '- previous strains of this ELI_ig_point were: '
!~    write(*,'(1X, A)', advance = "no") 'e_xx = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_previous(1)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 'e_yy = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_previous(2)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 'g_xy = ' 
!~    write(*,display_fmt) this%sdv_previous(3)
!~	
!~	write(*,'(1X, A)') '- previous stresses of this ELI_ig_point were: '
!~    write(*,'(1X, A)', advance = "no") 's_xx = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_previous(4)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 's_yy = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_previous(5)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 't_xy = ' 
!~    write(*,display_fmt) this%sdv_previous(6)
!~	
!~	write(*,'(1X, A)') '- current strains of this ELI_ig_point are: '
!~    write(*,'(1X, A)', advance = "no") 'e_xx = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_current(1)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 'e_yy = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_current(2)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 'g_xy = ' 
!~    write(*,display_fmt) this%sdv_current(3)
!~	
!~	write(*,'(1X, A)') '- current stresses of this ELI_ig_point are: '
!~    write(*,'(1X, A)', advance = "no") 's_xx = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_current(4)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 's_yy = ' 
!~    write(*,display_fmt, advance = "no") this%sdv_current(5)
!~	write(*,'(1X, A)', advance = "no") ''
!~	write(*,'(1X, A)', advance = "no") 't_xy = ' 
!~    write(*,display_fmt) this%sdv_current(6)
!~
!~	
!~end subroutine display_ELI_ig_point




end module material_module


