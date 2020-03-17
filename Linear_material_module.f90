module Linear_material_module
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
!    10/03/20  X. P. Ai              Used in linear elastic material
!
!	Debugging Log:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
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



interface empty
 	module procedure empty_ELI
end interface empty

interface set
	module procedure set_ELI
end interface set

interface display
 	module procedure display_ELI
end interface display

interface extract
	module procedure extract_ELI
end interface extract



public :: set, empty, display, extract 






contains







pure subroutine empty_ELI (this)
  ! Purpose:
  ! to reset this ELI object's components into their default values (ZERO)

    type(ELI_material), intent(inout) :: this

    this%E         = ZERO
    this%nu        = ZERO

end subroutine empty_ELI






pure subroutine set_ELI (this, E, nu)
  ! Purpose:
  ! to set this ELI object's components during preprocessing before analysis

    type(ELI_material),       	 intent(inout) :: this
    real(DP),        			 intent(in)    :: E
    real(DP),                    intent(in)    :: nu

    this%E        = E
    this%nu       = nu

end subroutine set_ELI





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
  
  
  



end module Linear_material_module


