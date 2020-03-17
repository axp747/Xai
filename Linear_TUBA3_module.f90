module Linear_TUBA3_module
!
!   Purpose:
!
!
!   Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!	 11/03/20  X. P. Ai	             Original conde 
!
!
!
!
!	Debugging Log for Original Code:
!    Date      Programmer            Description
!    ========  ====================  ========================================
!
!
!
!

use parameter_module,		only: DP, ZERO, ONE, TWO, THREE, FOUR, TEN,&
& SIXTY, HALF, TWO_THIRD, FIVE_TWO, STAT_SUCCESS, STAT_FAILURE, MSGLENGTH

use material_module,	only: ELI_material, extract

implicit none
private


!	include list of public element-dependent parameters
integer, parameter, public	:: 	NDIM     = 3  ! number of dimension 
integer, parameter, public	:: 	NST      = 3  ! number of stress
integer, parameter, public	:: 	NNODES   = 3  ! number of nodes
integer, parameter, public	::  NDOF     = 18 ! number of total degree of freedom
integer, parameter, public	:: 	NDOF_NODE = 6 ! number of degree of freedom per node
integer, parameter, public	:: 	NIGPOINTS = 7 ! number of IG points
!------------------------------------------------------------------

type, public :: TUBA3_element
!	This time both element and components are public, to
!	be accessed by the uel and directly transferred to the svars()
	real(DP)	::	w 		= ZERO
	real(DP)	::	curv(3) = ZERO
	real(DP)	::	mom(3)	= ZERO

end type





interface Qform
	module procedure Qform
end interface Qform

interface Dplate_form
!	Dplate_form->ELI material
	module procedure Dplate_form
end interface Dplate_form

interface Rform
	module procedure Rform
end interface Rform

interface moments_point
	module procedure moments_point
end interface moments_point

interface sdvs_centroid
	module procedure sdvs_centroid
end interface sdvs_centroid

public	::	 Qform, Dplate_form, Rform, moments_point, sdvs_centroid





contains


	



pure subroutine Qform(coord_nod,area,QXX,QYY,QXY,istat,emsg) ! 
!	Purpose:
!	Forms the matrices Qxx, Qyy, Qxy for TUBA3 element, as from appendix
!	of the article 'A higher order triangular plate bending element
!	revisited' by S. Dasgupta and D. Sengupta
!
!	Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!	 28/11/18  G. Tosti Balducci     Conversion original .f->.f90
!	 01/12/18  G. Tosti Balducci	 Changed coord_nod(NDIM,NNODE)
!	 							     to coord_nod(NDIM-1,NNODE) in Qform
!    10/03/20  X. P. Ai	             Copy
!       
!	Record of debugging:
!	 Version						 Debug status
!	=====================			 ========================================
!	Conversion original				 no compile-time errors
!
!


!   Declaration
!	I/O
	real(DP), intent(in)	::	coord_nod(NDIM-1,NNODES) !coords on mid-surface
	real(DP), intent(out)	::  area
	real(DP), intent(out)	::  QXX(10,NDOF), QYY(10,NDOF), QXY(10,NDOF)
	integer,  intent(out)   ::  istat
	character(len=MSGLENGTH), intent(out)   :: emsg
!	LOCAL variables
	integer, parameter		::	NDEST(NDOF)=(/4,5,6,8,7,9,7,6,9,4,8,5,8,9,&
											&5,7,4,6/)
	real(DP)				::  XT(NNODES), YT(NNODES)
	real(DP)				::	B(NNODES), C(NNODES)
	real(DP)				::  TWOD, AR
	real(DP)				::  P(NNODES)
	real(DP)				::  R(NNODES,NNODES)
!	counters
	integer					::	I,J,M,IND,JC,JJ,MM,K,IPER,NR
	
!	Initialization
	area=ZERO
	QXX=ZERO
	QYY=ZERO
	QXY=ZERO
	istat=STAT_SUCCESS
	emsg=''
	XT=coord_nod(1,:)
	YT=coord_nod(2,:)
	B=ZERO
	C=ZERO
	TWOD=ZERO
	AR=ZERO
	P=ZERO
	R=ZERO
	
	I=0; J=0; M=0;
	IND=0;
	JC=0; JJ=0; MM=0;
	K=0;
	IPER=0;
	NR=0;
	
!	Preliminaries
	B(1)=YT(2)-YT(3)
	B(2)=YT(3)-YT(1)
	B(3)=YT(1)-YT(2)
	C(1)=XT(3)-XT(2)
	C(2)=XT(1)-XT(3)
	C(3)=XT(2)-XT(1)
	TWOD=B(1)*C(2)-B(2)*C(1)
	AR=HALF*TWOD
	area=AR
	P(1)=B(2)*C(3)+B(3)*C(2)
	P(2)=B(3)*C(1)+B(1)*C(3)
	P(3)=B(1)*C(2)+B(2)*C(1)
	do I=1,NNODES
		do J=1,NNODES
			R(I,J)=-(B(I)*B(J)+C(I)*C(J))/(B(I)*B(I)+C(I)*C(I))
		end do
	end do
	IND=0
	
!	Main calculations
	main_loop: do I=1,NNODES
		JC=(I-1)*6
		QXX(I,JC+4)=TWOD**2
	    QYY(I,JC+6)=TWOD**2
	    QXY(I,JC+5)=2._DP*TWOD**2
		J=I+1
	    IF(J>3) J=J-3
		M=J+1
	    IF(M>3) M=M-3
	    JJ=J
	    MM=M
		QXX(10,JC+1)=120._DP*(B(I)*B(I)+2._DP*B(I)*(R(J,I)*B(M)+R(M,I)*B(J)))
	    QXX(10,JC+2)=-48._DP*B(I)*TWOD-120._DP*B(I)*(R(J,I)*B(M)*C(J)-R(M,I)*&
        &B(J)*C(M))
		QXX(10,JC+3)=120._DP*B(I)*B(J)*B(M)*(R(J,I)-R(M,I))
	    QXX(10,JC+4)=-6._DP*B(I)*B(I)*C(J)*C(M)+8._DP*B(I)*C(I)*P(I)+20._DP*&
		&B(I)*(R(J,I)*B(M)*C(J)*C(J)+R(M,I)*B(J)*C(M)*C(M))
		QXX(10,JC+5)=-2._DP*B(I)*B(I)*P(I)-B(I)*B(J)*B(M)*(16._DP*C(I)+&
		&40._DP*(R(J,I)*C(J)+R(M,I)*C(M)))
		QXX(10,JC+6)=10._DP*B(I)*B(I)*B(J)*B(M)+20._DP*B(I)*B(M)*B(J)*&
		&(R(J,I)*B(J)+R(M,I)*B(M))
!
!
		QYY(10,JC+1)=120._DP*(C(I)*C(I)+2._DP*C(I)*(R(J,I)*C(M)+R(M,I)*C(J)))
		QYY(10,JC+2)=120._DP*C(I)*C(J)*C(M)*(R(M,I)-R(J,I))
		QYY(10,JC+3)=-48._DP*C(I)*TWOD+120._DP*C(I)*(R(J,I)*C(M)*B(J)-&
		&R(M,I)*C(J)*B(M))
		QYY(10,JC+4)=10._DP*C(I)*C(I)*C(J)*C(M)+20._DP*C(I)*C(J)*C(M)*&
		&(R(J,I)*C(J)+R(M,I)*C(M))
		QYY(10,JC+5)=-2._DP*C(I)*C(I)*P(I)-C(I)*C(J)*C(M)*(16._DP*B(I)+40._DP*&
		&(R(J,I)*B(J)+R(M,I)*B(M)))
		QYY(10,JC+6)=-6._DP*C(I)*C(I)*B(J)*B(M)+8._DP*C(I)*B(I)*P(I)+&
		&20._DP*C(I)*(R(J,I)*C(M)*B(J)**2+R(M,I)*C(J)*B(M)**2)
!
!
		QXY(10,JC+1)=240._DP*(B(I)*C(I)+R(J,I)*P(J)+R(M,I)*P(M))
		QXY(10,JC+2)=-48._DP*C(I)*TWOD-120._DP*(R(J,I)*C(J)*P(J)-R(M,I)*&
		&C(M)*P(M))
		QXY(10,JC+3)=-48._DP*B(I)*TWOD+120._DP*(R(J,I)*B(J)*P(J)-&
		&R(M,I)*B(M)*P(M))
		QXY(10,JC+4)=4._DP*B(I)*C(I)*C(J)*C(M)+8._DP*C(I)*C(I)*P(I)+&
		&20._DP*(R(J,I)*C(J)*C(J)*P(J)+R(M,I)*C(M)*C(M)*P(M))
		QXY(10,JC+5)=12._DP*B(I)*C(I)*P(I)-16._DP*P(J)*P(M)-40._DP*(R(J,I)*&
		&B(J)*C(J)*P(J)+R(M,I)*B(M)*C(M)*P(M))
		QXY(10,JC+6)=4._DP*C(I)*B(I)*B(J)*B(M)+8._DP*B(I)*B(I)*P(I)+20._DP*&
		&(R(J,I)*B(J)*B(J)*P(J)+R(M,I)*B(M)*B(M)*P(M))
		
!
!
		row_selec_loop: do K=1,3
			iper_loop :	do IPER=1,2
				select case (IPER)
				case(1)
					J=JJ
					M=MM
				case(2)
					J=MM
					M=JJ
				case default
					istat=STAT_FAILURE
					emsg='IPER is diff. from 1 or 2, Qform, &
					& element_module'
					QXX=ZERO; QYY=ZERO; QXY=ZERO
					return
				end select
				IND=IND+1
				NR=NDEST(IND)
				select case (K)
				case(1)
					QXX(NR,JC+1)=60._DP*(-B(I)*B(I)+R(J,I)*B(M)*B(M)+2._DP*R(M,I)*B(J)*&
					&B(M))
					QXX(NR,JC+2)=12._DP*B(I)*TWOD+(-1._DP)**IPER*(24._DP*B(I)*B(I)*C(M)+18._DP*&
					&B(M)*B(M)*C(I)+30._DP*R(J,I)*B(M)*B(M)*C(J)-60._DP*R(M,I)*B(J)*B(M)*&
					&C(M))
					QXX(NR,JC+3)=-(-1._DP)**IPER*(24._DP*B(I)*B(I)*B(M)+18._DP*B(M)*B(M)*B(I)&
					&+30._DP*R(J,I)*B(M)*B(M)*B(J)-60._DP*R(M,I)*B(M)*B(M)*B(J))
					QXX(NR,JC+4)=-(-1._DP)**IPER*6._DP*B(I)*C(M)*TWOD+2._DP*B(M)*B(M)*C(I)*&
					&C(J)+4._DP*B(J)*B(M)*C(M)*C(I)+5._DP*R(J,I)*B(M)*B(M)*C(J)*C(J)+&
					&10._DP*R(M,I)*B(J)*B(M)*C(M)*C(M)
					QXX(NR,JC+5)=(-1._DP)**IPER*10._DP*B(I)*B(M)*TWOD-6._DP*B(M)*B(M)*P(M)-&
					&10._DP*R(J,I)*B(M)*B(M)*B(J)*C(J)-20._DP*R(M,I)*B(M)*B(M)*B(J)*C(M)
					QXX(NR,JC+6)=6._DP*B(I)*B(J)*B(M)*B(M)+5._DP*R(J,I)*B(J)*B(J)*B(M)*&
					&B(M)+10._DP*R(M,I)*B(J)*B(M)*B(M)*B(M)
!
!
					QYY(NR,JC+1)=60._DP*(-C(I)*C(I)+R(J,I)*C(M)*C(M)+2._DP*R(M,I)*C(J)*&
					&C(M))
					QYY(NR,JC+2)=(-1._DP)**IPER*(24._DP*C(I)*C(I)*C(M)+18._DP*C(M)*C(M)*C(I)&
					&+30._DP*R(J,I)*C(M)*C(M)*C(J)-60._DP*R(M,I)*C(M)*C(M)*C(J))
					QYY(NR,JC+3)=12._DP*C(I)*TWOD-(-1._DP)**IPER*(24._DP*C(I)*C(I)*B(M)+&
					&18._DP*C(M)*C(M)*B(I)+30._DP*R(J,I)*C(M)*C(M)*B(J)-60._DP*R(M,I)*C(J)*&
					&C(M)*B(M))
					QYY(NR,JC+4)=6._DP*C(I)*C(J)*C(M)*C(M)+5._DP*R(J,I)*C(J)*C(J)*C(M)*&
					&C(M)+10._DP*R(M,I)*C(J)*C(M)**3
					QYY(NR,JC+5)=-(-1._DP)**IPER*10._DP*C(I)*C(M)*TWOD-6._DP*C(M)*C(M)*&
					&P(M)-10._DP*R(J,I)*C(M)*C(M)*C(J)*B(J)-20._DP*R(M,I)*C(M)*C(M)*&
					&C(J)*B(M)
					QYY(NR,JC+6)=(-1._DP)**IPER*6._DP*C(I)*B(M)*TWOD+2._DP*C(M)*C(M)*B(I)*&
					&B(J)+4._DP*C(J)*C(M)*B(M)*B(I)+5._DP*R(J,I)*C(M)*C(M)*B(J)*B(J)+&
					&10._DP*R(M,I)*C(J)*C(M)*B(M)*B(M)
!
!
					QXY(NR,JC+1)=120._DP*(-B(I)*C(I)+R(J,I)*B(M)*C(M)+R(M,I)*P(I))
					QXY(NR,JC+2)=12._DP*C(I)*TWOD+(-1._DP)**IPER*(48._DP*B(I)*C(I)*C(M)+&
					&36._DP*B(M)*C(M)*C(I)+60._DP*R(J,I)*B(M)*C(M)*C(J)-60._DP*R(M,I)*&
					&C(M)*P(I))
					QXY(NR,JC+3)=12._DP*B(I)*TWOD-(-1._DP)**IPER*(48._DP*C(I)*B(I)*B(M)+&
					&36._DP*C(M)*B(M)*B(I)+60._DP*R(J,I)*C(M)*B(M)*B(J)-60._DP*R(M,I)*&
					&B(M)*P(I))
					QXY(NR,JC+4)=10._DP*C(I)*B(J)*C(M)*C(M)+2._DP*C(I)*C(J)*B(M)*C(M)+&
					&10._DP*R(J,I)*C(J)*C(J)*B(M)*C(M)+10._DP*R(M,I)*C(M)*C(M)*P(I)
					QXY(NR,JC+5)=-16._DP*B(M)*C(M)*P(M)+2._DP*P(I)*P(J)-20._DP*R(J,I)*B(J)*&
					&C(J)*B(M)*C(M)-20._DP*R(M,I)*B(M)*C(M)*P(I)
					QXY(NR,JC+6)=10._DP*B(I)*C(J)*B(M)*B(M)+2._DP*B(I)*B(J)*C(M)*B(M)+&
					&10._DP*R(J,I)*B(J)*B(J)*C(M)*B(M)+10._DP*R(M,I)*B(M)*B(M)*P(I)
!
!
				case(2)
					QXX(NR,JC+1)=60._DP*(B(I)*B(I)+2._DP*R(M,I)*B(M)*B(I))
					QXX(NR,JC+2)=12._DP*B(I)*TWOD-(-1._DP)**IPER*(36._DP*B(I)*B(I)*C(M)+&
					&60._DP*R(M,I)*B(M)*B(I)*C(M))
					QXX(NR,JC+3)=(-1._DP)**IPER*(36._DP*B(I)*B(I)*B(M)+60._DP*R(M,I)*&
					&B(M)*B(M)*B(I))
					QXX(NR,JC+4)=3._DP*B(I)*B(I)*C(M)*C(M)+4._DP*B(I)*B(M)*C(I)*C(M)+&
					&10._DP*R(M,I)*B(I)*B(M)*C(M)*C(M)
					QXX(NR,JC+5)=-10._DP*B(I)*B(I)*C(M)*B(M)-4._DP*B(M)*B(M)*C(I)*B(I)&
					&-20._DP*R(M,I)*B(M)*B(M)*B(I)*C(M)
					QXX(NR,JC+6)=7._DP*B(I)*B(I)*B(M)*B(M)+10._DP*R(M,I)*B(I)*B(M)**3
!
!
					QYY(NR,JC+1)=60._DP*(C(I)*C(I)+2._DP*R(M,I)*C(M)*C(I))
					QYY(NR,JC+2)=-(-1._DP)**IPER*(36._DP*C(I)*C(I)*C(M)+60._DP*R(M,I)*&
					&C(M)*C(M)*C(I))
					QYY(NR,JC+3)=12._DP*C(I)*TWOD+(-1._DP)**IPER*(36._DP*C(I)*C(I)*B(M)+60._DP&
					&*R(M,I)*C(M)*C(I)*B(M))
					QYY(NR,JC+4)=7._DP*C(I)*C(I)*C(M)*C(M)+10._DP*R(M,I)*C(I)*C(M)**3
					QYY(NR,JC+5)=-10._DP*C(I)*C(I)*B(M)*C(M)-4._DP*C(M)*C(M)*B(I)*C(I)-&
					&20._DP*R(M,I)*C(M)*C(M)*C(I)*B(M)
					QYY(NR,JC+6)=3._DP*C(I)*C(I)*B(M)*B(M)+4._DP*C(I)*C(M)*B(I)*B(M)+&
					&10._DP*R(M,I)*C(I)*C(M)*B(M)*B(M)
!
!
					QXY(NR,JC+1)=120._DP*(B(I)*C(I)+R(M,I)*P(J))
					QXY(NR,JC+2)=12._DP*C(I)*TWOD-(-1._DP)**IPER*(72._DP*B(I)*C(I)*C(M)+&
					&60._DP*R(M,I)*C(M)*P(J))
					QXY(NR,JC+3)=12._DP*B(I)*TWOD+(-1._DP)**IPER*(72._DP*C(I)*B(I)*B(M)+60._DP*&
					&R(M,I)*B(M)*P(J))
					QXY(NR,JC+4)=10._DP*C(I)*B(I)*C(M)*C(M)+4._DP*C(I)*C(I)*B(M)*C(M)+&
					&10._DP*R(M,I)*C(M)*C(M)*P(J)
					QXY(NR,JC+5)=-12._DP*B(I)*C(I)*B(M)*C(M)-4._DP*P(J)*P(J)-20._DP*R(M,I)*&
					&B(M)*C(M)*P(J)
					QXY(NR,JC+6)=10._DP*B(I)*C(I)*B(M)*B(M)+4._DP*B(I)*B(I)*C(M)*B(M)+&
					&10._DP*R(M,I)*B(M)*B(M)*P(J)
!
!
				case(3)
					QXX(NR,JC+1)=60._DP*R(M,I)*B(I)*B(I)
					QXX(NR,JC+2)=-(-1._DP)**IPER*(6._DP*B(I)*B(I)*C(I)+30._DP*R(M,I)*B(I)*&
					&B(I)*C(M))
					QXX(NR,JC+3)=(-1._DP)**IPER*(6._DP*B(I)*B(I)*B(I)+30._DP*R(M,I)*B(I)*&
					&B(I)*B(M))
					QXX(NR,JC+4)=2._DP*B(I)*B(I)*C(I)*C(M)+5._DP*R(M,I)*B(I)*B(I)*&
					&C(M)*C(M)
					QXX(NR,JC+5)=-2._DP*B(I)*B(I)*(P(J)+5._DP*R(M,I)*B(M)*C(M))
					QXX(NR,JC+6)=2._DP*B(I)*B(I)*B(I)*B(M)+5._DP*R(M,I)*B(I)*B(I)*B(M)*B(M)
!
!
					QYY(NR,JC+1)=60._DP*R(M,I)*C(I)*C(I)
					QYY(NR,JC+2)=-(-1._DP)**IPER*(6._DP*C(I)**3+30._DP*R(M,I)*C(I)*C(I)*C(M))
					QYY(NR,JC+3)=(-1._DP)**IPER*(6._DP*C(I)*C(I)*B(I)+30._DP*R(M,I)*C(I)*&
					&C(I)*B(M))
					QYY(NR,JC+4)=2._DP*C(I)*C(I)*C(I)*C(M)+5._DP*R(M,I)*C(I)*C(I)*C(M)*&
					&C(M)
					QYY(NR,JC+5)=-2._DP*C(I)*C(I)*(P(J)+5._DP*R(M,I)*C(M)*B(M))
					QYY(NR,JC+6)=2._DP*C(I)*C(I)*B(I)*B(M)+5._DP*R(M,I)*C(I)*C(I)*B(M)*&
					&B(M)
!
!
					QXY(NR,JC+1)=120._DP*R(M,I)*B(I)*C(I)
					QXY(NR,JC+2)=-(-1._DP)**IPER*(12._DP*B(I)*C(I)*C(I)+60._DP*R(M,I)*B(I)*&
					&C(I)*C(M))
					QXY(NR,JC+3)=(-1._DP)**IPER*(12._DP*C(I)*B(I)*B(I)+60._DP*R(M,I)*C(I)*&
					&B(I)*B(M))
					QXY(NR,JC+4)=4._DP*C(I)*C(I)*B(I)*C(M)+10._DP*R(M,I)*C(M)*C(M)*B(I)*&
					&C(I)
					QXY(NR,JC+5)=-4._DP*B(I)*C(I)*P(J)-20._DP*R(M,I)*B(I)*C(I)*B(M)*C(M)
					QXY(NR,JC+6)=4._DP*B(I)*B(I)*C(I)*B(M)+10._DP*R(M,I)*B(M)*B(M)*C(I)*&
					&B(I)
				case default
					istat=STAT_FAILURE
					emsg='K is diff. from 1,2 or 3, Qform, &
					& element_module'
					QXX=ZERO; QYY=ZERO; QXY=ZERO
					return
				end select
			end do iper_loop
		end do row_selec_loop
		
	end do main_loop
end subroutine Qform




pure subroutine Dplate_form(material,thickness,Dmat) ! Only used in linear part
!	Purpose:
!	Forms the D-matrix for a plate (relation between moments and curvatures)

!   Declaration
!	I/O
	type(ELI_material), intent(in)	:: material
	real(DP)	      , intent(in)	:: thickness
	real(DP)		  , intent(inout)	:: Dmat(NST,NST)
!	local variables
	real(DP)		  	::	E
	real(DP)		  	::  nu
	
!	Initialize and extract
	Dmat	= ZERO
	call extract(E, nu, material)
	
!	Build D-matrix
	
	Dmat(1,1) = (E*thickness**3)/(12._DP*(ONE-nu**2))
	Dmat(1,2) = nu*Dmat(1,1)
	Dmat(2,1) = Dmat(1,2)
	Dmat(2,2) = Dmat(1,1)
	Dmat(3,3) = (ONE-nu)/TWO*Dmat(1,1) 

end subroutine Dplate_form






pure subroutine Rform(Rmat)
!	Purpose:
!	Builds the R-matrix, part of the final K

	real(DP), intent(out)	::	Rmat(10,10)
	
!	Declare local variables
	real(DP)	::	mult
	real(DP)	::	R_T(size(Rmat,2),size(Rmat,1))
	real(DP)	::	diag(size(Rmat,1))
	
	integer		::  i,j !counters
	
!	Initialize
	RMAT	=	ZERO
	mult	=	ONE/3360._DP
	R_T		=	ZERO
	diag	=	ZERO
	
	i=0; j=0
	
! 	Build R
	Rmat(1,2:) = (/THREE, THREE, TEN, TEN, FOUR, ONE, FOUR, ONE, TWO/)
	Rmat(2,3:) = (/THREE, FOUR, ONE, TEN, TEN, ONE, FOUR, TWO/)
	Rmat(3,4:) = (/ONE, FOUR, ONE, FOUR, TEN, TEN, TWO/)
	Rmat(4,5:) = (/TWO, THREE, ONE, ONE, TWO_THIRD, ONE/)
	Rmat(5,6:) = (/ONE, TWO_THIRD, THREE, ONE, ONE/)
	Rmat(6,7:) = (/TWO, TWO_THIRD, ONE, ONE/)
	Rmat(7,8:) = (/ONE, THREE, ONE/)
	Rmat(8,9:) = (/TWO, ONE/)
	Rmat(9,10) = ONE
	
	diag 	   = (/SIXTY, SIXTY, SIXTY, FOUR, FOUR, FOUR, FOUR, FOUR, FOUR,&
			     & TWO_THIRD/)
			  
	R_T		   = transpose(Rmat)
	
	Rmat	   = Rmat+R_T
	
	do i=1, size(Rmat,1)
		Rmat(i,i)	=	diag(i)
	end do
	
	Rmat	   = mult*Rmat
	
end subroutine Rform






pure subroutine moments_point(coord_point,coord_nod,u,area,Qmat,Dmat,moments_vec)
!	Purpose:
!	Computes the three moments (Mxx,Myy,Mxy) at the given point.


!	Declare
!	I/O
	real(DP), intent(in)	::	coord_point(NDIM-1)
	real(DP), intent(in)	::	coord_nod(NDIM-1,NNODES)
	real(DP), intent(in)	::	u(NDOF) !nodal disp. vector
	real(DP), intent(in)	::	area
	real(DP), intent(in)	::	Qmat(30,NDOF)
	real(DP), intent(in)	::	Dmat(NST,NST)
	real(DP), intent(inout)	::	moments_vec(3)
	
!	locals
	real(DP)				:: curv(NST)
	real(DP)				:: xp !point x
	real(DP)				:: yp !point y
	real(DP)				:: Fmat(3,30)
	real(DP)				:: Bmat(NST,NDOF)
	real(DP)				:: f_shape(NDOF)
	real(DP)				:: Lvec(3) !contains the area coordinates
	
	
!	Initialize
	curv=ZERO
	xp=ZERO
	yp=ZERO
	Fmat=ZERO
	Bmat=ZERO
	f_shape=ZERO
	Lvec=ZERO

!	Local copy of the query point coordinates
	xp=coord_point(1)
	yp=coord_point(2)
	
!	Compute the the shape functions and the area coordinates
	call init_shape_areacoords((/xp,yp/), coord_nod, f_shape, Lvec)
		
!	Compute the F-matrix 
	call Fform(Lvec(1),Lvec(2),Lvec(3),Fmat)
	
!	Use F and Q to compute the B matrix
	Bmat=ONE/(FOUR*area**2)*matmul(Fmat,Qmat)
	
!	Compute the curatures at the point
	curv=matmul(Bmat,u)
	
!	Compute the moments at the point (output)
	moments_vec=-matmul(Dmat,curv)
	

end subroutine moments_point



pure subroutine sdvs_centroid(coord_nod,u,area,Qmat,Dmat,element)
!	Purpose:
!	The components of element (displacement, curvatures and moments) are
!	computed and are then transferred in the svars() in the uel


!	Declare
!	I/O
	real(DP), intent(in)	::	coord_nod(NDIM-1,NNODES)
	real(DP), intent(in)	::	u(NDOF) !nodal disp. vector
	real(DP), intent(in)	::	area
	real(DP), intent(in)	::	Qmat(30,NDOF)
	real(DP), intent(in)	::	Dmat(NST,NST)
	type(TUBA3_element), intent(inout)	::	element
	
!	locals
	real(DP)				:: w
	real(DP)				:: curv(NST)
	real(DP)				:: mom(NST)
	real(DP)				:: xn(NNODES) !nodal x
	real(DP)				:: yn(NNODES) !nodal y
	real(DP)				:: xG	!centroid x
	real(DP)				:: yG	!centroid y
	real(DP)				:: Fmat(3,30)
	real(DP)				:: Bmat(NST,NDOF)
	real(DP)				:: f_shape(NDOF)
	real(DP)				:: Lvec(3) !contains the area coordinates
	
	
!	Initialize
	w=ZERO
	curv=ZERO
	mom=ZERO
	xn=ZERO
	yn=ZERO
	Fmat=ZERO
	Bmat=ZERO
	f_shape=ZERO
	Lvec=ZERO
	
	xn=coord_nod(1,:)
	yn=coord_nod(2,:)
	
	
!	Compute the centroid coordinates
	xG=ONE/THREE*(xn(1)+xn(2)+xn(3))
	yG=ONE/THREE*(yn(1)+yn(2)+yn(3))
	
!	Compute the the shape functions and the area coordinates
	call init_shape_areacoords((/xG,yG/), coord_nod, f_shape, Lvec)
	
!	Displacement at the centroid
	w=dot_product(f_shape,u)
	
!	Compute the F-matrix 
	call Fform(Lvec(1),Lvec(2),Lvec(3),Fmat)
	
!	Use F and Q to compute the B matrix
	Bmat=ONE/(FOUR*area**2)*matmul(Fmat,Qmat)
	
!	Compute the curatures at the centroid
	curv=matmul(Bmat,u)
	
!	Compute the moments at the centroid
	mom=-matmul(Dmat,curv)
	
!	Substute local variables into the components of 'element'
	element%w=w
	element%curv=curv
	element%mom=mom

end subroutine sdvs_centroid




! list of private subroutines


pure subroutine init_shape_areacoords(coord_point,coord_nod,f,L) ! used in  both  linear  and  nonlinear parts
!	Purpose:
!	evaluates the shape functions for TUBA3 and the area coordinates at the given point

!	Declararion
!	I/O
	real(DP),	intent(in)		::	coord_point(NDIM-1) !coords on mid-surface
	real(DP),	intent(in)		::	coord_nod(NDIM-1,NNODES) !coords on mid-surface
	real(DP),	intent(out)		::	f(NDOF)
	real(DP),	intent(inout)	::  L(3)
	
!	locals
	real(DP)				::  x, y
	real(DP)				::  xn(NNODES), yn(NNODES)
	real(DP)				::	a(NNODES), b(NNODES), c(NNODES)
	real(DP)				::  Delta
	real(DP)				::  r(NNODES,NNODES)
	
	integer					::  iter, i, j, k
	integer					::	temp1, temp2
	
	
!	Initialize
	x	=	coord_point(1)
	y	=	coord_point(2)
	
	xn	=	coord_nod(1,:)
	yn	=	coord_nod(2,:)
	
	a	=	ZERO
	b	=	ZERO
	c	=	ZERO
	r	=	ZERO
	
	iter  = 0
	i 	  =	0
	j 	  =	0
	k 	  =	0
	temp1 =	0
	temp2 =	0
	
	a(1)=xn(2)*yn(3)-xn(3)*yn(2)
	a(2)=xn(3)*yn(1)-xn(1)*yn(3)
	a(3)=xn(1)*yn(2)-xn(2)*yn(1)
	b(1)=yn(2)-yn(3)
	b(2)=yn(3)-yn(1)
	b(3)=yn(1)-yn(2)
	c(1)=xn(3)-xn(2)
	c(2)=xn(1)-xn(3)
	c(3)=xn(2)-xn(1)
	do i=1,NNODES
		do j=1,NNODES
			r(i,j)=-(b(i)*b(j)+c(i)*c(j))/(b(i)*b(i)+c(i)*c(i))
		end do
	end do
	
	
	Delta=HALF*(b(1)*c(2)-b(2)*c(1))
	
	L(1)=ONE/(TWO*Delta)*(a(1)+b(1)*x+c(1)*y)
	L(2)=ONE/(TWO*Delta)*(a(2)+b(2)*x+c(2)*y)
	L(3)=ONE/(TWO*Delta)*(a(3)+b(3)*x+c(3)*y)
	
	
!	Build the shape functions
	i	=	1
	j	=	2
	k	=	3
	
	do iter=1,3
	
		f(NDOF_NODE*(iter-1)+1)	=	L(i)**5 + 5._DP*L(i)**4*L(j) + 5._DP*L(i)**4*L(k)&
		& + 10._DP*L(i)**3*L(j)**2 + 10._DP*L(i)**3*L(k)**2 + 20._DP*L(i)**3*L(j)*L(k)&
		& + 30._DP*r(j,i)*L(i)**2*L(j)*L(k)**2 + 30._DP*r(k,i)*L(i)**2*L(j)**2*L(k)
		
		f(NDOF_NODE*(iter-1)+2)	=	c(k)*L(i)**4*L(j) - c(j)*L(i)**4*L(k)&
		& + 4._DP*c(k)*L(i)**3*L(j)**2 - 4._DP*c(j)*L(i)**3*L(k)**2&
		& + 4._DP*(c(k)-c(j))*L(i)**3*L(j)*L(k) - (3._DP*c(i)+15._DP*r(j,i)*c(j))*L(i)**2*L(j)*L(k)**2&
		& + (3._DP*c(i)+15._DP*r(k,i)*c(k))*L(i)**2*L(j)**2*L(k)
		
		f(NDOF_NODE*(iter-1)+3)	=	-b(k)*L(i)**4*L(j) + b(j)*L(i)**4*L(k)&
		& - 4._DP*b(k)*L(i)**3*L(j)**2 + 4._DP*b(j)*L(i)**3*L(k)**2&
		& + 4._DP*(b(j)-b(k))*L(i)**3*L(j)*L(k) + (3._DP*b(i)+15._DP*r(j,i)*b(j))*L(i)**2*L(j)*L(k)**2&
		& - (3._DP*b(i)+15._DP*r(k,i)*b(k))*L(i)**2*L(j)**2*L(k)
		
		f(NDOF_NODE*(iter-1)+4)	=	HALF*c(k)**2*L(i)**3*L(j)**2 + HALF*c(j)**2*L(i)**3*L(k)**2&
		& - c(j)*c(k)*L(i)**3*L(j)*L(k) + (c(i)*c(j)+FIVE_TWO*r(j,i)*c(j)**2)*L(i)**2*L(j)*L(k)**2&
		& + (c(i)*c(k)+FIVE_TWO*r(k,i)*c(k)**2)*L(i)**2*L(j)**2*L(k)
		
		f(NDOF_NODE*(iter-1)+5)	=	-b(k)*c(k)*L(i)**3*L(j)**2 - b(j)*c(j)*L(i)**3*L(k)**2&
		& + (b(j)*c(k)+b(k)*c(j))*L(i)**3*L(j)*L(k)&
		& - (b(i)*c(j)+b(j)*c(i)+5._DP*r(j,i)*b(j)*c(j))*L(i)**2*L(j)*L(k)**2&
		& - (b(i)*c(k)+b(k)*c(i)+5._DP*r(k,i)*b(k)*c(k))*L(i)**2*L(j)**2*L(k)
		
		f(NDOF_NODE*(iter-1)+6)	=	HALF*b(k)**2*L(i)**3*L(j)**2 + HALF*b(j)**2*L(i)**3*L(k)**2&
		& - b(j)*b(k)*L(i)**3*L(j)*L(k) + (b(i)*b(j)+FIVE_TWO*r(j,i)*b(j)**2)*L(i)**2*L(j)*L(k)**2&
		& + (b(i)*b(k)+FIVE_TWO*r(k,i)*b(k)**2)*L(i)**2*L(j)**2*L(k)
	
		temp1=j
		temp2=k
		k=i
		i=temp1
		j=temp2
		
		temp1=0
		temp2=0
	end do
	

end subroutine init_shape_areacoords




pure subroutine Fform(L1,L2,L3,Fmat) ! used to calculate [B]Liner
!	Purpose:
!	evaluates the F-matrix for TUBA3, given the area coords evaluated in the query
!	point

!	Declaration
!	I/O
	real(DP), intent(in)	::	L1, L2, L3
	real(DP), intent(out)	::	Fmat(3,30)
	
!	locals
	real(DP)				::	Lvec(10)
	

!	Initialize
	Fmat = ZERO
	Lvec = ZERO
	
	
!	Build Lvec
	Lvec = (/L1**3, L2**3, L3**3, L1**2*L2, L1**2*L3, L2**2*L1, L2**2*L3,&
	& L3**2*L1, L3**2*L2, L1*L2*L3/)

!	Build F-matrix
	Fmat(1,:10)=Lvec
	Fmat(2,11:20)=Lvec
	Fmat(3,21:30)=Lvec

end subroutine Fform






end module Linear_TUBA3_module