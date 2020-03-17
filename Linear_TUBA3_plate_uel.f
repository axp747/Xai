!***************************************!
!   Abaqus UEL interface including      !    		!
!   TUBA3 linear plate element     		!
!   created by X.P.Ai, TU Delft         !
!   Date: 17-Mar-2020                   !
!                                       !
!***************************************!
!
!   Purpose:
!
!
!   Record of revision:
!    Date      Programmer            Description of change
!    ========  ====================  ========================================
!	 17/03/20  X. P. Ai		         Original code 
!
!
!
!
!
!
!------ include external modules --------------------------
include 'globals/parameter_module.f90'
include 'Linear_material_module.f90'
include 'Linear_TUBA3_module.f90'

subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
 &       props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
 &       kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
 &       lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
 
	
	
	
!   0-load modules (with parameters, functions, subroutines)	

	
	use parameter_module,           			 only:	DP, ZERO, ONE, TWO, THREE, EIGHT,&
	& MSGLENGTH, STAT_SUCCESS, STAT_FAILURE
	
	use material_module,		  			 	 only:  set, ELI_material

	use Linear_TUBA3_module,					 only:	Qform, Dplate_form,&
	& Rform, TUBA3_element, NDIM, PNDOF => NDOF, PNDOF_NODE => NDOF_NODE,&
	& PNNODES => NNODES, NST 
	

	
	
	
	include 'aba_param.inc'
	
	dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*)
    dimension svars(*),energy(8),coords(mcrd,nnode),u(ndofel)
    dimension du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*)
    dimension jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*)
    dimension predef(2,npredf,nnode),lflags(*),jprops(*)
	
	
!	1-Declaration	
	
!	Declare local copies of the uel inputs/outputs used
    real(DP) 		::	coordinates(NDIM-1,PNNODES) ! coordinates(2*3)
    real(DP) 		::	disp(PNDOF) ! 
	real(DP) 		::  Kmat_L(PNDOF,PNDOF)
	real(DP) 		::	Ddisp(PNDOF)    !local copy of the du (dof vector)
	real(DP) 		::	disp_old(PNDOF)    !local copy of the du (dof vector)
	real(DP)		::	Fvec(PNDOF) ! Fvec = Kmat * disp 
	
	
!	Declare other variables local to the subroutine
!   TUBA3 GLOBAL
!   elem : 
!   p_material: plate material
!   E :
!   nu :
!   thickness :
	type(TUBA3_element)			::	elem
	type(ELI_material)			::	p_material
	real(DP)					::	E,nu
	real(DP)					::	thickness
	integer						::	p_jelem_list(4) !number of entries depends on how many elements
									!one wants to query
									
									
!   Linear_TUBA3 PART	: used to calculate [K]L matrix
!   Kmat_L = f(AreaT,Dmat,Rmat,QXX,QYY,QXY)
!   QXX QYY QXY
!   Rmat
!   Dmat
!   AreaT
	real(DP)					::  QXX(10,PNDOF),QYY(10,PNDOF),QXY(10,PNDOF)
	real(DP)					::  Qmat(30,PNDOF)
	real(DP)					::  Rmat(10,10)
	real(DP)					::  AreaT ! Triangle aera
	real(DP)					::  Dmat(NST,NST) ! (3*3)
	
!	COMMON used 
!   x_query:
!
	real(DP)				    ::  x_query, y_query
    integer						::	p_clock
	integer						::	istatus
	character(len=MSGLENGTH)	::	error_msg
	
	integer						::	i !counter
	integer						::	j !counter	
!   Fold:	
	real(DP)					::	Fold(PNDOF)
	
!	2-Initialize amatrx and rhs
!   uel required outputs
	amatrx      = ZERO
	rhs(:,1)    = ZERO
	
!	3-Initialize all the local variables. Extremely high number is used to avoid printing
    p_jelem_list=(/10000000,10000000,10000000,10000000/)

!   3 props
	thickness=props(3) !mm
	E=props(1)	 !MPa
	nu=props(2)

    coordinates=ZERO
	disp=ZERO
	Ddisp=ZERO
	disp_old=ZERO
	
	AreaT=ZERO
	QXX=ZERO; QYY=ZERO; QXY=ZERO
	Qmat=ZERO
	Dmat=ZERO
	Rmat=ZERO
	Kmat=ZERO
	Fvec=ZERO
	Fold=ZERO

	p_clock=0
	
	istatus=STAT_SUCCESS
	error_msg=''
	
	i=0
	j=0	
	
!	4-Local copy of svars 
!	The local variables 2:end used as sdvs are the fields of the element type
!	The first svars is the analysis clock
	p_clock=svars(1)
	elem%w=svars(2)
	elem%curv=svars(3:5) ! 3,4,5
	elem%mom=svars(6:8)	! 6,7,8

!	5-Material setting
!   E=props(1)	 !MPa
!	nu=props(2)
	call set(p_material,E,nu)
	
!	6-Copy coords, u and du(:,1) (uel inputs) into local variables and set disp_old
	coordinates=coords(:NDIM-1,:)
	disp=u
	Ddisp=du(:,1)
	disp_old=disp-Ddisp
	
!---------------------START OF COMPUTATIONS--------------------
!	7-Build and print element area, K and the internal force vector
	call Qform(coordinates,AreaT,QXX,QYY,QXY,istatus,error_msg)
	call Dplate_form(p_material,thickness,Dmat)
	call Rform(Rmat)	
!   calculate Kmat_L
	Kmat=ONE/(EIGHT*AreaT**3)*(Dmat(1,1)*matmul(transpose(QXX),matmul(Rmat,QXX))&
	&+Dmat(1,2)*matmul(transpose(QXX),matmul(Rmat,QYY))&
	&+Dmat(1,2)*matmul(transpose(QYY),matmul(Rmat,QXX))&
	&+Dmat(2,2)*matmul(transpose(QYY),matmul(Rmat,QYY))&
	&+Dmat(3,3)*matmul(transpose(QXY),matmul(Rmat,QXY)))
	
	Fvec=matmul(Kmat,disp)
	
	
!	8-Compute w, k, M at the element's centroid
	Qmat(:10,:)  =QXX
	Qmat(11:20,:)=QYY
	Qmat(21:30,:)=QXY
	
	call sdvs_centroid(coordinates,disp,AreaT,Qmat,Dmat,elem)
		
!	9-The element is linear, hence K doesn't change with u. I can use u to compute the forces
!	at the last converged increment ([Fold]=[K]{u_old})
	Fold=matmul(Kmat,disp_old
	
!	10-Extract the forces at the nodes and print them in the std output
	if (((kstep==1) .and. (kinc==1)).or.(p_clock/=kinc*kstep)) then
		p_clock=kinc*kstep
		if (any(p_jelem_list==jelem)) then
			write(*,'(1X, A)') ''
			write(*,'(1X, A)') '*****New Increment*****'
			write(*,'(1X, A)') ''
			write(*,'(1X, A, I2, A, I4, A, I3)') 'Nstep: ', kstep,' Ninc: ',kinc,&
			&' plate elem no: ',jelem
			write(*,'(1X, A, F7.4)') 'Time: ',time(2)
			do i=1,PNNODES
				x_query=coordinates(1,i)
				y_query=coordinates(2,i)
		
				write(*,'(1X, A)') ''
				write(*,'(1X, A, 2F9.4)') 'Forces at LCI, for the node of coordinates: ', x_query, y_query
				write(*,'(1X, A, ES12.5)') 'F= ', Fold(1+PNDOF_NODE*(i-1))
			end do
		endif
	endif
!---------------------END OF COMPUTATIONS--------------------

!	11-Update svars from TUBA3_element components
	svars(1)=p_clock
	svars(2)=elem%w
	svars(3:5)=elem%curv
	svars(6:8)=elem%mom
	
!	12-Update amatrx and rhs(:,1)
	amatrx=Kmat
	rhs(:,1)=-Fvec

 end subroutine uel