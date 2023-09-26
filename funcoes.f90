module funcoes

use omp_lib
use variaveis

!#################################################################
!#                             					 #
!#     			CREATING THE MESH               	 #
!#                             					 #
!#################################################################

contains

subroutine mesh

write(*,*) 'Criando a malha'
write(*,*) ''
deltax=(xmax-xmin)/(n-1)

deltay=(ymax-ymin)/(m-1)

x(1)=xmin
do i=2,n
x(i)=x(i-1)+deltax
end do

y(1)=ymin
do i=2,m
y(i)=y(i-1)+deltay
end do

end subroutine mesh
 
!#################################################################
!#                             					 #
!#     		ASSEMBLING THE IDENTIFICATION VECTOR           	 #
!#                             					 #
!#################################################################
 

subroutine id_vector

do k=1,nt

do i=1,n
do j=1,m
!******************CIRCULAR TUMOR**********************************
if(CIRCULAR == .TRUE.) then

if(x(i).ge.(xc(k)-(raio(k))).and.x(i).le.(xc(k)+(raio(k)))) then

! Quandrantes inferiores
if(y(j).le.yc(k)) then
rmin = yc(k) - ((raio(k)**2.0) - ((x(i) - xc(k))**2.0))**0.5
if(y(j).gt.rmin) then
ID(i+((j-1)*n))=1 ! Dentro do tumor
end if
else
! Quandrantes inferiores
rmin = yc(k) + (raio(k)**2.0 - (x(i) - xc(k))**2.0)**0.5
if(y(j).lt.rmin) then
ID(i+((j-1)*n))=1 ! Dentro do tumor
end if
end if
end if



!**********************ELIPTICAL TUMOR******************************

else

if(x(i).ge.a_xmin.and.x(i).le.a_xmax) then

! Quadrantes Superiores
if(y(j).ge.yc(k).and.y(j).le.b_ymax) then
if((y(j) - yc(k)).le.((b_ymax-yc(k))/(a_xmax-xc(k)))*(((a_xmax-xc(k))**2.0)-((x(i) - xc(k))**2.0))**0.5) then
ID(i+((j-1)*n)) = 1
end if
end if
! Quadrantes inferiores
if(y(j).le.yc(k).and.y(j).ge.b_ymin) then
if((yc(k) - y(j)).le.((b_ymax-yc(k))/(a_xmax-xc(k)))*(((a_xmax-xc(k))**2.0)-((x(i) - xc(k))**2.0))**0.5) then
ID(i+((j-1)*n)) = 1
end if
end if
else
ID(i+((j-1)*n)) = 0
end if
end if
end do
end do

end do

end subroutine id_vector
 
!#########################################################################
!#                                                                       #
!#     		   CALCULATING THE GOVERNING EQUATIONS                   #
!#                                                                       #
!#########################################################################
subroutine principal 

525 FORMAT(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)
T=T_ini


!#########################################################################
!#                                                                       #
!#     	      SOLVING THE MAGNETIC FIELD INSIDE THE TUMOUR               #
!#                                                                       #
!#########################################################################

!if(CALCCAMPO == .TRUE.) then

!write(*,*)'#########################################################################'
!write(*,*)'#                                                                       #'
!write(*,*)'#     	    SOLVING THE MAGNETIC FIELD INSIDE THE TUMOUR           #'
!write(*,*)'#                                                                       #'
!write(*,*)'#########################################################################'

!call camponaouniforme
!end if

if(.not.STEADY) then

write(*,*) 'BEGGINING OF THE SIMULATION'
write(*,*) ''

lreal=0.0
write(*,*) 'SIMULATION PROGRESS'
write(*,*) lreal,'%'

!#########################################################################
!#  		      BEGGINING OF THE TEMPORAL LOOP                     #
!#########################################################################


write(*,*)'#########################################################################'
write(*,*)'#  		      BEGGINING OF THE TEMPORAL LOOP                       #'
write(*,*)'#########################################################################'


! CREATING A FILE TO RECORD THE TEMPERATURE EVOLUTION AT THE CENTER OF 
! THE TUMOUR

open(2112,file='evolucao_temporal.plt')
write(2112,*)'Variables="time","Temperature"'
contagem=0

do l=1,npast !TIME MARCHING

if(EXPLICIT) then

!#########################################################################
!#                                                                       #
!#     	    SOLVING THE EQUATIONS THROUGH AN EXPLICIT SCHEME             #
!#                                                                       #
!#########################################################################



! CALCULATING THE BLOOD PERFUSION LOCALLY 

call perfusao_sanguinea 
!blood perfusion w(T) is updated each time-step as T(t) evolves

! NOW WE EVOLVE WITH AN EXPLICIT SCHEME THE TEMPERATURE IN EACH REGION

tc_antes=T((m*n-m)/2)
call solucao_explicita
tc_agora=T((m*n-m)/2)

if(abs(tc_antes-tc_agora).lt.1.0E-08)then
if(contagem.lt.1)then
timesteady=l*dt
contagem=contagem+1
end if
end if

!#########################################################################
!#                                                                       #
!#     		   	     WRITTING DATA FILES     	                 #
!#  	            (FOR THE EXPLICIT SOLUTION SCHEME)                   #
!#                                                                       #
!#########################################################################

AUXINT=l/(t_int/dt)
AUXREAL=l/(t_dec/dt)

lreal=l
write(*,*) (lreal/npastreal)*100 ,'%'
if(AUXINT.EQ.AUXREAL) then
write(rea_char, '(I3)') AUXINT
open (AUXINT,file='saida'//rea_char//'.plt')
!write(AUXINT,*) 'Variables="x","y","ID","H","Hx","Hy","W","T"'
write(AUXINT,*) 'Variables="x","y","W","T"'
write(AUXINT,*) 'ZONE F=POINT,I='
write(AUXINT,*) n
write(AUXINT,*) ',J='
write(AUXINT,*) m
do i=1,n
do j=1,m
!write(AUXINT,525)x(i),y(j),ID(i+((j-1)*n)),H(i+((j-1)*n)),HX(i+((j-1)*n)),HY(i+((j-1)*n)),W(i+((j-1)*n)),T(i+((j-1)*n))
write(AUXINT,'(F12.4,F12.4,F12.4,F12.4)') x(i),y(j),W(i+((j-1)*n)),T(i+((j-1)*n))
end do
end do
end if

auxiliarinteiro=l/(t_tempi/dt)
auxiliarreal=lreal/(t_tempr/dt)

if(auxiliarinteiro.eq.auxiliarreal) then
write(2112,'(F12.4,F12.4)') l*dt, T((m*n-m)/2)
end if

! END OF THE EXPLCIT SCHEME IMPLEMENTATION

else 

!#########################################################################
!#                                                                       #
!#     		   SOLUTION THROUGH THE IMPLICIT SCHEME                  #
!#                                                                       #
!#########################################################################
 
  write(*,*) 'SECTION UNDER CONSTRUCTION'

end if 

!#########################################################################
!#  		      WE ARE STILL IN A TEMPORAL LOOP                    #
!#########################################################################


! CLOSING TEMPORAL LOOP
end do

else

!#########################################################################
!#                                                                       #
!#     	   		   STEADY STATE SOLUTION            		 #
!#                                                                       #
!#########################################################################

 call solucao_permanente

end if 

end subroutine principal
!#########################################################################

!#########################################################################
!#                                                                       #
!#     	   	          AUXILIARY SUBROUTINES                          #
!#                                                                       #
!#########################################################################

!subroutine camponaouniforme
!#########################################################################
!#     	        NON-UNIFORM FIELD SOLUTION SUBROUTINE 			 #
!#########################################################################
! do k=2,m-1
!   do j=2,n-1
! Quadrante 1
!if(y(k).le.((yc/xc)*x(j)).and.y(k).le.((yc/(xc-xmax))*(x(j)-xmax))) then
!ylinha=0.0
!xlinha = xc + ((ylinha-yc)*(x(j)-xc)/(y(k)-yc))
!ylinha= yc - sqrt(raiomaximo**2.0 - (xlinha-xc)**2.0)
!H(k+((j-1)*n)) = HZERO*sin(omega*l*dt)
!HY(k+((j-1)*n)) = H(k+((j-1)*n))*(xlinha-xc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!HX(k+((j-1)*n)) = H(k+((j-1)*n))*(ylinha-yc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!end if
! Quadrante 2
!if(y(k).gt.((yc/xc)*x(j)).and.y(k).le.((yc/(xc-xmax))*(x(j)-xmax))) then
!xlinha=0.0
!ylinha = yc + ((xlinha-xc)*(y(k)-yc)/(x(j)-xc))
!xlinha= xc - sqrt(raiomaximo**2.0 - (ylinha-yc)**2.0)
!H(k+((j-1)*n)) = HZERO*sin(omega*l*dt)
!HY(k+((j-1)*n)) = H(k+((j-1)*n))*(xlinha-xc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!HX(k+((j-1)*n)) = H(k+((j-1)*n))*(ylinha-yc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!end if
! Quadrante 3
!if(y(k).gt.((yc/xc)*x(j)).and.y(k).gt.((yc/(xc-xmax))*(x(j)-xmax))) then
!ylinha=ymax
!xlinha = xc + ((ylinha-yc)*(x(j)-xc)/(y(k)-yc))
!ylinha= yc + sqrt(raiomaximo**2.0 - (xlinha-xc)**2.0)
!H(k+((j-1)*n)) = HZERO*sin(omega*l*dt)
!HY(k+((j-1)*n)) = H(k+((j-1)*n))*(xlinha-xc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!HX(k+((j-1)*n)) = H(k+((j-1)*n))*(ylinha-yc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!end if
! Quadrante 4
!if(y(k).le.((yc/xc)*x(j)).and.y(k).gt.((yc/(xc-xmax))*(x(j)-xmax))) then
!xlinha=xmax
!ylinha = yc + ((xlinha-xc)*(y(k)-yc)/(x(j)-xc))
!xlinha= xc + sqrt(raiomaximo**2.0 - (ylinha-yc)**2.0)
!H(k+((j-1)*n)) = HZERO*sin(omega*l*dt)
!HY(k+((j-1)*n)) = H(k+((j-1)*n))*(xlinha-xc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!HX(k+((j-1)*n)) = H(k+((j-1)*n))*(ylinha-yc)/sqrt((xc-xlinha)**2.0 + (yc-ylinha)**2.0)
!end if
!end do
!end do
!end subroutine camponaouniforme
!#########################################################################

subroutine perfusao_sanguinea
!#########################################################################
!#     	        BLOOD PERFUSION RATE CALCULATION SUBROUTINE              #
!#########################################################################

do k=2,m-1 ! SWEEPS IN "Y" DIRECTION
do j=2,n-1 !SWEEPS IN "X" DIRECTION

! HEALTHY TISSUE REGION
if(ID(k+((j-1)*n)).eq.0) then
! FOR T < 37°C
if(T(k+((j-1)*n)).lt.37.0) then
W(k+((j-1)*n)) = 8.33E-04
end if
! FOR 37°C < T < 42°C
if(T(k+((j-1)*n)).ge.37.0.and.T(k+((j-1)*n)).le.42.0) then
W(k+((j-1)*n)) = 8.33E-04 - ((T(k+((j-1)*n)) - 37.0)**4.8)/(5.438E+06)
end if
! FOR T > 42°C
if(T(k+((j-1)*n)).lt.37.0) then
W(k+((j-1)*n)) = 4.16E-04
end if
end if

! TUMOUR TISSUE REGION
if(ID(k+((j-1)*n)).eq.1) then
! FOR T < 45°C
if(T(k+((j-1)*n)).lt.45.0) then
W(k+((j-1)*n)) = 4.5E-04 + 3.55E-03*EXP(-((T(k+((j-1)*n))-45.0)**2.0)/12.0)
end if
! FOR T >= 45°C
if(T(k+((j-1)*n)).ge.45.0) then
W(k+((j-1)*n)) = 4.00E-03
end if
end if
end do
end do

! WORING THE BOUNDARIES

! LOWER BOUNDARY
do j=1,n
i=j
W(i) = 8.33E-04
end do
! LEFT BOUNDARY
do j=1,m-2
i=j*n + 1
W(i) = 8.33E-04
end do
! RIGHT BOUNDARY 
do j=1,m-2
i=(j+1)*n 
W(i) = 8.33E-04
end do
! UPPER BOUNDARY 
do j=1,n
i=(n*(m-1)) + j
W(i) = 8.33E-04
end do

end subroutine perfusao_sanguinea
!#########################################################################

subroutine solucao_explicita

!#########################################################################
!#     	                 EXPLICIT SOLUTION SUBROUTINE                    #
!#########################################################################

do k=2,m-1 ! SWEEPS IN "Y" DIRECTION
do j=2,n-1 !SWEEPS IN "X" DIRECTION

! HEALTHY TISSUE REGION

if(ID(k+((j-1)*n)).eq.0) then 
DIFUSAO = (K1/(deltax**2.0))*(T(k+((j-1)*n)+1) + T(k+((j-1)*n)-1)-2.0*T(k+((j-1)*n))) + (K1/(deltay**2.0))*(T(k+((j-1)*n)+n) + T(k+((j-1)*n)-n)-2.0*T(k+((j-1)*n)))
PERFUSAO = RHOB*CB*W(k+((j-1)*n))*(TB-T(k+((j-1)*n)))
T(k+((j-1)*n)) = T(k+((j-1)*n)) + (dt/(RHO1*C1))*(DIFUSAO + PERFUSAO + Q1)  
end if
TIU= omp_get_num_threads() !Comando para contar a quantidade de threads em processo

! INSIDE THE TUMOR

if(ID(k+((j-1)*n)).eq.1) then 
DIFUSAO = (K2/(deltax**2.0))*(T(k+((j-1)*n)+1) + T(k+((j-1)*n)-1)-2.0*T(k+((j-1)*n))) + (K2/(deltay**2.0))*(T(k+((j-1)*n)+n) + T(k+((j-1)*n)-n)-2.0*T(k+((j-1)*n)))
PERFUSAO = RHOB*CB*W(k+((j-1)*n))*(TB-T(k+((j-1)*n))) !BLOOD PERFUSION TERM
!cor=exp(-((y(k)-yc)**2.0 + (x(j)-xc)**2.0)/((raio)**2.0))
if(QUIDUASLINHAS) then
GERACAO = muzero*acos(-1.0)*PHI*MD*HZERO*omega*ximag !MAGNETIC GENERATION TERM
else
GERACAO = pm !MAGNETIC GENERATION TERM
end if
T(k+((j-1)*n)) = T(k+((j-1)*n)) + (dt/(RHO2*C2))*(DIFUSAO + PERFUSAO + Q2 + cor*GERACAO)     
end if

end do ! CLOSING THE "J" LOOP
end do ! CLOSING THE "K" LOOP

! WORING THE BOUNDARIES

! LOWER BOUNDARY
do j=1,n
i=j
T(i) = 2.0*T(i+n) - T(i+(2*n))
end do
! LEFT BOUNDARY
do j=1,m-2
i=j*n + 1
T(i) = 2.0*T(i+1) - T(i+2) 
end do
! RIGHT BOUNDARY 
do j=1,m-2
i=(j+1)*n 
T(i) = 2.0*T(i-1) - T(i-2)
end do
! UPPER BOUNDARY 
do j=1,n
i=(n*(m-1)) + j
T(i) = 2.0*T(i-n) - T(i-(2*n))
end do

end subroutine solucao_explicita
!#########################################################################

subroutine solucao_permanente
real C1,C2,C3, tol
integer iter

!#########################################################################
!#     	        STEADY STATE SOLUTION SUBROUTINE 			 #
!#########################################################################

iter =1
tol =1.0E-3
525 FORMAT(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)
if(iter.eq.1) then
T_ANTES=TZERO
end if

! SOLVE THE BLOOD PERFUSION FIELD
180 call perfusao_sanguinea

T_ANTES=T

! BUILDING OUR LINEAR SYSTEM

do k=1,m-2
do j=1,n-2

i=(n*k) + 1 + j

! CALCULATING THE COEFFICIENTS FOR DIFERENT REGIONS

! HEALTHY TISSUE REGION

if(ID(k+((j-1)*n)).eq.0) then 
 C1 = K1/(deltax**2.0)
 C2 = K1/(deltay**2.0)
 C3 = 2.0/(deltax**2.0) + 2.0/(deltay**2.0) + RHOB*CB*W(k+((j-1)*n)) 
 
! ONLY METABOLIC HEAT GENERATION IN THE SOURCE TERM FOR HEALTHY REGION

 B(i) = Q1/C3
end if

! TUMOUR REGION

if(ID(k+((j-1)*n)).eq.1) then 
 C1 = K2/(deltax**2.0)
 C2 = K2/(deltay**2.0)
 C3 = 2.0/(deltax**2.0) + 2.0/(deltay**2.0) + RHOB*CB*W(k+((j-1)*n))  

! METABOLIC AND MAGNETIC GENERATION IN THE SOURCE TERM FOR TUMOUR REGION

if (QUIDUASLINHAS) then
 B(i) = (Q2 + muzero*acos(-1.0)*PHI*MD*HZERO*omega*ximag*cor)/C3  
else
 B(i) = (Q2 + pm*cor)/C3  
end if
end if

! ASSEMBLYING THE COEFFICIENT'S MATRIX FOR INTERNAL NODES	

A(i,i) = 1.0
A(i,i+1) = -C1/C3
A(i,i-1) = -C1/C3
A(i,i+n) = -C2/C3
A(i,i-n) = -C2/C3	
end do
end do

! DEALING WITH THE BOUNDARIES

! LOWER BOUNDARY
do j=1,n
i=j
A(i,i) = 1.0
B(i) = TZERO
end do

! LEFT BOUNDARY
do j=1,m-2
i=j*n + 1
A(i,i) = 1.0
B(i) = TZERO
end do

! RIGHT BOUNDARY
do j=1,m-2
i=(j+1)*n 
A(i,i) = 1.0
B(i) = TZERO
end do

! UPPER BOUNDARY 
do j=1,n
i=(n*(m-1)) + j
A(i,i) = 1.0
B(i) = TZERO
end do

! SOLVING THE LINEAR SYSTEM

 call jacobi(A,T,B,n*m,500)
! call seidel(A,T,B,n*m,500)
 
! call la_gesv(A, B) !LAPACK
! EVOLING ITERATION PARAMETER

iter=iter+1
T_AGORA=T

write(*,*) 'Iteration/Error:',iter, maxval(abs(T_AGORA-T_ANTES))

! CHECKING FOR TEMPERATURE CONVERGENCE

if(maxval(abs(T_AGORA-T_ANTES)).gt.tol) then
go to 180
end if

! WRITING DATA FILE
open (512,file='saida.plt')
!write(512,*) 'Variables="x","y","ID","H","Hx","Hy","W","T"'
write(512,*) 'Variables="x","y","ID","W","T"'
write(512,*) 'ZONE F=POINT,I='
write(512,*) n
write(512,*) ',J='
write(512,*) m

do i=1,n
do j=1,m
!write(512,525)x(k),y(j),ID(k+((j-1)*n)),H(k+((j-1)*n)),HX(k+((j-1)*n)),HY(k+((j-1)*n)),W(k+((j-1)*n)),T(k+((j-1)*n))
write(512,'(F12.4,F12.4,F12.4,F12.4,F12.4)') x(i),y(j),ID(i+((j-1)*n)),W(i+((j-1)*n)),T(i+((j-1)*n))
end do
end do


! NOTE: HERE WE MUST DO SOME ITERATIONS TO UPDATE W(T) UNTIL T CONVERGES

end subroutine solucao_permanente
!#########################################################################

subroutine jacobi(A,x,b,n,k)
!#########################################################################
!#      LINEAR SYSTEM SOLUTION BY THE JACOBI METHOD SUBROUTINE 		 #
!#########################################################################
implicit none
integer i,j,n,iter,k
real A(n,n), C(n,n), h(n), l(n)
real x(n),b(n),f(n),g

iter=0

do i=1,n
x(i)=1.0
end do

100 do j=1,n
	do i=1,n
	if(j.ne.i) then
	f(i)=(-A(i,j)*x(j))
	end if
	end do
	end do

	g=sum(f)

	do i=1,n
	x(i)=(b(i)+g)/A(i,i)
	end do

	iter=iter+1
	
	do i=1,n
	f(i)=0
	end do

	g=0

	if(iter.ne.k) then
	go to 100
	end if

	iter=0

	do i=1,n
	h(i)=b(i)/A(i,i)
	C(i,i)=0
	end do

	do i=1,n
	do j=1,n
	if(i.ne.j) then
	C(i,j)=-A(i,j)/A(i,i)
	end if
	end do
	end do

200	l=matmul(C,x)

	do i=1,n
	x(i)=l(i)+h(i)
	end do

	iter=iter+1

	if(iter.ne.k) then
	go to 200
	end if
end subroutine jacobi
!#########################################################################

subroutine seidel(A,x,b,n,k)
!#########################################################################
!#     LINEAR SYSTEM SOLUTION BY THE GAUSS-SEIDEL METHOD SUBROUTINE 	 #
!#########################################################################
implicit none
integer i,j,n,iter,k,saida
real*16 A(n,n), C(n,n), D(n,n)
real*16 x(n),b(n),f(n),g(n),y(n)

iter=0

do i=1,n
x(i)=b(i)/a(i,i)
g(i)=0
f(i)=0
end do

500 do i=2,n
    do j=1,i-1
       C(i,j)=A(i,j)*x(j)
    end do
    end do

f=sum(C,dim=2)


do i=1,n
x(i)=(b(i)-f(i)-g(i))/A(i,i)
end do


do i=1,n
do j=1,n
f(i)=0
g(i)=0
 C(i,j)=0
 D(i,j)=0
end do
end do


do i=1,n-1
do j=i+1,n
       D(i,j)=A(i,j)*x(j)
end do
end do

g=sum(D,dim=2)


iter=iter+1

if(iter.ne.k) then
go to 500
end if

y=matmul(x,A)

end subroutine seidel
!#########################################################################

end module funcoes 
