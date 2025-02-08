module functions

use omp_lib
use variables
contains


!##########################################
!#                             					  #
!#     			READING INPUTS		            #
!#                             					  #
!##########################################

 subroutine inputs 
!529 FORMAT(E9.2E2,1x,E9.2E2,1x,E9.2E2,1x,E9.2E2)  
 open(15,file='inputs.dat')
 do i=1,Nrea
 read(15,*) PHI(i), HZERO(i), omega(i), raio_part(i) 
 end do
 close(15)
 end subroutine inputs

!##########################################
!#                             			      #
!#     			CREATING THE MESH             #
!#                             					  #
!##########################################

subroutine mesh

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
 
!#########################################
!#                             					 #
!#     		ASSEMBLING THE ID VECTOR     	 #
!#                             					 #
!#########################################
 

subroutine id_vector

do i=1,n
do j=1,m

if(x(i).ge.(xc-(raio)).and.x(i).le.(xc+(raio))) then

! Quandrantes inferiores
if(y(j).le.yc) then
rmin = yc - ((raio**2.0) - ((x(i) - xc)**2.0))**0.5
if(y(j).gt.rmin) then
ID(i+((j-1)*n))=1 ! Dentro do tumor
end if
else
! Quandrantes inferiores
rmin = yc + (raio**2.0 - (x(i) - xc)**2.0)**0.5
if(y(j).lt.rmin) then
ID(i+((j-1)*n))=1 ! Dentro do tumor
end if
end if

end if

end do
end do

end subroutine id_vector
 
!##########################################
!#                                        #
!#   SOLVING GOVERNING EQUATIONS          #
!#                                        #
!##########################################
subroutine main

525 FORMAT(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)
T=T_ini


lreal=0.0


!##########################################
!#    BEGGINING OF THE TIME LOOP          #
!##########################################


write(*,*)'#########################################'
write(*,*)'#                             				   #'
write(*,*)'#   		  BEGGINING OF THE SIMULATIONS   #'
write(*,*)'#                             				   #'
write(*,*)'#########################################'

! CREATING A FILE TO RECORD THE TEMPERATURE EVOLUTION AT THE CENTER OF 
! THE TUMOUR

open(2112,file='outputs.dat')
write(2112,*)'Variables="SIMULATION NUMBER","PHI","H","W","PARTICLE_RADIUS","TC","TS"'

do p=1,Nrea
contagem(p)=0
call imag
write(*,*) 'PROGRESS OF SIMULATION NUMBER',p
do l=1,npast !TIME MARCHING

!##########################################
!#                                        #
!#    EXPLICIT SCHEME SOLUTION            #
!#                                        #
!##########################################

call perfusao_sanguinea 
T_ANTES=T
tc_antes=T((m*n-m)/2)
call solucao_explicita
tc_agora=T((m*n-m)/2)
T_AGORA=T

lreal=l
!(lreal/npastreal)*100,'%'
write(*,'(F12.4,a,F12.2,a)') l*dt,'s', T((m*n-m)/2),'°C'
!pause

if(T((m*n-m)/2).gt.55.0) then
contagem(p)=contagem(p)+1
T=TZERO
go to 102
end if

if(maxval(abs(T_ANTES-T_AGORA)).lt.1.0E-06)then
!write(*,*) abs(tc_antes-tc_agora) 
!pause
if(contagem(p).lt.1)then
timesteady(p)=l*dt
115 contagem(p)=contagem(p)+1
if(T((m*n-m)/2).le.80.0) then
write(2112,'(I12,F12.3,F12.1,F12.1,E12.2E2,F12.2,F12.2)') p, PHI(p), HZERO(p), omega(p), raio_part(p), T((m*n-m)/2), timesteady(p)
end if
T=TZERO
go to 102
end if
end if

!##########################################
!# 	     WRITTING DATA FILES              #
!#        (TECPLOT FORMAT)                #
!##########################################

!AUXINT=l/(t_int/dt)
!AUXREAL=l/(t_dec/dt)
!if(AUXINT.EQ.AUXREAL) then
!write(rea_char, '(I3)') AUXINT
!open (AUXINT,file='saida'//rea_char//'.plt')
!write(AUXINT,*) 'Variables="x","y","ID","H","Hx","Hy","W","T"'
!write(AUXINT,*) 'Variables="x","y","W","T"'
!write(AUXINT,*) 'ZONE F=POINT,I='
!write(AUXINT,*) n
!write(AUXINT,*) ',J='
!write(AUXINT,*) m
!do i=1,n
!do j=1,m
!write(AUXINT,525)x(i),y(j),ID(i+((j-1)*n)),H(i+((j-1)*n)),HX(i+((j-1)*n)),HY(i+((j-1)*n)),W(i+((j-1)*n)),T(i+((j-1)*n))
!write(AUXINT,'(F12.4,F12.4,F12.4,F12.4)') x(i),y(j),W(i+((j-1)*n)),T(i+((j-1)*n))
!end do
!end do
!end if

auxiliarinteiro=l/(t_tempi/dt)
auxiliarreal=lreal/(t_tempr/dt)

if(auxiliarinteiro.eq.auxiliarreal) then
call imag
end if

! CLOSING TEMPORAL LOOP
end do 

102 write(*,*) 'END OF SIMULATION NUMBER',p
end do

end subroutine main
!#########################################

!#########################################
!#                                       #
!#        AUXILIARY SUBROUTINES          #
!#                                       #
!#########################################


!#########################################
!#   BLOOD PERFUSION RATE CALCULATION    #
!#########################################

subroutine perfusao_sanguinea

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

! WRITING THE BOUNDARIES

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
!#########################################

subroutine solucao_explicita

!#########################################
!#    EXPLICIT SOLUTION SUBROUTINE       #
!#########################################

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
GERACAO = muzero*acos(-1.0)*PHI(p)*MD*HZERO(p)*omega(p)*(-ximag) !MAGNETIC GENERATION TERM
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

!#########################################

subroutine randomica(a,b,c,n,d)
real a,b       ! a,b = random number range
integer n, m 
real c(n)      ! c = generated random sequence
integer d,i,e
integer f(8)
integer, allocatable :: seed(:)

 call random_seed(size = m)
allocate (seed(m))

 CALL DATE_AND_TIME(values=f)
 CALL SYSTEM_CLOCK(count=e)

do i = 1,m
seed(i) =  47*d + f(8)*i*d*12000 + e*(3*d+i)
end do

 call random_seed(put = seed)

 call random_number(c)

 c = a+(b-a)*c

end subroutine randomica


!#########################################
!# IMAGINARY PART OF THE CHI CALCULATION #
!#########################################


subroutine imag
implicit none
double precision :: a1, b1, h1
integer :: M1,i1,n1,j1
double precision, allocatable :: y1(:),t1(:)
double precision :: k11,k22,k33,k44
double precision :: omega1, lambda1, alpha1

!**** dados de entrada1

omega1=(2.0d0*pi*omega(p)*6.0d0*visc*pi*raio_part(p)**3)/((1.38E-23)*(T((m*n-m)/2)+273))

alpha1=(4.0d0*pi*(raio_part(p)**3)*muzero*MD*HZERO(p))/(3.0d0*(1.38E-23)*(T((m*n-m)/2)+273))

lambda1=(muzero*(MD**2)*pi*raio_part(p)**3)/(18.0d0*(1.38E-23)*(T((m*n-m)/2)+273))

!**** dados de entrada2
pi = dacos(-1.d0)
a1=0.0d0
b1=(2.0d0*pi/omega1)*100
M1=10000.d0
!*** passo
h1=(b1-a1)/(1.0d0*M1)

!aloca variaveis na memoria
allocate(y1(0:M1),t1(0:M1))

!*** condições de contorno
y1(0)=0.0000001d0  
!****
t1(0)=a1
do i1=1,M1
t1(i1)=a1+i1*h1
end do

!**** RUNGE-KUTTA 4 ORDEM
do i1=0,M1-1
k11 = h1*f1(t1(i1),y1(i1))
k22 = h1*f1(t1(i1)+0.5d0*h1,y1(i1)+0.5d0*k11)
k33 = h1*f1(t1(i1)+0.5d0*h1,y1(i1)+0.5d0*k22)
k44 = h1*f1(t1(i1)+h1,y1(i1)+k33)
y1(i1+1) = y1(i1)+(1.0d0/6.0d0)*(k11+2.0d0*k22+2.0d0*k33+k44) 
end do

!#########################################

!Cálculo da parte imaginária X"(w)

ximag = 0.0

do j1=1,M1-1
  h1 = t1(j1+1)-t1(j1)
  ximag = ximag + (h1/2.0d0)*(omega1/(2.0d0*pi*100))*(((1.0d0/dtanh(y1(j1+1))-1.0d0/y1(j1+1)+&
  8.0d0*PHI(p)*lambda1*(1.0d0/(dtanh(y1(j1+1)))-1.0d0/(y1(j1+1)))*((1.0d0/y1(j1+1))**2-(1.0d0/dsinh(y1(j1+1)))**2)))*dcos(omega1*t1(j1+1))+&
  ((1.0d0/dtanh(y1(j1))-1.0d0/y1(j1)+8.0d0*PHI(p)*lambda1*(1.0d0/(dtanh(y1(j1)))-1.0d0/(y1(j1)))*((1.0d0/y1(j1))**2-&
  (1.0d0/dsinh(y1(j1)))**2)))*dcos(omega1*t1(j1)))
end do

!#########################################

!Função f = d(alpha_e)/dt para o cálculo de alpha_e

contains 
function f1(var1,var2) ! Aqui f = dy/dt !f(x,y)
implicit none
double precision :: f1, var1,var2
f1=((-2.0d0)*(0.75/var2)*(+1.0*var2-alpha1*dsin(omega1*var1))*(1.0d0/DTANH(var2)-1.0d0/var2))/((1.00)*(1.d0/var2**2-&
               1.0d0/DSINH(var2)**2+ & 
               8.0d0*PHI(p)*lambda1*((1.0d0/var2**2-1.0d0/DSINH(var2)**2)**2+(1.0d0/DSINH(var2)-1.0d0/var2)*(-2.0d0/var2**3+ &
                    2.0d0*(1.0d0/DTANH(var2))*(1.0d0/DSINH(var2)**2)))))
end function f1
!#########################################

end subroutine imag

!#########################################

end module functions
