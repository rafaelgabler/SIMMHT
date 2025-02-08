program simmht

use omp_lib


!##############################
!#                            #
!#    PROGRAM simmht          #
!#                            #
!##############################

! Wrote by Rafael Gabler at 04/03/2021
 
use variables
use functions

call input
call cpu_time(ti)

allocate(ID(m*n))
allocate(T(m*n))
allocate(T_ANTES(m*n))
allocate(T_AGORA(m*n))
allocate(W(m*n))
allocate(x(n))
allocate(y(m))
allocate(HZERO(Nrea))
allocate(raio_part(Nrea))
allocate(omega(Nrea))
allocate(PHI(Nrea))
allocate(timesteady(Nrea))
allocate(Tcenter(Nrea))
allocate(contagem(Nrea))


xc = (xmax-xmin)/2.0
yc = (ymax-ymin)/2.0

pi = acos(-1.0)
muzero = 4.0*pi*1.0E-07
HZERO=0.0
omega=0.0
PHI=0.0
raio_part=0.0
dt = 1.0E-03 !min(1.0E-03,1.0/omega)

npast = sim_time/dt
npastreal = npast

!**************************************************************************!
write(*,*)'#################################################'
write(*,*)'#                             				   #'
write(*,*)'#     		  READING INPUTS	               #'
write(*,*)'#                             				   #'
write(*,*)'################################################# '

if(fileinput)then
 call inputs
else
call randomica(1.5E+03,4.0E+03,HZERO,Nrea,9)
call randomica(1.0E+05,3.0E+05,omega,Nrea,8)
call randomica(3.0E-02,5.00E-02,PHI,Nrea,1)
call randomica(5.0E-09,8.0E-09,raio_part,Nrea,4)
open(unit=19,file='generated_inputs.dat')
do p=1,Nrea
write(19,'(I12,F12.3,F12.1,F12.1,E12.2E2)') p, PHI(p), HZERO(p), omega(p), raio_part(p)
end do
end if


write(*,*)'#################################################'
write(*,*)'#                             				   #'
write(*,*)'#     		  CREATING THE MESH                #'
write(*,*)'#                             				   #'
write(*,*)'#################################################'

 call mesh
 
 call id_vector

write(*,*)'#################################################'
write(*,*)'#                             				   #'
write(*,*)'#             ASSEMBLING ID VECTOR              #'
write(*,*)'#                             				   #'
write(*,*)'#################################################' 
 
 call main
 call cpu_time(tf)

tpro=tf-ti

write(*,*) ''
write(*,*) 'TOTAL SIMULATION TIME:',tpro,'seconds'
write(*,*) 'FINAL TEMPERATURE AT THE CENTER OF THE TUMOUR:', T((m*n-m)/2)
write(*,*) 'TIME TO REACH STEADY-STATE TEMPERATURE:', timesteady

deallocate(ID)
deallocate(T)
deallocate(T_ANTES)
deallocate(T_AGORA)
deallocate(W)
deallocate(x)
deallocate(y)
deallocate(HZERO)
deallocate(raio_part)
deallocate(omega)
deallocate(PHI)
deallocate(timesteady)
deallocate(Tcenter)
deallocate(contagem)
end