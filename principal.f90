program mht2d

use omp_lib
use lapack95
use f95_precision

!##############################
!#                            #
!#    PROGRAM mht2d           #
!#                            #
!##############################

! Wrote by Rafael Gabler at 04/03/2021
 
use variaveis
use funcoes


write(*,*) 'Type the lowest X coordinate of the calculation domain'
read*,xmin
write(*,*) 'Type the lowest Y coordinate of the calculation domain'
read*,ymin
write(*,*) 'Type the highest X coordinate of the calculation domain'
read*,xmax
write(*,*) 'Type the highest Y coordinate of the calculation domain'
read*,ymax

write(*,*) 'How many tumours you want to simulate?'
read*,nt

allocate(xc(nt))
allocate(yc(nt))
allocate(raio(nt))

do i=1,nt
write(*,*) 'Type the X coordinate of the center of tumour number:',i
read*,xc(i)
write(*,*) 'Type the Y coordinate of the center of tumour number:',i
read*,yc(i)
write(*,*) 'Type the radius of tumour number:',i
read*,raio(i)
end do

call entrada
call cpu_time(ti)

allocate(ID(m*n))
allocate(T(m*n))
allocate(T_ANTES(m*n))
allocate(T_AGORA(m*n))
!if(CALCCAMPO)then
!allocate(HX(m*n))
!allocate(HY(m*n))
!allocate(H(m*n))
!end if
allocate(W(m*n))
allocate(x(n))
allocate(y(m))
if(STEADY)then
allocate(A(m*n,m*n))
allocate(B(m*n))
end if

a_xmin = 0.003958
a_xmax = 0.016042							! Ponto máximo do tumor em x
b_ymin = 0.008113
b_ymax = 0.011887 							! Ponto máximo do tumor em y

if( ELIPTICO == .TRUE. ) then
 CIRCULAR = .FALSE.
else
  CIRCULAR = .TRUE.
end if

!xc = (xmax-xmin)/2.0
!yc = (ymax-ymin)/2.0
!raiomaximo = ((ymax-yc)**2.0 + (xmax-xc)**2.0)**0.5
pi = acos(-1.0)
muzero = 4.0*pi*1.0E-07

dt = 1.0E-03 !min(1.0E-03,1.0/omega)


npast = sim_time/dt
npastreal = npast

!**************************************************************************!

write(*,*)'#################################################################'
write(*,*)'#                             				   #'
write(*,*)'#     		  CREATING THE MESH               	   #'
write(*,*)'#                             				   #'
write(*,*)'#################################################################'

 call mesh
 
 call id_vector

write(*,*)'#################################################################'
write(*,*)'#                             				   #'
write(*,*)'#         ASSEMBLING THE IDENTIFICATION VECTOR              	   #'
write(*,*)'#                             				   #'
write(*,*)'#################################################################' 
 
 call principal
 call cpu_time(tf)


tpro=tf-ti

write(*,*) ''
write(*,*) 'TOTAL SIMULATION TIME:',tpro,'seconds'
write(*,*) 'FINAL TEMPERATURE AT THE CENTER OF THE TUMOUR:', T((m*n-m)/2)
write(*,*) 'TIME TO REACH STEADY-STATE TEMPERATURE:', timesteady

! ##### WRITING A TEXT FILE WITH SOME IMPORTANT INFORMATIONS #####
open (out, file='info.txt',status='unknown')
write(out,'(A25)')'INFORMACOES DA SIMULACAO'
write(out,'(A18,I3,A2)')'Tempo simulado = ',sim_time,'s'
write(out,'(A28,BZ,F7.3,A2)')'Tempo computacional gasto =',tpro,'s'
write(out,'(A8,BZ,I3,A3,I3)')'Malha: ',m,'x ',n
write(out,'(A22,BZ,I2,A2,A22)')'Plotando dados a cada',t_int,'s','do tempo da simulacao'
write(out,'(A47,BZ,I2)')'Quantidade de nucleos atualmente sendo usados:', TIU
write(out,'(A7,BZ,F7.4)')'alpha:', alphazero
write(out,'(A4,BZ,F7.4)')'w*:', omegaestrela
write(out,'(A6,BZ,F5.3)')'phi: ', PHI
write(out,'(A4,BZ,F8.5)')'X":', ximag
write(out,'(A8,BZ,F6.3)')'T_ini: ', T_ini


end
