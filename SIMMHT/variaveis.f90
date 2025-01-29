module variaveis

!############################################################
!############ VARIABLES DEFINITION MODULE ##################
!############################################################
implicit none

integer n,m,i,j,k,l,s, npast, Nrea, p
real, dimension(:), allocatable:: T_ANTES
real, dimension(:), allocatable:: T_AGORA
real, dimension(:), allocatable:: T
real, dimension(:), allocatable:: ID
real, dimension(:), allocatable:: x
real, dimension(:), allocatable:: y
real*16, dimension(:), allocatable:: W 
real, dimension(:), allocatable:: HZERO
real, dimension(:), allocatable:: raio_part
real, dimension(:), allocatable:: omega
real, dimension(:), allocatable::  PHI
real, dimension(:), allocatable:: timesteady
real, dimension(:), allocatable:: Tcenter
real tc_antes, tc_agora
integer, dimension(:), allocatable:: contagem
real xmin,xmax,ymin,ymax, ylinha, xlinha
real deltax,deltay
real dt, raio, xc, yc, rmin, rmax, eps, ax, by
character(3) rea_char
real*16 alphazero, raiomaximo, TZERO, T_ini
real RHO1, RHO2, RHO3, RHOB, C1, C2, C3, CB, Q1, Q2, Q3
real K1,K2,K3, TB
real DIFUSAO, PERFUSAO, GERACAO
real AUXREAL, auxiliarreal
integer AUXINT
real muzero, MD
real npastreal, lreal, VP, KB
logical decaicampo, ELIPTICO, CIRCULAR, fileinput
logical geracaodireta, STEADY, CALCCAMPO, EXPLICIT
real delta, VH, visc, gama, tauneel, taubrow, tauzero, taueff, ANISOTROPY
real pi, qsi, chii, chizero, cor, ximag,pimag, dip
real omegaestrela
real*16 aux
real a_xmin, a_xmax, b_ymin, b_ymax
integer t_int, TIU, auxiliarinteiro, t_tempi
real t_dec, tf, ti, tpro, sim_time, t_tempr


!############################################################
!################ VARIABLES’ DEFINITION ######################
!############################################################

! n - número de nós em y
! m - número de nós em x
! i,j,k,l - números inteiros utilizados em loops
! npast - numero de passos de tempo
! T - vetor de temperatura
! ID - vetor identificador da região (0,1,2)
! HX, HY - campo magnético (componentes em X e Y)
! H - módulo local do campo magnético
! ALPHA - campo alpha
! MX, MY - campo de magnetização (componentes em X e Y)
! ID - vetor identificador da região (0,1,2)
! 0 - fora do tumor, 1 - superfície do tumor, 2 - dentro do tumor
! x - vetor posição horizontal x
! y - vetor posição vertical y
! xmin, xmax, ymin, ymax - coordenadas dos nós de quina da malha
! ax, by - semi-eixo da elipse, a -x, b -y
! e - excentricidade da elipse
! ylinha, xlinha - posições do círculo externo hipotético para imposição de campo radial
! deltax, deltay - espaçamento entre nós em x e y respectivamente 
! dt - passo de tempo
! raio - raio do tumor
! xc, yc - posições x e y do centro do tumor
! rmin, rmax - coordenadas y do raio interno e do raio externo do tumor
! eps - tolerancia (camada de partículas na superficie do tumor)
! omega - frequência do campo
! alphazero - alpha base
! HZERO - campo magnético externo (valor máximo)
! TZERO - Temperatura de referencia
!
! PROPRIEDADES TERMOFÍSICAS:
!
! Os sub-indices 1,2,3 representam:    1 - tecido saudável
!				       2 - tecido cancerígeno
!				       3 - superfície do tumor com as partículas
!	
! W - taxa de perfusão sanguínea (nesse caso é um vetor, pois varia em função da temperatura)
! RHO - massa específica (RHOB seria a do sangue)
! C - calor específico (CB seria o do sangue)
! Q - taxa de geração metabólica
! PHI - fração volumétrica de partículas
! K - condutividade térmica
! TB - temperatura do sangue arterial
! T_ini - temperatura inicial
!
! DIFUSAO, PERFUSAO, GERACAO - termos locais computando cada um desses mecanismos 
! AUXREAL, AUXINT - variáveis reais e inteiras usadas para gravar arquivos de saída
! muzero - permeabilidade magnetica do vacuo
! MD - magnetização das partículas sólidas
! VP - volume médio das partículas
! raio_part - raio médio das partículas
! KB - constante de Boltzmann
! delta - espessura da camada de surfactante 
! VH - volume hidrodinâmico da partícula
! visc - viscosidade do fluido base 
! gama - variável usada na expressão do cálculo do tempo de Néel
! tauneel - tempo de Néel
! taubrow - tempo Browniano
! tauzero - attempt time (artigo de 2019 de Tang et al.) 
! taueff - escala de tempo de relaxação magnética efetiva 
! ANISOTROPY - constante de anisotropia 
! pi - é o número pi mesmo
! qsi, chii, chizero - parâmetros para implementação do termo de geração mag. diretamente
! cor - fator de correção para o termo de geração
! ximag - X''
! dip - magnitude do momento de dipolo da partícula
! aux - variável auxiliar para implementação do uso de X'' no termo de geração magnética
! t_real, t_int - variáveis usadas para determinar intervalo de extração dos dados

end module variaveis