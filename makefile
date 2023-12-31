
# Arquivos fonte padrão (linux) e outros sistemas
#
SRC-LNX = variaveis.f90 funcoes.f90 entrada.f90 principal.f90  \
          
OBJ-LNX = $(SRC-LNX:.f90=.o)

#
# Definição do compilador específico e das flags

#ENGINE-LNX = gfortran
ENGINE-LNX = ifort
FLAGS-LNX = -m64 -O2 -qopenmp -mkl
#FLAGS-LNX = -m64 -O2 -ldislin 
#FLAGS-LNX = -m64 -O2 -openmp -traceback -fpe0 -g

#
# Objetivo principal: geração de executável em linux


sims.ex : $(SRC-LNX) 
	$(ENGINE-LNX) $(FLAGS-LNX) -o simmht.ex $(SRC-LNX)


#
# Limpeza
clean :
	rm simmht.ex
