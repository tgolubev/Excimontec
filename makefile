CC = mpicxx
FLAGS = -Wall -Wextra -O3 -std=c++11
OBJS = main.o OSC_Sim.o Exciton.o Polaron.o Event.o Lattice.o Object.o Simulation.o Site.o Utils.o Setup_Poisson_bilayer.o Tridiag_solver.o Set_diagonals.o

Excimontec.exe : $(OBJS)
	$(CC) $(FLAGS) $(OBJS) -o Excimontec.exe

main.o : main.cpp OSC_Sim.h Exciton.h Polaron.h KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h 
	$(CC) $(FLAGS) -c main.cpp
	
OSC_Sim.o : OSC_Sim.h OSC_Sim.cpp Exciton.h Polaron.h KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h KMC_Lattice/Poisson/Setup_Poisson_bilayer.h 
	$(CC) $(FLAGS) -c OSC_Sim.cpp

Exciton.o : Exciton.h Exciton.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	$(CC) $(FLAGS) -c Exciton.cpp

Polaron.o : Polaron.h Polaron.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	$(CC) $(FLAGS) -c Polaron.cpp

Event.o : KMC_Lattice/Event.h KMC_Lattice/Event.cpp KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	$(CC) $(FLAGS) -c KMC_Lattice/Event.cpp

Lattice.o : KMC_Lattice/Lattice.h KMC_Lattice/Lattice.cpp KMC_Lattice/Site.h KMC_Lattice/Utils.h
	$(CC) $(FLAGS) -c KMC_Lattice/Lattice.cpp

Object.o : KMC_Lattice/Object.h KMC_Lattice/Object.cpp KMC_Lattice/Utils.h
	$(CC) $(FLAGS) -c KMC_Lattice/Object.cpp

Simulation.o : KMC_Lattice/Simulation.h KMC_Lattice/Simulation.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	$(CC) $(FLAGS) -c KMC_Lattice/Simulation.cpp
	
Site.o : KMC_Lattice/Site.h KMC_Lattice/Site.cpp
	$(CC) $(FLAGS) -c KMC_Lattice/Site.cpp
	
Utils.o : KMC_Lattice/Utils.h KMC_Lattice/Utils.cpp
	$(CC) $(FLAGS) -c KMC_Lattice/Utils.cpp
	
Setup_Poisson_bilayer.o: KMC_Lattice/Poisson/Setup_Poisson_bilayer.h KMC_Lattice/Poisson/Setup_Poisson_bilayer.cpp KMC_Lattice/Poisson/Set_diagonals.h  KMC_Lattice/Poisson/Tridiag_solver.h   
	$(CC) $(FLAGS) -c KMC_Lattice/Poisson/Setup_Poisson_bilayer.cpp
	
Tridiag_solver.o: KMC_Lattice/Poisson/Tridiag_solver.h  KMC_Lattice/Poisson/Tridiag_solver.cpp KMC_Lattice/Poisson/Set_diagonals.h KMC_Lattice/Poisson/Set_diagonals.cpp 
	$(CC) $(FLAGS) -c KMC_Lattice/Poisson/Tridiag_solver.cpp
	
Set_diagonals.o: KMC_Lattice/Poisson/Set_diagonals.h KMC_Lattice/Poisson/Set_diagonals.cpp
	$(CC) $(FLAGS) -c KMC_Lattice/Poisson/Set_diagonals.cpp
	
clean:
	\rm *.o *~ Excimontec.exe
