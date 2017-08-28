all: compilar-secuencial compilar-paralelo

compilar-secuencial:
	gcc NCSecuencialMPI.c -lm -o NCSecuencialMPI -DLEERVARIABLES=0

compilar-secuencial-leyendoParametros:
	gcc NCSecuencialMPI.c -lm -o NCSecuencialMPI -DLEERVARIABLES=1

ejecutar-secuencial:
	./NCSecuencialMPI

compilar-paralelo:
	mpicc -g NCParaleloMPI.c -lm -o NCParaleloMPI

compilar-paralelo-sinResultados:
	mpicc -g NCParaleloMPI.c -lm -DNO_SAL -o NCParaleloMPI

ejecutar-paralelo:
	@read -p "Numero de procesos para lanzar la ejecucion: " procesos; \
	mpirun -n $$procesos NCParaleloMPI

limpiar-directorio:
	rm NCSecuencialMPI
	rm NCParaleloMPI
