/*
Asignatura: Computadores avanzados 2016
Autor: Cristian Trapero Mora
Trabajo: Problema de los n cuerpos paralelo rapido
Compilar: mpicc -g NCParaleloMPI.c -lm -o NCParaleloMPI
Ejecutar: mpirun -n 2 NCParaleloMPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>

#define G 1
#define FICHERO "datos.dat"

// Estructura con las variables del problema
struct variablesProblema {
	int n;
	int k;
	int tp;
	double delta;
	double u;
};

// Estructura que define una coordenada (X,Y) de un cuerpo
struct coordenadas {
	int id;
	double x;
	double y;
};

// Estructuras MPI para especificar nuevos tipos de datos
MPI_Datatype MPI_VARIABLES;
MPI_Datatype MPI_COORDENADAS;

// Estructura de apoyo para CNCR
MPI_Datatype MPI_CNC;
MPI_Datatype MPI_CNCR;

// Variables del problema
struct variablesProblema variables;

// Todas las masas de los cuerpos
double *masas;

// Todas las posiciones de los cuerpos
struct coordenadas *posiciones;

// Todas las velocidades de los cuerpos
struct coordenadas *velocidades;

// Todas las velocidades de los cuerpos
struct coordenadas *aceleraciones;

// Estructuras necesarias para algoritmo en anillo
struct coordenadas *p_local;
struct coordenadas *p_anillo;
struct coordenadas *v_local;
struct coordenadas *a_local;
struct coordenadas *a_anillo;

// Identificador de proceso y numero de procesadores
int rank, npr;

// Variables para el calculo de cuerpos por proceso
int n_cuerpos, n_restantes;

// Tiempo de ejecucion maximo de entre todos los procesos
double tiempo_global;

int leerVariables();
int leerCuerpos();
void crearEstructuraMPIVariables();
void crearEstucturaMPICoordenadas();
void crearEstucturaMPICNCR(int n_cuerpos, int npr);
void calcularAceleraciones();
void actualizarAceleraciones(int fase, struct coordenadas *p_cuerpos);
void imprimirDatos(int cuerpos_totales);
void liberarMemoria();

int main(int argc, char **argv) {

	// Variables para medir el tiempo de ejecucion del programa
	double inicio_local, final_local, tiempo_local;

	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	inicio_local = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npr);

	// Creamos un MPI_Datatype con la estructura del struct variablesProblema
	crearEstructuraMPIVariables();

	//Creamos un MPI_Datatype con la estructura del struct coordenadas
	crearEstucturaMPICoordenadas();

	// Solo el proceso 0 lee el fichero
	if(rank == 0){
		//Leemos las variables del problema
		if(leerVariables()) return 1;
	}

	// Envio/recibo las variables del problema a todos los procesos
	MPI_Bcast(&variables, 1, MPI_VARIABLES, 0, MPI_COMM_WORLD);

	n_cuerpos = variables.n / npr;
	n_restantes = variables.n % npr;

	// Si hay cuerpos sin repartir, todos los procesos tendrán uno más
	if (n_restantes > 0) n_cuerpos++;

	// Variable usada para crear los cuerpos vacios
	int cuerpos_totales = n_cuerpos * npr;

	if(rank == 0){
		// Reservamos memoria para leer todos los cuerpos del fichero
		masas = malloc(sizeof(double) * cuerpos_totales);
		posiciones = malloc(sizeof(struct coordenadas) * cuerpos_totales);
		velocidades = malloc(sizeof(struct coordenadas) * cuerpos_totales);
		aceleraciones = malloc(sizeof(struct coordenadas) * cuerpos_totales);

		// Obtengo todos los datos de los cuerpos
		if(leerCuerpos(variables.n, cuerpos_totales)) return 1;
	}

	crearEstucturaMPICNCR(n_cuerpos, npr);

	p_local = malloc(sizeof(struct coordenadas) * n_cuerpos);
	p_anillo = malloc(sizeof(struct coordenadas) * n_cuerpos);
	v_local = malloc(sizeof(struct coordenadas) * n_cuerpos);
	a_local = malloc(sizeof(struct coordenadas) * n_cuerpos);
	a_anillo = malloc(sizeof(struct coordenadas) * n_cuerpos);

	// Todos los procesos deben tener una copia de las masas de todos los cuerpos
	if(rank != 0) masas = malloc(sizeof(double) * cuerpos_totales);

	// Enviamos/recibimos las masas a todos los procesos
	MPI_Bcast(masas, variables.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Enviamos/recibimos las posiciones de solo los cuerpos que le corresponde a cada proceso
	MPI_Scatter(posiciones, 1, MPI_CNCR, p_local, n_cuerpos, MPI_COORDENADAS, 0, MPI_COMM_WORLD);

	// Enviamos/recibimos las velocidades de solo los cuerpos que le corresponde a cada proceso
	MPI_Scatter(velocidades, 1, MPI_CNCR, v_local, n_cuerpos, MPI_COORDENADAS, 0, MPI_COMM_WORLD);

	// ----------COMENZAMOS LOS CALCULOS DE LAS ACELERACIONES----------
	// Incremento de tiempo
	double t = 0.0;

	//Calculamos la aceleracion inicial
	calcularAceleraciones();

	// Actualizamos los datos del proceso 0 con los datos locales de cada proceso (inversa a MPI_Scatter)
	MPI_Gather(a_local, n_cuerpos, MPI_COORDENADAS, aceleraciones, n_cuerpos, MPI_COORDENADAS, 0, MPI_COMM_WORLD);

	//Imprimimos aceleracion inicial
	#ifndef NO_SAL
	if(rank == 0) {
		printf("Por cada instante de tiempo y cada cuerpo, aparecen:\n");
		printf("	cuerpo(id)	posicion(x),	posicion(y),	velocidad(x),	velocidad(y),	aceleracion(x),	aceleracion(y)\n");
		printf("%.2f\n",t);
		imprimirDatos(cuerpos_totales);
	}
	#endif

	// Bucle de calculos
	for(int paso=1; paso<=variables.tp; paso++) {

		//Calculamos la posicion y velocidad de los cuerpos
		for(int q=0; q<n_cuerpos; q++) {
			p_local[q].x += v_local[q].x * variables.delta;
			p_local[q].y += v_local[q].y * variables.delta;

			v_local[q].x += a_local[q].x * variables.delta;
			v_local[q].y += a_local[q].y * variables.delta;
		}

		//Calculamos la aceleracion d los cuerpos
		calcularAceleraciones();

		// Actualizamos los datos del proceso 0 con los datos locales de cada proceso (inversa a MPI_Scatter)
		MPI_Gather(p_local, n_cuerpos, MPI_COORDENADAS, posiciones, n_cuerpos, MPI_COORDENADAS, 0, MPI_COMM_WORLD);
		MPI_Gather(v_local, n_cuerpos, MPI_COORDENADAS, velocidades, n_cuerpos, MPI_COORDENADAS, 0, MPI_COMM_WORLD);
		MPI_Gather(a_local, n_cuerpos, MPI_COORDENADAS, aceleraciones, n_cuerpos, MPI_COORDENADAS, 0, MPI_COMM_WORLD);

		// Imprimimos los cuerpos cada k pasos
		#ifndef NO_SAL
		if(paso%variables.k==0 || paso==variables.tp)	{
			if(rank == 0){
				printf("\n%.2f\n",t);
				imprimirDatos(cuerpos_totales);
			}
		}
		#endif

		t+=variables.delta;
	}

	final_local = MPI_Wtime();
	tiempo_local = final_local - inicio_local;

	// Obtenemos el mayor tiempo global de todos los procesos
	MPI_Reduce(&tiempo_local, &tiempo_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank==0)	printf("Tiempo de ejecucion maximo: %lf segundos\n", tiempo_global);

	MPI_Finalize();
	liberarMemoria();
	return 0;
}

int leerVariables() {
	int resultado = 0;
	FILE *fichero = fopen(FICHERO, "r");

	if (fichero == NULL) {
		printf("No se ha podido abrir el fichero.\n");
		resultado = 1;
	} else{
		// Leemos las variables del problema
		fscanf(fichero, "%d, %lf, %d, %lf, %d", &variables.n,  &variables.delta, &variables.tp, &variables.u, &variables.k);
	}
	fclose(fichero);
	return resultado;
}

int leerCuerpos(int cuerpos, int cuerpos_totales) {
	int resultado = 0;
	FILE *fichero = fopen(FICHERO, "r");

	if (fichero == NULL) {
		printf("No se ha podido abrir el fichero.\n");
		resultado = 1;
	} else{
		//No saltamos la primera linea de variables
		fscanf(fichero, "%*[^\n]");

		// Leemos los datos de cada cuerpo
		for(int i=0; i<cuerpos_totales; i++){
			// Si es un cuerpo valido
			if(i<cuerpos){
				//Recuperamos los datos del cuerpo: masa, posicionX, posicionY, velocidadX, velocidadY
				fscanf(fichero, "%lf, %lf, %lf, %lf, %lf", &(masas[i]), &(posiciones[i]).x, &(posiciones[i]).y, &(velocidades[i]).x, &(velocidades[i]).y);
				posiciones[i].id=i;
				velocidades[i].id=i;
				aceleraciones[i].id=i;
				aceleraciones[i].x = 0;
				aceleraciones[i].y = 0;
			}else{
				// Creamos un cuerpo vacio
				masas[i] = 0;
				posiciones[i].id = -1;
				posiciones[i].x = 0;
				posiciones[i].y = 0;
				velocidades[i].id = -1;
				velocidades[i].x = 0;
				velocidades[i].y = 0;
				aceleraciones[i].id = -1;
				aceleraciones[i].x = 0;
				aceleraciones[i].y = 0;
			}
		}
	}
	fclose(fichero);
	return resultado;
}

// Ejemplo extraido de http://mpi.deino.net/mpi_functions/MPI_Type_struct.html
void crearEstructuraMPIVariables() {
	int count = 2;
	int array_blocklens[2] = {3, 2};
	MPI_Datatype old_types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Aint intextsn;
	MPI_Type_extent(MPI_INT, &intextsn);
	MPI_Aint offsets[2] = {(MPI_Aint)0,3*intextsn};
	MPI_Type_struct(count, array_blocklens, offsets, old_types, &MPI_VARIABLES);
	MPI_Type_commit(&MPI_VARIABLES);
}

void crearEstucturaMPICoordenadas() {
	int count = 2;
	int blocklens[2] = {1, 2};
	MPI_Datatype old_types[2] = {MPI_LONG_LONG, MPI_DOUBLE};
	MPI_Aint intextsn;
	MPI_Type_extent(MPI_LONG_LONG, &intextsn);
	MPI_Aint offsets[2] = {(MPI_Aint)0, intextsn};
	MPI_Type_struct(count, blocklens, offsets, old_types, &MPI_COORDENADAS);
	MPI_Type_commit(&MPI_COORDENADAS);
}

void crearEstucturaMPICNCR(int n_cuerpos, int npr) {
	MPI_Type_vector(n_cuerpos, 1, npr, MPI_COORDENADAS, &MPI_CNC);
	MPI_Type_commit(&MPI_CNC);
	MPI_Aint lb, extent;
	MPI_Type_get_extent(MPI_COORDENADAS, &lb, &extent);
	MPI_Type_create_resized(MPI_CNC, lb, extent, &MPI_CNCR);
	MPI_Type_commit(&MPI_CNCR);
}

void calcularAceleraciones(){
	// Algoritmo en anillo pp. 38-39
	int fte = (rank+1) % npr;
	int dest = (rank-1+npr) % npr;

	// p_anillo <- p_local
	for (int i=0; i<n_cuerpos; i++){
		p_anillo[i].id = p_local[i].id;
		p_anillo[i].x = p_local[i].x;
		p_anillo[i].y = p_local[i].y;
	}

	// a_local <- a_anillo <- 0
	for (int i=0; i<n_cuerpos; i++){
		a_anillo[i].id = p_local[i].id;
		a_local[i].id = a_anillo[i].id;
		a_anillo[i].x = 0;
		a_local[i].x = 0;
		a_anillo[i].y = 0;
		a_local[i].y = 0;
	}

	// Actualizamos a_local y a_anillo sumando las aceleraciones debidas a las fuerzas
	actualizarAceleraciones(0, p_local);

	// Obligatorio para la directiva MPI_Sendrecv_replace
	MPI_Status status;

	for(int fase=1; fase<npr; fase++) {
		// Envio/recibo p_anillo, implementado como en la pp. 43
		MPI_Sendrecv_replace(p_anillo, n_cuerpos, MPI_COORDENADAS, dest, 1, fte, 1, MPI_COMM_WORLD, &status);

		// Envio/recibo a_anillo, implementado como en la pp. 43
		MPI_Sendrecv_replace(a_anillo, n_cuerpos, MPI_COORDENADAS, dest, 1, fte, 1, MPI_COMM_WORLD, &status);

		// Actualizamos a_local y a_anillo sumando las aceleraciones debidas a las fuerzas
		actualizarAceleraciones(fase, p_anillo);
	}

	//Enviamos/recibimos a_anillo, implementado como en la pagina 43
	MPI_Sendrecv_replace(a_anillo, n_cuerpos, MPI_COORDENADAS, dest, 1, fte, 1, MPI_COMM_WORLD, &status);

	// a_local = a_local + a_anillo
	for(int i=0; i<n_cuerpos; i++){
		a_local[i].x += a_anillo[i].x;
		a_local[i].y += a_anillo[i].y;
	}
}

void imprimirDatos(int cuerpos_totales) {
	for(int i=0; i<cuerpos_totales; i++){
		if(posiciones[i].id != -1)	printf("		%d	%lf	%lf	%lf	%lf	%lf	%lf\n", posiciones[i].id, posiciones[i].x, posiciones[i].y, velocidades[i].x, velocidades[i].y, aceleraciones[i].x, aceleraciones[i].y);
	}
}

void actualizarAceleraciones(int fase, struct coordenadas *p_cuerpos) {
	// Variables auxiliares para el calculo de las aceleraciones
	double modulo, denominador, distanciaX, distanciaY, aceleracionAuxX, aceleracionAuxY;

	// Identificador del cuerpo
	int id_cuerpo = (rank+fase) % npr;

	for(int i = 0, p = id_cuerpo; i < n_cuerpos, p < n_cuerpos*npr; i++, p += npr) {
		for(int j = 0, q = rank; j < n_cuerpos, q < n_cuerpos*npr; j++, q += npr) {

			// Calculamos las aceleraciones si el cuerpo p es mayor que q
			if(p>q) {
				// Calculamos la distancia (modulo) entre dos cuerpos:
				//1. Distancia entre las posiciones x
				//2. Distancia entre las posiciones y
				//3. Raiz cuadrada de la suma de dichas distancias
				distanciaX = p_local[j].x - p_cuerpos[i].x;
				distanciaY = p_local[j].y - p_cuerpos[i].y;
				modulo = sqrt(distanciaX * distanciaX + distanciaY * distanciaY);

				//Solo si la distancia entre los cuerpos es mayor que el umbral
				if(modulo >= variables.u) {
					// Denominador de la funcion de aceleracion
					denominador = pow(modulo, 3);

					//Calculamos la aceleracion de q con respecto p
					aceleracionAuxX = (G * -distanciaX) / denominador;
					aceleracionAuxY = (G * -distanciaY) / denominador;

					// Aceleracion Apq
					a_local[j].x += aceleracionAuxX * masas[p_cuerpos[i].id];
					a_local[j].y += aceleracionAuxY * masas[p_cuerpos[i].id];

					// Aceleracion Aqp ahorrando calculos
					a_anillo[i].x += -masas[p_local[j].id] * aceleracionAuxX;
					a_anillo[i].y += -masas[p_local[j].id] * aceleracionAuxY;
				}
			}
		}
	}
}

void liberarMemoria() {
	free(masas);
	free(posiciones);
	free(velocidades);
	free(aceleraciones);
	free(p_local);
	free(p_anillo);
	free(v_local);
	free(a_local);
	free(a_anillo);
}
