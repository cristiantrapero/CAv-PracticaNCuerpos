/*
	Asignatura: Computadores avanzados 2016
	Autor: Cristian Trapero Mora
	Trabajo: Problema de los n cuerpos secuencial
	Compilar: gcc NCSecuencialMPI.c -lm -o NCSecuencialMPI -DLEERVARIABLES=0
	Ejecutar: ./NCSecuencialMPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"

#define G 1
#define FICHERO "datos.dat"
//#define LEERVARIABLES 0

struct cuerpo {
	double masa;
	double posicionX;
	double posicionY;
	double velocidadX;
	double velocidadY;
	double aceleracionX;
	double aceleracionY;
};

struct cuerpo *cuerpos;
int n, tp, k;
double delta, u;

int leerFichero();
void leerEntradas();
void reservarEspacio();
void calcularAceleraciones();
void imprimirCuerpos();

int main(int argc, char *argv[]) {
	// Leemos los datos del problema
	leerEntradas();

	// Reservamos el espacio para los cuerpos
	cuerpos = malloc(sizeof(struct cuerpo)*n);

	// Variables para obtener el tiempo de ejecucion
	double inicio, fin;

	// t: tiempo simulacion; paso: paso en el que estamos; q: cuerpo;
	double t;
	int paso, q;

	// Si la lectura del fichero es correcta, continuamos
	if(leerFichero()) return 1;

	printf("\nPor cada instante de tiempo y cada cuerpo, aparecen:\n");
	printf("	posicion(x),	posicion(y),	velocidad(x),	velocidad(y),	aceleracion(x),	aceleracion(y)");

	GET_TIME(inicio);
	t = 0.0;

	calcularAceleraciones();
	imprimirCuerpos();

	for (paso=1; paso<=tp; paso++){

		//Calculamos la posicion y la velocidad de los cuerpos
		for(q=0; q<n; q++){
			cuerpos[q].posicionX += cuerpos[q].velocidadX * delta;
			cuerpos[q].posicionY += cuerpos[q].velocidadY * delta;

			cuerpos[q].velocidadX += cuerpos[q].aceleracionX * delta;
			cuerpos[q].velocidadY += cuerpos[q].aceleracionY * delta;
		}

		calcularAceleraciones();

		if( paso%k==0 || paso==tp )	{
			printf("\n\n%.2f", t);
			imprimirCuerpos();
		}
		t+=delta;
	}

	GET_TIME(fin);

	printf("\nTiempo de ejecución: %lf segundos.\n", fin-inicio);

	free(cuerpos);
	return 0;
}

void leerEntradas() {

	if(LEERVARIABLES)	{
		printf("Introduce el numero de cuerpos (n): ");
		scanf("%d", &n);

		printf("Introduce el incremento del tiempo en cada paso (delta): ");
		scanf("%lf", &delta);

		printf("Introduce el numero total de pasos (Tp): ");
		scanf("%d", &tp);

		printf("Introduce la distancia umbral (U): ");
		scanf("%lf", &u);

		printf("Introduce cada cuantos pasos imprimir los resultados (k): ");
		scanf("%d", &k);
	}else{
		int resultado = 0;
		FILE *fichero = fopen(FICHERO, "r");

		if (fichero == NULL) {
			printf("No se ha podido abrir el fichero.\n");
			resultado = 1;
		} else{
			fscanf(fichero, "%d, %lf, %d, %lf, %d", &n,  &delta, &tp, &u, &k);
		}
		fclose(fichero);
	}
}

int leerFichero() {

	int resultado = 0;
	FILE *fichero = fopen(FICHERO, "r");

	if (fichero == NULL) {
		printf("No se ha podido abrir el fichero.\n");
		resultado = 1;
	} else{
		//No saltamos la primera linea de variables
		fscanf(fichero, "%*[^\n]");
		for(int i=0; i<n; i++){
			//Recuperamos los datos del cuerpo: masa, posicionX, posicionY, velocidadX, velocidadY
			fscanf(fichero, "%lf, %lf, %lf, %lf, %lf", &(cuerpos[i]).masa, &(cuerpos[i]).posicionX, &(cuerpos[i]).posicionY, &(cuerpos[i]).velocidadX, &(cuerpos[i]).velocidadY);
		}
	}
	fclose(fichero);
	return resultado;
}

void calcularAceleraciones(){
	int p, q;

	// Variables auxiliares para el calculo de la aceleracion entre dos cuerpos
	double modulo, denominador, distanciaX, distanciaY, aceleracionAuxX, aceleracionAuxY;

	for(q=0; q<n; q++){
		cuerpos[q].aceleracionX = 0.0;
		cuerpos[q].aceleracionY = 0.0;
	}

	//Calculamos la aceleracion inicial de los cuerpos
	for (q = 0; q<n; q++) {
		for(p = q+1; p<n; p++) {

				// Calculamos la distancia (modulo) entre dos cuerpos:
					//1. Distancia entre las posiciones x
					//2. Distancia entre las posiciones y
					//3. Raiz cuadrada de la suma de dichas distancias
				distanciaX=cuerpos[q].posicionX - cuerpos[p].posicionX;
				distanciaY=cuerpos[q].posicionY - cuerpos[p].posicionY;
				modulo = sqrt(distanciaX * distanciaX + distanciaY * distanciaY);

				//Solo si la distancia entre los cuerpos es mayor que el umbral
				if(modulo >= u) {
					// Denominador de la funcion de aceleracion
					denominador = pow(modulo, 3);

					//Calculamos la aceleracion de q con respecto p
					aceleracionAuxX = (G * -distanciaX) / denominador;
					aceleracionAuxY = (G * -distanciaY) / denominador;

					//Aceleracion del cuerpo q
					cuerpos[q].aceleracionX += aceleracionAuxX * cuerpos[p].masa;
					cuerpos[q].aceleracionY += aceleracionAuxY * cuerpos[p].masa;

					//Aceleracion del cuerpo p ahorrando calculos
					cuerpos[p].aceleracionX += -cuerpos[q].masa * aceleracionAuxX;
					cuerpos[p].aceleracionY += -cuerpos[q].masa * aceleracionAuxY;
				}
		}
	}
}

void imprimirCuerpos() {
	for(int i=0 ; i<n ; i++){
		printf("\n     %d	%lf	%lf	%lf	%lf	%lf	%lf", i, cuerpos[i].posicionX, cuerpos[i].posicionY, cuerpos[i].velocidadX, cuerpos[i].velocidadY, cuerpos[i].aceleracionX, cuerpos[i].aceleracionY);
	}
}
