#include <stdio.h>
#include <time.h>
#include <iostream>
#include "index.h"
#include "basics.h"
using namespace std;
extern int MaxHeight;
typedef struct
{
	Obj elem;
	Tdist dist;
	int who;
} qelem;

typedef struct
{
	qelem *data;
	int ndata;
} arrdata;

typedef struct
{
	int *num;
	int nnum;
} arrnum;

//node
typedef struct
{
	//int height;
	Obj obj;
	Tdist maxd;
	arrnum vec;
	arrdata queue;
} nodo;

typedef struct
{
	int id;
	Tdist dist;  /* dist to q */
	Tdist lbound;  /* lower bound */
	Tdist mind; /* best distsance */
} heapElem;

typedef struct
{
	nodo *nodos;
	int nnodos;
	int np;
	char *descr;
} grafo;

static double avgprof (grafo G, int n, int height, int *num);

void printTree(Index *S);
void traverse(nodo *root, grafo *G);
