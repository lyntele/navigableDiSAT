#ifndef OBJINCLUDED
#define OBJINCLUDED

	/* object database. for now, only one active database at a time will
	   be permitted. all the operations fail if there is no loaded DB */

typedef int Obj;	/* object id */

#define NullObj (-1)   /* null object */

#define NewObj 0      /* new object */

#include "basics.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
class Objvector
{
	vector<float> value;
	//int dim;
public:
	Objvector() {}
	Objvector(const Objvector &obj);
	Objvector(vector<float> v);
	Objvector & operator = (const Objvector & obj);
	vector<float> & getValue();
	void setValue(vector<float> v);
	~Objvector() {}
};
	/* loads a DB given its name, and returns its size. if there is
	   another one already open, it will be closed first */
int openDB (char *name);

	/* frees the currently open DB, if any */
void closeDB (void);

	/*(internal) computes the distance between two objects */
#define CONT
#ifdef CONT
double distanceInter (Obj o1, Obj o2);
double qdistanceInter(Objvector o1, Obj o2);
#else
int distanceInter (Objvector o1, Objvector o2);
#endif

      /* exported, computes distance and does the accounting */
#define mydistance(o1,o2)(distanceInter(o1,o2))
#define myqdistance(o1,o2)(qdistanceInter(o1,o2))

	/* returns an object identifier from description s. if the object
	   does not belong to the DB, it must be identified as NewObj */
Obj parseobj (char *str);

	/* prints obj in a user-friendly way, terminated with a newline */
void printobj (Obj obj);

void printDBelements();
#endif
