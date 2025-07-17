
#include "obj.h"

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

extern double numDistances;

Objvector::Objvector(const Objvector & obj)
{
	value = obj.value;
}

Objvector::Objvector(vector<float> v)
{
	value = v;
}

Objvector & Objvector::operator=(const Objvector & obj)
{

	if (this == &obj) {
		return *this;
	}
	else {
		value = obj.value;
		return *this;
	}
}

vector<float> & Objvector::getValue()
{
	return value;
}



void Objvector::setValue(vector<float> v)
{
	value = v;
}

typedef struct sEuclDB
{
	float *nums;  /* coords all together */
    int nnums;	  /* number of vectors (with space for one more) */
    int coords;  /* coordinates */
    double (*df) (float *p1, float *p2, int k); /* distance to use */
	float *nobj;
} EuclDB;

static int never = 1;
static EuclDB DB;

#define db(p) (DB.nums + DB.coords*(int)p)

	/* L2 distance */
static double L2D (float *p1, float *p2, int k)
{
	register int i;
    double tot = 0,dif;
    for (i=0;i<k;i++)
	{
		tot += pow(p1[i] - p2[i],2);
	}
    return sqrt(tot);
}

static double L1D (float *p1, float *p2, int k)
{
	register int i;
    double tot = 0,dif;
    for (i=0;i<k;i++)
	{
		dif = (p1[i]-p2[i]);
		if (dif < 0) dif = -dif;
		tot += dif;
	}
    return tot;
}

static double LiD (float *p1, float *p2, int k)
{
	register int i;
    double max = 0,dif;
    for (i=0;i<k;i++)
	{
		dif = (p1[i]-p2[i]);
		if (dif < 0) dif = -dif;
		if (dif > max) max = dif;
	}
    return max;
}

double distanceInter (Obj o1, Obj o2)
{
	numDistances++;
	return DB.df (db(o1),db(o2),DB.coords);
}

double qdistanceInter(Objvector o1,Obj o2)
{
	numDistances++;
	for (int d = 0; d < DB.coords; ++d) {
		DB.nobj[d] = o1.getValue()[d];
	}
	return DB.df(DB.nobj, db(o2), DB.coords);
}
Obj parseobj (char *p)
{
	float *d = db(NewObj);
    int i,step;
    for (i=0;i<DB.coords-1;i++)
	{
		sscanf (p,"%f,%n",d+i,&step);
		p += step;
	}
    sscanf (p,"%f",d+i);
    return NewObj;
}

void printobj (Obj obj)
{
	int i;
    float *p = db(obj);
    for (i=0;i<DB.coords-1;i++) printf ("%f,",p[i]);
    printf ("%f\n",p[i]);
}

int openDB(char *name)
{
	FILE *f = fopen(name, "r");
	if (!f) {
		printf("Error: Cannot open file %s\n", name);
		exit(1);
	}
	int dim, num, func;
	fscanf(f, "%d %d %d\n", &dim, &num, &func);
	printf("Reading database: dim=%d, num=%d, func=%d\n", dim, num, func);
	if (func == 1) DB.df = L1D;
	else if (func == 2) DB.df = L2D;
	else DB.df = LiD;
	DB.coords = dim;
	DB.nnums = num;
	printf("Allocating memory for %d vectors of %d dimensions\n", DB.nnums, DB.coords);
	DB.nums = (float *)mymalloc((DB.nnums + 1) * sizeof(float)* DB.coords);
	DB.nobj = (float *)mymalloc(DB.coords * sizeof(float));

	for (int i = 0; i < DB.nnums; ++i)
	{
		for (int j = 0; j < DB.coords; ++j)
		{
			fscanf(f, "%f", &DB.nums[i*DB.coords + j]);
		}
		fscanf(f, "\n");

	}
	fclose(f);
	return DB.nnums;
}

void closeDB (void)
{
	if (never) { DB.nums = NULL; never=0;}
    if (DB.nums == NULL) return;
    myfree (DB.nums);
    DB.nums = NULL;
}

void printDBelements() {
	for (int i = 0; i < DB.nnums; ++i)
	{
		for (int j = 0; j < DB.coords; ++j)
		{
			cout << DB.nums[i*DB.coords + j] << " ";
		}
		cout << endl;
	}
}
