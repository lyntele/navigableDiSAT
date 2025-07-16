#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include "obj.h"
#include "index.h"
#include "sat.h"
#include <string>
#include <vector>
using namespace std;

double numDistances;
int MaxHeight;
int main(int argc, char **argv)
{
	Index S;
	int n, np;
	struct stat sdata;
	char *fname, *finame;
	clock_t begin, buildEnd, queryEnd;
	double buildComp, queryComp;


	 if (argc < 4)
	{
		fprintf (stderr,"Usage: %s db-name size idxfile ... (extra index args)\n",argv[0]);
		system("pause");
		exit(1);
	}

	fname = argv[2];
	finame = argv[1];


	int qcount = 100;
	FILE * f = fopen(argv[3], "w");


	double radius[7];
	int kvalues[] = { 1, 5, 10, 20, 50, 100 };
	int mdim=2;
	char * querydata;
	querydata = argv[4];
	MaxHeight = atoi(argv[5]);

	if (string(finame).find("LA") != -1) {
		double r[] = { 473, 692, 989, 1409, 1875, 2314, 3096 };
		memcpy(radius, r, sizeof(r));
		mdim = 2;
	}
	else if (string(finame).find("integer") != -1) {
		double r[] = { 2321, 2733, 3229, 3843, 4614, 5613, 7090 };
		memcpy(radius, r, sizeof(r));
		mdim = 20;
	}
	else if (string(finame).find("mpeg_1M") != -1) {
		double r[] = { 3838, 4092, 4399, 4773, 5241, 5904, 7104 };
		memcpy(radius, r, sizeof(r));
		mdim = 282;
	}
	else if (string(finame).find("sf") != -1) {
		double r[] = { 100, 200, 300, 400, 500, 600, 700 };
		memcpy(radius, r, sizeof(r));
	}



    begin = clock();
    cout << "pn=" << pn << endl;
    n = openDB(finame);
    cout << "obj number:" << n << endl;
    buildEnd = clock() - begin;
    np = pn;
    if (!np || (np > n)) np = n;
    printf("indexing %li objects out of %li...\n", np, n);
    fprintf(f, "pivotnum: %d\n", np);


    cout<<"start building SAT......"<<endl;
    begin = clock();
    numDistances = 0;
    S = build(finame, np, &argc, &argv); //build SAT
    buildEnd += clock() - begin;
    buildComp = numDistances;
    fprintf(f, "building...\n");
    fprintf(f, "finished... %f build time\n", (double)buildEnd / CLOCKS_PER_SEC);
    fprintf(f, "finished... %f distances computed\n", buildComp);
    printf("\n\nsaving...\n");
    saveIndex(S, fname);
    stat(fname, &sdata);
    fprintf(f, "saved... %lli bytes\n", (long long)sdata.st_size);
    fprintf(f, "\nquerying...\n");



    cout << "start knnSearching......" << endl;
    FILE * fp;
    double rad;
    for (int k = 0; k < 6; k++) {
        ifstream fquery(querydata, ios::in);
        vector<float> obj(mdim, 0);
        begin = clock();
        numDistances = 0;
        rad = 0;
        for (int i = 0; i < qcount; i++) {
            for (int j = 0;j < mdim;j++)
            {
                fquery >> obj[j];
            }
            Objvector q(obj);
            rad += searchNN(S, q, kvalues[k], false); //knnSearch
        }
        queryEnd = clock() - begin;
        queryComp = numDistances;
        fprintf(f, "k: %d\n", kvalues[k]);
        fprintf(f, "finished... %f query time\n", (double)queryEnd / CLOCKS_PER_SEC / qcount);
        fprintf(f, "finished... %f distances computed\n", queryComp / qcount);
        fprintf(f, "finished... %f radius\n", rad / qcount);
        fprintf(f, "\n");
        fflush(f);
        fquery.close();
    }

    cout << "start rangeSearching......" << endl;
    for (int k = 0; k < 7; ++k) {
        ifstream fquery(querydata, ios::in);
        vector<float> obj(mdim, 0);
        begin = clock();
        numDistances = 0;
        rad = 0;
        for (int i = 0; i < qcount; i++) {
            for (int j = 0;j < mdim;j++)
            {
                fquery >> obj[j];
            }
            Objvector q(obj);
            rad += search(S, q, radius[k], false); //rangeSearch
        }
        queryEnd = clock() - begin;
        queryComp = numDistances;
        fprintf(f, "r: %f\n", radius[k]);
        fprintf(f, "finished... %f query time\n", (double)queryEnd / CLOCKS_PER_SEC / qcount);
        fprintf(f, "finished... %f distances computed\n", queryComp / qcount);
        fprintf(f, "finished... %f objs\n", rad / qcount);
        fprintf(f, "\n");
        fflush(f);
        fquery.close();
    }

	return 0;
}

