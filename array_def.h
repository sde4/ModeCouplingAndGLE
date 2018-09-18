#ifndef arr_def
#define arr_def

int*		CreateIntVector	(int 	dim1);
double* 	CreateVector	(int 	dim1);
int** 		CreateIntMatrix	(int 	dim1, int dim2);
double**	CreateMatrix	(int 	dim1, int dim2);
double*** 	CreateMatrix3D	(int 	dim1, int dim2, int dim3);
void 		FreeIntMatrix	(int 	**m);
void 		FreeMatrix	(double 	**m);

#endif
