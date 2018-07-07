#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int* CreateIntVector(int dim1)
{
  return calloc((unsigned int) dim1, sizeof(unsigned int));
}

float* CreateVector(int dim1)
{
  return calloc((unsigned int) dim1, sizeof(float));
}

int** CreateIntMatrix(int dim1, int dim2)
{
  int      i;
  int    **m;
  m = calloc((unsigned int) dim1, sizeof(unsigned int *));
  for (i=0; i<dim1; i++)
  {
    m[i] = calloc((unsigned int) dim2, sizeof(unsigned int));
  }
  return m;
}

float**	CreateMatrix(int dim1, int dim2)
{
  int      i;
  float    **m;
  m = calloc((unsigned int) dim1, sizeof(float *));
  for (i=0; i<dim1; i++)
  {
    m[i] = calloc((unsigned int) dim2, sizeof(float));
  }
  return m;
}

float*** CreateMatrix3D(int dim1, int dim2, int dim3)
{
  int      i, j;
  float    ***m;
  m = calloc((unsigned int) dim1, sizeof(float **));
  for (i=0; i<dim1; i++)
  {
    m[i] = calloc((unsigned int) dim2, sizeof(float *));
    for (j=0; j<dim2; j++)
    {
      m[i][j] = calloc((unsigned int) dim3, sizeof(float));
    }
  }
  return m;
}

void FreeIntMatrix(int **m)
{
  int      i;
  for (i=0; m[i]!=NULL; i++)
  {
    free(m[i]);
  }
  free(m);
}

void FreeMatrix(float **m)
{
  int      i;
  for (i=0; m[i]!=NULL; i++)
  {
    free(m[i]);
  }
  free(m);
}
