#ifndef GRAPH_H_
#define GRAPH_H_

#include "spmat.h"
/*
 * Takes a adjacency matrix A (sparse matrix),
 * number of nodes n as input and represents a graph
 */

typedef struct _graph {

	/* indexes vector - g length vector */
	int *index_vector;

	/* belongings vector  +-1 */
	/*int	*s_vector;*/

	/* current number of vertices */
	int current_size;

	/* ki vector */
	int *k_vector;

	/* Sum of ranks -> sum of rank vector */
	int M;

	/* A[g] */
	spmat *A_g;

	int *K_g_vector;

	double *f_vector;

	double shift;

	void (*free_graph)(struct _graph *G);

} graph;

/*graph* graph_create(int n, spmat* mother_A, int *s_vector, int *k_vector, int M);*/
graph* graph_create(spmat* mother_A, int *index_vector, int index_vector_size, int *k_vector, int M);

#endif /* GRAPH_H_ */
