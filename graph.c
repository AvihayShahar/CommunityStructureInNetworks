/*
graph* create_mother_graph(int n, spmat *mother_A, int *k_vector, int M) {
	graph *G;
	int *index_vector, *s_vector, i, *ptr;

	G = (graph*)malloc(sizeof(graph));
	index_vector = (int*)malloc(n*sizeof(int));
	ptr = index_vector;
	for (i = 0; i < n; i++) {
		*index_vector = i;
		index_vector++;
	}
	index_vector = ptr;

	s_vector = (int*)malloc(n*sizeof(int));
	ptr = s_vector;
	for (i = 0; i < n; i++) {
		*s_vector = 1;
		s_vector++;
	}
	s_vector = ptr;

	 G->prev = mother_A;
	G->current_size = n;
	 G->previous_size = n;
	G->A_g = mother_A;
	G->index_vector = index_vector;
	G->s_vector = s_vector;
	G->k_vector = k_vector;
	G->M = M;
	 G->free = &free;

	return G;
}
*/

void printVecInt(int *vector){
	int i;
	for(i = 0; i < 20; i++){
		printf("%d\n", *vector);
		vector++;
	}
	printf("%s", "\n");
}

void printVecDouble(double *vector){
	int i;
	for(i = 0; i < 20; i++){
		printf("%f\n", *vector);
		vector++;
	}
	printf("%s", "\n");
}
void printKmahmid(int *k_vector, double *b0) {
	int i;
	double result = 0;
	double k0 = (double)*k_vector;

	for (i = 0; i < 20; i++) {
		result += ((k0 * (double)*k_vector)/82)*(*b0);
		b0++;
		k_vector++;
	}
	printf("%lf", result);
}

#define IS_POSITIVE(X) ((X) > 0.00001)

double epsilon = 0.00001;

#include "glist.h"
#include "graph.h"
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

//graph* graph_create(spmat *mother_A, int *s_vector, int *k_vector, int M) {
//	graph *G;
//	int i, *index_vector, index_vector_size, *K_g_vector;
//	spmat *A_g;
//	double *f_vector;
//
//	G = (graph*)malloc(sizeof(graph));
//
//	index_vector_size = calculate_index_vector_size(s_vector, mother_A->n);
//	index_vector = (int*)malloc(index_vector_size*sizeof(int));
//	calculate_index_vector(s_vector, mother_A->n, index_vector);
//
//	/* CALCULATING Ki of G vector */
//	K_g_vector = (int*)malloc(index_vector_size*sizeof(int));
//	extract_k_vector(k_vector, K_g_vector, index_vector_size, index_vector);
//
//	f_vector = (double*)malloc(index_vector_size*sizeof(double));
//
//	A_g = spmat_allocate_list(index_vector_size);
//	A_g->n = index_vector_size;
//
//	G->index_vector = index_vector;
//	G->s_vector = s_vector;
//	G->current_size = index_vector_size;
//	G->k_vector = k_vector;
//	G->M = M;
//
//	/* Actually create A_g */
//	create_Ag(G , mother_A, A_g, index_vector, index_vector_size, K_g_vector , f_vector);
//
//	G->K_g_vector = K_g_vector;
//	G->f_vector = f_vector;
//	G->A_g = A_g;
//
//	return G;
//}

graph* graph_create(spmat *mother_A, int *index_vector, int index_vector_size, int *k_vector, int M);
int calculate_index_vector_size(int *s_vector, int n);
void calculate_index_vector(int *s_vector, int n, int *index_vector);
void create_Ag(graph *G , spmat *mother_A, spmat *A_g, int *index_vector, int g_size, int *K_g_vector, double *f_vector);
double calculate_fi__and_shift(int *K_g_vector, int i, int g_size, int *row, double sum_of_row_i_A, double *f_vector, int M);
double scalar_product(double *vector1, double *vector2, int n);
void algorithm2(graph *G, double *division_vector);
double* power_iteration(graph *G, double *result_vector);
double calculate_modularity(double *division_vector, graph *G);
void modularity_maximization(double *division_vector, graph *G);
double optimized_score_i(graph *G, double *division_vector, int i, int d_i, node *node);
double calculate_eigen_value(double *bk, graph *G);
void sign_transformation(double *division_vector, double *result_vector, int g_size);
void calculate_shift_fi_part(double *result_vector, double *b0, double *f_vector, double shift, int g_size);
void calculate_fi_part(double *result_vector, double *b0, double *f_vector, int g_size);
void calculate_K_g_part(double *b0 , int *K_g_vector, double *result_vector, int M, int g_size);
void extract_k_vector(int *k_vector, int *K_g_vector, int g_size, int *index_vector);
void divideVectorByK(double *vector, double K, int N);
int isLegal(double *vector1, double *vector2, int N);
void randomizeVector(int n, double* vector);

graph* graph_create(spmat *mother_A, int *index_vector, int index_vector_size, int *k_vector, int M) {
	graph *G;
	int *K_g_vector;
	spmat *A_g;
	double *f_vector;

	G = (graph*)malloc(sizeof(graph));

	/* CALCULATING Ki of G vector */
	K_g_vector = (int*)malloc(index_vector_size*sizeof(int));
	extract_k_vector(k_vector, K_g_vector, index_vector_size, index_vector);

	f_vector = (double*)malloc(index_vector_size*sizeof(double));

	A_g = spmat_allocate_list(index_vector_size);
	A_g->n = index_vector_size;

	G->index_vector = index_vector;
	/*G->s_vector = s_vector;*/
	G->current_size = index_vector_size;
	G->k_vector = k_vector;
	G->M = M;

	/* Actually create A_g */
	create_Ag(G , mother_A, A_g, index_vector, index_vector_size, K_g_vector, f_vector);
	/*printf("\n");
	printMat(A_g);
	printf("\n");
	printVecDouble(f_vector);
	printf("\n");
	printf("%lf", G->shift);*/
	G->K_g_vector = K_g_vector;
	G->f_vector = f_vector;
	/*printVecDouble(f_vector);*/
	G->A_g = A_g;

	return G;
}

int calculate_index_vector_size(int *s_vector, int n) {
	int i, size;
	size = 0;

	for (i = 0; i < n; i++) {
		if (*s_vector == 1) {
			size++;
		}
		s_vector++;
	}
	return size;
}

void calculate_index_vector(int *s_vector, int n, int *index_vector) {
	int i;

	for (i = 0; i < n; i++) {
		if (*s_vector == 1) {
			*index_vector = i;
			index_vector++;
		}
		s_vector++;
	}
}

void create_Ag(graph *G , spmat *mother_A, spmat *A_g, int *index_vector, int g_size, int *K_g_vector, double *f_vector) {
	int i, j, *index_vector_ptr, *row, *row_start, current_row, prev_row, *index_vector_start, row_i_of_A;
	linked_List_Struct *mother_list_struct, *A_g_struct;
	node *curr_node, **list_mother_A, **start_list_A_g;
	double max_shift, current_shift;

	/* reads the private into the list struct */
	mother_list_struct = (linked_List_Struct*)mother_A->private;
	list_mother_A = mother_list_struct->node_list;

	A_g_struct = (linked_List_Struct*)A_g->private;
	start_list_A_g = A_g_struct->node_list;

	row = (int*)calloc(g_size, sizeof(int));
	index_vector_start = index_vector;
	row_i_of_A = 0;
	current_shift = 0.0;
	max_shift = 0.0;
	prev_row = 0;
	row_start = row;
	index_vector_ptr = index_vector;

	for (i = 0; i < g_size; i++) {
		current_row = *index_vector_ptr;
		list_mother_A += current_row - prev_row;
		curr_node = *list_mother_A;
		j = 0;

		while (curr_node != NULL && j < g_size) {
			if (curr_node->column == *index_vector) {
				row_i_of_A++; /* this is used to calculate fi. This is the sum of row i of A: Sum[j] Aij */
				*row = 1;
				curr_node = curr_node->next;
				index_vector++;
				j++;
			}
			else if (curr_node->column < *index_vector) {
				curr_node = curr_node->next;
			}
			else {
				index_vector++;
				j++;
			}
			row++;
		}
		/* THIS IS CALCULATING fi coordinate of f_vector */
		/* THIS IS CALCULATING Ki*Kj/M and subtracting it from sum of row i of A[g]
		*  */

		row = row_start;
		current_shift = calculate_fi__and_shift(K_g_vector, i, g_size, row, (double)row_i_of_A, f_vector, G->M);
		if (current_shift > max_shift) {
			max_shift = current_shift;
		}

		/* NEED TO PROMOTE F POINTER */
		f_vector++;
		row_i_of_A = 0;

		/* NOW WE NEED TO ADD ROW TO MATRIX
		 * NEED TO RESET POINTER ON LIST STRUCT OF A_g
		 * */
		A_g->add_row(A_g, row, i);

		/* resetting vector to 0's */
		memset(row, 0, g_size*sizeof(int));
		index_vector = index_vector_start;
		index_vector_ptr++;
		prev_row = current_row;
	}

	/*printVecDouble(f_vector);*/

	G->shift = max_shift;

	/* reset pointer */
	A_g_struct->node_list = start_list_A_g;
}

double calculate_fi__and_shift(int *K_g_vector, int i, int g_size, int *row, double sum_of_row_i_A, double *f_vector, int M) {
	int j;
	double ki, total_K_g, shift_i, current_K_j, fi, /*diagonal_i, */double_M;

	double_M = (double)M;

	shift_i = 0.0;
	total_K_g = 0.0;
	ki = (double)K_g_vector[i];

	for (j = 0; j < g_size; j++) {
		current_K_j = (ki* ((double)*K_g_vector)) / double_M;
		total_K_g += current_K_j; /* for the fi calculation */
		if (j != i) {
			shift_i += fabs((double)*row - current_K_j);
		}
		/*else {
			diagonal_i = (double)*row;
		}*/
		K_g_vector++;
		row++;
	}
	fi = sum_of_row_i_A - total_K_g;
	*f_vector = fi;
	/*printf("sum_of_row_i_A: %lf\n", sum_of_row_i_A);
	printf("total_K_g: %lf\n", total_K_g);

	printf("f: %lf\n", fi);*/

	shift_i += fabs(-((ki*ki)/double_M) - fi);
	/*printf("shift: %lf\n", shift_i);*/

	return shift_i;
}

double scalar_product(double *vector1, double *vector2, int n) {
	int i;
	double product;
	product = 0.0;

	for (i = 0; i < n; i++) {
		product += (*vector1)*(*vector2);
		vector1++;
		vector2++;
	}
	return product;
}

void algorithm2(graph *G, double *division_vector) {
	double *power_iteration_result_vector, modularity;
	int g_size, i;

	g_size = G->current_size;
	power_iteration_result_vector = (double*)malloc(sizeof(double)*g_size);

	power_iteration_result_vector = power_iteration(G, power_iteration_result_vector);
	/*printVecDouble(power_iteration_result_vector);*/
	sign_transformation(division_vector, power_iteration_result_vector, g_size);
	/*printVecDouble(division_vector);*/
	modularity_maximization(division_vector, G);
	printVecDouble(division_vector);
	printf("\n");
	exit(34232);
	free(power_iteration_result_vector);

	modularity = calculate_modularity(division_vector, G);
	if (!IS_POSITIVE(modularity)) {
		for (i = 0; i < G->current_size; i++) {
			*division_vector = 1;
		}
		division_vector++;
	}

	/*return IS_POSITIVE(modularity);*/
}

double* power_iteration(graph *G, double *result_vector){
	double *b0, *tmp, *f_vector, eigen_value;
	int legal, *K_g_vector, g_size, i;
	/*spmat *A_g;*/


	g_size = G->current_size;
	K_g_vector = G->K_g_vector;
	f_vector = G->f_vector;

	b0 = (double*)malloc(g_size*sizeof(double));
	randomizeVector(g_size, b0);
	/*printVecDouble(b0);*/

	/*result_vector = (double*)malloc(g_size*sizeof(double));*/
	/*fi_shift_vector = (double*)malloc(G->current_size*sizeof(double));*/

	legal = 0;
	while (legal == 0) {

		G->A_g->mult(G->A_g, b0, result_vector);
		calculate_K_g_part(b0, K_g_vector, result_vector,G->M, g_size);
		calculate_shift_fi_part(result_vector, b0, f_vector, G->shift, g_size);

		/* HERE COMES THE NORMALIZATION
		 * NEED TO CHECK WHAT IS GOING ON WITH LIFE AND THIS NORMALIZION
		 * MAYBE WE NEED TO BRIBE MOSHE
		 *
		 *
		 *  */
		divideVectorByK(result_vector, sqrt(scalar_product(result_vector, result_vector, g_size)), g_size);
		legal = isLegal(result_vector, b0, g_size);

		/* swap pointers */
		tmp = b0;
		b0 = result_vector;
		result_vector = tmp;
	}

	/* MOSHE SAID TO MAXIMIZE MODULARITY EVEN IF THE EIGENVALUE
	 * IS NEGATIVE!!!!
	 * IT MEANS THIS WHOLE SECTION IS IRRELAVENT!!
	 *  */


	eigen_value = calculate_eigen_value(result_vector, G);
	if (!IS_POSITIVE(eigen_value)) {
		for (i = 0; i < g_size; i++) {
			*result_vector = 1;
			result_vector++;
		}
	}
	/*s_vector = (int*)malloc(sizeof(int)*g_size);*/
	/*sign_transformation(s_vector, result_vector, g_size);*/ /* creates s partition vector */
	free(b0);
	return result_vector;
}

double calculate_modularity(double *division_vector, graph *G) {
	int g_size, i;
	double *B_g_s, *ptr, s_B_g_s;

	g_size = G->current_size;
	B_g_s = (double*)malloc(sizeof(double)*g_size);
	ptr = B_g_s;

	G->A_g->mult(G->A_g, division_vector, B_g_s);
	calculate_K_g_part(division_vector, G->K_g_vector, B_g_s, G->M, g_size);
	calculate_shift_fi_part(B_g_s, division_vector, G->f_vector, 0, g_size);

	s_B_g_s = 0.0;
	for (i = 0; i < g_size; i++) {
		s_B_g_s += (*ptr) * (*division_vector);
		ptr++;
		division_vector++;
	}
	free(B_g_s);
	return s_B_g_s;
}

void modularity_maximization(double *division_vector, graph *G) {
	int isFirst, maxImproveIndex, *indices_ptr, *unmoved_ptr, i, *unmoved, g_size, j_score_index, out_index, *indices;
	double deltaQ, maxScore, currScore, *improve, *improve_ptr, maxImprove, *s_ptr, Q0;
	node *node, **arr, **arr_ptr;

	/*Q0 = calculate_initial_modularity(s_vector, G);*/
	g_size = G->current_size;

	/* WE CALL THIS HERE SO WE CAN PROMOTE POINTER ON A_g
	 * SO WE CAN AVOID JUMPING ACCROSS ROWS IN EACH ITERATION
	 *
	 *  */
	arr = ((linked_List_Struct*)G->A_g->private)->node_list;
	node = *arr;
	arr_ptr = arr;

	/* calloc reset to 0
	 * cell i represents vertex i moved if unmoved[i] == 1
	 * else
	 * vertex i unmoved (and moved[i] == 0)
	 *
	 * s = (1,-1,1,1,1,-1,.....)
	 *
	 * */
	unmoved = (int*)calloc(g_size, sizeof(int));
	unmoved_ptr = unmoved;
	indices = (int*)malloc(g_size*sizeof(int));
	indices_ptr = indices;
	improve = (double*)malloc(g_size*sizeof(double));
	improve_ptr = improve;
	/* calculates modularity (Q) difference
	 * for all vertices i that unmoved[i] == 0
	 *
	 *  WE ENDED UP TAKING MAX SCORE AT EVERY ITERATION
	 *  */


	/*unmoved = unmoved_ptr;
	memset(unmoved, 0, g_size*sizeof(int));
	unmoved = unmoved_ptr;
	memset(unmoved, 0, g_size*sizeof(int));*/

	do {
		for (out_index = 0; out_index < g_size; out_index++) {
				/*j_score_index = 0;*/
				/*maxScore = 0.0;*/
				isFirst = 1; /* BOOLEAN */
				s_ptr = division_vector;
				Q0 = calculate_modularity(division_vector, G);
				for (i = 0; i < g_size; i++) {
					/**s_ptr *= -1;*/
					if (*unmoved != 1) {
						*s_ptr *= -1;
						currScore = calculate_modularity(division_vector, G) - Q0;
						/*currScore = optimized_score_i(G, division_vector, i, ((int)*s_ptr), node);*/
						*s_ptr *= -1;
						/*arr++;
						node = *arr;*/
						if (isFirst == 1) {
							j_score_index = i;
							maxScore = currScore;
							isFirst = 0;
						}
						else if (currScore > maxScore) {
							j_score_index = i;
							maxScore = currScore;
						}
						/*s_ptr += 1;*/
					}
					/**s_ptr *= -1;*/
					s_ptr += 1;
					unmoved++;
					arr++;
					node = *arr;
				}
				/* RESET POINTERS */
				s_ptr = division_vector;
				arr = arr_ptr;
				node = *arr;

				unmoved = unmoved_ptr;
				unmoved[j_score_index] = 1;
				division_vector[j_score_index] *= -1;
				*indices = j_score_index;
				indices++;

				if (out_index == 0) {
					*improve = maxScore;
					maxImprove = maxScore;
					maxImproveIndex = out_index;
				}
				else {
					*improve = *(improve-1) + maxScore;
					if (maxImprove < *improve) {
						maxImprove = *improve;
						maxImproveIndex = out_index;
					}
				}
				improve++;
			}

			indices -= 1;
			for (i = g_size-1; i > maxImproveIndex; i--) {
				division_vector[*indices] *= -1;
				indices--;
			}

			if (maxImproveIndex == g_size-1) {
				/* DELTA Q == 0!!
				 * CANT PARTITION ANY LONGER!!
				 *
				 *  */
				deltaQ = 0.0;
			}
			else {
				/* CHECK WHAT THE HELL IS GOING ON!!!?!!??*/
				deltaQ = maxImprove;
			}
			/* ??? */
			unmoved = unmoved_ptr;
			memset(unmoved, 0, g_size*sizeof(int));
			improve = improve_ptr;
			memset(improve, 0.0, g_size*sizeof(double));
			indices = indices_ptr;
			memset(indices, 0, g_size*sizeof(int));
	}
	while (IS_POSITIVE(deltaQ));
}

/* MOHAMAD CALCULATION */
double optimized_score_i(graph *G, double *division_vector, int i, int d_i, node *node) {
	int *k_ptr, index, g_size, node_col, count;
	double row_i_A_minus_D_times_dj_sum, k_i, M, *s_ptr;

	row_i_A_minus_D_times_dj_sum = 0.0;
	k_ptr = G->k_vector;
	k_i = (double)k_ptr[i];
	index = 0;
	g_size = G->current_size;
	M = (double)G->M;
	s_ptr = division_vector;

	node_col = node->column;
	while (node != NULL && index < g_size) {
		if (node_col > index) {
			while (node_col > index) {
				row_i_A_minus_D_times_dj_sum -= (k_i*((double)*k_ptr)/M)*((double)*s_ptr);
				k_ptr++;
				index++;
				s_ptr++;
			}
			row_i_A_minus_D_times_dj_sum += (1-k_i*((double)*k_ptr)/M)*((double)*s_ptr);
		}
		else {
			row_i_A_minus_D_times_dj_sum += (1-k_i*((double)*k_ptr)/M)*((double)*s_ptr);
		}
		k_ptr++;
		index++;
		s_ptr++;
		node = node->next;
		if (node != NULL) {
			node_col = node->column;
		}
	}

	for (count = index; count < g_size; count++) {
		row_i_A_minus_D_times_dj_sum -= (k_i*((double)*k_ptr)/M)*((double)*s_ptr);
		k_ptr++;
		index++;
		s_ptr++;
	}

	return 4*(((double)d_i)*row_i_A_minus_D_times_dj_sum + ((double)k_i)*((double)k_i)/M);
}

double calculate_eigen_value(double *bk, graph *G) {
	double eigen_value, *result;
	int g_size;

	g_size = G->current_size;

	result = (double*)malloc(g_size*sizeof(double));

	/*for (i = 0; i < g_size; i++) {*/
	G->A_g->mult(G->A_g, bk, result);
	calculate_K_g_part(bk, G->K_g_vector, result, G->M, g_size);
	calculate_shift_fi_part(result, bk, G->f_vector, G->shift, g_size);
		/* HERE COMES THE NORMALIZATION
				 * NEED TO CHECK WHAT IS GOING ON WITH LIFE AND THIS NORMALIZION
				 * MAYBE WE NEED TO BRIBE MOSHE
				 *
				 *
				 *  */
	eigen_value = scalar_product(result, bk, g_size) / scalar_product(bk, bk, g_size);
	/*}*/

	free(result);
	return eigen_value - (G->shift);
}

/* according to s definition, the eigen vector should go through a sign transformation */
void sign_transformation(double *division_vector, double *result_vector, int g_size) {
	int i;

	for (i = 0; i < g_size; i++) {
		if (IS_POSITIVE(*result_vector)) {
			*division_vector = 1;
		}
		else {
			*division_vector = -1;
		}
		division_vector++;
		result_vector++;
	}
}

void calculate_shift_fi_part(double *result_vector, double *b0, double *f_vector, double shift, int g_size) {
	int i;

	for (i = 0; i < g_size; i++) {
		*result_vector += ((-(*f_vector) + shift) * (*b0));
		result_vector++;
		f_vector++;
		b0++;
	}
}

void calculate_fi_part(double *result_vector, double *b0, double *f_vector, int g_size) {
	int i;

	for (i = 0; i < g_size; i++) {
		*result_vector += (-(*f_vector)) * (*b0);
		result_vector++;
		f_vector++;
		b0++;
	}
}

void calculate_K_g_part(double *b0 , int *K_g_vector, double *result_vector, int M, int g_size) {
	int i, *k_ptr;
	double product;

	product = 0.0;
	k_ptr = K_g_vector;

	for (i = 0; i < g_size; i++) {
		product += ((double)*K_g_vector) * (*b0);
		K_g_vector++;
		b0++;
	}

	for (i = 0; i < g_size; i++) {
		*result_vector -= ((((double)*k_ptr) / (double)M) * product);
		k_ptr++;
		result_vector++;
	}
}

void extract_k_vector(int *k_vector, int *K_g_vector, int g_size, int *index_vector) {
	int i, current_index, prev_index;

	prev_index = 0;
	for (i = 0; i < g_size; i++) {
		current_index = *index_vector;
		k_vector += current_index - prev_index;
		*K_g_vector = *k_vector;
		index_vector++;
		K_g_vector++;
		prev_index = current_index;
	}
}

void divideVectorByK(double *vector, double K, int N){
	int i;

	for(i = 0; i < N; i++){
		*vector = (double)(*vector / K);
		vector++;
	}
}

int isLegal(double *vector1, double *vector2, int N){
	int i;
	double result;

	for(i = 0; i < N; i++){
		result = fabs(*vector1 - *vector2);
		if(result > epsilon){
			return 0;
		}
		vector1 = vector1 + 1;
		vector2= vector2 + 1;
	}
	return 1;
}

void randomizeVector(int n, double* vector){
	int i;

	for(i = 0; i < n; i++){
		*vector = (double)rand();
		vector++;
	}
}
