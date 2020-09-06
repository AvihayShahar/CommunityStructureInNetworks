/*
 * main.c
 *
 *  Created on: 3 Sep 2020
 *      Author: avikef
 */


#include "spmat.h"
#include <stddef.h>
#include "graph.h"
#include "glist.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

void printMat(spmat* A) {
	int i, col;
	node **arr, *node;

	arr = (((linked_List_Struct*)A->private)->node_list);
	node = *arr;

	for (i = 0; i < A->n; i++) {
		col = 0;
		while (col < A->n && node != NULL) {
			if (col != node->column) {
				printf("%d, ", 0);
			}
			else {
				printf("%d, ", 1);
				node = node->next;
			}
			col++;
		}
		printf("\n");
		arr++;
		node = *arr;
	}
}

void readMotherAIntoSpMat(spmat *graph, FILE *input_graph_file, int n, int* k_vector, int *M);
/*void initialize_mother_A(FILE *input_graph_file, spmat* mother_A, int n, /*linked_List_Struct *list_struct,*/ /*int* k_vector, int *M);*/
void divide_graph(graph *G, double *division_vector, int *index_vector1, int *index_vector2);
void algorithm2(graph *G, double *division_vector);

/* function takes stuff and repersents a graph using sp_mat */
void readMotherAIntoSpMat(spmat* mother_A, FILE *input_graph_file, int n, int *k_vector, int *M) {
	int *vertex_vector, *insert_vector, i, j, curr_edges_count, prev_index, counter, current_index,
	*start_vertex, *start_insert;
	linked_List_Struct *mother_list_struct;
	node **tmp_nodes_list;

	mother_list_struct = (linked_List_Struct*)mother_A->private;
	tmp_nodes_list = mother_list_struct->node_list;

	/* calloc initializes to 0 */
	vertex_vector = (int*)malloc(n*sizeof(int));
	insert_vector = (int*)calloc(n, sizeof(int));

	start_vertex = vertex_vector;
	start_insert = insert_vector;

	for (i = 0; i < n; i++) {
		prev_index = 0;
		current_index = 0;

		/* get number of edges for vertex i */
		fread(&curr_edges_count, sizeof(int), 1, input_graph_file);
		*k_vector = curr_edges_count;
		k_vector++;
		*M += curr_edges_count;

		/* read connected vertexes */
		fread(vertex_vector, sizeof(int), curr_edges_count, input_graph_file);


		for (j = 0; j < curr_edges_count; j++) {
			current_index = *vertex_vector;
			insert_vector += current_index - prev_index;
			*insert_vector = 1;
			vertex_vector++;
			prev_index = current_index;
		}

		insert_vector = start_insert;
		vertex_vector = start_vertex;

		/* NEED TO EDIT SPMAT TO GET VALUES OF 1 AND 0.
		 * ADD ROW TAKES A DOUBLE* VECTOR */
		mother_A->add_row(mother_A, insert_vector, i);

		/* resetting vector to 0's */
		memset(insert_vector, 0, n*sizeof(int));
	}
	mother_list_struct->node_list = tmp_nodes_list;
	free(vertex_vector);
	free(insert_vector);
}

/* initiating graph_spmat */
//void initialize_mother_A(FILE *input_graph_file, spmat* mother_A, int n/*, linked_List_Struct *list_struct*/, int* k_vector, int *M) {
//	node **tmp_nodes_list;
//	linked_List_Struct *list_struct;
//	mother_A = spmat_allocate_list(n);
//
//	list_struct = (linked_List_Struct*)mother_A->private;
//
//	tmp_nodes_list = list_struct->node_list;
//
//	readMotherAIntoSpMat(mother_A, input_graph_file, n, k_vector, M);
//
//	list_struct->node_list = tmp_nodes_list;
//}

/*
 * argv[0] = program name
 * argv[1] = input file path
 * argv[2] = output file path
 */
int main(int argc, char* argv[]) {
	FILE *input_graph_file, *output_file;
	spmat *mother_A;
	char *graph_file_path, *output_file_path;
	double *division_vector, *division_vector_ptr;
	int n, *k_vector, M, /*isDividible, */*initial_index_vector, i, *ptr, count1, count2 ,*index_vector1,*index_vector2;
	/*linked_List_Struct *list_struct;*/
	glist *O, *P;
	graph *mother_G, *g1, *g2;
	graph_node *mother_node, *g1_node, *g2_node, *current_graph_node;
	node **tmp_nodes_list;
	linked_List_Struct *mother_list_struct;

	/* for the randomizing thingy*/
	srand(time(NULL));

	/* Starting the run with creating the "home" graph A */
	graph_file_path = argv[1];
	input_graph_file = fopen(graph_file_path, "rb");

	fread(&n, sizeof(int), 1, input_graph_file);

	/* rank vector */
	k_vector = (int*)malloc(n*sizeof(int));

	M = 0;

	/* SETTING MOTHER A */
	mother_A = spmat_allocate_list(n);
	readMotherAIntoSpMat(mother_A, input_graph_file, n, k_vector, &M);
	fclose(input_graph_file);
	/*printMat(mother_A);*/

	/* Now we are finished with inserting the first A into an spmat - mother_A */

	/* WE HAVE M, n, mother_A, k_vector
	 * M we have here
	 * mother_A we have here
	 * k_vector we have here
	 * and n we have here AND it is on the mother_A field!!!
	 *
	 * WE NEED: ?!?!?!?
	 * TO CREATE A GRAPH WITH AN s_vector!!
	 *
	 * */
	/*initial_s_vector = (int*)malloc(sizeof(int)*n);
	memset(initial_s_vector, 1, sizeof(int)*n);*/

	initial_index_vector = (int*)malloc(sizeof(int)*n);
	ptr = initial_index_vector;

	for (i = 0; i < n; i++) {
		*ptr = i;
		ptr++;
	}

	mother_G = graph_create(mother_A, initial_index_vector, n, k_vector, M);

	/*printf("\n");
	printf("\n");*/
	/*printMat(mother_G->A_g);*/


	mother_node = (graph_node*)malloc(sizeof(graph_node));
	mother_node->G = mother_G;
	mother_node->next = NULL;

	/* ALGORITHM 3 STARTS HERE!! */
	/* CREATING O, P G-LISTS */

	P = glist_allocate_list();
	O = glist_allocate_list();
	P->enque(P, mother_node);

	while (P->current_graph_node != NULL) {
		current_graph_node = P->deque(P);
		/*printf("\n");
		printMat(current_graph_node->G->A_g);*/
		/* CAN USE N SIZED VECTOR AND SEND IT WITH A G SIZE INPUT */
		division_vector = (double*)malloc(sizeof(double)*(current_graph_node->G->current_size));

		algorithm2(current_graph_node->G, division_vector);
		division_vector_ptr = division_vector;

		count1 = 0;
		count2 = 0;

		for (i = 0; i < current_graph_node->G->current_size; i++) {
			if (*division_vector == 1) {
				count1++;
			}
			else {
				count2++;
			}
			division_vector++;
		}
		division_vector = division_vector_ptr;

		if (count1 == 0 || count2 == 0) {
			O->enque(O, current_graph_node);
		}
		else {
			index_vector1 = (int*)malloc(sizeof(int)*count1);
			index_vector2  = (int*)malloc(sizeof(int)*count2);

			divide_graph(current_graph_node->G, division_vector,index_vector1,index_vector2);

			g1 = graph_create(mother_A, index_vector1, count1, k_vector, M);
			g2 = graph_create(mother_A, index_vector2, count2, k_vector, M);

			g1_node = (graph_node*)malloc(sizeof(graph_node));
			g1_node->G = g1;
			g1_node->next = NULL;

			g2_node = (graph_node*)malloc(sizeof(graph_node));
			g2_node->G = g2;
			g2_node->next = NULL;

			if(count1 == 1){
				O->enque(O, g1_node);
			}else{
				P->enque(P, g1_node);
			}
			if(count2 == 1){
				O->enque(O, g2_node);
			}else{
				P->enque(P, g2_node);
			}
		}
		free(division_vector);
	}

	return 1;
}

void divide_graph(graph *G, double *division_vector, int *index_vector1, int *index_vector2) {
	int i;
	for (i = 0; i <G->current_size; i++) {
		if (*division_vector == 1) {
			*index_vector1 = *(G->index_vector);
			index_vector1++;
		}else{
			*index_vector2 = *(G->index_vector);
			index_vector2++;
		}
		division_vector++;
	}
}
