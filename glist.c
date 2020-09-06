/*
 * glist.c
 *
 *  Created on: 1 Sep 2020
 *      Author: avikef
 */
#include "glist.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* THIS IS ACTUALLY A GODDAMN STACK1!!!
 * NOT A LINKED LIST THAT REPRESENTS A QUEUE!!!
 *
 *  */

glist* glist_allocate_list();
void free_glist(glist *L);
graph_node* deque(glist *L);
void enque(glist *L, graph_node *new_graph);


glist* glist_allocate_list(/*graph *G*/) {
	glist *L;
	/*graph_node *start_node;*/

	L = (glist*)malloc(sizeof(glist));
	/*L->current_graph_node = (graph_node*)malloc(sizeof(graph_node));*/


	/*current_graph_node = (graph_node*)malloc(sizeof(current_graph_node));*/
/*
	*if (G != NULL) {
		*start_node = (graph_node*)malloc(sizeof(graph_node));
		*start_node->G = G;
		*start_node->next = NULL;
		 *
		 * SHOULD WE ASSIGN NULL!?@?@!?@!?@!?@!?@!?@!?!?!?!?!??!
		 *
		 * */
	/*	L->current_graph_node = start_node;
	}*/

	L->deque = &deque;
	L->enque = &enque;
	L->free_glist = &free_glist;

	return L;
}

void free_glist(glist *L) {
	graph_node *next, *current;

	current = L->current_graph_node;
	next = current->next;

	while (current != NULL) {
		current->G->free_graph(current->G);
		free(current);
		current = next;
		next = current->next;
	}
	free(current);
}

graph_node* deque(glist *L) {
	graph_node *removed;

	removed = L->current_graph_node;
	L->current_graph_node = L->current_graph_node->next;

	removed->next = NULL;

	return removed;
}

void enque(glist *L, graph_node *new_graph) {
	new_graph->next = L->current_graph_node;
	L->current_graph_node = new_graph;
}
