/*
 * glist.h
 *
 *  Created on: 1 Sep 2020
 *      Author: avikef
 */

#include "graph.h"

#ifndef GLIST_H_
#define GLIST_H_

/* to represent a graph that sits
 * as a node on glist
 * */
typedef struct graph_node {
	graph* G;
    struct graph_node *next;
} graph_node;

typedef struct _glist {
	/* list size */
	/*int		n;*/

	/* pointer to the linked list */
	graph_node *current_graph_node;

	/* removes graph to list beginning */
	graph_node*	(*deque)(struct _glist *L);

	void	(*free_glist)(struct _glist *L);

	/* add graph to list beginning */
	void	(*enque)(struct _glist *L, graph_node *new_graph);



} glist;

glist* glist_allocate_list();

#endif /* GLIST_H_ */
