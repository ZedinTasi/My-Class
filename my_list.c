#include "my_list.h"
#include <stdlib.h>
#include "stdafx.h"


struct my_list {
	int list_size;
	struct my_list_node *head;
	struct my_list_node *tail;
};
struct my_list_node {
	double val;
	struct my_list_node *next;
};

int get_list_size(struct my_list *lst) {
	if (lst != NULL) {
		return lst->list_size;
	}
	else { return 0; }
}

struct my_list * create_list()
{
	struct my_list *obj = (struct my_list*)malloc(sizeof(struct my_list));
	obj->head = NULL;
	obj->tail = NULL;
	obj->list_size = 0;
	return obj;
}

void destroy_list(struct my_list * lst)
{
	if ((lst == NULL) || (lst->list_size == 0) || (lst->head == NULL)) return;
	struct my_list_node *current = lst->head;
	struct my_list_node *nxt = NULL;
	while (current) {
		nxt = current->next;
		free(current);
		current = nxt;
	}
}

struct my_list * change_list_head(struct my_list *lst, struct my_list_node * obj)
{
	if (lst != NULL && obj != NULL) {
		lst->head = obj;
	}
	return lst;
}

struct my_list * change_list_tail(struct my_list * lst, struct my_list_node * obj)
{
	if (lst != NULL && obj != NULL) {
		lst->tail = obj;
	}
	return lst;
}

static struct my_list_node *create_list_node(double val) {
	struct my_list_node *o = (struct my_list_node *)malloc(sizeof(struct my_list_node));
	o->val = val;
	o->next = NULL;
	return o;
}

struct my_list_node *get_list_head(struct my_list * lst)
{
	if (lst != NULL) {
		return lst->head;
	}
	else { return NULL; }
}

struct my_list_node *get_list_tail(struct my_list * lst)
{
	if (lst != NULL) {
		return lst->tail;
	}
	else { return NULL; }
}

struct my_list_node *ith_list_node(struct my_list *lst, int idx) {
	
	if (lst == NULL || (idx - 1) > lst->list_size) return NULL;
	 
	struct my_list_node *curr = lst->head;
	if (idx == 1) { return curr; }
	for (int i = 0; i < idx-1; i++) {
		curr = curr->next;
	}
	return curr;
}


struct my_list_node * insert_to_list_index(struct my_list * lst, double val, int idx)
{
	//(void)val;
	if (lst == NULL) return NULL;
	if ((idx-1) > lst->list_size) return NULL;

	struct my_list_node *obj = create_list_node(val);

	struct my_list_node *nxt = lst->head;
	struct my_list_node *curr = lst->head;
	
	int c_idx = 0;
	while (c_idx < (idx - 1)) { 
		curr = nxt;
		nxt = nxt->next;
		c_idx++; 
	}
	if (nxt == curr) {
		obj->next = lst->head;
		lst->head = obj;
	} else {
		obj->next = nxt;
		curr->next = obj;
	}
	lst->list_size++;
	if ((idx - 1) > lst->list_size) lst->tail = obj;
	return obj;
}

struct my_list_node *exchange_list_node(struct my_list *lst, int idxi, int idxj) {

	if (lst == NULL || (idxj - 1) > lst->list_size || (idxj - 1) > lst->list_size) return NULL;
	if (idxi >= idxj) return NULL;
	/*
	double tmp = ith_list_node(lst, idxi)->val;
	ith_list_node(lst, idxi)->val = ith_list_node(lst, idxj)->val;
	ith_list_node(lst, idxj)->val = tmp;
	*/

	struct my_list_node *pre_i, *node_i, *pre_j = NULL, *node_j, *tmp;
	if (idxj - idxi > 1) {
		if (idxi == 1) {
			node_i = lst->head;			
			pre_j = ith_list_node(lst, idxj - 1);
			node_j = pre_j->next;
			tmp = node_j->next;
			lst->head = node_j;
		}
		else {
			pre_i = ith_list_node(lst, idxi - 1);
			node_i = pre_i->next;
			pre_j = ith_list_node(lst, idxj - 1);
			node_j = pre_j->next;
			tmp = node_j->next;
			pre_i->next = node_j;
		}

		node_j->next = node_i->next;
		pre_j->next = node_i;
		node_i->next = tmp;
	}
	else{
		if (idxi == 1) {
			node_i = lst->head;
			node_j = node_i->next;
			tmp = node_j->next;
			lst->head = node_j;

		}
		else {
			pre_i = ith_list_node(lst, idxi - 1);
			node_i = pre_i->next;
			node_j = node_i->next;
			tmp = node_j->next;
			pre_i->next = node_j;
		}
		node_j->next = node_i;
		node_i->next = tmp;
	}



	if (idxj == lst->list_size) {
		lst->tail = node_i;
	}
	return lst->head;
}

struct my_list_node *append_list_node(struct my_list *lst, double val) {
	struct my_list_node *obj = create_list_node(val);
	if (lst->list_size == 0) lst->head = obj;
	else lst->tail->next = obj;
	lst->tail = obj;
	lst->list_size++;
	return obj;
}

double get_node_value(struct my_list_node *obj) {
	if (obj != NULL) {
		return obj->val;
	}
	else return 0;
}

struct my_list_node * get_next_node(struct my_list_node * obj)
{
	if (obj != NULL) {
		return obj->next;
	}
	else return NULL;
}



struct my_list_node * change_node_next(struct my_list_node * main, struct my_list_node * obj)
{
	if (main != NULL) {
		main->next = obj;
		return main;
	}
	else return NULL;
}

double exchange_node_value(struct my_list_node * obj1, struct my_list_node * obj2)
{

	if (obj1 == NULL || obj2 == NULL || (obj1->val == obj2->val)) return 1;
	double tmp = obj1->val;
	obj1->val = obj2->val;
	obj2->val = tmp;
	return obj1->val;
}
