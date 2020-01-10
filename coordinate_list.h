#pragma once

struct clist {
	int list_size;
	struct clist_node *head;
	struct clist_node *tail;

	double ith_node_value(int idx, int j);

};
struct clist_node {
	double coo[4];
	struct clist_node *next;
};


int clist_get_size(clist *lst);
clist_node * clist_get_next_node(clist_node *node);
void clist_new(clist *lst);
void clist_delete(clist *lst);
void clist_node_new(clist_node *node, double sx, double sy, double px, double py);

clist_node *clist_get_head(clist *lst);
clist_node *clist_get_tail(clist *lst);
clist_node *clist_get_ith_node(clist *lst, int idx);

void clist_insert_list_node(clist *lst, clist_node *node, int idx);
void clist_exchange_list_node(clist *lst, int idxi, int idxj);
void clist_append_list_node(clist *lst, clist_node *node);
double clist_get_node_value(clist_node *node, int j);
