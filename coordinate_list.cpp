#include "coordinate_list.h"
#include <stdlib.h>

#define MIN(x,y) ((x>y)?y:x)
#define MAX(x,y) ((x<y)?y:x)

int clist_get_size(clist * lst)
{
	if (!lst) {
		return lst->list_size;
	}
	else { return 0; }
}

clist_node * clist_get_next_node(clist_node * node)
{
	return node->next;
}

void clist_new(clist * lst)
{
	lst = (clist*)malloc(sizeof(clist));
	lst->list_size = 0;
	lst->head = NULL;
	lst->tail = NULL;
}

void clist_delete(clist * lst)
{
	if ((!lst) || (lst->list_size == 0) || (!lst->head))return;

	clist_node *curr = lst->head;
	clist_node *nxt = NULL;
	while (curr) {
		nxt = curr->next;
		free(curr);
		curr = nxt;
	}
}

void clist_node_new(clist_node * node, double sx, double sy, double px, double py)
{
	node = (clist_node*)malloc(sizeof(clist_node));
	node->coo[0] = sx;
	node->coo[1] = sy;
	node->coo[2] = px;
	node->coo[3] = py;
	node->next = NULL;
}

clist_node * clist_get_head(clist * lst)
{
	if (!lst)return nullptr;
	return lst->head;
}

clist_node * clist_get_tail(clist * lst)
{
	if (!lst)return nullptr;
	return lst->tail;
}

clist_node *clist_get_ith_node(clist *lst, int idx)
{
	if (!lst||lst->list_size<=idx)return nullptr;
	if (idx = 0)return lst->head;
	clist_node *curr = lst->head;
	for (int i = 0; i < idx; i++) {
		curr = curr->next;
	}
	return curr;
}

void clist_insert_list_node(clist * lst, clist_node * node, int idx)
{
	if (!lst || !node)return;
	if (idx<0 || idx>lst->list_size)return;

	if (idx == 0) {
		node->next = lst->head;
		lst->head = node;
		lst->list_size++;
		return;
	}

	clist_node* curr = lst->head;
	clist_node* nxt = curr->next;
	int c_idx = 1;
	while (c_idx < idx) {
		curr = nxt;
		nxt = curr->next;
		c_idx++;
	}
	curr->next = node;
	node->next = nxt;
	if (idx == lst->list_size - 1)lst->tail = node;
	lst->list_size++;
}

void clist_exchange_list_node(clist * lst, int idxi, int idxj)
{
	if (!lst)return;
	if (idxi<0 || idxi>lst->list_size)return;
	if (idxj<0 || idxj>lst->list_size)return;

	int c_idx = 0;

	int i = MIN(idxi, idxj);
	int j = MAX(idxi, idxj);
	clist_node *A = lst->head, *B = lst->head;

	while (c_idx < i) {
		A = A->next;
		B = B->next;
		c_idx++;
	}
	while (c_idx < j) {
		B = B->next;
		c_idx++;
	}
	double tmp[4] = { A->coo[0],A->coo[1],A->coo[2],A->coo[3] };
	for (int i = 0; i < 4; i++) {
		A->coo[i] = B->coo[i];
		B->coo[i] = tmp[i];
	}
}

void clist_append_list_node(clist * lst, clist_node * node)
{
	if (!lst)return;
	if (lst->list_size == 0) {
		lst->head = node;
	}
	else {
		lst->tail->next = node;
	}
	lst->tail = node;
	lst->list_size++;

}

double clist_get_node_value(clist_node * node, int j)
{
	return node->coo[j];
}

double clist::ith_node_value(int idx, int j)
{
	return clist_get_node_value(clist_get_ith_node(this, idx),j);
}
