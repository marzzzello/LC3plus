/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

/*
	Tool for dynamic memory estimation
    Anisse Taleb, November 2003
  
*/

/* turn off stdlib function warnings in visual studio */
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "dynmem.h"

#define MAX_FUNC_NAME_LENGTH	128



struct C_Path{
               char		func_name[MAX_FUNC_NAME_LENGTH];
			   long		mem_usage;
			   long		mem_usage_acc;
               struct	C_Path *next;
			   struct	C_Path *parent;
};


typedef struct C_Path C_Path;

static C_Path *Max_Mem_path; /* Snap shot of the maximum memory usage path */
static C_Path *Head_path;
static C_Path *Curr_path;	  /* Current code path							*/


static long Max_mem_usage;

static C_Path *Static_mem_head;
static C_Path *Static_mem_curr;

/*
 * initialisation. Must be called at start, sets up all the
 * structures and resets counters.
 *
 * Called externally.
 */

static void path_free(C_Path * path)
{
	C_Path *path_end,*tmp;
	path_end = path;

	if(path == NULL) return;
	/* find last node */
	while(path_end->next != NULL) path_end = path_end->next;

	/* walk back */
	while(path_end != NULL) {
		tmp = path_end;
		path_end = path_end->parent;
		free(tmp);
	}

}
static void path_save(void)
{
	C_Path * node;
	C_Path * node_save;
	C_Path * previous;
	/* free the old path */
	path_free(Max_Mem_path);
	Max_Mem_path = NULL;
	/* */

	Max_Mem_path = (C_Path *) malloc(sizeof(C_Path));

	node = Head_path;
	node_save = Max_Mem_path;
	previous = NULL;
	while(node != NULL){

		strcpy(node_save->func_name,node->func_name);
		node_save->mem_usage = node->mem_usage;
		node_save->mem_usage_acc = node->mem_usage_acc;

		node_save->parent = previous;
		node_save->next	  =	NULL;

		if(node->next != NULL) {
			node_save->next = (C_Path *) malloc(sizeof(C_Path));
		}
		node = node->next;
		previous = node_save;
		node_save = node_save->next;
	} while(node != NULL);

}

void P_Dyn_Mem_Init(void)
{

	/* initialize the path */
	Head_path = (C_Path *) malloc(sizeof(C_Path));
	strcpy(Head_path->func_name,"---Top---");
	Head_path->mem_usage = 0;
	Head_path->mem_usage_acc = 0; /*memory usage accumulator */
	Head_path->parent = NULL;
	Head_path->next	  =	NULL;

	Curr_path = Head_path;
	/* Initilaize the memory consumptions */
	Max_mem_usage = 0;
	Max_Mem_path = NULL;
}

void P_Dyn_Mem_In(const char *func_name,
				long mem_usage)
{
#ifdef DEBUG_DYNMEM
	printf( "Entering: %s\n", func_name );
#endif

	/* Enter a new function, push on the stack  */
	if(Curr_path->next != NULL) {
		fprintf(stderr,"\n Something went wrong !!");
		exit(0);

	}
	Curr_path->next = (C_Path*) malloc (sizeof(C_Path));

	if(Curr_path->next == NULL) {
		fprintf(stderr,"\n Can't allocate memory");
		exit(0);
	}

	Curr_path->next->parent = Curr_path;
	Curr_path = Curr_path->next;
	Curr_path->next = NULL;
	/* save function name */
	strcpy(Curr_path->func_name,func_name);
	/* save memory usage */
	Curr_path->mem_usage = mem_usage;
	/* update memory usage accumulator */
	Curr_path->mem_usage_acc = Curr_path->parent->mem_usage_acc + mem_usage;
}

/* Add memory size to current function */
void P_Dyn_Mem_Add(long mem_usage)
{
#ifdef DEBUG_DYNMEM
	printf( "Staying in: %s\n", Curr_path->func_name);
#endif

	/* Staying in the current function, nothing to push  */

	/* save memory usage */
	Curr_path->mem_usage = Curr_path->mem_usage + mem_usage;
	/* update memory usage accumulator */
	Curr_path->mem_usage_acc = Curr_path->mem_usage_acc + mem_usage;
}


/* before any return */

void P_Dyn_Mem_Out(void)
{
	C_Path *tmp;
	if(Curr_path->mem_usage_acc > Max_mem_usage) {
		/* set new memory usage record */
		Max_mem_usage = Curr_path->mem_usage_acc;
		/* save snap shot */
		path_save();
	}

#ifdef DEBUG_DYNMEM
	printf( "Exiting: %s\n", Curr_path->func_name );
#endif

	/* delete last node */
	tmp = Curr_path;

	Curr_path = Curr_path->parent;
	Curr_path->next = NULL;
	free(tmp);
}


/* Write data and exit */

void P_Dyn_Mem_Exit(void)
{
	/*FILE *f_mem_statistics;*/
	C_Path *node;
	/*f_mem_statistics = fopen("mem_stat.txt","wt");
	fprintf(f_mem_statistics,"\n Maximum dynamic memory usage = %ld b, %ld kW \n Critical Memory Usage Path", Max_mem_usage, Max_mem_usage/(1000*2));*/
	/* fprintf(stderr,"\n Maximum dynamic memory usage = %ld kWords (%ld bytes)\n Critical Memory Usage Path", Max_mem_usage/(1000*2),Max_mem_usage); */
	fprintf(stderr,"\n Maximum dynamic memory usage = %.2f kWords (%ld bytes)\n Critical Memory Usage Path", (ceil((float)Max_mem_usage/20.f)/100.f), Max_mem_usage);

	node = Max_Mem_path;
        while (node != NULL) {
          /*fprintf(f_mem_statistics,"\n %-30s %10ld bytes  (+ %10ld bytes)",node->func_name, node->mem_usage_acc, node->mem_usage);*/
          fprintf(stderr,"\n %-30s %10ld bytes (+ %10ld bytes)",node->func_name, node->mem_usage_acc, node->mem_usage);
          node=node->next;
	}

	/*fclose(f_mem_statistics);*/
	fprintf(stderr,"\n");

#ifdef CTEST_MEASUREMENTS
  printf("<DartMeasurement name=\"%s\" type=\"%s\">%.3f</DartMeasurement>\n", "dynMem [kWords]", "numeric/float", Max_mem_usage/(1000.0*2));
#endif

	path_free(Max_Mem_path);
	path_free(Head_path);
}

void P_Dyn_Mem_Exit_noprint(void)
{
	path_free(Max_Mem_path);
	path_free(Head_path);
}
void P_Sta_Mem_Init(void)
{
	/* initialize the head object */
	Static_mem_head = (C_Path *) malloc(sizeof(C_Path));
	strcpy(Static_mem_head->func_name,"---Top---");
	Static_mem_head->mem_usage = 0;
	Static_mem_head->mem_usage_acc = 0;
	Static_mem_head->parent = NULL;
	Static_mem_head->next = NULL;
  Static_mem_curr = Static_mem_head;
}

void P_Sta_Mem_Add(const char *func_name,
                  long mem_usage)
{
	/* Enter a new function, push on the stack  */
	Static_mem_curr->next = (C_Path*) malloc (sizeof(C_Path));
	Static_mem_curr->next->parent = Static_mem_curr;
	Static_mem_curr = Static_mem_curr->next;
	Static_mem_curr->next = NULL;
	/* save function name */
	strcpy(Static_mem_curr->func_name,func_name);
	/* save memory usage */
	Static_mem_curr->mem_usage = mem_usage;
	/* update memory usage accumulator */
	Static_mem_curr->mem_usage_acc = Static_mem_curr->parent->mem_usage_acc + mem_usage;
}

void P_Sta_Mem_Exit(void)
{
	C_Path *node;
	/* fprintf(stderr,"\n Static memory usage = %ld kWords (%ld bytes)\n Detailed Memory Usage", Static_mem_curr->mem_usage_acc/(1000*2),Static_mem_curr->mem_usage_acc); */
	fprintf(stderr,"\n Static memory usage = %.2f kWords (%ld bytes)\n Detailed Memory Usage", (ceil((float)Static_mem_curr->mem_usage_acc/20.f)/100.f),Static_mem_curr->mem_usage_acc);

	node = Static_mem_head->next;
	while (node != NULL) {
		fprintf(stderr,"\n %-30s %10ld bytes (+ %10ld bytes)",node->func_name, node->mem_usage_acc, node->mem_usage);
		node=node->next;
	}

	fprintf(stderr,"\n");

#ifdef CTEST_MEASUREMENTS
  printf("<DartMeasurement name=\"%s\" type=\"%s\">%.3f</DartMeasurement>\n", "statMem [kWords]", "numeric/float", Static_mem_curr->mem_usage_acc/(1000.0*2));
#endif

	path_free(Static_mem_head);
}

void P_Sta_Mem_Exit_noprint(void)
{
	path_free(Static_mem_head);
}
