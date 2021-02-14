/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef DYNMEM_H
#define DYNMEM_H

#ifdef _MSC_VER
/* This disables warnings about anonymous temporary struct declarations */
#pragma warning(disable:4115)
#pragma warning(disable:4116)
#endif

#ifdef __cplusplus
extern "C" {
#endif

void P_Dyn_Mem_Init(void);
void P_Dyn_Mem_In(const char *func_name,
                  long mem_usage);
void P_Dyn_Mem_Add(long mem_usage);
void P_Dyn_Mem_Out(void);
void P_Dyn_Mem_Exit(void);
void P_Dyn_Mem_Exit_noprint(void);

void P_Sta_Mem_Init(void);
void P_Sta_Mem_Add(const char *func_name,
                  long mem_usage);
void P_Sta_Mem_Exit(void);
void P_Sta_Mem_Exit_noprint(void);

#ifdef DONT_COUNT_MEM

#define Dyn_Mem_Init()	
#define Dyn_Mem_Exit()	
#define Dyn_Mem_Exit_noprint()	
#define Dyn_Mem_In(a,b)	
#define Dyn_Mem_Add(b)	
#define Dyn_Mem_Out()	

#define DYN_MEM_IN		Dyn_Mem_In
#define DYN_MEM_ADD		Dyn_Mem_Add
#define DYN_MEM_OUT		Dyn_Mem_Out	
	
#define Sta_Mem_Init()
#define Sta_Mem_Exit()
#define Sta_Mem_Exit_noprint()
#define Sta_Mem_Add(a,b)

#else /* DONT_COUNT_MEM */

#define Dyn_Mem_Init	P_Dyn_Mem_Init
#define Dyn_Mem_Exit	P_Dyn_Mem_Exit 
#define Dyn_Mem_Exit_noprint	P_Dyn_Mem_Exit_noprint
#define Dyn_Mem_In		P_Dyn_Mem_In
#define Dyn_Mem_Add		P_Dyn_Mem_Add
#define Dyn_Mem_Out		P_Dyn_Mem_Out

#define DYN_MEM_IN		P_Dyn_Mem_In
#define DYN_MEM_ADD		P_Dyn_Mem_Add
#define DYN_MEM_OUT		P_Dyn_Mem_Out

#define Sta_Mem_Init  P_Sta_Mem_Init
#define Sta_Mem_Exit  P_Sta_Mem_Exit
#define Sta_Mem_Exit_noprint  P_Sta_Mem_Exit_noprint
#define Sta_Mem_Add   P_Sta_Mem_Add

#endif /* DONT_COUNT_MEM */

#ifdef __cplusplus
}
#endif

#endif
