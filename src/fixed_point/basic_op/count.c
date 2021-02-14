/*
  ===========================================================================
   File: COUNT.C                                         v.2.3 - 30.Nov.2009
  ===========================================================================

            ITU-T   STL   BASIC   OPERATORS

            COMPLEXITY EVALUATION FUNCTIONS

   History:
   03 Nov 04   v2.0     Incorporation of new 32-bit / 40-bit / control
                        operators for the ITU-T Standard Tool Library as
                        described in Geneva, 20-30 January 2004 WP 3/16 Q10/16
                        TD 11 document and subsequent discussions on the
                        wp3audio@yahoogroups.com email reflector.
                        norm_s()      weight reduced from 15 to 1.
                        norm_l()      weight reduced from 30 to 1.
                        L_abs()       weight reduced from  2 to 1.
                        L_add()       weight reduced from  2 to 1.
                        L_negate()    weight reduced from  2 to 1.
                        L_shl()       weight reduced from  2 to 1.
                        L_shr()       weight reduced from  2 to 1.
                        L_sub()       weight reduced from  2 to 1.
                        mac_r()       weight reduced from  2 to 1.
                        msu_r()       weight reduced from  2 to 1.
                        mult_r()      weight reduced from  2 to 1.
                        L_deposit_h() weight reduced from  2 to 1.
                        L_deposit_l() weight reduced from  2 to 1.
   March 06    v2.1     Changed to improve portability.
   Dec 06      v2.2     Changed to specify frame rate using setFrameRate()
                        Adding WMOPS_output_avg() for global average computation
                        L_mls() weight of 5.
                        div_l() weight of 32.
                        i_mult() weight of 3.
  Nov 09       v2.3     STL2009

  Jan 12       v2.4     Add call stack functionality, best case and global operation counters.
  ============================================================================
*/


/*****************************************************************************
 *
 * This file contains functions for the automatic complexity calculation
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "stl.h"

#if WMOPS
/* Million frames per second */
static double frameRate = FRAME_RATE; /* default value : 10 ms */
#endif /* if WMOPS */

#if WMOPS
/* Global counter variable for calculation of complexity weight */
BASIC_OP multiCounter[MAXCOUNTERS];
BASIC_OP glob_multiCounter;
int currCounter=0; /* Zero equals global counter */
#endif /* if WMOPS */

void setFrameRate(int samplingFreq, int frameLength)
{
#if WMOPS
  if(frameLength > 0)
    {
        frameRate = samplingFreq / 1000000.0 / frameLength;
    }
    return;
#else
    (void)samplingFreq;
    (void)frameLength;
#endif /* if WMOPS */
}

#if WMOPS
/*
 * Below list is used for displaying the code profiling information in
 * the file which name is defined by CODE_PROFILE_FILENAME.
 * For further details see generic_WMOPS_output() function.
 * Attention, the ordering in this table must be kept in synchronisation
 * with the structure definition BASIC_OP.
 */
char* BasicOperationList[] =
{
   "add",           "sub",             "abs_s",         "shl",             "shr",
   "extract_h",     "extract_l",       "mult",          "L_mult",          "negate",
   "round",         "L_mac",           "L_msu",         "L_macNs",         "L_msuNs",
   "L_add",         "L_sub",           "L_add_c",       "L_sub_c",         "L_negate",
   "L_shl",         "L_shr",           "mult_r",        "shr_r",           "mac_r",

   "msu_r",         "L_deposit_h",     "L_deposit_l",   "L_shr_r",         "L_abs",
   "L_sat",         "norm_s",          "div_s",         "norm_l",          "move16",
   "move32",        "Logic16",         "Logic32",       "Test",            "s_max",
   "s_min",         "L_max",           "L_min",         "L40_max",         "L40_min",
   "shl_r",         "L_shl_r",         "L40_shr_r",     "L40_shl_r",       "norm_L40",

   "L40_shl",       "L40_shr",         "L40_negate",    "L40_add",         "L40_sub",
   "L40_abs",       "L40_mult",        "L40_mac",       "mac_r40",
   "L40_msu",       "msu_r40",         "Mpy_32_16_ss",  "Mpy_32_32_ss",    "L_mult0",
   "L_mac0",        "L_msu0",          "lshl",          "lshr",            "L_lshl",
   "L_lshr",        "L40_lshl",        "L40_lshr",      "s_and",           "s_or",

   "s_xor",         "L_and",           "L_or",          "L_xor",           "rotl",
   "rotr",          "L_rotl",          "L_rotr",        "L40_set",         "L40_deposit_h",
   "L40_deposit_l", "L40_deposit32",   "Extract40_H",   "Extract40_L",     "L_Extract40",
   "L40_round",     "L_saturate40",    "round40",       "IF",              "GOTO",
   "BREAK",         "SWITCH",          "FOR",           "WHILE",           "CONTINUE"

 , "L_mls",         "div_l",           "i_mult"
};
#endif /* if WMOPS */


#if WMOPS
const BASIC_OP op_weight =
{
    1,     1,     1,     1,     1,
    1,     1,     1,     1,     1,
    1,     1,     1,     1,     1,
    1,     1,     2,     2,     1,
    1,     1,     1,     3,     1,

    1,     1,     1,     3,     1,
    4,     1,     18,    1,     1,
    2,     1,     2,     2,     1,
    1,     1,     1,     1,     1,
    3,     3,     3,     3,     1,

    1,     1,     1,     1,     1,
    1,     1,     1,     2,
    1,     2,     2,     4,     1,
    1,     1,     1,     1,     1,
    1,     1,     1,     1,     1,

    1,     1,     1,     1,     3,
    3,     3,     3,     3,     1,
    1,     1,     1,     1,     1,
    1,     1,     1,     4,     4,
    4,     8,     3,     4,     4

  , 5,     32,    3
};
#endif /* if WMOPS */


Word32 TotalWeightedOperation (void);
Word32 DeltaWeightedOperation (void);


#if WMOPS
/* Counters for separating counting for different objects */

/**
 maxCounter: current number of counters. Each scope initialized with BASOP_sub_start() gets a own counter assigned.
 objectName: Name of each counter passed to BASOP_sub_start().
 fwc_corr:
 nbTimeObjectIsCalled: number of times a counter (object) is referenced in the current frame.
 */

static int maxCounter=0;
static char* objectName[MAXCOUNTERS+1];

static Word16 fwc_corr[MAXCOUNTERS+1];
static long int nbTimeObjectIsCalled[MAXCOUNTERS+1];

#define NbFuncMax  (5000L*20L)

/**
 funcid: current function call for each counter
 bc : best case for each counter and function call
 wc : worst case for each counter and function call
 nbframe: number of frames for each counter.
 glob_bc: global best case self time for each counter for current frame.
 glob_wc: global worst case self time for each counter for current frame.
 glob_sum_curr: global cummulative time for each counter for current frame.
 glob_sum_bc: global best case cummulative time for each counter for current frame.
 glob_sum_wc: global worst case cummulative time for each counter for current frame.

 total_wmops: total wmops self time for each counter for current frame.
 total_sum: total wmops cummulative time for each counter for current frame.
 LastWOper: values used for WMOPS deltas
 */

static Word32 funcid[MAXCOUNTERS+1]={0}, nbframe[MAXCOUNTERS+1] ={0}, nbcalls[MAXCOUNTERS+1]={0}, nbcalls_wf[MAXCOUNTERS+1]={0}, nbcalls_wf_tmp[MAXCOUNTERS+1]={0};
static Word32 glob_bc[MAXCOUNTERS+1]={0}, glob_wc[MAXCOUNTERS+1]={0}, glob_wf[MAXCOUNTERS+1]={0}, glob_wf_tmp[MAXCOUNTERS+1]={0};
static Word32 *wc[MAXCOUNTERS+1] = {NULL};
static float total_wmops[MAXCOUNTERS+1];
static float total_sum[MAXCOUNTERS+1];

static Word32 LastWOper[MAXCOUNTERS+1];

static Word16 call_tree[MAXCOUNTERS+1][MAXCOUNTERS+1];
static int sum_curr[MAXCOUNTERS+1];
static int sum_bc[MAXCOUNTERS+1];
static int sum_wc[MAXCOUNTERS+1];
static int sum_wf[MAXCOUNTERS+1];
static int sum_wf_tmp[MAXCOUNTERS+1];
static int glob_sum_curr[MAXCOUNTERS+1];
static int glob_sum_bc[MAXCOUNTERS+1];
static int glob_sum_wc[MAXCOUNTERS+1];
#if MAX_CALLERS_SAVED_FRAMES
#define MAX_CALLERS_PRINT 20
static float callers_frames[MAX_CALLERS_SAVED_FRAMES+1][MAXCOUNTERS+1];
static int callers_frames_nos[MAX_CALLERS_SAVED_FRAMES];
static float callers_totals[MAX_CALLERS_SAVED_FRAMES];
#endif

#endif /* if WMOPS */


#if WMOPS
static char* my_strdup(const char *s) {
/*
 * duplicates UNIX function strdup() which is not ANSI standard:
 * -- malloc() memory area big enough to hold the string s
 * -- copy string into new area
 * -- return pointer to new area
 *
 * returns NULL if either s==NULL or malloc() fails
 */
    char *dup;

    if (s == NULL)
        return NULL;

    /* allocate memory for copy of ID string (including string terminator) */
    /* NOTE: the ID strings will never be deallocated because there is no
             way to "destroy" a counter that is not longer needed          */
    if ((dup = (char *) malloc(strlen(s)+1)) == NULL)
        return NULL;

    return strcpy(dup, s);
}
#endif /* if WMOPS */

int getCounterId(const char *objectNameArg) {
#if WMOPS

    assert(maxCounter < MAXCOUNTERS-1);

  if(maxCounter>=MAXCOUNTERS-1) return 0;
  objectName[++maxCounter]=my_strdup(objectNameArg);
  return maxCounter;

#else /* if WMOPS */
  (void)objectNameArg;
  return 0; /* Dummy */
#endif /* if WMOPS */
}

#if WMOPS
int readCounterId(void) {
   return currCounter;
}
#endif /* if WMOPS */


#if WMOPS
char * readCounterIdName(void) {
   return objectName[currCounter];
}
#endif /* if WMOPS */

void setCounter( int counterId) {
#if WMOPS

    assert(counterId < MAXCOUNTERS);

   if( (counterId > maxCounter)
    || (counterId < 0)) {
      currCounter=0;
      return;
   }
   currCounter=counterId;
   call_occurred = 1;
#else
   (void)counterId;
#endif /* if WMOPS */
}


void incrementNbTimeObjectIsCalled( int counterId) {
#if WMOPS

    assert(counterId < MAXCOUNTERS);

  if( (counterId > maxCounter)
   || (counterId < 0)) {
      nbTimeObjectIsCalled[0]++;
      return;
    }
  nbTimeObjectIsCalled[counterId]++;
#else
  (void)counterId;
#endif /* if WMOPS */
}


#if WMOPS
static Word32 WMOPS_frameStat(void) {
/* calculate the WMOPS seen so far and update the global
   per-frame maximum (glob_wc)
 */
    Word32 tot;

    tot = TotalWeightedOperation ();
    if (tot > glob_wc[currCounter])
        glob_wc[currCounter] = tot;

    if (tot < glob_bc[currCounter])
        glob_bc[currCounter] = tot;

    /* check if fwc() was forgotten at end of last frame */
    if (tot > LastWOper[currCounter]) {
        if (!fwc_corr[currCounter]) {
            fprintf(stderr,
                    "count: operations counted after last fwc() for '%s'; "
                    "-> fwc() called\n",
                    objectName[currCounter]?objectName[currCounter]:"");
        }
        fwc();
    }

    return tot;
}
#endif /* if WMOPS */


#if WMOPS
static void WMOPS_clearMultiCounter(void) {
    int i;

    Word32 *ptr = (Word32 *) &multiCounter[currCounter];
    for( i = 0; i < (int)(sizeof (multiCounter[currCounter])/ sizeof (Word32)); i++) {
        *ptr++ = 0;
    }
}
#endif /* if WMOPS */


void ClearNbTimeObjectsAreCalled(void) {
#if WMOPS
    int i;

    for (i = 0; i < (int)(sizeof (multiCounter[currCounter])/ sizeof (Word32)); i++) {
        nbTimeObjectIsCalled[i] = 0;
    }
#endif /* if WMOPS */
}

Word32 TotalWeightedOperation (void) {
#if WMOPS
    int i;
    Word32 tot, *ptr;
    const Word32 *ptr2;

    tot = 0;
    ptr = (Word32 *) &multiCounter[currCounter];
    ptr2 = (const Word32 *) &op_weight;
    for (i = 0; i < (int)(sizeof (multiCounter[currCounter])/ sizeof (Word32)); i++) {
        tot += ((*ptr++) * (*ptr2++));
    }

    return ((Word32) tot);

#else /* if WMOPS */
    return 0; /* Dummy */

#endif /* if WMOPS */

}

Word32 DeltaWeightedOperation (void) {
#if WMOPS
    Word32 NewWOper, delta;

    NewWOper = TotalWeightedOperation ();
    delta = NewWOper - LastWOper[currCounter];
    LastWOper[currCounter] = NewWOper;
    return (delta);

#else /* if WMOPS */
    return 0; /* Dummy */

#endif /* if WMOPS */
}


void Init_WMOPS_counter (void) {
#if WMOPS
    int i;

    if(wc[currCounter] == NULL){
        wc[currCounter] = malloc(NbFuncMax * sizeof(Word32));
    }
    /* reset function weight operation counter variable */

    for (i = 0; i < NbFuncMax; i++) {
        wc[currCounter][i] = (Word32) 0;
    }
    glob_wc[currCounter] = 0;
    glob_bc[currCounter] = MAX_32;
    glob_wf[currCounter] = 0;
    glob_wf_tmp[currCounter] = 0;

    nbframe[currCounter] = 0;
    nbcalls[currCounter] = 0;
    nbcalls_wf[currCounter] = 0;
    nbcalls_wf_tmp[currCounter] = 0;
    total_wmops[currCounter] = 0.0;
    total_sum[currCounter] = 0.0;

    /* initially clear all counters */
    WMOPS_clearMultiCounter();
    LastWOper[currCounter] = 0;
    funcid[currCounter] = 0;

    /* Following line is useful for incrIf(), see control.h */
    call_occurred = 1;
    funcId_where_last_call_to_else_occurred=MAXCOUNTERS;

    sum_bc[currCounter] = MAX_32;
    sum_wc[currCounter] = 0;
    sum_curr[currCounter] = 0;
    sum_wf[currCounter] = 0;
    sum_wf_tmp[currCounter] = 0;

    for (i = 0; i < MAXCOUNTERS; i++)
    {
      call_tree[currCounter][i] = -1;
    }

    glob_sum_curr[currCounter] = 0;
    glob_sum_wc[currCounter] = 0;
    glob_sum_bc[currCounter] = MAX_32;

#endif /* if WMOPS */
}


Word32 Reset_WMOPS_counter (void) {
#if WMOPS
    Word32 tot = WMOPS_frameStat();

    /* increase the frame counter --> a frame is counted WHEN IT BEGINS */
    nbframe[currCounter]++;
    /* Call counter */
    nbcalls[currCounter] += funcid[currCounter];
    /* add wmops used in last frame to count, then reset counter */
    /* (in first frame, this is a no-op                          */
    total_wmops[currCounter] += (float)( tot * frameRate );
    total_sum[currCounter] += (float)( (glob_sum_curr[currCounter]) * frameRate );

    /* clear counter before new frame starts */
    WMOPS_clearMultiCounter();
    LastWOper[currCounter] = 0;
    funcid[currCounter] = 0;           /* new frame, set function id to zero */

    glob_sum_curr[currCounter] = 0;
    sum_curr[currCounter] = 0;

    return tot;
#else /* if WMOPS */
    return 0;
#endif /* if WMOPS */
}


Word32 fwc (void) {
/* function worst case */
#if WMOPS
    Word32 tot;

    tot = DeltaWeightedOperation ();

    if(wc[currCounter]){
        if (tot > wc[currCounter][funcid[currCounter]]){
            wc[currCounter][funcid[currCounter]] = tot;
        }
    }

    funcid[currCounter]++;
    /*make sure that BASOP_frame_update(); is put at the end of the main loop*/
    if (funcid[currCounter] >= NbFuncMax) {
      printf("to many function calls\n");
    }
    assert(funcid[currCounter] < NbFuncMax);

    return (tot);

#else /* if WMOPS */
    return 0; /* Dummy */

#endif /* if WMOPS */
}

void WMOPS_output (Word16 dtx_mode) {
#if WMOPS
    int i;
    Word32 tot, tot_wm, tot_wc;

    /* get operations since last reset (or init),
       but do not update the counters (except the glob_wc[] maximum)
       so output CAN be called in each frame without problems.
       The frame counter is NOT updated!
     */
    tot = WMOPS_frameStat();
    tot_wm = (Word32)(total_wmops[currCounter] + ((float) tot) * frameRate);

    fprintf (stdout, "%10s:WMOPS=%.3f",
         objectName[currCounter]?objectName[currCounter]:"",
         ((float) tot) * frameRate);

    if (nbframe[currCounter] != 0)
    {
        fprintf (stdout, "  Average=%.3f",
                 tot_wm / (float) nbframe[currCounter]);
    }
    fprintf (stdout, "  WorstCase=%.3f",
             ((float) glob_wc[currCounter]) * frameRate);

    /* Worst worst case printed only when not in DTX mode */
    if (dtx_mode == 0)
    {
        tot_wc = 0L;
        for (i = 0; i < funcid[currCounter]; i++)
            tot_wc += wc[currCounter][i];
        fprintf (stdout, "  WorstWC=%.3f", ((float) tot_wc) * frameRate);
    }
    fprintf (stdout, " (%d frames)\n", nbframe[currCounter]);
#else
    (void)dtx_mode;
#endif /* if WMOPS */
}


void WMOPS_output_avg (Word16 dtx_mode, Word32 *tot_wm, Word16 *num_frames) {
#if WMOPS
    int i;
    Word32 tot, tot_wc;

    /* get operations since last reset (or init),
       but do not update the counters (except the glob_wc[] maximum)
       so output CAN be called in each frame without problems.
       The frame counter is NOT updated!
     */
    tot = WMOPS_frameStat();
    *tot_wm = (Word32)(total_wmops[currCounter] + ((float) tot) * frameRate);
    *num_frames = (Word16)nbframe[currCounter];

    fprintf (stdout, "%10s:WMOPS=%.3f",
         objectName[currCounter]?objectName[currCounter]:"",
         ((float) tot) * frameRate);

    if (nbframe[currCounter] != 0)
    {
        fprintf (stdout, "  Average=%.3f",
                 *tot_wm / (float) nbframe[currCounter]);
    }
    fprintf (stdout, "  WorstCase=%.3f",
             ((float) glob_wc[currCounter]) * frameRate);

    /* Worst worst case printed only when not in DTX mode */
    if (dtx_mode == 0)
    {
        tot_wc = 0L;
        for (i = 0; i < funcid[currCounter]; i++)
            tot_wc += wc[currCounter][i];
        fprintf (stdout, "  WorstWC=%.3f", ((float) tot_wc) * frameRate);
    }
    fprintf (stdout, " (%d frames)\n", nbframe[currCounter]);
#else
    (void)dtx_mode;
    (void)tot_wm;
    (void)num_frames;
#endif /* if WMOPS */
}

void generic_WMOPS_output (Word16 dtx_mode, char *test_file_name)
{
#if WMOPS
   int		saved_value;
   int i;
   Word32	tot, tot_wm, tot_wc, *ptr;
   const Word32 *ptr2;
   Word40   grand_total;
   FILE	*WMOPS_file;

   saved_value = currCounter;

   /*Count the grand_total WMOPS so that % ratio per function group
   can be displayed. */
   grand_total = 0;
   for( currCounter = 0; currCounter <= maxCounter; currCounter++) {
      tot = WMOPS_frameStat();
      grand_total += tot;

   }


   if( (WMOPS_file=fopen(WMOPS_DATA_FILENAME,"a"))!=NULL) {

      printf( "opened file %s in order to print WMOPS for each function group.\n", WMOPS_DATA_FILENAME);

      /* Print the file header line. */
      fprintf (WMOPS_file, "Test file name\tFunction Name \tFrame\tNb Times Called\tWMOPS\t%% versus grand total");

      if (nbframe[saved_value] != 0)
         fprintf (WMOPS_file, "\tAverage");

      fprintf (WMOPS_file, "\tWorstCase");

      /* Worst worst case printed only when not in DTX mode */
      if (dtx_mode == 0)
         fprintf (WMOPS_file, "\tWorstWC");

      fprintf (WMOPS_file, "\n");

      /* Print the WMOPS for each Function Group by scanning
         all the function groups with currCounter index.*/
      for( currCounter = 0; currCounter <= maxCounter; currCounter++) {

         fprintf (WMOPS_file, "%s", test_file_name);
         fprintf (WMOPS_file, "\t%s",
                  objectName[currCounter] ? objectName[currCounter] : "");
         fprintf (WMOPS_file, "\t%d", nbframe[currCounter]);

         tot = WMOPS_frameStat();
         tot_wm = (Word32)(total_wmops[currCounter] + ((float) tot) * frameRate);

         fprintf (WMOPS_file, "\t\t%ld", nbTimeObjectIsCalled[currCounter]);
         fprintf (WMOPS_file, "\t%.6f", ((float) tot) * frameRate);
         fprintf (WMOPS_file, "\t%.3f", ((float) tot) / grand_total * 100);

         if (nbframe[currCounter] != 0)
            fprintf (WMOPS_file, "\t%.3f", tot_wm / (float) nbframe[currCounter]);

         fprintf (WMOPS_file, "\t%.3f", ((float) glob_wc[currCounter]) * frameRate);

         /* Worst worst case printed only when not in DTX mode */
         if (dtx_mode == 0) {
            tot_wc = 0L;
            for (i = 0; i < funcid[currCounter]; i++)
               tot_wc += wc[currCounter][i];
            fprintf (WMOPS_file, "\t%.3f", ((float) tot_wc) * frameRate);
         }
         fprintf (WMOPS_file, "\n");

      }

      /* Print the file Grand Total line */
      fprintf (WMOPS_file, "%s", test_file_name);
      fprintf (WMOPS_file, "\tGrand Total");
      fprintf (WMOPS_file, "\t%d", nbframe[saved_value]);
      fprintf (WMOPS_file, "\t\t%.6f", ((float) grand_total) * frameRate);
      fprintf (WMOPS_file, "\t100.000");
      fprintf (WMOPS_file, "\n");
      fclose(WMOPS_file);

   } else
      printf( "Can not open file %s for WMOPS editing\n", WMOPS_DATA_FILENAME);


   if( (WMOPS_file=fopen(WMOPS_TOTAL_FILENAME,"a"))!=NULL) {
      printf( "opened file %s in order to print application's total WMOPS.\n", WMOPS_TOTAL_FILENAME);
      fprintf (WMOPS_file, "%s", test_file_name);
      fprintf (WMOPS_file, "\tframe=%d", nbframe[currCounter]);
      fprintf (WMOPS_file, "\tWMOPS=%.6f", ((float) grand_total) * frameRate);
      fprintf (WMOPS_file, "\n");
      fclose(WMOPS_file);

   } else
      printf( "Can not open file %s for WMOPS editing.\n", WMOPS_TOTAL_FILENAME);


   if( (WMOPS_file=fopen(CODE_PROFILE_FILENAME,"a"))!=NULL) {

      printf( "opened file %s in order to print basic operation distribution statistics.\n", CODE_PROFILE_FILENAME);

      /* Print the file header line. */
      fprintf (WMOPS_file, "Test file name\tBasic Operation Name\tframe\tWMOPS\t\t%% versus grand total\n");

      /* Print the WMOPS for each Basic Operation across all the defined */
      /* Function Groups. */
      for( i = 0; i < (int)(sizeof(op_weight) / sizeof(Word32)); i++) {
         fprintf (WMOPS_file, "%-16s", test_file_name);
         fprintf (WMOPS_file, "\t%s", BasicOperationList[i]);
         fprintf (WMOPS_file, "\t%d", nbframe[0]);

         tot = 0;
         ptr = (Word32 *) &multiCounter[0] + i;
         ptr2 = (const Word32 *) &op_weight + i;
         for( currCounter = 0; currCounter <= maxCounter; currCounter++) {
            tot += ((*ptr) * (*ptr2));
            ptr += (sizeof(op_weight) / sizeof(Word32));
         }

         fprintf (WMOPS_file, "\t%.6f", ((float) tot) * frameRate);
         fprintf (WMOPS_file, "\t%.3f", ((float) tot) / grand_total * 100);
         fprintf (WMOPS_file, "\n");
      }

      /* Print the file Grand Total line */
      fprintf (WMOPS_file, "%s", test_file_name);
      fprintf (WMOPS_file, "\tGrand Total");
      fprintf (WMOPS_file, "\t%d", nbframe[saved_value]);
      fprintf (WMOPS_file, "\t%.6f", ((float) grand_total) * frameRate);
      fprintf (WMOPS_file, "\t100.000");
      fprintf (WMOPS_file, "\n");
      fclose(WMOPS_file);

   } else
      printf( "Can not open file %s for basic operations distribution statistic editing\n", CODE_PROFILE_FILENAME);

   currCounter = saved_value;
#else
    (void)dtx_mode;
    (void)test_file_name;
#endif /* if WMOPS */
}

#if WMOPS
#define MAX_STACK       715
static int stack[MAX_STACK];
static int sptr;
static int sum_stack[MAX_STACK];
#endif

#if WMOPS
/* jdr 20120117: add FLC similar functions */
void BASOP_frame_update(void)
{
  int i, current;
#if MAX_CALLERS_SAVED_FRAMES
  int k;
#endif
  float total = 0.0f;

  {
    static int sptr_target=-2;

    if (sptr_target == -2) {
      sptr_target = sptr;
    } else {
      if (sptr_target != sptr) {
        fprintf(stderr, "BASOP_sub_start/BASOP_sub_end imbalance detected!!!\n");
        sptr_target = sptr;
      }
    }
  }

  /* Get current counter */
  current = readCounterId();

  /* Update global operation counters */
  for (i=1; i<=maxCounter; i++)
  {
    int j;

    for (j=0; j< (int)(sizeof(BASIC_OP)/sizeof(UWord32)); j++)
    {
      ((UWord32*)&glob_multiCounter)[j] += ((UWord32*)&multiCounter[i])[j];
    }
  }

#if MAX_CALLERS_SAVED_FRAMES
  /* Reset all counters */
  for (i=1; i<=maxCounter; i++)
  {
      callers_frames[0][i] = 0.0f;
  }
#endif
  /* Reset all counters */
  for (i=1; i<=maxCounter; i++)
  {
    if (current != i && funcid[i] > 0) {
      setCounter(i);

      glob_sum_curr[currCounter] += sum_curr[currCounter];

      if (glob_sum_curr[currCounter] > glob_sum_wc[currCounter]) {
        glob_sum_wc[currCounter] = glob_sum_curr[currCounter];
      }
      if (glob_sum_curr[currCounter] < glob_sum_bc[currCounter]) {
        glob_sum_bc[currCounter] = glob_sum_curr[currCounter];
      }
#if MAX_CALLERS_SAVED_FRAMES
      /* Keep a Copy before it is Reset */
      callers_frames[0][currCounter] = (float)Reset_WMOPS_counter();
      total += callers_frames[0][currCounter];
#else
      total += (float)Reset_WMOPS_counter();
#endif
    }
  }

#if MAX_CALLERS_SAVED_FRAMES
  /* Keep Callers for this Worst Case Frame */
  /* Select Slot to Use (Slot 0 is the Current) */
  k = 0;
  for (i=k+1; i<MAX_CALLERS_SAVED_FRAMES; i++)
  {
    /* Is it the Min? */
    if (callers_totals[i] < callers_totals[k])
    { /* Yes */
       k = i;
    }
  }
  /* Current Greater than the Min? */
  if (total > callers_totals[k])
  {
    k+=1;
    /* Save Info of Callers */
    for (i=1; i<=maxCounter; i++)
    {
      callers_frames[k][i] = callers_frames[0][i];
    }
    if (i < MAXCOUNTERS)
      callers_frames[k][i] = -1;
    /* Save Total */
    callers_totals[k-1] = total;
    /* Save Frame Number */
    callers_frames_nos[k-1] = nbframe[0];
  }
#endif
  if (total < glob_bc[0]) {
    glob_bc[0] = (Word32) total;
  }
  if (total > glob_wc[0]) {
    glob_wc[0] = (Word32) total;

    for (i=0; i < MAXCOUNTERS+1; i++)
    {
      glob_wf[i] = glob_wf_tmp[i];
      sum_wf[i] = sum_wf_tmp[i];
      nbcalls_wf[i] = nbcalls_wf_tmp[i];
    }
  }

  for (i=0; i < MAXCOUNTERS+1; i++)
  {
    glob_wf_tmp[i] = 0;
    sum_wf_tmp[i] = 0;
    nbcalls_wf_tmp[i] = 0;
  }

  /* Restore current counter */
  setCounter(current);

#ifdef BASOP_OVERFLOW2
  /* make sure overflow warnings are always re-enabled (if this fails, there's most probably
     a BASOP_SATURATE_WARNING_OFF without subsequent BASOP_SATURATE_WARNING_ON) */
  /* assert(overflow_warning_disable_counter == 0); */ /* Disable this santiy check to enable inverse scope strategy */
#endif

    nbframe[0]++;
}
#endif /* WMOPS */

#if WMOPS
void printStack(char *text, char* Id)
{
   int i;
   if (!Id) return;
   if (!strcmp("*", Id) || (!strcmp(Id,objectName[currCounter])))
   {
     printf ("%s %s", text, objectName[currCounter]);
     for (i=sptr-1; i>0; i--) {
       printf(" <- %s", objectName[stack[i]]);
     }
     printf("\n");
   }
}
#endif /* WMOPS */

#if WMOPS
void BASOP_push_wmops (const char *label)
{

  int new_flag, prev_counter;
  int i, j;

  /* Check if new counter label */
  new_flag = 1;
  for (i = 1; i <= maxCounter; i++)
    {
      if (strcmp(objectName[i], label) == 0)
        {
          new_flag = 0;
          break;
        }
    }

  prev_counter = readCounterId();

  /* Configure new record */
  if (new_flag)
    {
      i = (int)getCounterId(label);
      setCounter(i);
      Init_WMOPS_counter();
    }
  else
    {
      setCounter(i);
    }


  /* Push current context onto stack */
  if (currCounter >= 0)
    {
      if (sptr >= MAX_STACK)
        {
          fprintf (stderr, "\r push_wmops(): stack exceeded, try inreasing MAX_STACK\n");
          exit (-1);
        }
      stack[sptr++] = prev_counter;

      /* Reset accumulated WMOPS */
      sum_stack[sptr] = 0;

      /* update call tree */
      for (j = 0; j < MAXCOUNTERS; j++)
        {
          if (call_tree[i][j] == prev_counter)
            {
              break;
            }
          else if (call_tree[i][j] == -1)
            {
              call_tree[i][j] = prev_counter;
              break;
            }
        }
    }

  /*wmops[currCounter].start_selfcnt = ops_cnt;
  wmops[currCounter].start_cnt = ops_cnt;
  nbTimeObjectIsCalled[currCounter]++;*/

  incrementNbTimeObjectIsCalled(currCounter);

  sum_curr[currCounter] = 0;

#ifdef DEBUG_COUNTER
  printf("Entering: %s\n", readCounterIdName());
#endif
}
#endif /* WMOPS */

#if WMOPS
Word32 BASOP_pop_wmops (void)
{
  Word32 ops_cnt;

#ifdef DEBUG_COUNTER
  printf("Exiting: %s\n", readCounterIdName());
#endif

  ops_cnt = fwc();
  glob_wf_tmp[currCounter] += ops_cnt;
  nbcalls_wf_tmp[currCounter]++;

  /* Get back previous context from stack */
  if (sptr > 0)
    {
      int prevCounter;
      sum_curr[currCounter] += ops_cnt;
      sum_wf_tmp[currCounter] += ops_cnt;
      prevCounter = currCounter;
      setCounter(stack[--sptr]);
      sum_curr[currCounter] += sum_curr[prevCounter];
      sum_wf_tmp[currCounter] += sum_curr[prevCounter];
    }
  else
    {
      /* current_record = -1; */
      setCounter(0);
    }

  if (sum_curr[currCounter] > sum_wc[currCounter]) {
    sum_wc[currCounter] = sum_curr[currCounter];
  }
  if (sum_curr[currCounter] < sum_bc[currCounter]) {
    sum_bc[currCounter] = sum_curr[currCounter];
  }

  return ops_cnt;
}
#endif /* WMOPS */

#ifdef EVS_WMOPS_COUNT
Word32 BASOP_get_wops (void)
{
  return BASOP_pop_wmops();
}
#endif

#if WMOPS
static Word32 prom_cnt = 0;
#endif

void WMOPS_destroy(void)
{
#if WMOPS
  int i;

  /* release the memory allocated for wc[], bc[] */
  for (i = 0; i < MAXCOUNTERS; i++)
  {

    if(wc[i] != NULL){
        free(wc[i]);
        wc[i] = NULL;
    }
  }

  /* release the memory allocated for the objectName array */
  for (i = 0; i < MAXCOUNTERS+1; i++)
  {
    if (NULL != objectName[i])
      {
        free(objectName[i]);
        objectName[i] = NULL;
      }
  }

  maxCounter = 0;
#endif

  return;
}

void WMOPS_output_all(Word16 dtx_mode)
{
#if WMOPS
  float ops_cnt = 0.0f;
  int i;

  char *sfmts = "%-40s %8s %8s %7s %7s\n";
  char *dfmts = "%-40s %8.2f %8.3f %7.3f %7.3f\n";
  char *sfmt =  "%-40s %8s %8s %7s %7s  %7s %7s %7s\n";
  char *dfmt =  "%-40s %8.2f %8.3f %7.3f %7.3f  %7.3f %7.3f %7.3f\n";

  fprintf (stderr, "\nProgram Memory Analysis: %12.0f words\n", (float)prom_cnt);
  /*fprintf (stderr, "\nInstruction Type Analysis (for worst case frame):\n\n");*/
  fprintf (stderr, "\nInstruction Type Analysis (for worst case frame number %ld):\n\n", (long int)nbframe[0]);   /* added -- JPA */
  for (i = 0; i < (int)(sizeof(BasicOperationList)/sizeof(char*)) ; i++)
    {
      fprintf (stderr, "\t%16s:          %12d\n", BasicOperationList[i], ((UWord32*)&glob_multiCounter)[i]);
    }

  fprintf (stderr, "\n\nWeighted MOPS Analysis:\n");
  fprintf (stderr, "%74s  %23s\n", "|------  SELF  ------|"
                                  ,"|---  CUMULATIVE  ---|");
  fprintf (stderr, sfmt, "                  routine", " calls", "  min ", "  max ", "  avg ", "  min ", "  max ", "  avg ");
  fprintf (stderr, sfmt, " ------------------------", "  ------", "------", "------", "------", "------", "------", "------");
  for (i = 0; i <= maxCounter; i++)
    {
      fprintf (stderr, dfmt,
               objectName[i],
               (nbframe[i] == 0) ? 0 : (float)nbcalls[i]/(float)nbframe[0],
               (glob_bc[i] == 0) ? 0 : frameRate*(float)glob_bc[i],
               (glob_wc[i] == 0) ? 0 : frameRate*(float)glob_wc[i],
               (nbframe[i] == 0) ? 0 : (float)total_wmops[i]/(float)nbframe[i],
               frameRate*(glob_sum_bc[i]),
               frameRate*(glob_sum_wc[i]),
               (nbframe[i] == 0) ? 0 : (float)(total_sum[i] /(float)nbframe[i]));
               /* frameRate*(glob_bc[i]+wmops_children_bc[i]), */
               /* frameRate*(glob_wc[i]+wmops_children_wc[i]), */
               /* (nbframe[i] == 0) ? 0 : (float)((total_wmops[i] + total_wmops_children[i]) /(float)nbframe[i])); */

      ops_cnt += total_wmops[i];
    }

  fprintf (stderr, sfmts, " -----------------", "  ------", "------", "------", "------");
  fprintf (stderr, dfmts,
           "total",
           (double)nbframe[0],
           frameRate*glob_bc[0],
           frameRate*glob_wc[0],
           (nbframe[0] == 0) ? 0 : ops_cnt/nbframe[0]);


  (void)dtx_mode;

#ifdef BASOP_OVERFLOW2
#ifndef HIDE_UNUSED_BASOP
  fprintf(stderr, "Saturation warnings counted: %d\n", overflow_count);
#endif
#ifdef CTEST_MEASUREMENTS
#ifndef HIDE_UNUSED_BASOP
  printf("<DartMeasurement name=\"%s\" type=\"%s\">%d</DartMeasurement>\n", "overflows", "numeric/integer", overflow_count);
#endif
#endif
#endif

#if MAX_CALLERS_SAVED_FRAMES
  for (i = 1; i <= MAX_CALLERS_SAVED_FRAMES; i++)
  {
      int j, k, l, m;
      const char *frame_rank[] = { "st", "nd", "rd", "th" };
      float current;

      k = 0;
      for (j=k+1; j<MAX_CALLERS_SAVED_FRAMES; j++)
      {
          /* Is it the Max? */
          if (callers_totals[j] > callers_totals[k])
          { /* Yes */
              k = j;
          }
      }
      k+=1;

      fprintf(stderr, "\nActive Callers Report for %i%s Worst Case Frame #: %i\n",
                      i, i<=3?frame_rank[i-1]:frame_rank[3],
                      callers_frames_nos[k-1]);
      /* Print up to 'MAX_CALLERS_PRINT' Callers */
      current = 0.0f;
      for (l = 0; l < MAX_CALLERS_PRINT; l++)
      {
          /* Find Highest Complexity */
          m = 1;
          for (j = m+1; j <= maxCounter; j++)
          {
            if (callers_frames[k][j] < 0)
                break;
            if (callers_frames[k][j] > callers_frames[k][m])
                m = j;
          }
          /* Done ? */
          if (callers_frames[k][m] == 0)
              break;
          fprintf(stderr, "   %-52s %10.3f\n", objectName[m], callers_frames[k][m]*frameRate);
          /* Count it */
          current += callers_frames[k][m];
          /* Mark as Done */
          callers_frames[k][m] = 0.0f;
      }
      /* Check if All Printed */
      if (current+0.001f < callers_totals[k-1])
      {
        fprintf(stderr, " Only first %i Callers have been Printed!\n", MAX_CALLERS_PRINT);
        fprintf(stderr, "   %-52s %10.3f\n", "Total for non Printed", (callers_totals[k-1]-current)*frameRate);
      }
      fprintf(stderr, "   %-52s %10.3f\n", "Total", callers_totals[k-1]*frameRate);
      /* Mark as Done */
      callers_totals[k-1] = 0.0f;
  }
#endif
  WMOPS_destroy();
#else
  (void)dtx_mode;
#endif /* if WMOPS */
}

typedef struct result_t {
  char * name;
  float calls;
  float wf_calls;
  float bc;
  float wc;
  float avg;
  float wf;
  float bc_sum;
  float wc_sum;
  float avg_sum;
  float wf_sum;
} result_t;

#if WMOPS
int cmpfunc (const void * a, const void * b);
int cmpfunc (const void * a, const void * b)
{
  const result_t * A = a;
  const result_t * B = b;

  if( (1.0f*A->avg_sum) > (1.0f*B->avg_sum))
    return -1;
  else if( (1.0f*B->avg_sum) > (1.0f*A->avg_sum))
    return 1;
  else
    return 0;
}
#endif

void WMOPS_output_all_std(Word16 dtx_mode)
{
#if WMOPS
  float ops_cnt = 0.0f, wf_ops_cnt=0.0f;
  int i;

  char *sfmts = "%-40s %8s          %8s %7s %7s %7s\n";
  char *dfmts = "%-40s %8.2f          %8.3f %7.3f %7.3f %7.3f\n";
  char *sfmt =  "%-40s %8s %8s %8s %7s %7s %7s   %7s %7s %7s %7s\n";
  char *dfmt =  "%-40.40s %8.2f %8.0f %8.3f %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f %7.3f\n";

  result_t result[MAXCOUNTERS+1];

  fprintf (stdout, "\nProgram Memory Analysis: %12.0f words\n", (float)prom_cnt);
  /*fprintf (stdout, "\nInstruction Type Analysis (for worst case frame):\n\n");*/
  fprintf (stdout, "\nInstruction Type Analysis (for worst case frame number %ld):\n\n", (long int)nbframe[0]);   /* added -- JPA */
  for (i = 0; i < (int)(sizeof(BasicOperationList)/sizeof(char*)) ; i++)
  {
    fprintf (stdout, "\t%16s:          %12d\n", BasicOperationList[i], ((UWord32*)&glob_multiCounter)[i]);
  }

  for (i = 1; i <= maxCounter; i++)
  {
    if(glob_bc[i] > glob_wc[i])
      glob_bc[i] = -1;
    if(glob_sum_bc[i] > glob_sum_wc[i])
      glob_sum_bc[i] = -1;
    if(nbcalls_wf[i] == 0)
      sum_wf[i]=-1;

    result[i-1].name = (objectName[i] == NULL) ? "\0" : objectName[i];
    result[i-1].calls =(nbframe[i] == 0) ? 0 : (float)nbcalls[i]/(float)nbframe[0];
    result[i-1].wf_calls =(nbframe[i] == 0) ? 0 : (float)nbcalls_wf[i];
    result[i-1].bc =(glob_bc[i] == 0) ? 0 : (float)frameRate*(float)glob_bc[i];
    result[i-1].wc =(glob_wc[i] == 0) ? 0 : (float)frameRate*(float)glob_wc[i];
    result[i-1].avg = (nbframe[i] == 0) ? 0 : (float)total_wmops[i]/(float)nbframe[i];
    result[i-1].wf =(glob_wf[i] == 0) ? 0 : (float)frameRate*(float)glob_wf[i];
    result[i-1].bc_sum =(float)frameRate*(1.0f*glob_sum_bc[i]);
    result[i-1].wc_sum =(float)frameRate*(1.0f*glob_sum_wc[i]);
    result[i-1].avg_sum =(nbframe[i] == 0) ? 0 : (float)((1.0f* total_sum[i]) /(float)nbframe[i]);
    result[i-1].wf_sum =(float)frameRate*(1.0f*sum_wf[i]);
    ops_cnt += total_wmops[i];
    wf_ops_cnt += (float)(frameRate*glob_wf[i]);
  }
  fprintf (stdout, "\n\nWeighted MOPS Analysis:\n");
  fprintf (stdout, "%59s %22s  %22s\n","|---- CALLS ----|" , "|-----------  SELF  -----------|"
                                    ,"|--------  CUMULATIVE  --------|");
  fprintf (stdout, sfmt, "                  routine", "    avg ", "  wf  ", "  min ", "  max ", "  avg ", "  wf  ", "  min ", "  max ", "  avg ", "  wf  ");
  fprintf (stdout, sfmt, " ------------------------", "  ------", "------", "------", "------", "------", "------", "------", "------", "------", "------");

  qsort(result, maxCounter, sizeof(result_t), cmpfunc);

  for (i = 0; i < maxCounter; i++)
    {
      fprintf (stdout, dfmt,
          result[i].name,
          result[i].calls,
          result[i].wf_calls,
          result[i].bc,
          result[i].wc,
          result[i].avg,
          result[i].wf,
          result[i].bc_sum,
          result[i].wc_sum,
          result[i].avg_sum,
          result[i].wf_sum
               );
    }

  fprintf (stdout, sfmts, " -----------------", "  ------", "------", "------", "------", "------");
  fprintf (stdout, dfmts,
           "total",
           (double)nbframe[0],
           frameRate*glob_bc[0],
           frameRate*glob_wc[0],
           (nbframe[0] == 0) ? 0 : ops_cnt/nbframe[0],
           (nbframe[0] == 0) ? 0 : wf_ops_cnt);


  (void)dtx_mode;

#ifdef BASOP_OVERFLOW2
#ifndef HIDE_UNUSED_BASOP
  fprintf(stdout, "Saturation warnings counted: %d\n", overflow_count);
#endif
#ifdef CTEST_MEASUREMENTS
  printf("<DartMeasurement name=\"%s\" type=\"%s\">%d</DartMeasurement>\n", "overflows", "numeric/integer", overflow_count);
#endif
#endif

  WMOPS_destroy();
#else
  (void)dtx_mode;
#endif /* if WMOPS */
}

#if WMOPS
void Reset_all_WMOPS_counter (void)
{
  int i;
  int currCounterSave;

  currCounterSave = currCounter;

  for (i=2; i <= maxCounter; i++)
  {
    setCounter(i);
    Init_WMOPS_counter();
    objectName[i] = 0;
  }

  currCounter = currCounterSave;
  maxCounter = 1;
}
#endif /* WMOPS */

/* Returns the total min/max/avg WMOPS values like printed in BASOP_end(). */
#if WMOPS
void BASOP_get_total_wmops(double *min, double *max, double *avg)
{
  if(min != NULL)
    *min = frameRate * glob_bc[0];
  if(max != NULL)
    *max = frameRate * glob_wc[0];
  if(avg != NULL) {
    int i;
    double ops_cnt = 0;
    for(i = 1; i <= maxCounter; i++)
      ops_cnt += total_wmops[i];
    *avg = (nbframe[0] == 0) ? 0 : ops_cnt / nbframe[0];
  }
}
#endif /* WMOPS */

/* end of file */

#if WMOPS
/* helper functions for matlab version */
void BASOP_sub_sub_start(char *msg)
{
     BASOP_sub_start(msg);
}
#endif

#if WMOPS
void BASOP_sub_sub_end(void)
{
     BASOP_sub_end();
}
#endif
