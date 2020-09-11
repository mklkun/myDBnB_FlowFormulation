#ifndef FUNCION_HEADER
#define FUNCION_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip> // for precision


#include "global_variables.h"


using namespace std;


/*************************************************************************************************/
void SORT_NON_INCR_INT(int *item,int *score,int n);
/*************************************************************************************************/

/*************************************************************************************************/
void SORT_NON_INCR_DOUBLE(int *item,double *score,int n);
/*************************************************************************************************/

/*************************************************************************************************/
template <typename VALMAP, typename CAPMAP>
bool feasableAnyL(unsigned const pos, VALMAP const &VARVM, CAPMAP &VARCM);
/*************************************************************************************************/

/*************************************************************************************************/
int kHNDP_first_fit();
/*************************************************************************************************/

/*************************************************************************************************/
bool extract_number(string &str, int &number);
/*************************************************************************************************/

/*************************************************************************************************/
void kHNDP_instance_read(const char *probname);
/*************************************************************************************************/

/*************************************************************************************************/
void transform_graph(int const s_i, int const t_i, unsigned const pos);
/*************************************************************************************************/

/*************************************************************************************************/
void kHNDP_instance_free();
/*************************************************************************************************/

#endif