#ifndef BP_FormB_HEADER
#define BP_FormB_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <errno.h>
#include <sstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <set>
#include <iomanip> // for precision

//#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include <lemon/network_simplex.h>

#include "global_variables.h"
#include "global_functions.h"
#include "control_panel.h"


using namespace std;


typedef struct
{
    int id;
    int level;
    int br_var;
    int side;
//    int val;
//    int cur_LB;
//    int cur_UB
} BB_node_info;

typedef struct BB_node
{
    struct BB_node *father;

    BB_node_info info;

    double dual_bound;
    int primal_bound;
    int status;

    struct BB_node *left;
    struct BB_node *right;
} BB_node, *PTR_BB_node;


enum Task_status { EVALUATED, CALCULATING, WAITING };
enum Processor_status { PROC_FREE, PROC_COMPUTING, PROC_RETIRED };
enum LP_status { INFEASIBLE, INTEGER, FRACTIONAL, PRUNED, LP_ERROR };
enum node_sides { LEFT, RIGHT, ROOT };


/*************************************************************************************************/
void MPI_MyBBInit();
/*************************************************************************************************/

/*************************************************************************************************/
void MPI_MyBBfree();
/*************************************************************************************************/

/*************************************************************************************************/
bool myOperator (PTR_BB_node i_node, PTR_BB_node j_node);
/*************************************************************************************************/

/*************************************************************************************************/
bool myPairsOperator (pair<int,double> i_fracVar, pair<int,double> j_fracVar);
/*************************************************************************************************/

/*************************************************************************************************/
int my_floor(double a);
/*************************************************************************************************/

/*************************************************************************************************/
int my_ceil(double a);
/*************************************************************************************************/

/*************************************************************************************************/
void display_BB_node(int i_node);
/*************************************************************************************************/

/*************************************************************************************************/
bool branch_kHNDP_select_branching(IloNumArray const &solX, IloArray<IloNumArray> const &solF, int *item_i);
/*************************************************************************************************/

/*************************************************************************************************/
void branch_kHNDP_right(int item_i);
/*************************************************************************************************/

/*************************************************************************************************/
void branch_kHNDP_left(int item_i);
/*************************************************************************************************/

/*************************************************************************************************/
void branch_kHNDP_free(int item_i);
/*************************************************************************************************/

/*************************************************************************************************/
void add_children(int i_node, int i_var);
/*************************************************************************************************/

/*************************************************************************************************/
void master_flow_build();
/*************************************************************************************************/

/*************************************************************************************************/
void master_kHNDP_initialize();
/*************************************************************************************************/

/*************************************************************************************************/
double master_solve_lp();
/*************************************************************************************************/

/*************************************************************************************************/
bool master_check_integrality(IloNumArray const &solX, IloArray<IloNumArray> const &solF);
/*************************************************************************************************/

/*************************************************************************************************/
void master_display_status();
/*************************************************************************************************/

/*************************************************************************************************/
void master_kHNDP_free();
/*************************************************************************************************/

/*************************************************************************************************/
void master_kHNDP_Output();
/*************************************************************************************************/

/*************************************************************************************************/
int master_kHNDP_heur_solve(IloNumArray const &solX);
/*************************************************************************************************/

/*************************************************************************************************/
int branch_and_bound(int level);
/*************************************************************************************************/

/*************************************************************************************************/
int evaluate_node(int &primal_bound, int &i_item);
/*************************************************************************************************/

/*************************************************************************************************/
void affect_task();
/*************************************************************************************************/

/*************************************************************************************************/
void recieve_results();
/*************************************************************************************************/

/*************************************************************************************************/
void return_results(int stat, double current_LB, int primal_bound, int i_item);
/*************************************************************************************************/

/*************************************************************************************************/
bool update_global_LB(PTR_BB_node _node);
/*************************************************************************************************/

/*************************************************************************************************/
void recieve_task();
/*************************************************************************************************/

/*************************************************************************************************/
void get_ready_for_solving();
/*************************************************************************************************/

/*************************************************************************************************/
void restore_master_LP();
/*************************************************************************************************/

/*************************************************************************************************/
void master();
/*************************************************************************************************/

/*************************************************************************************************/
void slave();
/*************************************************************************************************/

#endif
