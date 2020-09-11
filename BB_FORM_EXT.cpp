

#include "BB_FORM_EXT.h"


///////////////////////////////////CPLEX///////////////////////////////////////
IloEnv env_MAS;
IloModel mod_MAS(env_MAS);
IloCplex cplex_MAS(mod_MAS);

IloNumVarArray *X;
IloArray<IloNumVarArray> *F;
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////BRANCH & BOUND///////////////////////////////////
int nb_evaluated,nb_calculating,nb_waiting,nb_pruned; // nodes

int nb_f_frac;

int lpstat,nodecount,numrows,numcols;
double bestobjval, objval;

bool solved;

double first_LB_time;
int best_UB;
double best_LB;
int nb_BB_nodes;
bool STOP_TIME_LIMIT = false;

deque < PTR_BB_node > waiting_line;
vector < pair <int,double> > frac_vars;
vector < PTR_BB_node > BB_nodes;
///////////////////////////////////////////////////////////////////////////////


/////////////////////////////////// MPI ///////////////////////////////////////
int *slaves_status;
int *task_affected;
int nb_proc_calculating,nb_proc_waiting;

vector < pair <int,int> > branching_chain;

MPI_Status status;

#define WORKTAG         1
#define COMPLETIONTAG   2
#define DIETAG          3
///////////////////////////////////////////////////////////////////////////////


#define EPSILON_BB  0.000001
#define INFINITE    999999999


///////////////////////////////////OPTIONS/////////////////////////////////////
#define use_heuristic_master 1 //compute an heuristic solution using the master at each node using the value as timelimit

#define BEST_FIRST
//#define FIFO

#define branch_min_reduced_cost
//#define branch_first_frac_var
//#define branch_closer_0_5

//#define cplex_LP_write

#define cplex_alg_auto
//#define cplex_alg_Primal_Simplex
//#define cplex_alg_Dual_Simplex
//#define cplex_alg_Network_Simplex
//#define cplex_alg_Barrier
//#define cplex_alg_Sifting
//#define cplex_alg_Concurrent

//#define desable_cplex_presolve_and_heur
//#define disable_cplex_cuts
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////PRINTINGS////////////////////////////////////
//#define outputFile_by_rank_MPI
//#define print_column
//#define print_pi
//#define print_MASTER
//#define writer_PRICER
//#define verbose_BP
//#define verbose_PRUNING
//#define verbose_HEUR
//#define write_lp_FORM_B
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////MPI OUTPUT////////////////////////////////////
#ifdef outputFile_by_rank_MPI
    string proc_OutFile_name;
    ofstream procOutFile;
#endif // outputFile_by_rank_MPI
///////////////////////////////////////////////////////////////////////////////


/*************************************************************************************************/
void MPI_MyBBInit()
/*************************************************************************************************/
{
    slaves_status = new int [world_size];

    for (unsigned i(1); i < world_size; i++)
        slaves_status[i] = PROC_FREE;

    task_affected = new int [world_size];

    nb_proc_calculating = 1;
    nb_proc_waiting = world_size - 1;

#ifdef outputFile_by_rank_MPI
    proc_OutFile_name = "outMPI/P"+ to_string(world_rank) + "_" + name + "_" + "_L" + to_string(L) + "_k" + to_string(k) + ".txt";

    procOutFile.open(proc_OutFile_name, ios::out);
#endif // outputFile_by_rank_MPI
}


/*************************************************************************************************/
void MPI_MyBBfree()
/*************************************************************************************************/
{
    delete [] slaves_status;
    delete [] task_affected;

#ifdef outputFile_by_rank_MPI
    procOutFile.close();
#endif // outputFile_by_rank_MPI
}


/*************************************************************************************************/
bool myOperator (PTR_BB_node i_node, PTR_BB_node j_node)
/*************************************************************************************************/
{
    return (i_node->dual_bound < j_node->dual_bound);
}


/*************************************************************************************************/
bool myPairsOperator (pair<int,double> i_fracVar, pair<int,double> j_fracVar)
/*************************************************************************************************/
{
    return (i_fracVar.second < j_fracVar.second);
}


/*************************************************************************************************/
int my_floor(double a)
/*************************************************************************************************/
{
	int q;
	q = (int) (a + EPSILON_BB);
	return (q);
}


/*************************************************************************************************/
int my_ceil(double a)
/*************************************************************************************************/
{
	int q;
	double b;
	b = a + 1.0 - EPSILON_BB;
	q = (int) (b);
	return (q);
}


/*************************************************************************************************/
void display_BB_node(int i_node)
/*************************************************************************************************/
{
    cout << "Node " << i_node << " from between " << nb_BB_nodes << " :" << endl;
    if (i_node < nb_BB_nodes)
    {
        if (BB_nodes[i_node]->father != NULL)
            cout << "   father id " << BB_nodes[i_node]->father->info.id << ";" << endl;
        cout << "   id " << BB_nodes[i_node]->info.id << ";" << endl;
        cout << "   level " << BB_nodes[i_node]->info.level << ";" << endl;
        cout << "   branching node " << BB_nodes[i_node]->info.br_var << ";" << endl;
        switch (BB_nodes[i_node]->info.side)
        {
        case ROOT:
            cout << "   ROOT node;" << endl;
            break;
        case LEFT:
            cout << "   side LEFT;" << endl;
            break;
        case RIGHT:
            cout << "   side RIGHT;" << endl;
        }

        cout << "   status " << BB_nodes[i_node]->status << ";" << endl;
        cout << "   dual bound " << BB_nodes[i_node]->dual_bound << ";" << endl;
        cout << "   primal bound " << BB_nodes[i_node]->primal_bound << ";" << endl;
        cout << "   left pointer " << BB_nodes[i_node]->left << ";" << endl;
        cout << "   right pointer " << BB_nodes[i_node]->right << ";" << endl;
    }
}


/*************************************************************************************************/
bool branch_kHNDP_select_branching(IloNumArray const &solX, IloArray<IloNumArray> const &solF, int *item_i)
/*************************************************************************************************/
{
    int i_var(0);

#ifdef branch_min_reduced_cost
    frac_vars.clear();

    for (unsigned i(0); i < G.m; i++,i_var++)
    {
        //if ((fabs(solX[i] - floor(solX[i])) > EPSILON_BB) && (fabs(solX[i] - floor(solX[i]) - 1.0) > EPSILON_BB))
		if (solX[i] < 1 - EPSILON_BB && solX[i] > EPSILON_BB)
		{
			frac_vars.push_back( make_pair(i_var , cplex_MAS.getReducedCost((*X)[i])) );
		}
    }

    sort(frac_vars.begin(), frac_vars.end(), myPairsOperator);

    if (frac_vars.size() >= 1)
    {
//        cout << "Selecting variable " << frac_vars[0].first << " that has reduced cost " << frac_vars[0].second << endl;
        *item_i = frac_vars[0].first;
        return true;
    }

    frac_vars.clear();

	for(unsigned i(0); i < demandsNumber; i++)
    {
        for(unsigned j(0); j < DGS[i]->m; j++,i_var++)
        {
            //if ((fabs(solF[i][j] - floor(solF[i][j])) > EPSILON_BB) && (fabs(solF[i][j] - floor(solF[i][j]) - 1.0) > EPSILON_BB))
            if(solF[i][j] < 1 - EPSILON_BB && solF[i][j] > EPSILON_BB)
            {
                frac_vars.push_back( make_pair(i_var , cplex_MAS.getReducedCost((*F)[i][j])) );
            }
        }
    }

    sort(frac_vars.begin(), frac_vars.end(), myPairsOperator);

    if (frac_vars.size() >= 1)
    {
//        cout << "Selecting variable " << frac_vars[0].first << " that has reduced cost " << frac_vars[0].second << endl;
        *item_i = frac_vars[0].first;
        return true;
    }
#endif // branch_min_reduced_cost

#ifdef branch_closer_0_5
    frac_vars.clear();

    for (unsigned i(0); i < G.m; i++,i_var++)
    {
        //if ((fabs(solX[i] - floor(solX[i])) > EPSILON_BB) && (fabs(solX[i] - floor(solX[i]) - 1.0) > EPSILON_BB))
		if (solX[i] < 1 - EPSILON_BB && solX[i] > EPSILON_BB)
		{
			frac_vars.push_back( make_pair(i_var , fabs(solX[i] - 0.5)) );
		}
    }

    sort(frac_vars.begin(), frac_vars.end(), myPairsOperator);

    if (frac_vars.size() >= 1)
    {
//        cout << "Selecting variable " << frac_vars[0].first << " that has value " << frac_vars[0].second << endl;
        *item_i = frac_vars[0].first;
        return true;
    }

    frac_vars.clear();

	for(unsigned i(0); i < demandsNumber; i++)
    {
        for(unsigned j(0); j < DGS[i]->m; j++,i_var++)
        {
            //if ((fabs(solF[i][j] - floor(solF[i][j])) > EPSILON_BB) && (fabs(solF[i][j] - floor(solF[i][j]) - 1.0) > EPSILON_BB))
            if(solF[i][j] < 1 - EPSILON_BB && solF[i][j] > EPSILON_BB)
            {
                frac_vars.push_back( make_pair(i_var , fabs(solF[i][j] - 0.5)) );
            }
        }
    }

    sort(frac_vars.begin(), frac_vars.end(), myPairsOperator);

    if (frac_vars.size() >= 1)
    {
//        cout << "Selecting variable " << frac_vars[0].first << " that has value " << frac_vars[0].second << endl;
        *item_i = frac_vars[0].first;
        return true;
    }

#endif // branch_closer_0_5

#ifdef branch_first_frac_var
	for (unsigned ii = 0; ii < G.m; ii++,i_var++)
	{
	    if(solX[ii] < 1 - EPSILON_BB && solX[ii] > EPSILON_BB)
        {
            *item_i = i_var;

//            cout << "Selecting X " << i_var << " ! ii = " << ii << " value = " << solX[ii] << endl;

            return true;
        }
    }

    for (unsigned ii = 0; ii < demandsNumber; ii++)
	{
	    for (unsigned jj(0); jj < DGS[ii]->m; jj++,i_var++)
            if(solF[ii][jj] < 1 - EPSILON_BB && solF[ii][jj] > EPSILON_BB)
            {
                *item_i = i_var;

//                cout << "Selecting F " << i_var << " ! ii = " << ii << "  jj = " << jj << " value = " << solX[ii] << endl;

                return true;
            }

    }
#endif // branch_first_frac_var

	return false;
}


/*************************************************************************************************/
void branch_kHNDP_right(int item_i)
/*************************************************************************************************/
{
	////////////////////////////////////////////////////////////////////////////////////////

	double d = 1.0;

    if (item_i < G.m)
        (*X)[item_i].setLB(d);
    else
    {
        int _item_i = item_i - G.m;
        unsigned i;
        for(i=0; i < demandsNumber && _item_i >= DGS[i]->m; _item_i-=DGS[i]->m, i++);

//        cout << "branching right on F[" << i << "][" << _item_i << "]!" << endl;

        (*F)[i][_item_i].setLB(d);
    }


//    if (status != 0)
//    {
//        cout  << "error in CPXchgobj \n";
//        exit(-1);
//    }

	////////////////////////////////////////////////////////////////////////////////////////
}


/*************************************************************************************************/
void branch_kHNDP_left(int item_i)
/*************************************************************************************************/
{
	////////////////////////////////////////////////////////////////////////////////////////

	double d = 0.0;

    if (item_i < G.m)
        (*X)[item_i].setUB(d);
    else
    {
        int _item_i = item_i - G.m;
        unsigned i;
        for(i=0; i < demandsNumber && _item_i >= DGS[i]->m; _item_i-=DGS[i]->m, i++);

//        cout << "branching left on F[" << i << "][" << _item_i << "]!" << endl;

        (*F)[i][_item_i].setUB(d);
    }

//    if (status != 0)
//    {
//        cout  << "error in CPXchgobj \n";
//        exit(-1);
//    }

	////////////////////////////////////////////////////////////////////////////////////////
}


/*************************************************************************************************/
void branch_kHNDP_free(int item_i)
/*************************************************************************************************/
{
	////////////////////////////////////////////////////////////////////////////////////////

	double ld, ud;

	ld = 0.0;
    ud = 1.0;

    if (item_i < G.m)
        (*X)[item_i].setBounds(ld,ud);
    else
    {
        int _item_i = item_i - G.m;
        unsigned i;
        for(i=0; i < demandsNumber && _item_i >= DGS[i]->m; i++, _item_i-=DGS[i]->m);

        (*F)[i][_item_i].setBounds(ld,ud);
    }

//    if (status != 0)
//    {
//        cout  << "error in CPXchgobj \n";
//        exit(-1);
//    }

	////////////////////////////////////////////////////////////////////////////////////////
}


/*************************************************************************************************/
void add_children(int i_node, int i_var)
/*************************************************************************************************/
{
//#ifdef outputFile_by_rank_MPI
//    ofstream procOutFile(proc_OutFile_name, ios::out | ios::trunc);
//#endif // outputFile_by_rank_MPI

    PTR_BB_node _node;

    // creating left son node
    _node = new BB_node;

    _node->father = BB_nodes[i_node];
    _node->info.id = nb_BB_nodes;
    _node->info.level = BB_nodes[i_node]->info.level + 1;
    _node->info.side = LEFT;
    _node->info.br_var = i_var;
    _node->left = NULL;
    _node->right = NULL;
    _node->dual_bound = 0.0;
    _node->primal_bound = 999999999;
    _node->status = WAITING;

    BB_nodes.push_back(_node);
    nb_BB_nodes++;
    nb_waiting++;

    waiting_line.push_back(_node);

    // adjusting link with father
    BB_nodes[i_node]->left = _node;

    //display_BB_node(nb_BB_nodes-1);

    // creating right son node
    _node = new BB_node;

    _node->father = BB_nodes[i_node];
    _node->info.id = nb_BB_nodes;
    _node->info.level = BB_nodes[i_node]->info.level + 1;
    _node->info.side = RIGHT;
    _node->info.br_var = i_var;
    _node->left = NULL;
    _node->right = NULL;
    _node->dual_bound = 0.0;
    _node->primal_bound = 999999999;
    _node->status = WAITING;

    BB_nodes.push_back(_node);
    nb_BB_nodes++;
    nb_waiting++;

    waiting_line.push_back(_node);

    // adjusting link with father
    BB_nodes[i_node]->right = _node;

    //display_BB_node(nb_BB_nodes-1);
}


/*************************************************************************************************/
void master_flow_build()
/*************************************************************************************************/
{
    char var[50];

    /******************************************* Creating & Setting variables ***********************************************/

    X = new IloNumVarArray (env_MAS, G.m, 0, 1, ILOFLOAT);

    F = new IloArray<IloNumVarArray> (env_MAS, demandsNumber);

    for(unsigned i = 0; i < demandsNumber ; i++)
        (*F)[i] = IloNumVarArray(env_MAS, DGS[i]->m, 0, 1, ILOFLOAT);

    for (unsigned i(0); i < G.m ; i++)
    {
        sprintf(var,"x_(%d_%d)", G.idFromNode(G.u(G.edgeFromId(i))), G.idFromNode(G.v(G.edgeFromId(i))));   // ou i directement     g_str.EIM[G.edgeFromId(i)]
        (*X)[i].setName(var);
    }

    for (unsigned i = 0; i < demands.size(); i++)
    {
        for (unsigned j = 0; j < DGS[i]->m; j++)
        {
            /// Faut il corriger en utilisant les myIdFromDiNode???
            sprintf(var, "f_d%d_(%d_%d)", i , DGS[i]->idFromDiNode(DGS[i]->source(DGS[i]->arcFromId(j))), DGS[i]->idFromDiNode(DGS[i]->target(DGS[i]->arcFromId(j))));  //j, dg_strs[i]->AIM[DGS[i]->arcFromId(j)]);
            (*F)[i][j].setName(var);
        }
    }

    /************************************************************************************************************/

    /******************************************* Objective function ***********************************************/

    IloExpr obj(env_MAS);

    for (unsigned i(0); i < G.m; i++)
        obj += G.weightFromEdge(G.edgeFromId(i)) * (*X)[i];

    mod_MAS.add(IloMinimize(env_MAS, obj));

    /************************************************************************************************************/

    /**************************************** Creating constraints ********************************************/

    int param_form(2);

    // creation contraintes de flot
    for(unsigned i(0); i < demands.size(); i++)
    {
        for(unsigned j(0); j < DGS[i]->n; j++)
        {
            IloExpr expr(env_MAS);

            lemonDiNode nod = DGS[i]->nodeFromId(j);
            for (ListDigraph::OutArcIt a(DGS[i]->m_dg, nod); a != INVALID; ++a)
                expr += (*F)[i][DGS[i]->m_dg.id(a)];

            for (ListDigraph::InArcIt a(DGS[i]->m_dg, nod); a != INVALID; ++a)
                expr -= (*F)[i][DGS[i]->m_dg.id(a)];

            IloConstraint cons;

            if (nod == DGS[i]->m_s)
                cons = (expr == k);                                                                             /// demander exactement k flots ou >= k ?!?!
            else if (nod == DGS[i]->m_t)
                cons = (expr == -(int)k);
            else
                cons = (expr == 0);

            sprintf(var, "cf_d%d_n%d",i,j);
            cons.setName(var);
            mod_MAS.add(cons);
        }
    }

        //int param_form(3);

        // Contraintes de couplage
        switch(param_form)
        {
        // Forme désagrégée
        case 1:
            for (unsigned i(0); i < demands.size(); i++)
            {
                for (unsigned j = 0; j < DGS[i]->m; j++)
                {
                    int indice = DGS[i]->idFromArc(DGS[i]->arcFromId(j));
                    if(indice < 276447231)
                    {
                        IloConstraint cons;
                        cons = ((*F)[i][j] <= (*X)[indice]);
                        sprintf(var, "cc_d%d_a%d",i ,j);
                        cons.setName(var);
                        mod_MAS.add(cons);
                    }
                }
            }
            break;
        // Forme semi agrégée
        case 2:
            for (unsigned i(0); i < demands.size(); i++)
            {
                IloExprArray exprs(env_MAS,G.m);
                for (unsigned j(0); j < G.m; j++)
                    exprs[j] = IloExpr(env_MAS);
                for (unsigned j = 0; j < DGS[i]->m; j++)
                {
                    int indice = DGS[i]->idFromArc(DGS[i]->arcFromId(j));
                    if(indice < 276447231)
                        exprs[indice] += (*F)[i][j];
                }


                for (unsigned j(0); j < G.m; j++)
                {
                    IloConstraint cons;
                    IloExpr Rs;
                    Rs = (*X)[j];
                    cons = (exprs[j] <= Rs);
                    sprintf(var, "cc_d%d_e%d", i, j);
                    cons.setName(var);

                    mod_MAS.add(cons);
                }
            }
            break;
        // Forme agrégée
        case 3:
            IloExprArray exprs(env_MAS,G.m);
            for (unsigned i(0); i < G.m; i++)
                exprs[i] = IloExpr(env_MAS);

            for (unsigned i(0); i < demands.size(); i++)
            {
                for (unsigned j = 0; j < DGS[i]->m; j++)
                {
                    int indice = DGS[i]->idFromArc(DGS[i]->arcFromId(j));
                    if(indice < 276447231)
                        exprs[indice] += (*F)[i][j];
                }

            }

            for (unsigned i(0); i < G.m; i++)
            {
                IloConstraint cons;
                IloExpr Rs;
                Rs = (IloInt)demands.size() * (*X)[i];
                cons = (exprs[i] <= Rs);
                sprintf(var, "cc_e%d",i);
                cons.setName(var);

                mod_MAS.add(cons);
            }
            break;
        }

    /************************************************************************************************************/

#ifdef cplex_LP_write
    string OutputName = "mytest.lp";
    cplex_MAS.exportModel(OutputName.c_str());
#endif // cplex_LP_write


TERMINATE:
    ;  // deleting variables?
}


/*************************************************************************************************/
double master_solve_lp()
/*************************************************************************************************/
{
	//cout << "SOLVE MASTER\n";

	bool optimality=false;
//	double *pi=new double[itemTypeNumber];
//	double *column=new double[itemNumber];
	double current_LP;
//	double lagrangian_bound;
//	double val_KP;

	IloBool current_status;

	////////////////////////////////////////////////////////
//	int var_feasibility=var_feasibility_insert();
	////////////////////////////////////////////////////////

	do{
		current_status = cplex_MAS.solve();

		//cout << "current status = " << current_status << endl;

		if(current_status == IloTrue)
        {
			current_LP = cplex_MAS.getObjValue();
			optimality = true;
		}
		else
            break;

#ifdef	print_MASTER
		cout << "current_LP\t" << current_LP << "\t";
#endif

//		lagrangian_bound=current_LP;


//		status=CPXgetpi(env_MASTER,lp_MASTER,pi,0,itemTypeNumber-1);
//		if(status!=0) {
//			printf("error in CPXgetpi\n");
//			exit(-1);
//		}

#ifdef print_pi
		for(int j=0;j<itemTypeNumber;j++)
		{
			printf("%.3f\t",pi[j]);
		}
		cout << endl;
#endif

//		//solve the KP
//		val_KP=pricer_solve_cplex(itemTypeNumber,pi,column);
//
//#ifdef print_column
//		cout << "val_KP\t" << val_KP << endl;
//		for(int j=0;j<itemTypeNumber;j++)
//		{
//			printf("%f\n",column[j]);
//		}
//		cout << endl;
//		cin.get();
//#endif
//
//		//positive reduced profit
//		if(val_KP-1>EPSILON_BB){
//
//			optimality=false;
//
//			master_column_add(column);
//			master_heur_column_add(column);
//
//			BP_cols++;
//
//#ifdef print_LP_round
//			status=CPXwriteprob(env_MASTER,lp_MASTER,"current_master.lp",NULL);
//			if(status!=0) {
//				printf("error in CPXwriteprob\n");
//				exit(-1);
//			}
//#endif
//
//		}

//		lagrangian_bound=(lagrangian_bound/val_KP);    /// val_KP is 0

		time_finish     =   clock();
		current_time    =   (double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;

#ifdef	print_MASTER
		cout <<  "lagrangian_bound\t" << lagrangian_bound << "\t time \t" <<  current_time <<  endl;
#endif

	}while( optimality == false );




//	delete[] column;
//	delete[] pi;

	//if status == -1 the current LP is infeasible
//	int current_status=1;
//	current_status=var_feasibility_remove(var_feasibility);

#ifdef	print_MASTER
	master_display_status();
#endif

	if( current_status == IloFalse )
    {
		cout << "Current INFEASIBLE\n";
		//cin.get();
		return INFINITE;
	}

    return current_LP;
}


/*************************************************************************************************/
void master_display_status()   ///à redéfinir
/*************************************************************************************************/
{
    //Declaring current solution
    IloNumArray solX(env_MAS,G.m);
    IloArray<IloNumArray> solF(env_MAS,demandsNumber);

    try
    {
        // getting current X values
        cplex_MAS.getValues(solX,*X);
        // getting the current F values
        for(unsigned i(0); i < demandsNumber; i++)
        {
            solF[i] = IloNumArray(env_MAS,DGS[i]->m);
            cplex_MAS.getValues(solF[i],(*F)[i]);
        }
    }
    catch (const IloException& e)
    {
        cerr << "Exception caught: " << e << endl;
        printf("error in getting values master_display_status\n");
    }
    catch (...)
    {
        cerr << "Unknown exception caught!" << endl;
    }

	cout << "\nCurrent MASTER STATUS\n";
    for(int j=0; j<G.m; j++)
    {

        if(solX[j]<EPSILON_BB)
            continue;

        printf("%d\t%.3f\t",j,solX[j]);
    }

    double _ub;
    double _lb;

    _ub = cplex_MAS.getObjValue();
    _lb = cplex_MAS.getBestObjValue();

    cout << "\tbounds\t"<< _lb << "\t" << _ub << endl;


    solX.end();
    for(unsigned i(0); i < demandsNumber; i++)
        solF[i].end();
    solF.end();

	cin.get();
}


/*************************************************************************************************/
bool master_check_integrality(IloNumArray const &solX, IloArray<IloNumArray> const &solF)
/*************************************************************************************************/
{
	bool integer=true;

    for (unsigned i(0); i < G.m; i++)
    {
//        if ((fabs(solX[i] - floor(solX[i])) > EPSILON_BB)
//				&& (fabs(solX[i] - floor(solX[i]) - 1.0)
//                    > EPSILON_BB))
        if(solX[i] < 1 - EPSILON_BB && solX[i] > EPSILON_BB)
		{
#ifdef outputFile_by_rank_MPI
		    procOutFile << "Variable X[" << i << "] = " << solX[i] << " is not integer!" << endl;
#endif // outputFile_by_rank_MPI
			integer = false;
			break;
		}
    }

	for(unsigned i(0); i < demandsNumber && integer; i++)
    {
        for(unsigned j(0); j < DGS[i]->m; j++)
        {
//            if ((fabs(solF[i][j] - floor(solF[i][j])) > EPSILON_BB)
//                && (fabs(solF[i][j] - floor(solF[i][j]) - 1.0)
//                    > EPSILON_BB))
            if(solF[i][j] < 1 - EPSILON_BB && solF[i][j] > EPSILON_BB)
            {
#ifdef outputFile_by_rank_MPI
                procOutFile << "Variable F[" << i << "][" << j << "] = " << solF[i][j] << " is not integer!" << endl;
#endif // outputFile_by_rank_MPI
                integer = false;
                nb_f_frac++;
                break;
            }
        }
    }

	return integer;
}


/*************************************************************************************************/
void master_kHNDP_initialize()
/*************************************************************************************************/
{
//	status_item_couple=new int*[itemTypeNumber];
//	for(int i=0;i<itemTypeNumber;i++){
//		status_item_couple[i]=new int[itemTypeNumber];
//	}
//
//	for(int i=0;i<itemTypeNumber;i++){
//		for(int j=0;j<itemTypeNumber;j++){
//			status_item_couple[i][j]=0;
//		}
//	}

	master_flow_build();
//	master_heur_build();

//	//initialize KP01 pricer
//	double *itemWeight_double=new double[itemTypeNumber];
//	for(int i=0;i<itemTypeNumber;i++){
//		itemWeight_double[i]=itemWeight[i];
//	}
//
//	pricer_load_cplex(itemTypeNumber,capacity,itemWeight_double);
//
//	delete[]itemWeight_double;
//
//	//initialize the master
//	double *column=new double[itemTypeNumber];
//
//	//bins of one item
//	for(int i=0;i< itemTypeNumber;i++){
//		for(int j=0;j< itemTypeNumber;j++){
//			if(j!=i){
//				column[j]=0;
//			}
//			else{
//				column[j]=1;
//			}
//		}
//
//		master_column_add(column);
//		master_heur_column_add(column);
//	}
//
//	//first fit solution
//	for(int i=0;i< primal_bound;i++){
//		for(int j=0;j< itemTypeNumber;j++){
//			column[j]=0;
//		}
//		for(int j=0;j< itemTypeNumber;j++){
//			if(vector_item_bin[j]==i){
//				column[j]=1;
//			}
//		}
//		master_column_add(column);
//		master_heur_column_add(column);
//	}
//	delete[] column;

#ifdef cplex_alg_auto //0
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::AutoAlg);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::AutoAlg);
#endif

#ifdef cplex_alg_Primal_Simplex //1
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::Primal);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::Primal);
#endif

#ifdef cplex_alg_Dual_Simplex //2
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::Dual);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::Dual);
#endif

#ifdef cplex_alg_Network_Simplex //3
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::Network);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::Network);
#endif

#ifdef cplex_alg_Barrier //4
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::Barrier);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::Barrier);
#endif

#ifdef cplex_alg_Sifting //5
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::Sifting);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::Sifting);
#endif

#ifdef cplex_alg_Concurrent //6
	cplex_MAS.setParam(IloCplex::RootAlg, IloCplex::Concurrent);
    cplex_MAS.setParam(IloCplex::NodeAlg, IloCplex::Concurrent);
#endif

	cplex_MAS.setParam(IloCplex::Threads, number_of_CPU);
    cplex_MAS.setParam(IloCplex::ParallelMode, 1);

    cplex_MAS.setParam(IloCplex::MIPInterval,100); // display interval

    cplex_MAS.setParam(IloCplex::TiLim, TIME_LIMIT);

    cplex_MAS.setOut(env_MAS.getNullStream());

#ifdef desable_cplex_presolve_and_heur
    // turn off presolve and heuristics
    cplex_MAS.setParam(IloCplex::PreInd, false);
    cplex_MAS.setParam(IloCplex::HeurFreq, -1);

    cplex_MAS.setParam(IloCplex::DataCheck, CPX_ON);
#endif // desable_cplex_presolve_and_heur

#ifdef disable_cplex_cuts
    cplex_MAS.setParam(IloCplex::Cliques, -1);
    cplex_MAS.setParam(IloCplex::Covers, -1);
    cplex_MAS.setParam(IloCplex::DisjCuts, -1);
    cplex_MAS.setParam(IloCplex::FlowCovers, -1);
    cplex_MAS.setParam(IloCplex::FlowPaths, -1);
    cplex_MAS.setParam(IloCplex::FracCuts, -1);
    cplex_MAS.setParam(IloCplex::GUBCovers, -1);
    cplex_MAS.setParam(IloCplex::ImplBd, -1);
    cplex_MAS.setParam(IloCplex::MIRCuts, -1);
    cplex_MAS.setParam(IloCplex::ZeroHalfCuts, -1);
#endif // disable_cplex_cuts
}


/*************************************************************************************************/
void master_kHNDP_free()
/*************************************************************************************************/
{

//	pricer_free_cplex();
//
//	master_heur_free(istname);

	objval=best_UB;
	bestobjval=best_UB;
	if(STOP_TIME_LIMIT){
		bestobjval=best_LB;
	}
	nodecount=nb_BB_nodes;

	/////////////////////////////////////////////////////////////
	time_finish=clock();
	computation_time=(double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
	/////////////////////////////////////////////////////////////

	numcols=cplex_MAS.getNcols();
	numrows=cplex_MAS.getNrows();

	lpstat=101;
	if(STOP_TIME_LIMIT)
    {
		lpstat=107;
	}

	for (unsigned i(0); i < BB_nodes.size(); i++)
        delete BB_nodes[i];

	X->end();

	for(unsigned i = 0; i < demandsNumber ; i++)
        (*F)[i].end();

    delete [] F;

    env_MAS.end();

//	for(int i=0;i<itemTypeNumber;i++){
//		delete[] status_item_couple[i];
//	}
//	delete[] status_item_couple;
}


/*************************************************************************************************/
void master_kHNDP_Output()
/*************************************************************************************************/
{
    /////////////////////////////////////////////////////////////////////////////
    cout << endl << endl << endl << "Number of weird solutions = " << nb_f_frac << endl << endl;

    cout << "\nOptimal Value\t" << objval << endl;
	cout << "LP\t" << best_LB << endl;
	cout << "Time\t" << computation_time << endl;
	cout << "Nodes\t" << nodecount << endl;
	cout << "STATUS\t" << lpstat << endl;

	cout << endl << endl << endl;
	cout << G.n << " & " << demandsNumber << " & ";
	cout << best_UB << " & " << best_LB << " & ";
	cout << nb_BB_nodes << " & " << nb_evaluated << " & " << nb_waiting << " & ";
	cout << nb_pruned << " & " << numcols << " & " << numrows << " & ";
	cout << first_LB_time << " & " << computation_time;
	cout << endl << endl;
	//cout << heures << ":" << setw(2) << setfill('0') << minutes << ":";
	//cout << setw(2) << setfill('0') << secondes << endl;
	/////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_Exensive.txt", ios::app);
	compact_file << fixed
			<< name << "\t"
			<< G.n << "\t"
			<< G.m << "\t"
			<< demandsNumber << "\t"
			<< best_UB << "\t"
			<< best_LB << "\t"
//			<< objval << "\t"
//			<<  bestobjval<< "\t"
//			<< nodecount << "\t"
			<< nb_BB_nodes << "\t"
			<< nb_evaluated << "\t"
			<< nb_pruned << "\t"
			<< nb_waiting << "\t"
			<< numcols << "\t"
			<< numrows << "\t"
			<< lpstat << "\t"
			<< first_LB_time << "\t"
			<< computation_time << "\t"
			<< nb_f_frac << "\t"
			<< endl;
	compact_file.close();
	/////////////////////////////////////////////////////////////////////////////
}


/*************************************************************************************************/
int master_kHNDP_heur_solve(IloNumArray const &solX)
/*************************************************************************************************/
{
    for (unsigned i(0); i < G.m; i++)
        G.setValue(G.edgeFromId(i), 0);

//    for (unsigned i(0); i < demandsNumber; i++)
//        G.setValue(G.findEdge(G.nodeFromId(demands[i][0]), G.nodeFromId(demands[i][1])), 1);

    /// setting in 1 all fixed edges
    for (unsigned i(0); i < G.m; i++)
        if (solX[i] > 0.75)
            G.setValue(G.edgeFromId(i), 1);

    //#pragma omp parallel for
    for (unsigned i(0); i < demandsNumber; i++)
    {
        lemonArcWeightMap _AWMs(DGS[i]->m_dg, 0);

        // Adjusting costs!!!
        for (unsigned j(0); j < DGS[i]->m; j++)
        {
            lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

            _AWMs[a] = DGS[i]->weightFromArc(a);

            if(DGS[i]->idFromArc(a) < 276447231)
            {
                lemonEdge e = G.m_Edges[DGS[i]->idFromArc(a)];

                if(G.valueFromEdge(e) == 1)
                    _AWMs[a] = 0;
            }
        }

        //cout << " k = " << k << endl;

        //Calculating minimum cost max flow using NetworkSimplex from LEMON
        NetworkSimplex<lemonDigraph, int> myNetworkSimplex(DGS[i]->m_dg);
        myNetworkSimplex.upperMap(DGS[i]->m_ACM).costMap(_AWMs);
        myNetworkSimplex.stSupply(DGS[i]->m_s, DGS[i]->m_t, k);
        myNetworkSimplex.run();

        // cout << "Min " << k << "-flow cost: " << myNetworkSimplex.totalCost() << endl << endl;

        myNetworkSimplex.flowMap(DGS[i]->m_AVM);

        /// Vérification pour L >= 4
        lemonArcCapacityMap cap(DGS[i]->m_dg, 1);



        while (!feasableAnyL(i,DGS[i]->m_AVM,cap))
        {
            //Calculating minimum cost max flow using NetworkSimplex from LEMON
//                NetworkSimplex<lemonDigraph, int> myNetworkSimplex(DGS[i]->m_dg);
            myNetworkSimplex.upperMap(cap).costMap(_AWMs);
            myNetworkSimplex.stSupply(DGS[i]->m_s, DGS[i]->m_t, k);
            myNetworkSimplex.run();

            myNetworkSimplex.flowMap(DGS[i]->m_AVM);
        }


        // updating superposing values in e
        for(unsigned j(0); j < DGS[i]->m; j++)
        {
            lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

            //cout << "updating edge " << DGS[i]->idFromArc(a) << " verified index " << EIM[LEdges[DGS[i]->idFromArc(a)]] << endl;
            if(DGS[i]->idFromArc(a) < 276447231)
            {
                lemonEdge e = G.m_Edges[DGS[i]->idFromArc(a)];

                if(G.valueFromEdge(e) < 1 && DGS[i]->valueFromArc(a) == 1)
                    G.setValue(e, 1);
            }
        }
    }


    //for (unsigned i = 0; i < demands.size(); i++)
    //write_graph(*DGS[i], *dg_strs[i]);

    //HERE CALCULATING THE OBJECTIVE FORMULA
    int Z(0);



    // Calculating realisable solution and Zopt
    for(unsigned i=0; i < G.m; i++)
        if (G.valueFromEdge(G.m_Edges[i]) == 1)
            Z += G.weightFromEdge(G.m_Edges[i]);

//	//cout << "HEUR MASTER SOLVE\n";
//	status=CPXmipopt(env_MASTER_heur,lp_MASTER_heur);
//	if(status!=0) {
//		printf("error in CPXmipopt\n");
//	}
//
//	double primal_bound;
//	status=CPXgetobjval(env_MASTER_heur,lp_MASTER_heur,&primal_bound);
//	if(status!=0) {
//		primal_bound=binNumber;
//		printf("error in CPXgetmipobjval\n");
//	}

	return Z;

}


/*************************************************************************************************/
int branch_and_bound(int level)
/*************************************************************************************************/
{
	////////////////////////////////////
    if (!STOP_TIME_LIMIT)
    {
        time_finish=clock();
        current_time=(double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
        if(current_time>TIME_LIMIT)
            STOP_TIME_LIMIT=true;

        if (STOP_TIME_LIMIT)
            return -1;
    }
    else
        return -1;
	////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//solve the current LP
	nb_BB_nodes++;

	//if infeasible it returns +infinite
	double current_LP=master_solve_lp();

	cout << "current LP! " << current_LP << endl;

	if(level==0){
		best_LB=current_LP;
		time_finish=clock();
		current_time=(double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
		first_LB_time=current_time;
	}

	int dual_bound=my_ceil(current_LP);

#ifdef	verbose_BP
	time_finish=clock();
	current_time=(double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
	cout << "Level\t" << level << "\tbest_incumbent\t" << best_incumbent << "\tcurrent bound\t"<< current_LP << "\ttime\t"<< current_time << endl;
#endif

	///////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////

    //Declaring current solution
    IloNumArray solX(env_MAS,G.m);
    IloArray<IloNumArray> solF(env_MAS,demandsNumber);

    ///////////////////////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////////////////////
    //try to prune and run heuristics
    if(dual_bound>=best_UB)
    {

#ifdef	verbose_PRUNING
        cout << "PRUNE!\t" << "level\t" << level << "\tcols\t" <<BP_cols << "\tbound\t" << current_LP << "\tbest_incumbent\t" << best_incumbent <<  endl;
#endif

        return 1;

    }
    else
    {

#ifdef	use_heuristic_master
        // getting current X values
        cplex_MAS.getValues(solX,*X);

        //compute an heuristic solution solving to integratility the current master
        int primal_bound=my_floor(master_kHNDP_heur_solve(solX));  /// adjusting SH to the fixed goal
        cout << "primal current value from SH = " << primal_bound << endl;   /// Primal! not initial!
        if(best_UB>primal_bound)
        {
            best_UB=primal_bound;
        }

        if(dual_bound>=best_UB)
        {

#ifdef	verbose_PRUNING
            cout << "PRUNE with heur!\t" << "level\t" << level << "\tcols\t" <<BP_cols << "\tbound\t" << current_LP << "\tbest_incumbent\t" << best_incumbent <<  endl;
#endif

            return 1;
        }
#endif
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////
    // getting the current F values
    for(unsigned i(0); i < demandsNumber; i++)
    {
        solF[i] = IloNumArray(env_MAS,DGS[i]->m);
        cplex_MAS.getValues(solF[i],(*F)[i]);
    }


//	if(status!=0)
//	{
//		printf("error in CPXgetmipx check_integrality\n");
//		exit(-1);
//	}

    //check if the solution is integer then return
    bool integer = master_check_integrality(solX,solF);
    if(integer)
    {

        cout << "integer!!\t" << integer << "\tvalue\t" << dual_bound <<  endl;

        if(best_UB>dual_bound)
        {
            best_UB=dual_bound;
        }

#ifdef	verbose_BP
        cout << "return..." << endl;
#endif

        return 1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
	int item_i;

	cout << "branching phase!" << endl;

	bool _found = branch_kHNDP_select_branching(solX,solF,&item_i);

#ifdef	verbose_BP
	cout << "Couple Branching\t" << item_i << "\t" << item_j << endl;
#endif
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef	verbose_BP
	cout << "together branching" << endl;
#endif

	//together branching
	//set master
	branch_kHNDP_left(item_i);
	branch_and_bound(level+1);

	//reset master
	branch_kHNDP_free(item_i);
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//check bound then return
	if(dual_bound>=best_UB){

#ifdef	verbose_BP
		cout << "return..." << endl;
#endif
		return 1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef	verbose_BP
	cout << "separate branching\t"  << endl;
#endif

	//separate branching
	//set master
	branch_kHNDP_right(item_i);
	branch_and_bound(level+1);

	//reset master
	branch_kHNDP_free(item_i);
	///////////////////////////////////////////////////////////////////////////////////////////////////

    solX.end();
    for(unsigned i(0); i < demandsNumber && integer; i++)
        solF[i].end();
    solF.end();

#ifdef	verbose_BP
	cout << "return..." << endl;
#endif

	return 1;
}


/*************************************************************************************************/
int evaluate_node(double &dual_bound, int &primal_bound, int &i_item)
/*************************************************************************************************/
{
    /***************************************************************
        Solve the current LP and return the status and branching
        item if the solution is fractional
        status can be:
            INFEASIBLE
            INTEGER
            FRACTIONAL
            ERROR
    ****************************************************************/

///////////////////////////////////////////////////////////////////////////////////////////////////
	//solve the current LP

	//if infeasible it returns INFEASIBLE
	dual_bound = master_solve_lp();

	if ( dual_bound == INFINITE )
    {
        dual_bound = 0.0;
        primal_bound = INFINITE;
        return INFEASIBLE;
    }

#ifdef	verbose_BP
	time_finish=clock();
	current_time=(double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
	cout << "Level\t" << level << "\tbest_incumbent\t" << best_incumbent << "\tcurrent bound\t"<< dual_bound << "\ttime\t"<< current_time << endl;
#endif

	///////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////

    //Declaring current solution
    IloNumArray solX(env_MAS,G.m);
    IloArray<IloNumArray> solF(env_MAS,demandsNumber);

    ///////////////////////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////////////////////
    //try to prune and run heuristics
//    if(dual_bound>=best_UB)
//    {
//
//#ifdef	verbose_PRUNING
//        cout << "PRUNE!\t" <<<< "level\t" << level << "\tcols\t" <<BP_cols << "\tbound\t" << dual_bound << "\tbest_incumbent\t" << best_incumbent <<  endl;
//#endif
//        i_item == -1;
//        solX.end();
//        for(unsigned i(0); i < demandsNumber; i++)
//            solF[i].end();
//        solF.end();
//        return PRUNED;
//
//    }
//    else
//    {

        // getting current X values
        cplex_MAS.getValues(solX,*X);
//	}
	///////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////
    // getting the current F values
    for(unsigned i(0); i < demandsNumber; i++)
    {
        solF[i] = IloNumArray(env_MAS,DGS[i]->m);
        cplex_MAS.getValues(solF[i],(*F)[i]);
    }


//	if(status!=0)
//	{
//		printf("error in CPXgetmipx check_integrality\n");
//		exit(-1);
//	}

    //check if the solution is integer then return
    bool integer = master_check_integrality(solX,solF);
    if(integer)
    {

        //cout << "integer!!\t" << integer << "\tvalue\t" << dual_bound <<  endl;

        i_item = -1;
        primal_bound = dual_bound;

#ifdef	verbose_BP
        cout << "return..." << endl;
#endif
        solX.end();
        for(unsigned i(0); i < demandsNumber; i++)
            solF[i].end();
        solF.end();
        return INTEGER;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////

//	cout << "branching phase!" << endl;

#ifdef	use_heuristic_master
        //compute an heuristic solution solving to integratility the current master
        primal_bound=my_floor(master_kHNDP_heur_solve(solX));  /// adjusting SH to the fixed goal
        //cout << "primal current value from SH = " << primal_bound << endl;   /// Primal! not initial!

        /// why pruning here?! A bit strange

//        if(dual_bound>=best_UB)
//        {
//
//#ifdef	verbose_PRUNING
//            cout << "PRUNE with heur!\t" << "level\t" << level << "\tcols\t" <<BP_cols << "\tbound\t" << dual_bound << "\tbest_incumbent\t" << best_incumbent <<  endl;
//#endif
//
//            return PRUNED;
//        }
#endif

	bool _found = branch_kHNDP_select_branching(solX,solF,&i_item);

#ifdef	verbose_BP
	cout << "Couple Branching\t" << item_i << "\t" << item_j << endl;
#endif

    solX.end();
    for(unsigned i(0); i < demandsNumber; i++)
        solF[i].end();
    solF.end();

    return FRACTIONAL;
	///////////////////////////////////////////////////////////////////////////////////////////////////
}


/*************************************************************************************************/
void affect_task()
/*************************************************************************************************/
{
//#ifdef outputFile_by_rank_MPI
//    ofstream procOutFile(proc_OutFile_name, ios::out | ios::trunc);
//#endif // outputFile_by_rank_MPI

    int i_proc, i_task;

    // looking for the first free processor
    for (unsigned i(1); i < world_size; i++)
    {
        if ( slaves_status[i] == PROC_FREE )
        {
            i_proc = i;
            break;
        }
    }

    // looking for the first waiting node
#ifdef BEST_FIRST
    sort(waiting_line.begin(), waiting_line.end(), myOperator);
    i_task = waiting_line[0]->info.id;
#endif // BEST_FIRST

#ifdef FIFO
    for (unsigned i(1); i < BB_nodes.size(); i++)
    {
        if ( BB_nodes[i]->status == WAITING )
        {
            i_task = i;
            break;
        }
    }
#endif // FIFO

#ifdef outputFile_by_rank_MPI
    procOutFile << "Affecting task " << i_task << " to proc " << i_proc << endl;
//    cout << "Size of BB_nodes = " << nb_BB_nodes << endl;
#endif // outputFile_by_rank_MPI

    int *_node_info = new int [3];

    _node_info[0] = BB_nodes[i_task]->info.level;
    _node_info[1] = BB_nodes[i_task]->info.br_var;
    _node_info[2] = BB_nodes[i_task]->info.side;

    MPI_Send(_node_info, 3, MPI_INT, i_proc, WORKTAG, MPI_COMM_WORLD);

    delete [] _node_info;

    /// the rest of the branching chain
    int _level = BB_nodes[i_task]->info.level;
    BB_node *_BBnd = BB_nodes[i_task];
    while ( _level > 1 )
    {
        _BBnd = _BBnd->father;
        _node_info = new int [2];
        _node_info[0] = _BBnd->info.br_var;
        _node_info[1] = _BBnd->info.side;

        MPI_Send(_node_info, 2, MPI_INT, i_proc, COMPLETIONTAG, MPI_COMM_WORLD);

        delete [] _node_info;
        _level--;
    }

//    if ( _BBnd->father->info.side != ROOT )
//        cout << "INCOHERENCE in branching chain completion !  level = " << BB_nodes[i_task]->info.level << " and father id = " << _BBnd->father->info.id << endl;

//    cout << "node " << BB_nodes[i_task]->info.id << " son of node " << BB_nodes[i_task]->father->info.id << " obtained by branching ";
//    if (BB_nodes[i_task]->info.side == LEFT)
//        cout << "left";
//    else if (BB_nodes[i_task]->info.side == RIGHT)
//        cout << "right";
//    else if (BB_nodes[i_task]->info.side == ROOT)
//        cout << "root";
//    cout << " on variable " << BB_nodes[i_task]->info.br_var << endl;

    BB_nodes[i_task]->status = CALCULATING;
    task_affected[i_proc] = BB_nodes[i_task]->info.id;
    nb_waiting--;
    nb_calculating++;
    slaves_status[i_proc] = PROC_COMPUTING;
    nb_proc_waiting--;
    nb_proc_calculating++;

    waiting_line.pop_front();

//#ifdef outputFile_by_rank_MPI
//    procOutFile.close();
//#endif // outputFile_by_rank_MPI
}


/*************************************************************************************************/
void recieve_task()
/*************************************************************************************************/
{
    int *_node_info = new int [3];
    int _level;

    MPI_Recv(_node_info, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if ( status.MPI_TAG == WORKTAG )
    {
        branching_chain.clear();

        _level = _node_info[0];
        branching_chain.push_back(make_pair(_node_info[1],_node_info[2]));

#ifdef outputFile_by_rank_MPI
        procOutFile << "Work info:" << endl;
        procOutFile << "   level = " << _level << endl;
        procOutFile << "var" << branching_chain[branching_chain.size()-1].first << "_s" << branching_chain[branching_chain.size()-1].second;
#endif // outputFile_by_rank_MPI

        while ( _level > 1 )
        {
            MPI_Recv(_node_info, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if ( status.MPI_TAG != COMPLETIONTAG )
                cout << "Supposed to be a completion tag!!!! " << endl;

            branching_chain.push_back(make_pair(_node_info[0],_node_info[1]));

#ifdef outputFile_by_rank_MPI
        procOutFile << " --> var" << branching_chain[branching_chain.size()-1].first << "_s" << branching_chain[branching_chain.size()-1].second;
#endif // outputFile_by_rank_MPI

            _level--;
        }
#ifdef outputFile_by_rank_MPI
        procOutFile << endl;
#endif // outputFile_by_rank_MPI
    }

    delete [] _node_info;
}


/*************************************************************************************************/
void return_results(int stat, double dual_bound, int primal_bound, int i_item)
/*************************************************************************************************/
{
#ifdef outputFile_by_rank_MPI
    procOutFile << "Slave " << world_rank << " sending back: " << endl;
    procOutFile << "   status = " << stat << endl;
    procOutFile << "   dual_bound = " << dual_bound << endl;
    procOutFile << "   primal_bound = " << primal_bound << endl;
    procOutFile << "   i_item = " << i_item << endl;
#endif // outputFile_by_rank_MPI

    int *results = new int[2];

    results[0] = primal_bound;
    results[1] = i_item;

    //MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    if ( stat == FRACTIONAL )
    {
        MPI_Send(results, 2, MPI_INT, 0, stat, MPI_COMM_WORLD);
        MPI_Send(&dual_bound, 1, MPI_DOUBLE, 0, stat, MPI_COMM_WORLD);
    }
    if ( stat == INTEGER )
    {
        //cout << "Integer! :O" << endl;
        MPI_Send(results, 2, MPI_INT, 0, stat, MPI_COMM_WORLD);
    }
    else if ( stat == INFEASIBLE )
    {
        MPI_Send(results, 2, MPI_INT, 0, stat, MPI_COMM_WORLD);
        //cout << "Problem not feasible :/" << endl;
    }
    else if ( stat == LP_ERROR )
    {
        results[0] = -1;
        results[1] = -1;
        MPI_Send(results, 2, MPI_INT, 0, stat, MPI_COMM_WORLD);
        //cout << "Solving ERROR :(" << endl;
    }

    delete [] results;
}


/*************************************************************************************************/
bool update_global_LB(PTR_BB_node _node)
/*************************************************************************************************/
{
    if ( _node->info.side == RIGHT )
    {
        if ( _node->father->left->status == EVALUATED )
        {
            double _min = min(_node->dual_bound, _node->father->left->dual_bound);

            if (_min > _node->father->dual_bound)
            {
                _node->father->dual_bound = _min;
                return update_global_LB(_node->father);
            }
        }
    }
    else if ( _node->info.side == LEFT )
    {
        if ( _node->father->right->status == EVALUATED )
        {
            double _min = min(_node->dual_bound, _node->father->right->dual_bound);

            if (_min > _node->father->dual_bound)
            {
                _node->father->dual_bound = _min;
                return update_global_LB(_node->father);
            }
        }
    }
    else
    {
        // root node
        if ( _node->dual_bound > best_LB )  //eventually always true
        {
            best_LB = _node->dual_bound;
            return true;
        }
    }

    return false;
}


/*************************************************************************************************/
void recieve_results()
/*************************************************************************************************/
{
//#ifdef outputFile_by_rank_MPI
//    ofstream procOutFile(proc_OutFile_name, ios::out | ios::trunc);
//#endif // outputFile_by_rank_MPI

    double dual_bound;
    int *results = new int[2];

    bool flag = false;

    MPI_Recv(results, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

#ifdef outputFile_by_rank_MPI
    procOutFile << "Master recieved result with status: " << status.MPI_TAG << endl;
#endif // outputFile_by_rank_MPI

    if ( status.MPI_TAG == FRACTIONAL )
    {
        MPI_Recv(&dual_bound, 1, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

        //cout << "Master completing result with status: " << status.MPI_TAG << endl;

        int i_node = task_affected[status.MPI_SOURCE];

        task_affected[status.MPI_SOURCE] = -1;
        slaves_status[status.MPI_SOURCE] = PROC_FREE;
        nb_proc_calculating--;
        nb_proc_waiting++;

        //cout << "Freed proc " << status.MPI_SOURCE << endl;

        BB_nodes[i_node]->dual_bound = dual_bound;
        BB_nodes[i_node]->primal_bound = results[0];

        if(my_ceil(dual_bound) > best_UB)
        {
#ifdef	verbose_PRUNING
            cout << "PRUNE!\t" << "level\t" << level << "\tcols\t" <<BP_cols << "\tbound\t" << dual_bound << "\tbest_incumbent\t" << best_incumbent <<  endl;
#endif //verbose_PRUNING
#ifdef outputFile_by_rank_MPI
            procOutFile << "Pruned! Dual_bound = " << dual_bound << " best_UB = " << best_UB << endl;
#endif // outputFile_by_rank_MPI

            nb_pruned++;
        }
        else
        {
            add_children(i_node,results[1]);
            //cout << "New size of tree = " << BB_nodes.size() << endl;
        }

        flag = update_global_LB(BB_nodes[i_node]);

        if (BB_nodes[i_node]->primal_bound < best_UB)
        {
            best_UB = BB_nodes[i_node]->primal_bound;
            flag = true;
        }

        //cout << "updating task " << task_affected[status.MPI_SOURCE] << endl;

        BB_nodes[i_node]->status = EVALUATED;

        nb_evaluated++;
        nb_calculating--;
    }
    else if ( status.MPI_TAG == INTEGER )
    {
        int i_node = task_affected[status.MPI_SOURCE];

        task_affected[status.MPI_SOURCE] = -1;
        slaves_status[status.MPI_SOURCE] = PROC_FREE;
        nb_proc_calculating--;
        nb_proc_waiting++;

        BB_nodes[i_node]->dual_bound = results[0];
        BB_nodes[i_node]->primal_bound = results[0];

        /// mettre à jour car on prend le min des deux frères que quand les deux sont évalués
        flag = update_global_LB(BB_nodes[i_node]);

        if(BB_nodes[i_node]->primal_bound < best_UB)
        {
            best_UB = BB_nodes[i_node]->primal_bound;
            flag = true;
        }

#ifdef outputFile_by_rank_MPI
            procOutFile << "Integer! Dual_bound = " << BB_nodes[i_node]->dual_bound << " primal_bound = " << BB_nodes[i_node]->primal_bound << endl;
#endif // outputFile_by_rank_MPI

        BB_nodes[i_node]->status = EVALUATED;
        nb_evaluated++;
        nb_calculating--;
    }
    else if ( status.MPI_TAG == INFEASIBLE )
    {
        int i_node = task_affected[status.MPI_SOURCE];

        task_affected[status.MPI_SOURCE] = -1;
        slaves_status[status.MPI_SOURCE] = PROC_FREE;
        nb_proc_calculating--;
        nb_proc_waiting++;

        BB_nodes[i_node]->dual_bound = 0.0;
        BB_nodes[i_node]->primal_bound = INFINITE;

#ifdef outputFile_by_rank_MPI
            procOutFile << "Infeasible!" << endl;
#endif // outputFile_by_rank_MPI

        BB_nodes[i_node]->status = EVALUATED;
        nb_evaluated++;
        nb_calculating--;
    }
    else if ( status.MPI_TAG == LP_ERROR )
    {
        int i_node = task_affected[status.MPI_SOURCE];

        task_affected[status.MPI_SOURCE] = -1;
        slaves_status[status.MPI_SOURCE] = PROC_FREE;
        nb_proc_calculating--;
        nb_proc_waiting++;

        BB_nodes[i_node]->dual_bound = 0.0;
        BB_nodes[i_node]->primal_bound = INFINITE;

#ifdef outputFile_by_rank_MPI
            procOutFile << "Error while solving!" << endl;
#endif // outputFile_by_rank_MPI

        BB_nodes[i_node]->status = EVALUATED;
        nb_evaluated++;
        nb_calculating--;
    }

    delete [] results;

    if (flag)
        cout << "Primal     " << best_UB << "   dual    " << best_LB << endl;

    if (best_UB <= my_ceil(best_LB)) /// Or time limit!
    {
#ifdef outputFile_by_rank_MPI
        procOutFile << "PROBLEM SOLVED!" << endl;
        procOutFile << "   best_UB = " << best_UB << endl;
        procOutFile << "   best_LB = " << best_LB << endl;
#endif // outputFile_by_rank_MPI

        solved = true;
    }

//#ifdef outputFile_by_rank_MPI
//    procOutFile.close();
//#endif // outputFile_by_rank_MPI
}


/*************************************************************************************************/
void get_ready_for_solving()
/*************************************************************************************************/
{
    for (unsigned i(0); i < branching_chain.size(); i++)
    {
        if ( branching_chain[i].second == LEFT )
            branch_kHNDP_left(branching_chain[i].first);
        else
            branch_kHNDP_right(branching_chain[i].first);
    }
}


/*************************************************************************************************/
void restore_master_LP()
/*************************************************************************************************/
{
    for (unsigned i(0); i < branching_chain.size(); i++)
        branch_kHNDP_free(branching_chain[i].first);
}


/*************************************************************************************************/
void master()
/*************************************************************************************************/
{
//#ifdef outputFile_by_rank_MPI
//    ofstream procOutFile(proc_OutFile_name, ios::out | ios::trunc);
//#endif // outputFile_by_rank_MPI

    solved          =   false;
    best_UB         =   kHNDP_first_fit();
    nb_BB_nodes     =   0;
    nb_evaluated    =   0;
    nb_calculating  =   0;
    nb_waiting      =   0;
    nb_pruned       =   0;
    nb_f_frac       =   0;
    best_LB         =   -1;

    int	        i_rank;
    double      dual_bound;

    int i_item, primal_bound;

    // creating root node
    PTR_BB_node _node;

    _node = new BB_node;

    _node->father       =   NULL;
    _node->info.id      =   nb_BB_nodes;
    _node->info.level   =   0;
    _node->info.br_var  =   INFINITE;
    _node->info.side    =   ROOT;
    _node->left         =   NULL;
    _node->right        =   NULL;
    _node->status       =   EVALUATED;

    BB_nodes.push_back(_node);

    nb_BB_nodes++;

    time_start=clock();

    int stat = evaluate_node(dual_bound, primal_bound, i_item);
    nb_evaluated++;

    if ( stat == INTEGER )
    {
        best_UB = dual_bound;
        BB_nodes[0]->dual_bound     = dual_bound;
        BB_nodes[0]->primal_bound   = dual_bound;
        cout << "Optimum from root! :O" << endl;
        goto CLOSING;
    }
    else if ( stat == INFEASIBLE )
    {
        cout << "Problem not feasible :/" << endl;
        goto CLOSING;
    }
    else if ( stat == LP_ERROR )
    {
        cout << "Solving root ERROR :(" << endl;
        goto CLOSING;
    }
    else // stat == FRACTIONAL
        ;

    BB_nodes[0]->dual_bound     =   dual_bound;
    BB_nodes[0]->primal_bound   =   primal_bound;

    display_BB_node(0);

    if(best_UB>primal_bound)
        best_UB=primal_bound;

    best_LB = dual_bound;

    cout << "first branching node " << i_item << endl;

    // branching by creating and adding two nodes to the tree
    add_children(0, i_item);

    time_finish     =   clock();
	current_time    =   (double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
	first_LB_time   =   current_time;

	////////////////////////////////////
    if (!STOP_TIME_LIMIT)
    {
        if(current_time > TIME_LIMIT)
            STOP_TIME_LIMIT=true;

        if (STOP_TIME_LIMIT)
            goto CLOSING;
    }
    else
        goto CLOSING;
	////////////////////////////////////

    /*
    * Seed the slaves.PRUNED
    *
    * IF WE CAN SEED ALL SLAVES FROM THE BEGINNING,
    * THE COMMUNICATION CAN BE BETTER, BECAUSE WE
    * AFFECT TASK IF AND ONLY IF A PROCESSOR SEND
    * BACK RESULT  DETECT HIS ID THROUGHT:
    *              status.MPI_SOURCE
    */
//        MPI_Send(&work,         /* message buffer */ /* need to get_next_work_request */
//        1,              /* one data item */
//        MPI_INT,        /* data item is an integer */
//        i_rank,         /* destination process rank */
//        WORKTAG,        /* user chosen message tag */
//        MPI_COMM_WORLD);/* always use this */
/*    for (i_rank = 1; i_rank < world_size; ++i_rank)
    {
        work = i_rank * i_rank;

        MPI_Send(&work, 1, MPI_INT, i_rank, WORKTAG, MPI_COMM_WORLD);
    }*/

    /*
    * Receive a result from any slave and dispatch a new work
    * request work requests have been exhausted.
    */
//        MPI_Recv(&result,       /* message buffer */
//        1,              /* one data item */
//        MPI_DOUBLE,     /* of type double real */
//        MPI_ANY_SOURCE, /* receive from any sender */
//        MPI_ANY_TAG,    /* any type of message */
//        MPI_COMM_WORLD, /* always use this */
//        &status);       /* received message info */
    while ( (nb_waiting > 0 || nb_calculating > 0) && !solved ) /* valid new work request */
    {
        ////////////////////////////////////
        if (!STOP_TIME_LIMIT)
        {
            time_finish=clock();
            current_time=(double)(time_finish-time_start)/(double)CLOCKS_PER_SEC;
            if(current_time > TIME_LIMIT)
                STOP_TIME_LIMIT=true;

            if (STOP_TIME_LIMIT)
                goto PRECLOSING;
        }
        else
            goto PRECLOSING;
        ////////////////////////////////////

        while ( nb_waiting > 0 && nb_proc_waiting > 0 && !solved )
        {
#ifdef outputFile_by_rank_MPI
            procOutFile << "Master affecting task..." << endl;
#endif // outputFile_by_rank_MPI

            ////////////////////////////////////////////
            affect_task();
            ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
            procOutFile << "task sent." << endl;
#endif // outputFile_by_rank_MPI
        }
#ifdef outputFile_by_rank_MPI
        procOutFile << "Master recieving results..." << endl;
#endif // outputFile_by_rank_MPI

        ////////////////////////////////////////////
        recieve_results();
        ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
        procOutFile << "Results recieved." << endl;
#endif // outputFile_by_rank_MPI
    }
    /*
    * Receive results for outstanding work requests.
    */
//    for (i_rank = 1; i_rank < world_size; ++i_rank)
//    {
//        MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//        cout << "Master recieved " << result << " from process " << status.MPI_SOURCE << endl;
//    }

    /*
    *   If some slaves are still computing recieve wait for their results
    */
PRECLOSING:

    while( nb_proc_calculating > 1 )
    {
        recieve_results();
    }


    /*
    * Tell all the slaves to exit.
    */
CLOSING:

    for (i_rank = 1; i_rank < world_size; ++i_rank)
    {
        MPI_Send(0, 0, MPI_INT, i_rank, DIETAG, MPI_COMM_WORLD);
    }

//#ifdef outputFile_by_rank_MPI
//    procOutFile.close();
//#endif // outputFile_by_rank_MPI
}


/*************************************************************************************************/
void slave()
/*************************************************************************************************/
{
//#ifdef outputFile_by_rank_MPI
//    ofstream procOutFile(proc_OutFile_name, ios::out | ios::trunc);
//#endif // outputFile_by_rank_MPI

    int stat;
    double dual_bound;
    int i_item, primal_bound;
    int work;

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
    processor_name, world_rank, world_size);

    while (true)
    {
#ifdef outputFile_by_rank_MPI
        procOutFile << "Slave recieving task..." << endl;
#endif // outputFile_by_rank_MPI

        ////////////////////////////////////////////
        recieve_task();
        ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
        procOutFile << "Task recieved." << endl;
#endif // outputFile_by_rank_MPI

        /*
        * Check the tag of the received message.
        */
        if (status.MPI_TAG == DIETAG)
        {
#ifdef outputFile_by_rank_MPI
            procOutFile << "Slave to diet!" << endl;
#endif // outputFile_by_rank_MPI
            return;
        }

#ifdef outputFile_by_rank_MPI
        procOutFile << "Slave " << world_rank << " recieved " << status.MPI_TAG << " from process " << status.MPI_SOURCE << endl;
        procOutFile << "Slave preparing LP to be solved..." << endl;
#endif // outputFile_by_rank_MPI

        ////////////////////////////////////////////
        get_ready_for_solving();
        ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
        procOutFile << "Slave evaluating node..." << endl;
#endif // outputFile_by_rank_MPI

        ////////////////////////////////////////////
        stat = evaluate_node(dual_bound, primal_bound, i_item);
        ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
        procOutFile << "Node evaluated." << endl;
        procOutFile << "Slave retoring initial LP..." << endl;
#endif // outputFile_by_rank_MPI

        ////////////////////////////////////////////
        restore_master_LP();
        ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
        procOutFile << "Slave returning results..." << endl;
#endif // outputFile_by_rank_MPI

        ////////////////////////////////////////////
        return_results(stat, dual_bound, primal_bound, i_item);
        ////////////////////////////////////////////

#ifdef outputFile_by_rank_MPI
        procOutFile << "Results returned." << endl;
        procOutFile << "------------------------------------------" << endl;
#endif // outputFile_by_rank_MPI
    }
//#ifdef outputFile_by_rank_MPI
//    procOutFile.close();
//#endif // outputFile_by_rank_MPI
}
