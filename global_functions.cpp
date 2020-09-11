#include <lemon/network_simplex.h>

#include "global_functions.h"
#include "control_panel.h"


/*************************************************************************************************/
void SORT_NON_INCR_INT(int *item,int *score,int n)
/*************************************************************************************************/
{
	int salto,i,j,tempItem;
	int tempScore;

	for (salto = n/2; salto > 0; salto /=2)
		for (i = salto; i < n; i++)
			for (j = i-salto; j >= 0; j-=salto)
			{
				if (score[j] >= score[j+salto]) break;
				tempScore = score[j];
				score[j]=score[j+salto];
				score[j+salto]=tempScore;
				tempItem = item[j];
				item[j]=item[j+salto];
				item[j+salto]=tempItem;
			}
}


/*************************************************************************************************/
void SORT_NON_INCR_DOUBLE(int *item,int *score,int n)
/*************************************************************************************************/
{
	int salto,i,j,tempItem;
	int tempScore;

	for (salto = n/2; salto > 0; salto /=2)
		for (i = salto; i < n; i++)
			for (j = i-salto; j >= 0; j-=salto)
			{
				if (score[j] >= score[j+salto]) break;
				tempScore = score[j];
				score[j]=score[j+salto];
				score[j+salto]=tempScore;
				tempItem = item[j];
				item[j]=item[j+salto];
				item[j+salto]=tempItem;
			}
}


/*************************************************************************************************/
template <>
bool feasableAnyL(unsigned const pos, lemonArcValMap const &VARVM, lemonArcCapacityMap &VARCM)
/*************************************************************************************************/
{
    bool flag(true);
    lemonEdgeValMap occ(G.m_g, 0);

    for (unsigned i(0); i < DGS[pos]->m; i++)
    {
        lemonArc a = DGS[pos]->arcFromId(i);

        if(DGS[pos]->idFromArc(a) < 276447231)
        {
            lemonEdge e = G.edgeFromId(DGS[pos]->idFromArc(a));

//            cout << "occ[e] = " << occ[e] << endl;
            if (VARVM[a] > 0)
                occ[e]++;

            if (occ[e] > 1)
            {
                VARCM[a] = 0;
                flag = false;
            }
        }
    }

    return (flag);
}


/*************************************************************************************************/
int kHNDP_first_fit()
/*************************************************************************************************/
{
    /// setting in 1 all fixed edges
    for (unsigned i(0); i < demands.size(); i++)
        G.setValue(G.findEdge(G.nodeFromId(demands[i][0]), G.nodeFromId(demands[i][1])), 1);

    for (unsigned i(0); i < demandsNumber; i++)
    {
        for (unsigned j(0); j < DGS[i]->m; j++)
        {
            lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

            if(DGS[i]->idFromArc(a) < 276447231)
            {
                lemonEdge e = G.m_Edges[DGS[i]->idFromArc(a)];

                if(G.valueFromEdge(e) == 1)
                    DGS[i]->setWeight(a, 0);
            }
        }

        //cout << " k = " << k << endl;

        //Calculating minimum cost max flow using NetworkSimplex from LEMON
        NetworkSimplex<lemonDigraph, int> myNetworkSimplex(DGS[i]->m_dg);
        myNetworkSimplex.upperMap(DGS[i]->m_ACM).costMap(DGS[i]->m_AWM);
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
            myNetworkSimplex.upperMap(cap).costMap(DGS[i]->m_AWM);
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

    long Z(0);

    for(unsigned i=0; i < G.m; i++)
        if (G.valueFromEdge(G.m_Edges[i]) == 1)
            Z += G.weightFromEdge(G.m_Edges[i]);

    return(Z);
}


/*************************************************************************************************/
bool extract_number(string &str, int &number)
/*************************************************************************************************/
{
    string temp;

    unsigned i=0;

    while ( i < str.size() )
    {
        if(isdigit(str[i]))
        {
            while (i < str.size() && (isdigit(str[i]) || str[i] == '.'))
            {
                temp += str[i];
                i++;
            }
/*            else if (str[i] == '.')
            {
                str.erase(str.begin() + i);
                while(isdigit(str[i]) && str[i] == '0')
                    str.erase(str.begin() + i);
            }*/
            istringstream stream(temp);
            stream >> number;

            if(i == str.size())
                str.clear();
            else
                str.erase( 0, i);

            return(true);
        }

        i++;
    }

    str.clear();
    return(false);
}


/*************************************************************************************************/
void kHNDP_instance_read(const char *probname)
/*************************************************************************************************/
{
    string s(pathNameI + probname + ".tsp");

    cout << s << endl;

    ifstream fichier(s.c_str(), ios::in);

    int number;
    vector <unsigned> temp;
    vector < vector <unsigned> > v;
    string ligne;  // déclaration d'une chaîne qui contiendra la ligne lue
    unsigned flag(0);

    G.m_name = probname;

    while (getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
    {
        if(flag == 1)
        {
            while (!ligne.empty())
                if (extract_number(ligne, number))
                    temp.push_back(number);

            if (temp.size() == 3)
                v.push_back(vector <unsigned> (temp));
            else if(temp.size() == 0)
                flag = 0;
        }
        else if (flag == 2)
        {
            while (!ligne.empty())
                if (extract_number(ligne, number))
                    temp.push_back(number);

            if (temp.size() == 3 && temp[0]!=temp[1])
            {
                //k = temp[2];
                temp.erase(temp.begin() + 2);
                demands.push_back(vector <unsigned> (temp));
            }
        }
        else
        {
            if(!ligne.compare("NODE_COORD_SECTION"))
                flag = 1;
            else if (!ligne.compare("DEMAND_SECTION"))
                flag = 2;
        }

        temp.clear();
    }

    fichier.close();

    for (unsigned i=0; i < v.size(); i++)
    {
        G.m_Nodes.push_back(G.m_g.addNode());

        G.m_linking[v[i][1]] = i;
        G.m_NIM[G.m_Nodes[G.m_Nodes.size()-1]] = i;
        G.m_NDM[G.m_Nodes[G.m_Nodes.size()-1]] = v[i][1];
        G.m_NVM[G.m_Nodes[G.m_Nodes.size()-1]] = v[i][2];

        G.n++;
    }

    for (unsigned i=0; i < v.size(); i++)
    {
        for (unsigned j=i+1; j < v.size(); j++)
        {
            G.m_Edges.push_back(G.m_g.addEdge( G.m_Nodes[i], G.m_Nodes[j]));

            linkX.push_back(NBCOLS);
            NBCOLS++;

            G.m_EIM[G.m_Edges[G.m_Edges.size()-1]] = G.m;
            G.m_ECM[G.m_Edges[G.m_Edges.size()-1]] = 1;

            double dist = (G.m_NDM[G.m_Nodes[i]] - G.m_NDM[G.m_Nodes[j]]) * (G.m_NDM[G.m_Nodes[i]] - G.m_NDM[G.m_Nodes[j]]);
            dist += (G.m_NVM[G.m_Nodes[i]] - G.m_NVM[G.m_Nodes[j]]) * (G.m_NVM[G.m_Nodes[i]] - G.m_NVM[G.m_Nodes[j]]);
            dist = floor(sqrt(dist) + 0.5);

            G.m_EWM[G.m_Edges[G.m_Edges.size()-1]] = dist;

            G.m++;
        }
    }

    for(unsigned i=0; i < demands.size(); i++)
    {
        demands[i][0]--;
        demands[i][1]--;
    }

    G.n = countNodes(G.m_g);
    G.m = countEdges(G.m_g);

    demandsNumber = demands.size();

#ifdef print_network_requirements
    cout << " n = " << G.n << "       m = " << G.m << "    d = " << demands.size() << endl;
    G.write_graph();
    for(unsigned i(0); i < demandsNumber; i++)
        cout << "d" << i << ": " << demands[i][0] << "--" << demands[i][1] << endl;
#endif

    for (unsigned i(0); i < demandsNumber; i++)
    {
        DGS.push_back(new digraph);
        transform_graph(demands[i][0], demands[i][1], i);
    }
}


/*************************************************************************************************/
template <typename VERTEXLIST>
unsigned nodeFromMyId(VERTEXLIST vertices, int index)
/*************************************************************************************************/
{
    for(unsigned i=0; i < vertices.size(); i++)
    {
        if(G.m_NIM[vertices[i]]==index)
            return(i);
    }
    return(-1);
}


/*************************************************************************************************/
void transform_graph(int const s_i, int const t_i, unsigned const pos)
/*************************************************************************************************/
{
    DGS[pos]->m_s = DGS[pos]->m_dg.addNode();
    DGS[pos]->m_DNIM[DGS[pos]->m_s] = G.m_NIM[G.m_Nodes[s_i]];
    DGS[pos]->m_linking[G.m_NIM[G.m_Nodes[s_i]]] = 0;

    DGS[pos]->m_t = DGS[pos]->m_dg.addNode();
    DGS[pos]->m_DNIM[DGS[pos]->m_t] = G.m_NIM[G.m_Nodes[t_i]];
    DGS[pos]->m_linking[G.m_NIM[G.m_Nodes[t_i]]] = 1;

    vector<lemonNode> LNodes_st(G.m_Nodes);

    //cout << "erasing the " << find_vertex(LNodes_st, g_str, s_i) << " vertex" << endl;
    LNodes_st.erase(LNodes_st.begin() + nodeFromMyId(LNodes_st, s_i));
    //cout << "erasing the " << find_vertex(LNodes_st, g_str, t_i) << " vertex" << endl;
    LNodes_st.erase(LNodes_st.begin() + nodeFromMyId(LNodes_st, t_i));

    int supplement(0), cpt(1);

    for (unsigned i = 0; i < L-1; i++)
    {
        DGS[pos]->m_DNodes.push_back(vector <lemonDiNode> ());

        for (unsigned j = 0; j < LNodes_st.size(); j++)
        {
            DGS[pos]->m_DNodes[i].push_back(DGS[pos]->m_dg.addNode());
            cpt++;
            DGS[pos]->m_DNIM[DGS[pos]->m_DNodes[i][j]] = G.m_NIM[LNodes_st[j]] + supplement;
            DGS[pos]->m_linking[G.m_NIM[LNodes_st[j]] + supplement] = cpt;
        }
        supplement+=1000000;
    }

    //for(map<int,int>::const_iterator it = DGS[pos]->m_linkiDGS[pos]->begin(); it != DGS[pos]->m_linkiDGS[pos]->end(); ++it)
        //cout << it->first << " " << it->second << endl;
    linkF.push_back(vector<unsigned> ());

    for (unsigned i = 0; i < G.m_Edges.size(); i++)
    {
        int source_i, target_i;
        source_i = G.m_NIM[G.m_g.u(G.m_Edges[i])];
        target_i = G.m_NIM[G.m_g.v(G.m_Edges[i])];

        if(source_i == s_i) // case 1: arc sX
        {
            if(target_i == t_i) // case 1.1: arc st
            {
                int new_pos_x = DGS[pos]->m_linking[s_i];
                int new_pos_y = DGS[pos]->m_linking[t_i];
                DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x), DGS[pos]->m_dg.nodeFromId(new_pos_y)));

                linkF[linkF.size()-1].push_back(NBCOLS);
                NBCOLS++;

                DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
                DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
                DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
            }
            else // case 1.2: arc su
            {
                int new_pos_x = DGS[pos]->m_linking[s_i];
                int new_pos_y = DGS[pos]->m_linking[target_i];
                DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

                linkF[linkF.size()-1].push_back(NBCOLS);
                NBCOLS++;

                DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
                DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
                DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
            }
        }
        else if(target_i == s_i) // case 2: arc Xs
        {
            if(source_i == t_i) // case 2.1: arc ts
            {
                int new_pos_x = DGS[pos]->m_linking[s_i];
                int new_pos_y = DGS[pos]->m_linking[t_i];
                DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

                linkF[linkF.size()-1].push_back(NBCOLS);
                NBCOLS++;

                DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
                DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
                DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
            }
            else // case 2.2: arc us
            {
                int new_pos_x = DGS[pos]->m_linking[s_i];
                int new_pos_y = DGS[pos]->m_linking[source_i];
                DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

                linkF[linkF.size()-1].push_back(NBCOLS);
                NBCOLS++;

                DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
                DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
                DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
            }
        }
        else if(target_i == t_i) // case 3: arc ut
        {
            int new_pos_x = DGS[pos]->m_linking[source_i + 1000000 * (L-2)];
            int new_pos_y = DGS[pos]->m_linking[t_i];
            DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

            linkF[linkF.size()-1].push_back(NBCOLS);
            NBCOLS++;

            DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
            DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
            DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
        }
        else if(source_i == t_i) // case 4: arc tu
        {
            int new_pos_x = DGS[pos]->m_linking[1000000 * (L-2) + target_i];
            int new_pos_y = DGS[pos]->m_linking[t_i];
            DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

            linkF[linkF.size()-1].push_back(NBCOLS);
            NBCOLS++;

            DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
            DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
            DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
        }
        else // case 5: arc uv
        {
            supplement = 0;
            for (unsigned j = 0; j < L-2; j++)
            {
                int new_pos_x = DGS[pos]->m_linking[source_i + supplement];
                int new_pos_y = DGS[pos]->m_linking[target_i + 1000000 + supplement];
                DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

                linkF[linkF.size()-1].push_back(NBCOLS);
                NBCOLS++;

                DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
                DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
                DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];

                new_pos_x = DGS[pos]->m_linking[target_i + supplement];
                new_pos_y = DGS[pos]->m_linking[source_i + 1000000 + supplement];
                DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

                linkF[linkF.size()-1].push_back(NBCOLS);
                NBCOLS++;

                DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EIM[G.m_Edges[i]];
                DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_ECM[G.m_Edges[i]];
                DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = G.m_EWM[G.m_Edges[i]];
                supplement+=1000000;
            }
        }
    }

    /// THIS CAN BE DONE UP!!!! (TO VERIFY CAPACITY AND WEIGHT IN ARCS?!)
    for (unsigned i = 0; i < LNodes_st.size(); i++)
    {
        supplement = 0;
        for (unsigned j = 0; j < L-2; j++)
        {
            int new_pos_x = DGS[pos]->m_linking[G.m_NIM[LNodes_st[i]] + supplement];
            int new_pos_y = DGS[pos]->m_linking[G.m_NIM[LNodes_st[i]] + 1000000 + supplement];
            DGS[pos]->m_Arcs.push_back(DGS[pos]->m_dg.addArc(DGS[pos]->m_dg.nodeFromId(new_pos_x),DGS[pos]->m_dg.nodeFromId(new_pos_y)));

            linkF[linkF.size()-1].push_back(NBCOLS);
            NBCOLS++;

            DGS[pos]->m_AIM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = 276447231;
            DGS[pos]->m_ACM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = 1;
            DGS[pos]->m_AWM[DGS[pos]->m_Arcs[DGS[pos]->m_Arcs.size()-1]] = 0;
            supplement+=1000000;
        }
    }

    DGS[pos]->n = countNodes(DGS[pos]->m_dg);
    DGS[pos]->m = countArcs(DGS[pos]->m_dg);

    //DGS[pos]->write_graph();
}


/*************************************************************************************************/
void kHNDP_instance_free()
/*************************************************************************************************/
{
	for (unsigned i=0; i < demandsNumber; i++)
        delete DGS[i];

    delete sol;
}
