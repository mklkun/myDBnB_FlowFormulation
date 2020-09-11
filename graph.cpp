#include "graph.h"

#include <lemon/lgf_writer.h>


void graph::write_graph_viz(string s)
{
    ofstream fichier(s.c_str(), ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    if(fichier)
    {
        fichier << "graph G {" << endl;

        for(unsigned i(0); i < n; i++)
            fichier << i << " ;" << endl;                                                                                  /// ADDING NODE PROPERTIES

        for(unsigned i(0); i < m; i++)
            fichier << m_g.id(m_g.u(m_g.edgeFromId(i))) << "--" << m_g.id(m_g.v(m_g.edgeFromId(i))) << " ;" << endl;                   /// ADDING EDGE PROPERTIES

        fichier << "}" << endl;
        fichier.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
    //for ()
    //cout << " [index=\"" << g_str.NIM[v] << "\", degree=\"" << g_str.NDM[e] << "\"]";

    //cout << " [index=\"" << g_str.EIM[e] << "\", weight=\"" << g_str.EWM[e] << "\", capacity=\"" << g_str.ECM[e] << "\"]";
}


bool exists_in_position(int const indice, vector < vector <int> > const demands, int const pos)
{
    for(unsigned i(0); i < demands.size(); i++)
        if (demands[i][pos] == indice)
            return(true);
    return(false);
}


void graph::write_graph_viz_sol(vector < vector <int> > const demands, string const s)
{
    ofstream fichier(s.c_str(), ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    if(fichier)
    {
        fichier << "graph G {" << endl;

        for(unsigned i(0); i < n; i++)
        {
            if (exists_in_position(i,demands,0))
                fichier << i << " [color=blue];" << endl;
            else if (exists_in_position(i,demands,1))
                fichier << i << " [color=red];" << endl;
            else
                fichier << i << " ;" << endl;
        }

        for(unsigned i(0); i < m; i++)
            if(m_EVM[m_g.edgeFromId(i)] == 1)
                fichier << m_g.id(m_g.u(m_g.edgeFromId(i))) << "--" << m_g.id(m_g.v(m_g.edgeFromId(i))) << " ;" << endl;

        fichier << "}" << endl;
        fichier.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
}


void graph::write_graph()
{
    graphWriter(m_g, cout)
    .nodeMap("Index", m_NIM)
    .nodeMap("Degree", m_NDM)
    .edgeMap("Index", m_EIM)
    .edgeMap("Capacity", m_ECM)
    .edgeMap("Weight", m_EWM)
    .edgeMap("Value", m_EVM)
    .run();
}


void graph::write_graph(string s)
{
    ofstream filename(s.c_str());

    graphWriter(m_g, filename)
    .nodeMap("Index", m_NIM)
    .nodeMap("Degree", m_NDM)
    .edgeMap("Index", m_EIM)
    .edgeMap("Capacity", m_ECM)
    .edgeMap("Weight", m_EWM)
    .run();
}




graph::~graph()
{
    //dtor
}
