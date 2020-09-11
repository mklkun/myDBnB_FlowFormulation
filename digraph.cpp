#include "digraph.h"

#include <lemon/lgf_writer.h>


void digraph::write_graph_viz(string s)
{
    ofstream fichier(s.c_str(), ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    if(fichier)
    {
        fichier << "graph G {" << endl;

        for(unsigned i(0); i < n; i++)
            fichier << i << " ;" << endl;                                                                                  /// ADDING NODE PROPERTIES

        for(unsigned i(0); i < m; i++)
            fichier << m_dg.id(m_dg.source(m_dg.arcFromId(i))) << "->" << m_dg.id(m_dg.target(m_dg.arcFromId(i))) << " ;" << endl;                   /// ADDING EDGE PROPERTIES

        fichier << "}" << endl;
        fichier.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
    //for ()
    //cout << " [index=\"" << g_str.NIM[v] << "\", degree=\"" << g_str.NDM[e] << "\"]";

    //cout << " [index=\"" << g_str.EIM[e] << "\", weight=\"" << g_str.EWM[e] << "\", capacity=\"" << g_str.ECM[e] << "\"]";
}


void digraph::write_graph()
{
    digraphWriter(m_dg, cout)
    .nodeMap("Index", m_DNIM)
    .nodeMap("Degree", m_DNDM)
    .nodeMap("Value", m_DNVM)
    .arcMap("Index", m_AIM)
    .arcMap("Capacity", m_ACM)
    .arcMap("Weight", m_AWM)
    .arcMap("Value", m_AVM)
    //.arcMap("lambda", dg_str.ALM)
    .run();
}


void digraph::write_graph(string s)
{
    ofstream filename(s.c_str());

    digraphWriter(m_dg, filename)
    .nodeMap("Index", m_DNIM)
    .nodeMap("Degree", m_DNDM)
    .nodeMap("Value", m_DNVM)
    .arcMap("Index", m_AIM)
    .arcMap("Capacity", m_ACM)
    .arcMap("Weight", m_AWM)
    .arcMap("Value", m_AVM)
    .run();
}


digraph::~digraph()
{
    //dtor
}
