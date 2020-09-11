#ifndef DIGRAPH_H
#define DIGRAPH_H

#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/lgf_writer.h>


using namespace std;
using namespace lemon;


    typedef lemon::ListDigraph lemonDigraph;

    typedef lemonDigraph::ArcIt lemonArcIt;
    typedef lemonDigraph::NodeIt lemonDiNodeIt;

    typedef lemonDigraph::Node lemonDiNode;
    typedef lemonDigraph::Arc lemonArc;

    // une Map sert à donner des valeurs aux composants (noeuds, arêtes, ...) d'un graphe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    typedef lemonDigraph::NodeMap<unsigned> lemonDiNodeIndexMap;
    typedef lemonDigraph::NodeMap<int> lemonDiNodeDegreeMap;
    typedef lemonDigraph::NodeMap<int> lemonDiNodeValMap;

    typedef lemonDigraph::ArcMap<int> lemonArcIndexMap;
    typedef lemonDigraph::ArcMap<int> lemonArcCapacityMap;
    typedef lemonDigraph::ArcMap<double> lemonArcWeightMap;
    typedef lemonDigraph::ArcMap<int> lemonArcValMap;

    typedef lemonDigraph::ArcMap<float> lemonArcLambdaMap;


    using lemon::INVALID;


class digraph
{
    public:
        lemonDigraph m_dg;

        vector < vector<lemonDiNode> > m_DNodes;
        vector<lemonArc> m_Arcs;

        lemonDiNodeIndexMap m_DNIM;
        lemonDiNodeDegreeMap m_DNDM;
        lemonDiNodeValMap m_DNVM;

        lemonArcIndexMap m_AIM;
        lemonArcCapacityMap m_ACM;
        lemonArcWeightMap m_AWM;
        lemonArcValMap m_AVM;

        map <int, int> m_linking;

        lemonDiNode m_s;
        lemonDiNode m_t;

        unsigned n;
        unsigned m;

        digraph() : m_DNIM(m_dg,0),m_DNDM(m_dg,0),m_DNVM(m_dg,9999999),m_AIM(m_dg,0),m_ACM(m_dg,0),m_AWM(m_dg,0),m_AVM(m_dg,0),n(0),m(0) {};
        virtual ~digraph();

        lemonDiNode nodeFromId(unsigned j) { return (m_dg.nodeFromId(j)); };
        lemonArc arcFromId(unsigned j) { return (m_dg.arcFromId(j)); };
        lemonDiNode source(lemonArc a) { return(m_dg.source(a)); };
        lemonDiNode target(lemonArc a) { return(m_dg.target(a)); };

        int idFromDiNode(lemonDiNode u) { return (m_DNIM[u]); };
        int degreeFromDiNode(lemonDiNode u) { return (m_DNDM[u]); };
        int valueFromDiNode(lemonDiNode u) { return (m_DNVM[u]); };

        int idFromArc(lemonArc a) { return (m_AIM[a]); };
        int capacityFromArc(lemonArc a) { return (m_ACM[a]); };
        double weightFromArc(lemonArc a) { return (m_AWM[a]); };
        int valueFromArc(lemonArc a) { return (m_AVM[a]); };

        void setId(lemonDiNode u, int x) { m_DNIM[u] = x; };
        void setDegree(lemonDiNode u, int x) { m_DNDM[u] = x; };
        void setValue(lemonDiNode u, int x) { m_DNVM[u] = x; };

        void setId(lemonArc a, int x) { m_AIM[a] = x; };
        void setCapacity(lemonArc a, int x) { m_ACM[a] = x; };
        void setWeight(lemonArc a, double x) { m_AWM[a] = x; };
        void setValue(lemonArc a, int x) { m_AVM[a] = x; };

/******************added***********************/
        lemonDiNode nodeFromMyId(unsigned j)
        {
            if(m_DNIM[m_s] == j)
            {
                return m_s;
            }
            if(m_DNIM[m_t] == j)
            {
                return m_t;
            }
            unsigned pos = j % 1000000;
            for (unsigned i(0) ; i < m_DNodes[pos].size(); i++)
                if(m_DNIM[m_DNodes[pos][i]] == j)
                    return m_DNodes[pos][i];
            cout << "WEIRD for j=" << j << " !!!" << endl;
            lemonDiNode err;
            return err;
        };
        lemonArc findArc(lemonDiNode u, lemonDiNode v) { return (lemon::findArc(m_dg, u, v)); };

        int myIdFromDiNode(lemonDiNode u) { return (m_DNIM[u]); };
        int myIdFromArc(lemonArc a) { return (m_AIM[a]); };
/************************************************/

        void write_graph_viz(string const s);
        void write_graph();
        void write_graph(string s);

    protected:
    private:
};

#endif // DIGRAPH_H
