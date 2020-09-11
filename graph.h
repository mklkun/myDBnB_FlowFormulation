#ifndef GRAPH_H
#define GRAPH_H

#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>


using namespace std;
using namespace lemon;


enum colors { BLACK, GREY, RED, ORANGE, YELLOW, GREEN, BLUE, WHITE};


    typedef lemon::ListGraph lemonGraph;

    typedef lemonGraph::EdgeIt lemonEdgeIt;
    typedef lemonGraph::NodeIt lemonNodeIt;

    typedef lemonGraph::Node lemonNode;
    typedef lemonGraph::Edge lemonEdge;

    // une Map sert à donner des valeurs aux composants (noeuds, arêtes, ...) d'un graphe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    typedef lemonGraph::NodeMap<int> lemonNodeIndexMap;
    typedef lemonGraph::NodeMap<int> lemonNodeDegreeMap;
    typedef lemonGraph::NodeMap<int> lemonNodeValMap;

    typedef lemonGraph::EdgeMap<int> lemonEdgeIndexMap;
    typedef lemonGraph::EdgeMap<int> lemonEdgeCapacityMap;
    typedef lemonGraph::EdgeMap<double> lemonEdgeWeightMap;
    typedef lemonGraph::EdgeMap<double> lemonEdgeValMap;

    typedef lemonGraph::EdgeMap<float> lemonEdgeLambdaMap;

    using lemon::INVALID;

class graph
{
    public:
        lemonGraph m_g;

        string m_name;

        vector<lemonNode> m_Nodes;
        vector<lemonEdge> m_Edges;

        lemonNodeIndexMap m_NIM;
        lemonNodeDegreeMap m_NDM;
        lemonNodeValMap m_NVM;

        lemonEdgeIndexMap m_EIM;
        lemonEdgeCapacityMap m_ECM;
        lemonEdgeWeightMap m_EWM;
        lemonEdgeValMap m_EVM;

        map <int, int> m_linking;

        unsigned n;
        unsigned m;

        graph(string const name) : m_name(name),m_NIM(m_g,0),m_NDM(m_g,0),m_NVM(m_g,0),m_EIM(m_g,0),m_ECM(m_g,0),m_EWM(m_g,0),m_EVM(m_g,0),n(0),m(0) {};
        virtual ~graph();

        lemonNode nodeFromId(unsigned ind) { return (m_g.nodeFromId(ind)); };
        lemonEdge edgeFromId(unsigned ind) { return (m_g.edgeFromId(ind)); };
        lemonEdge findEdge(lemonNode u, lemonNode v) { return (lemon::findEdge(m_g, u, v)); };
        lemonNode u(lemonEdge e) { return(m_g.u(e)); };
        lemonNode v(lemonEdge e) { return(m_g.v(e)); };


        int idFromNode(lemonNode u) { return (m_NIM[u]); };
        int degreeFromNode(lemonNode u) { return (m_NDM[u]); };
        int valueFromNode(lemonNode u) { return (m_NVM[u]); };

        int idFromEdge(lemonEdge e) { return (m_EIM[e]); };
        int capacityFromEdge(lemonEdge e) { return (m_ECM[e]); };
        double weightFromEdge(lemonEdge e) { return (m_EWM[e]); };
        int valueFromEdge(lemonEdge e) { return (m_EVM[e]); };

        void setId(lemonNode u, int x) { m_NIM[u] = x; };
        void setDegree(lemonNode u, int x) { m_NDM[u] = x; };
        void setValue(lemonNode u, int x) { m_NVM[u] = x; };

        void setId(lemonEdge e, int x) { m_EIM[e] = x; };
        void setCapacity(lemonEdge e, int x) { m_ECM[e] = x; };
        void setWeight(lemonEdge e, double x) { m_EWM[e] = x; };
        void setValue(lemonEdge e, int x) { m_EVM[e] = x; };

        void write_graph_viz(string const s);
        void write_graph_viz_sol(vector < vector <int> > const demands, string const s);
        void write_graph();
        void write_graph(string s);

    protected:
    private:
};

#endif // GRAPH_H
