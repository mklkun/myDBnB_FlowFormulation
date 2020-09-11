#ifndef SOLUTION_H
#define SOLUTION_H

#include <map>

#include "graph.h"
#include "digraph.h"


using namespace std;
using namespace lemon;


class solution
{
    public:
        vector < bool > m_xe;                   // variables x(e)
        vector < vector < bool > >  m_fa;       // variables f(a)
        unsigned m_Z;                           // objective function value
        unsigned m_StateFlag;                   // feasability state flag 0:not feasable 1:feasable 2:not verified

        solution();
        solution(vector < bool > xe, vector < vector < bool > > fa, unsigned Z, int StateFlag) : m_xe(xe), m_Z(Z), m_StateFlag(StateFlag) { for (unsigned i=0; i < fa.size(); i++)  m_fa.push_back(vector < bool > (fa[i])); };
        solution(solution const& inSol) : m_xe(inSol.m_xe), m_Z(inSol.m_Z), m_StateFlag(inSol.m_StateFlag) { for (unsigned i=0; i < inSol.m_fa.size(); i++)  m_fa.push_back(vector < bool > (inSol.m_fa[i])); };
        solution(graph &G, vector < digraph* > const& DGS, lemonEdgeValMap const& EVM, vector<lemonArcValMap*> const& AVM, unsigned const Z, unsigned const StateFlag);
        virtual ~solution();
    protected:
    private:
};

#endif // GRAPH_H
