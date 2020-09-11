#include "solution.h"

#include <lemon/lgf_writer.h>

solution::solution()
{
    //ctor
}

solution::solution(graph &G, vector < digraph* > const& DGS, lemonEdgeValMap const& EVM, vector<lemonArcValMap*> const& AVM, unsigned const Z, unsigned const StateFlag)
{
    for(unsigned i(0); i < G.m; i++)
        if (EVM[G.edgeFromId(i)] == 1)
            m_xe.push_back(true);
        else
            m_xe.push_back(false);

    for(unsigned i(0); i < DGS.size(); i++)
    {
        m_fa.push_back(vector < bool > ());
        for(unsigned j(0); j < DGS[i]->m; j++)
            if ((*AVM[i])[DGS[i]->arcFromId(j)] == 1)
                m_fa[i].push_back(true);
            else
                m_fa[i].push_back(false);
    }

    m_Z = Z;
    m_StateFlag = StateFlag;
}

solution::~solution()
{
    //dtor
}

