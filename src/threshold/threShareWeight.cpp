#ifndef THRESHAREWEIGHT_CPP
#define THRESHAREWEIGHT_CPP

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <bitset>
#include <set>
#include <map>
#include <queue>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>

#include "threShareWeight.h"

#include "threshold.h"
#include "base/abc/abc.h"
#include "misc/vec/vec.h"
#include "bdd/extrab/extraBdd.h"
#include "bool/kit/kit.h"
#include "test.c"
#include "../../lib/lp_solve_5.5.2.5_dev/lp_lib.h"

using namespace std;

extern bool bdd2lp(DdManager * dd, DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& minus);
extern "C" DdNode * Th2Bdd(Thre_S* f, DdManager * dd);



#define PTIME
#define LPCHECK
//#define GREEDY

//////////////////////////////////////////////
//////////////////// Node ////////////////////
//////////////////////////////////////////////

unsigned long Node::set(int i) 
{
   bitset<VAR_NUM> b = _bits;
   b.set(i);
   return b.to_ulong();
}

unsigned long Node::flip(int i) 
{
   bitset<VAR_NUM> b = _bits;
   b.flip(i);
   return b.to_ulong();
}

void Node::add_edge(Edge* e, bool in)
{
   if (in)
   {
      _in.push_back(e);
      ++_inSize;
   }
   else
   {
      _out.push_back(e);
      ++_outSize;
   }
}

void Node::add_temp_edge(Edge* e, bool in)
{
   if (in)
      _in.push_back(e);
   else
      _out.push_back(e);
}

void Node::remove_temp_edges()
{
   _in.resize(_inSize);
   _out.resize(_outSize);
}

void Node::topoStep(vector< vector<Node*> >& topo, int& fromStep, int& fromJ)
{
   fromStep = -1;
   for (int s=0; s<topo.size() && fromStep==-1; ++s)
   {
      for (int j=0; j<topo[s].size(); ++j)
      {
         if (topo[s][j] == this)
         {
            fromStep = s;
            fromJ = j;
            break;
         }
      }
   }
}

string Node::str() 
{
   /*char str[] = "";
   int k=0;
   for (int i=0; i<VAR_NUM; ++i)
   {
      if (_bits[i]) 
      {
         str[k] = i+1+'0';
         str[k+1] = '+';
         k += 2;
      }
   }
   if (k!=0) str[k-1]=' ';
   return str; */
   return _bits.to_string(); 
}

//////////////////////////////////////////////
//////////////////// Edge ////////////////////
//////////////////////////////////////////////

Edge::Edge(Lattice* lat, Node* from, Node* to, Type type)
{
   //_lat = lat;
   _from = from;
   _to = to;
   _type = type;

   if (type != REQUIRED && type!=R_IMPLY)
   {
      from->add_edge(this, false);
      to->add_edge(this, true);
   }

   if (type==ORDER || (type==IMPLY && from->lvl()<VAR_NUM-1))
   {
      // add IMPLY edge
      for (int i=0; i<VAR_NUM; ++i)
      {
         if (from->at(i) || to->at(i)) continue;
         lat->add_implied_edge(lat->at(from->set(i)), lat->at(to->set(i)));
      }
   }
}

/////////////////////////////////////////////////
//////////////////// Lattice ////////////////////
/////////////////////////////////////////////////

Lattice::Lattice()
{
   vector< vector< Node* > > _nodes; // temporary
   // create all nodes
   _nodes.resize(VAR_NUM+1);
   for (unsigned long t=0; t<pow(2,VAR_NUM); ++t)
   {
      Node* n = new Node(t);
      _nodes[n->lvl()].push_back(n);
      _nodeList.push_back(n);
   }

   // add LATTICE edges
   for (int i=0; i<VAR_NUM; ++i)
      for (int j=0; j<_nodes[i].size(); ++j)
         for (int k=0; k<_nodes[i+1].size(); ++k)
            if (_nodes[i][j]->cover(_nodes[i+1][k]))
               _edges.insert(new Edge(this, _nodes[i][j], _nodes[i+1][k], LATTICE));

   // add ORDER edges
   for (int i=0; i<_nodes[1].size()-1; ++i) 
      _edges.insert(new Edge(this, _nodes[1][i], _nodes[1][i+1], ORDER)); // 001-->010-->100

   //for (set<Edge*, classcomp>::iterator ite=_edges.begin(); ite!=_edges.end(); ite++)
   //   cout<<(*ite)->from()->str()<<"-->"<<(*ite)->to()->str()<<endl;

   // build topological order
   set<size_t> sorted;
   set<size_t> * processing = new set<size_t>;
   set<size_t> * next = new set<size_t>;

   vector<Node*> tmp;
   next->insert((size_t)(_nodeList[0]));

   while(next->size()>0)
   {
      /*cout<<"next (processing): ";
      for (set<size_t>::iterator ite=next->begin(); ite!=next->end(); ite++)
         cout<< ((Node*)(*ite))->str() << " ";
      cout<<endl;

      cout<<"sorted: ";
      for (set<size_t>::iterator ite=sorted.begin(); ite!=sorted.end(); ite++)
         cout<<(*ite)<<" ";
      cout<<endl;*/

      delete processing;
      processing = next;
      next = new set<size_t>;

      _topological.push_back(tmp);
      vector<Node*>& thisLevel = _topological.back();

      for (set<size_t>::iterator ite=processing->begin(); ite!=processing->end(); ite++)
      {
         Node* n = (Node*)(*ite);

         bool noIn = true;
         for (int j=0; j<n->in().size(); ++j)
         {
            if (sorted.find((size_t)(n->in()[j]->from())) == sorted.end()) // one of its fanin is not sorted yet
            {
               noIn = false;
               break;
            }
         }

         if (noIn) // no fanin
         {
            thisLevel.push_back(n);
            for (int j=0; j<n->out().size(); ++j)
               next->insert((size_t)(n->out()[j]->to()));
         }
         else
            next->insert((size_t)n);
      }
      for (int i=0; i<thisLevel.size(); ++i)
         sorted.insert((size_t)(thisLevel[i]));
   }

   delete processing;
   delete next;

   /*for (int i=0; i<_topological.size(); ++i)
   {
      cout<<i<<" :";
      for (int j=0; j<_topological[i].size(); ++j)
         cout<<_topological[i][j]->str()<<" ";
      cout<<endl;
   }*/
}

void Lattice::add_implied_edge(Node* from, Node* to)
{
   Edge* e = new Edge(this, from, to);
   if (_edges.find(e) == _edges.end())
      _edges.insert(new Edge(this, from, to, IMPLY));
   delete e;
}

//////////////////////////////////////////////////
//////////////////// Function ////////////////////
//////////////////////////////////////////////////

Function::Function(Thre_S* g) 
{ 
   assert(Vec_IntSize(g->weights) == VAR_NUM);
   _gate = g;
   _threshold = g->thre;
   int Entry, i;
   Vec_IntForEachEntry(g->weights, Entry, i)
   {
      if (Entry<0)
      {
         _weights.push_back(-1*Entry);
         _threshold -= Entry;
      }
      else
         _weights.push_back(Entry);
   }
   sort(_weights.begin(), _weights.end());
   //reverse(_weights.begin(), _weights.end());

   //_shouldCheck = true;
}

bool Function::eval(unsigned long t)
{
   bitset<VAR_NUM> b(t);
   int sum = 0;
   for (int i=0; i<_weights.size() && i<b.size(); ++i)
      sum += b[i]*_weights[i];
   return (sum >= _threshold);
}

void Function::findRequiredEdges(Lattice* lat)
{
   map<size_t, bool> values;
   for (unsigned long t=0; t<pow(2,VAR_NUM); ++t)
      values.insert(pair<size_t,bool>((size_t)(lat->at(t)), eval(t)));
   
   vector<Node*> frontier0, frontier1;
   int frontier0Step=-1, frontier1Step=-1;
   bool findingStart = true;
   for (int i=0; i<lat->topoSize(); ++i)
   {
      vector<Node*>& v = lat->topoStep(i);
      bool allTheSame = true;
      for (int j=0; j<v.size(); ++j)
      {
         if (findingStart == values[(size_t)(v[j])]) //findingStart: looking for "not all 0" step; !findingStart: looking for "all 1" step
         {
            allTheSame = false;
            break;
         }
      }
      if (findingStart && !allTheSame)
      {
         frontier0Step = (i==0)? 0: i-1;
         findingStart = false;
         --i;
      }
      else if (!findingStart && allTheSame)
      {
         frontier1Step = (i==lat->topoSize()-1)? i: i+1; // seems don't need +1 ?
         break;
      }
   }

   if (frontier0Step==-1 || frontier1Step==-1) // const 0, or something wrong
   {
      cout<<"frontier not found, const 0 or something wrong "<<frontier0Step<<" "<<frontier1Step<<endl;
      printGate(_gate);
      cout<<endl;
      return;
   }

   for (int i=frontier0Step; i<=frontier1Step; ++i)
   {
      vector<Node*>& v = lat->topoStep(i);
      for (int j=0; j<v.size(); ++j)
      {
         bool value = values[(size_t)(v[j])];
         bool allDifferent = true; // all different from me
         if (!value) // finding frontier0: all out-edge are to 1
         {
            for (int k=0; k<v[j]->out().size(); ++k)
            {
               if (!values[(size_t)(v[j]->out()[k]->to())])
               {
                  allDifferent = false;
                  break;
               }
            }
            if (allDifferent)
               frontier0.push_back(v[j]);
         }
         else // finding frontier1: all in-edge are from 0
         {
            for (int k=0; k<v[j]->in().size(); ++k)
            {
               if (values[(size_t)(v[j]->in()[k]->from())])
               {
                  allDifferent = false;
                  break;
               }
            }
            if (allDifferent)
               frontier1.push_back(v[j]);
         }
         //cout<<"checking: "<<v[j]->str()<<" "<<value<<" "<<allDifferent<<endl;
      }
   }

   for (int i=0; i<frontier0.size(); ++i)
   {
      for (int j=0; j<frontier1.size(); ++j)
      {
         // find whether path frontier0[i]-->frontier1[j] exists
         // assume no need for reachability checking here, just check for single edge
         bool found=false;
         for (int k=0; k<frontier0[i]->out().size(); ++k)
         {
            if (frontier0[i]->out()[k]->to()==frontier1[j])
            {
               found=true;
               break;
            }
         }

         // if not, add this edge
         if (!found)
         {
            _required.push_back(new Edge(lat, frontier0[i], frontier1[j], REQUIRED));

            // find its imply edges
            findRImply(lat, frontier0[i], frontier1[j]);
         }
      }
   }

   cout<<"Finished finding required edges of this gate:"<<endl;
   printGate(_gate);
   cout<<"the required edges:"<<endl;
   for (int i=0; i<_required.size(); ++i)
      cout<<_required[i]->str()<<" ";
   cout<<endl<<endl;
}

void Function::findRImply(Lattice* lat, Node* from, Node* to)
{
   for (int i=0; i<VAR_NUM; ++i)
   {
      if (from->at(i) ^ to->at(i)) continue;
      else
      {
         Node* f = lat->at(from->flip(i));
         Node* t = lat->at(to->flip(i));
         Edge* e = new Edge(lat, f, t);
         if (*(lat->edges().find(e)) == *(lat->edges().end()))
         {
            bool added = false;
            for (int j=0; !added && j<_required.size(); ++j)
               if (*_required[j]==*e)
                  added = true;
            for (int j=0; !added && j<_R_imply.size(); ++j)
               if (*_R_imply[j]==*e)
                  added = true;
            if (!added)
            {
               _R_imply.push_back(new Edge(lat, f, t, R_IMPLY));
               findRImply(lat, f, t);
            } 
         }
         delete e;
      }
   }
}

void Function::checkAllEdges(vector< vector<Node*> >& topo, Function* f, Graph* g, bool secondCheck, bool cleanUp)
{
   if (cleanUp)
   {
      vector< vector<Node*> > topoCopy = topo;
      checkAndAddAllEdges(topoCopy, f, g, secondCheck);
      f->removeAllEdges(topo);
      f->addAllEdges(topo);
   }
   else
      checkAndAddAllEdges(topo, f, g, secondCheck);
}

void Function::checkAndAddAllEdges(vector< vector<Node*> >& topo, Function* f, Graph* g, bool secondCheck)
{
   /*cout<<"f1: "<<f->gate()->Id<<" f2: "<<_gate->Id<<endl;
   for (int i=0; i<topo.size(); ++i)
   {
      cout<<i<<" :";
      for (int j=0; j<topo[i].size(); ++j)
         cout<<topo[i][j]->str()<<" ";
      cout<<endl;
   }cout<<endl;*/

   bool noConflict = true, implication = true;
   for (int i=0; i<_required.size()+_R_imply.size(); ++i)
   {
      Edge* e = (i<_required.size())? _required[i]: _R_imply[i - _required.size()];
      int fromStep = -1, fromJ, toStep = -1, toJ;

      e->from()->topoStep(topo, fromStep, fromJ);
      e->to()->topoStep(topo, toStep, toJ);

      if (fromStep == toStep)
      {  
         implication = false;
         e->from()->add_temp_edge(e, false);
         e->to()->add_temp_edge(e, true);
         reorder(topo, toStep, toJ, fromStep+1);
      }
      else if (fromStep > toStep)
      {
         implication = false;
         // check for loop
         if (findPath(e->to(), e->from()))
         {
            noConflict = false;
            break;
         }
         e->from()->add_temp_edge(e, false);
         e->to()->add_temp_edge(e, true);
         reorder(topo, toStep, toJ, fromStep+1);
      }
      else 
      {
         // check implication
         if (!findPath(e->from(), e->to()))
            implication = false;
         e->from()->add_temp_edge(e, false);
         e->to()->add_temp_edge(e, true);
      }
   }
   if (!noConflict) g->setConflict(this, f);
   else if (implication) g->setImply(f, this); // f implies this
   else if (!secondCheck) g->setCompatible(this, f);
}

/*bool Function::findPath(Node* A, int Astep, Node* B, int Bstep, vector< vector<Node*> >& topo, set<int>& knobs)
{
   set<int>::iterator ite;
   ite = knobs.lower_bound(Astep);
   if (ite == knobs.end() || *ite >= Bstep)
      return findPath(A, B, topo, knobs);
   else
   {
      bool path1;
      if (*ite == Astep) path1 = true;
      else path1 = findPath(A, topo[*ite][0], topo, knobs);
      ite = knobs.upper_bound(Bstep);

   }
}*/

bool Function::findPath(Node* A, Node* B/*, vector< vector<Node*> >& topo, set<int>& knobs*/)
{
   queue<Node*> Q;
   set<Node*> visited;
   Q.push(A);
   visited.insert(A);
   while(!Q.empty())
   {
      Node* q = Q.front();
      Q.pop();
      for (int i=0; i<q->out().size(); ++i)
      {
         Node* n = q->out()[i]->to();
         if (n==B)
            return true;
         if (visited.find(n) == visited.end())
         {
            Q.push(n);
            visited.insert(n);
         }
      }
   }
   return false;
}

void Function::addAllEdges(vector< vector<Node*> >& topo)//, set<int>& knobs)
{
   for (int i=0; i<_required.size()+_R_imply.size(); ++i)
   {
      Edge* e = (i<_required.size())? _required[i]: _R_imply[i - _required.size()];
      int fromStep = -1, fromJ, toStep = -1, toJ;

      e->from()->topoStep(topo, fromStep, fromJ);
      e->to()->topoStep(topo, toStep, toJ);

      e->from()->add_temp_edge(e, false);
      e->to()->add_temp_edge(e, true);
      if (fromStep >= toStep)
         reorder(topo, toStep, toJ, fromStep+1);
   }
   /*for (int i=0; i<topo.size(); ++i)
      if (topo[i].size()==1)
         knobs.insert(i);*/
}

void Function::removeAllEdges(vector< vector<Node*> >& topo)
{
   for (int i=0; i<topo.size(); ++i)
      for (int j=0; j<topo[i].size(); ++j)
         topo[i][j]->remove_temp_edges();
}

void Function::reorder(vector< vector<Node*> >& topo, int s, int j, int targetStep)
{
   Node* n = topo[s][j];
   if (targetStep==topo.size()) 
   {
      vector<Node*> v;
      topo.push_back(v);
   }
   topo[targetStep].push_back(n);
   topo[s].erase(topo[s].begin()+j);
   for (int i=0; i<n->out().size(); ++i)
   {
      int s2, j2;
      n->out()[i]->to()->topoStep(topo, s2, j2);
      if (s2 <= targetStep)
         reorder(topo, s2, j2, targetStep+1);
   }
}

void Function::updateNeighborD(vector< vector<Function*> >& D, set<Function*>& _V)
{
   //cout<<"updating "<<gate()->Id<<"'s neighbor: "<<endl;
   for (set<Function*>::iterator ite=_adj.begin(); ite!=_adj.end(); ite++)
   {
      if (_V.find(*ite)==_V.end()) continue;
      //cout<<(*ite)->gate()->Id<<endl;
      for (int i=0; i<D[(*ite)->d()].size(); ++i)
         if (D[(*ite)->d()][i]==(*ite))
         {
            D[(*ite)->d()].erase(D[(*ite)->d()].begin()+i);
            D[(*ite)->d()-1].push_back(*ite);
            break;
         }
      (*ite)->Dminus();
   }
   /*for (int i=0; i<D.size(); ++i)
   {
      cout<<"D["<<i<<"]: ";
      for (int j=0; j<D[i].size(); ++j)
         cout<<D[i][j]->gate()->Id<<" ";
      cout<<endl;
   }*/
}


///////////////////////////////////////////////
//////////////////// Graph ////////////////////
///////////////////////////////////////////////

void Graph::setConflict(Function* f1, Function* f2)
{
   cout<<"These two gates conflict: "<<f1->gate()->Id<<" "<<f2->gate()->Id<<endl;
   //printGate(f1->gate());
   //printGate(f2->gate());
   //cout<<endl;
}

void Graph::setImply(Function* f1, Function* f2)
{
   cout<<"The first gate implies the second: "<<f1->gate()->Id<<" "<<f2->gate()->Id<<endl;
   //printGate(f1->gate());
   //printGate(f2->gate());
   //cout<<endl;

   f1->setImply(f2);
   _V.erase(f2);
}

void Graph::setCompatible(Function* f1, Function* f2)
{
   cout<<"These two gates are compatible: "<<f1->gate()->Id<<" "<<f2->gate()->Id<<endl;
   //printGate(f1->gate());
   //printGate(f2->gate());
   //cout<<endl;

   f1->setAdj(f2);
   f2->setAdj(f1);
}

void Graph::findCliques()
{
   #ifdef PTIME
      abctime clk = Abc_Clock();
   #endif
   vector<Function*> R, P, X;
   /*for (set<Function*>::iterator ite=_V.begin(); ite!=_V.end(); ite++)
      P.push_back(*ite);*/
   // making P in degeneracy order
   vector< vector<Function*> > D;
   for (set<Function*>::iterator ite=_V.begin(); ite!=_V.end(); ite++)
   {
      (*ite)->initD();
      if ((*ite)->d() >= D.size()) D.resize((*ite)->d()+1);
      D[(*ite)->d()].push_back(*ite);
   }
   /*for (int i=0; i<D.size(); ++i)
   {
      cout<<"D["<<i<<"]: ";
      for (int j=0; j<D[i].size(); ++j)
         cout<<D[i][j]->gate()->Id<<" ";
      cout<<endl;
   }*/

   for (int i=0; i<_V.size(); ++i)
   {
      int j=0;
      for (j=0; j<D.size(); ++j)
         if (D[j].size()>0)
            break;
      P.push_back(D[j].back()); // reverse degeneracy order
      Function* f = D[j].back();
      D[j].pop_back();
      f->updateNeighborD(D, _V);
   }

   #ifdef GREEDY
   _Vcopy=_V;
   while (_Vcopy.size()>2)
   {
   #endif

   while (P.size()>0)
   {
      int i=P.size()-1;
      vector<Function*> R1, P1, X1;
      R1 = R;
      R1.push_back(P[i]);
      for (int j=0; j<P.size(); ++j)
         if (P[j]->compatible(P[i]))
            P1.push_back(P[j]);
      for (int j=0; j<X.size(); ++j)
         if (X[j]->compatible(P[i]))
            X1.push_back(X[j]);

      /*#ifdef GREEDY
      if (findCliques_rec(R1, P1, X1))
      {
         vector<Function*>& c = _cliques.back();
         for (int j=0; j<c.size(); ++j)
         {
            X.push_back(c[j]);
            for (int k=0; k<P.size(); ++k)
            {
               if (P[k]==c[j])
               {
                  P.erase(P.begin()+k);
                  break;
               }
            }
         }
      }
      cout<<P.size()<<endl;
      #else*/
      findCliques_rec(R1, P1, X1);
      X.push_back(P[i]);
      P.pop_back();
      //#endif
   }

   #ifdef GREEDY
   cout<<_Vcopy.size()<<endl;
   for (set<Function*>::iterator it=_Vcopy.begin(); it!=_Vcopy.end(); ++it)
      P.push_back(*it);
   R = vector<Function*>();
   X = vector<Function*>();
   }
   #endif

   #ifdef PTIME
      Abc_PrintTime(1, "find cliques: ", Abc_Clock()-clk);
      clk = Abc_Clock();
   #endif

   for (int i=0; i<_cliques.size(); ++i)
   {
      vector<int> ans;
      #ifndef LPCHECK
      if (solveLp(_cliques[i], ans)==0)
      #else
      if (solveLpDd(_cliques[i], ans)==0)
      #endif 
         _implementations.push_back(ans);
      else
      {
         cout<<"clique solve failed: ";
         for (int j=0; j<_cliques[i].size(); ++j)
            cout<<_cliques[i][j]->gate()->Id<<" ";
         cout<<endl;

         vector<Function*> clique = _cliques[i];
         _cliques.erase(_cliques.begin()+i);
         --i;
         FindSubClique(clique, i);
      }
   }
   assert(_cliques.size()==_implementations.size());

   #ifdef PTIME
      Abc_PrintTime(1, "solve for implementation: ", Abc_Clock()-clk);
   #endif
}

void Graph::FindSubClique(vector<Function*>& clique, int& i)
{
   // keep i to be the last one that has been processed (i+1 is the next new one)
   
   if (clique.size()<2) return;
   for (int j=0; j<clique.size(); ++j)
   {
      vector<Function*> c = clique;
      c.erase(c.begin()+j);
      vector<int> ans;
      #ifndef LPCHECK
      if (solveLp(c, ans)==0)
      #else
      if (solveLpDd(c, ans)==0)
      #endif 
      {
         _implementations.push_back(ans);
         _cliques.insert(_cliques.begin()+i+1, c);
         ++i;
      }
      else
      {
         // keep recursive?
      }
   }

}

bool Graph::findCliques_rec(vector<Function*>& R, vector<Function*>& P, vector<Function*>& X)
{
   /*cout<<"R: "; for (int i=0; i<R.size(); ++i) cout<<R[i]->gate()->Id<<" "; cout<<endl;
   cout<<"P: "; for (int i=0; i<P.size(); ++i) cout<<P[i]->gate()->Id<<" "; cout<<endl;
   cout<<"X: "; for (int i=0; i<X.size(); ++i) cout<<X[i]->gate()->Id<<" "; cout<<endl;
   cout<<endl;*/

   if (P.size()==0 && X.size()==0)
   {
      if (R.size()>1)
      {
         _cliques.push_back(R);
         #ifdef GREEDY
         for (int i=0; i<R.size(); ++i)
            _Vcopy.erase(R[i]);
         #endif
         return true;
      }
      //for (int i=0; i<R.size(); ++i) cout<<R[i]->gate()->Id<<" "; cout<<endl;
      return false;
   }

   Function* u = (P.size()>0)? P[0]: X[0];
   for (int i=0; i<P.size(); ++i)
   {
      if (P[i]->compatible(u)) continue;

      vector<Function*> R1, P1, X1;
      R1 = R;
      R1.push_back(P[i]);
      for (int j=0; j<P.size(); ++j)
         if (P[j]->compatible(P[i]))
            P1.push_back(P[j]);
      for (int j=0; j<X.size(); ++j)
         if (X[j]->compatible(P[i]))
            X1.push_back(X[j]);

      #ifdef GREEDY
      if (findCliques_rec(R1, P1, X1))
      {
         //vector<Function*>& c = _cliques.back();
         return true;
         /*for (int j=0; j<c.size(); ++j)
         {
            X.push_back(c[j]);
            for (int k=0; k<P.size(); ++k)
            {
               if (P[k]==c[j])
               {
                  if (k<=i) --i;
                  P.erase(P.begin()+k);
                  break;
               }
            }
         }
         cout<<P.size()<<endl;
         if (P.size()==0) return true;*/
      }
      #else
      findCliques_rec(R1, P1, X1);
      X.push_back(P[i]);
      P.erase(P.begin()+i);
      --i;
      #endif
   }
   return false;
}

void Graph::reportCliques()
{
   set<Function*> Vcopy = _V;

   cout<<"Cliques found:"<<endl;
   for (int i=0; i<_cliques.size(); ++i)
   {
      cout<<"clique "<<i<<": ";
      for (int j=0; j<_cliques[i].size(); ++j)
      {
         cout<<_cliques[i][j]->gate()->Id<<" ";
         //printGate(_cliques[i][j]->gate());
         Vcopy.erase(_cliques[i][j]);
      }
      cout<<endl;
      cout<<"implementation: ";
      for (int j=0; j<_implementations[i].size(); ++j)
         cout<<_implementations[i][j]<<", ";
      cout<<endl<<endl;
   }

   for (set<Function*>::iterator ite=Vcopy.begin(); ite!=Vcopy.end(); ite++)
   {
      vector<Function*> vec;
      vec.push_back(*ite);
      vector<int> ans;
      #ifndef LPCHECK
      if (solveLp(vec, ans)==0)
      #else
      if (solveLpDd(vec, ans)==0)
      #endif 
      {
         _cliques.push_back(vec);
         _implementations.push_back(ans);
         cout<<"single gate "<<(*ite)->gate()->Id<<" ";
         cout<<"implementation: ";
         for (int j=0; j<ans.size(); ++j)
            cout<<ans[j]<<", ";
         cout<<endl;
      }
      else
         cout<<"single gate "<<(*ite)->gate()->Id<<"solve failed"<<endl;
   }
}

void Graph::printV()
{
   cout<<endl<<"_V: ";
   for (set<Function*>::iterator ite=_V.begin(); ite!=_V.end(); ite++)
      cout<<(*ite)->gate()->Id<<" ";
   cout<<endl<<endl;
}

void Graph::minCoverHeuristic()
{
   while (!_V.empty())
   {
      int max = 0, idx = 0;
      for (int i=0; i<_cliques.size(); ++i)
      {
         int c = 0;
         for (int j=0; j<_cliques[i].size(); ++j)
            if (_V.find(_cliques[i][j])!=_V.end())
               ++c;
         if (c > max)
         {
            max = c;
            idx = i;
         }
      }
      cout<<"chosen prime: "<<idx<<" [";
      for (int j=0; j<VAR_NUM; ++j)
         cout<<_implementations[idx][j]<<", ";
      cout<<"]"<<endl;

      for (int j=0; j<_cliques[idx].size(); ++j)
         _V.erase(_cliques[idx][j]);
      vector<Function*> vec;
      _cliques[idx] = vec;

      #if 0
      cout<<"size of V not covered yet: "<<_V.size()<<endl;
      #endif
   }
}

int Graph::solveLp(vector<Function*>& funcs, vector<int>& ans)
{
   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0;
   REAL *row = NULL;
   
   Ncol = VAR_NUM;
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=0; i<Ncol; ++i)
         set_int(lp,i+1,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
   }
   
   int k=0;
   while ((ret == 0)&&(k<funcs.size())) 
   {
      int k2=0;
      while ((ret==0)&&(k2<funcs[k]->required().size()))
      {
         Edge* e = funcs[k]->required()[k2];
         j = 0;
         for (int i=0; i<Ncol; ++i)
         {
            if (e->to()->at(i) && e->from()->at(i)) continue;
            else if (e->to()->at(i))
            {
               colno[j]=i+1;
               row[j++]=1;
            }
            else if (e->from()->at(i))
            {
               colno[j]=i+1;
               row[j++]=-1;
            }
         }
         
         /* add the row to lpsolve */
         if(!add_constraintex(lp, j, row, colno, GE, 1))
            ret = 3;
         k2++;
      }
      k++;
   }

   k=0;
   while ((ret==0)&&(k<Ncol-1))
   {
      //k+1 < k+2
      colno[0]=k+1;
      row[0]=-1;
      colno[1]=k+2;
      row[1]=1;
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, 2, row, colno, GE, 1))
         ret = 3;
      k++;
   }

   if (ret==0)
   {
      // the first (smallest) var >=1 
      colno[0]=1;
      row[0]=1;
      if(!add_constraintex(lp, 1, row, colno, GE, 1))
         ret = 3;
   }

   if(ret == 0) {
      set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
      
      /* set the objective function */
      j = 0;
      
      for (int i=0; i<VAR_NUM; ++i)
      {
         colno[j]=i+1;
         row[j++]=1;
      }
      
      /* set the objective in lpsolve */
      if(!set_obj_fnex(lp, j, row, colno))
         ret = 4;
   }
   
   if(ret == 0) {
      /* set the object direction to minimize */
      set_minim(lp);
      
      /* just out of curioucity, now show the model in lp format on screen */
      //write_LP(lp, stdout);
      
      /* I only want to see important messages on screen while solving */
      set_verbose(lp,NEUTRAL);// IMPORTANT);
      
      /* Now let lpsolve calculate a solution */
      ret = solve(lp);
      if(ret == OPTIMAL)
         ret = 0;
      else
         ret = 5;
   }

   if (ret==0) {
      /* objective value */
      //printf("Objective value: %f\n", get_objective(lp));
      
      /* variable values */
      get_variables(lp, row);
      for(j = 0; j < Ncol; j++)
      {
         //printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
         ans.push_back(row[j]);
         //cout<<row[j]<<" ";
      }
      //cout<<endl;
      /* we are done now */
   }
   
   /* free allocated memory */
   if(row != NULL)
      free(row);
   if(colno != NULL)
      free(colno);
   
   if(lp != NULL) {
      /* clean up such that all used memory by lpsolve is freed */
      delete_lp(lp);
   } 
   return(ret);
}

int Graph::solveLpDd(vector<Function*>& funcs, vector<int>& ans)
{
   DdManager * dd = Cudd_Init( VAR_NUM , 0 , CUDD_UNIQUE_SLOTS , CUDD_CACHE_SLOTS , 0 );

   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0;
   REAL *row = NULL;
   
   Ncol = VAR_NUM+funcs.size();
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=0; i<Ncol; ++i)
         set_int(lp,i+1,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
   } 

   int f = 0;
   while (ret==0 && f<funcs.size())
   {
      vector< vector<int> > A;
      vector<int> vec,B,minus;
      vec.resize(dd->size+1);
      A.push_back(vec); 
      minus.resize(dd->size,0); 
   
      DdNode* bFunc = Th2Bdd(funcs[f]->gate(), dd);
      assert(bdd2lp(dd,bFunc,A,B,minus));
   
      int k=0;
      while ((ret == 0)&&(k<A.size())) {
         j = 0;
         for (int i=0; i<A[k].size(); ++i)
         {
            if (A[k][i]==0) continue;
            colno[j]=(i==VAR_NUM)? i+1+f: i+1; // different T
            row[j++]=A[k][i];
         }
         
         /* add the row to lpsolve */
         if(!add_constraintex(lp, j, row, colno, GE, B[k]))
            ret = 3;
         k++;
      }
      f++;
   }

   if(ret == 0) {
      set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
      
      /* set the objective function */
      j = 0;
      
      for (int i=0; i<Ncol; ++i)
      {
         colno[j]=i+1;
         row[j++]=1;
      }
      
      /* set the objective in lpsolve */
      if(!set_obj_fnex(lp, j, row, colno))
         ret = 4;
   }
   
   if(ret == 0) {
      /* set the object direction to minimize */
      set_minim(lp);
      
      /* just out of curioucity, now show the model in lp format on screen */
      //write_LP(lp, stdout);
      
      /* I only want to see important messages on screen while solving */
      set_verbose(lp,NEUTRAL);// IMPORTANT);
      
      /* Now let lpsolve calculate a solution */
      ret = solve(lp);

      if(ret == OPTIMAL)
         ret = 0;
      else
         ret = 5;
   }

   if (ret==0) {
      // objective value
      //printf("Objective value: %f\n", get_objective(lp));
      
      // the answer
      get_variables(lp, row);
      for(j = 0; j < VAR_NUM; j++)
      {
         //printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
         ans.push_back(row[j]);
         //cout<<row[j]<<" ";
      }
      //cout<<endl;
   }

   
   /* free allocated memory */
   if(row != NULL)
      free(row);
   if(colno != NULL)
      free(colno);
   
   if(lp != NULL) {
      /* clean up such that all used memory by lpsolve is freed */
      delete_lp(lp);
   }

   Cudd_Quit(dd);
   return ret;
}

void Graph::checkCompatibilityLP(Function* f1, Function* f2)
{
   DdManager * dd = Cudd_Init( VAR_NUM , 0 , CUDD_UNIQUE_SLOTS , CUDD_CACHE_SLOTS , 0 );
   DdNode* bFunc1 = Th2Bdd(f1->gate(), dd);
   DdNode* bFunc2 = Th2Bdd(f2->gate(), dd);

   vector< vector<int> > A1,A2;
   vector<int> vec,B1,B2,minus;
   vec.resize(dd->size+1);
   A1.push_back(vec); 
   A2.push_back(vec); 
   minus.resize(dd->size,0); 

   assert(bdd2lp(dd,bFunc1,A1,B1,minus));
   assert(bdd2lp(dd,bFunc2,A2,B2,minus));

   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0;
   REAL *row = NULL;
   
   Ncol = VAR_NUM+2;
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=0; i<Ncol; ++i)
         set_int(lp,i+1,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
   }

   int k=0;
   while ((ret == 0)&&(k<A1.size())) {
      j = 0;
      for (int i=0; i<A1[k].size(); ++i)
      {
         if (A1[k][i]==0) continue;
         colno[j]=i+1;
         row[j++]=A1[k][i];
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, GE, B1[k]))
         ret = 3;
      k++;
   }
   k=0;
   while ((ret == 0)&&(k<A2.size())) {
      j = 0;
      for (int i=0; i<A2[k].size(); ++i)
      {
         if (A2[k][i]==0) continue;
         colno[j]=(i==VAR_NUM)? i+2: i+1; // the second T
         row[j++]=A2[k][i];
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, GE, B2[k]))
         ret = 3;
      k++;
   }

   k=0;
   if(ret == 0) {
      set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
      
      /* set the objective function */
      j = 0;
      
      for (int i=0; i<A1[0].size(); ++i)
      {
         colno[j]=i+1;
         row[j++]=1;
      }
      
      /* set the objective in lpsolve */
      if(!set_obj_fnex(lp, j, row, colno))
         ret = 4;
   }
   
   if(ret == 0) {
      /* set the object direction to minimize */
      set_minim(lp);
      
      /* just out of curioucity, now show the model in lp format on screen */
      //write_LP(lp, stdout);
      
      /* I only want to see important messages on screen while solving */
      set_verbose(lp,NEUTRAL);// IMPORTANT);
      
      /* Now let lpsolve calculate a solution */
      ret = solve(lp);

      if(ret == OPTIMAL)
         ret = 0;
      else
         ret = 5;
   }

   if (ret==0) {
      
      setCompatible(f1, f2);
      // objective value
      //printf("Objective value: %f\n", get_objective(lp));
      
      // the answer
      //get_variables(lp, row);

   }
   else
      setConflict(f1, f2);
   
   /* free allocated memory */
   if(row != NULL)
      free(row);
   if(colno != NULL)
      free(colno);
   
   if(lp != NULL) {
      /* clean up such that all used memory by lpsolve is freed */
      delete_lp(lp);
   }
   

   Cudd_Quit(dd);
}

extern "C" Vec_Ptr_t* Th_ShareWeightSyn(Vec_Ptr_t* vFuncs)
{
   return SharingWeightSynthesis(vFuncs);
}

Vec_Ptr_t* SharingWeightSynthesis(Vec_Ptr_t* vFuncs)
{
   Vec_Ptr_t* vRes = Vec_PtrAlloc(Vec_PtrSize(vFuncs));

   #ifdef PTIME
      abctime clk = Abc_Clock();
   #endif

   #ifndef LPCHECK
   Lattice* lat = new Lattice();
   
   #ifdef PTIME
      Abc_PrintTime(1, "build basic lattice time: ", Abc_Clock()-clk);
      clk = Abc_Clock();
   #endif
   #endif //not LPCHECK

   vector<Function*> functions;
   Thre_S* pFunc;
   int i;
   Vec_PtrForEachEntry(Thre_S*, vFuncs, pFunc, i)
   {
      functions.push_back(new Function(pFunc));
      #ifndef LPCHECK
      functions.back()->findRequiredEdges(lat);
      #endif //not LPCHECK
   }
   Graph* g = new Graph(functions);

   #ifdef PTIME
      Abc_PrintTime(1, "process all functions (find required edges): ", Abc_Clock()-clk);
      clk = Abc_Clock();
   #endif


   /////////////////////////////////////////////////////
   /////////////// compatibility check /////////////////
   /////////////////////////////////////////////////////
   #ifndef LPCHECK
   for (int i=0; i<functions.size(); ++i)
   {
      Function* f1 = functions[i];
      // (copy lattice) and add f1's R, R_imply
      vector< vector<Node*> > topo = lat->topo(); // copy
      //set<int> knobs;
      f1->addAllEdges(topo);

      for (int j=0 /*i+1*/; j<functions.size(); ++j)
      {
         Function* f2 = functions[j];
         if (j==i) continue;
         else if (j<i)
         {
            if (f2->compatible(f1))
               f2->checkAllEdges(topo, f1, g, true, true); 
         }
         else
            f2->checkAllEdges(topo, f1, g, false, true); 
      }

      f1->removeAllEdges(topo);
   }

   #else
   for (int i=0; i<functions.size()-1; ++i)
   {
      Function* f1 = functions[i];
      for (int j=i+1; j<functions.size(); ++j)
      {
         Function* f2 = functions[j];
         g->checkCompatibilityLP(f1, f2);
      }
   }
   #endif // LPCHECK

   g->printV();
   #ifdef PTIME
      Abc_PrintTime(1, "pair-wise compatibility check: ", Abc_Clock()-clk);
   #endif

   g->findCliques();
   g->reportCliques();

   #ifdef PTIME
      clk = Abc_Clock();
   #endif
   g->minCoverHeuristic();
   #ifdef PTIME
      Abc_PrintTime(1, "min cover: ", Abc_Clock()-clk);
   #endif

   #ifndef LPCHECK
   delete lat;
   #endif

   delete g;
   return vRes;
}

#endif
