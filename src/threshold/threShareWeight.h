#ifndef THRESHAREWEIGHT_H
#define THRESHAREWEIGHT_H

#include <vector>
#include <bitset>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>

#include "threshold.h"
#include "base/abc/abc.h"
#include "misc/vec/vec.h"

using namespace std;

#define VAR_NUM 5
//#define TRUTH_LEN 16 //2^VAR_NUM


Vec_Ptr_t* SharingWeightSynthesis(Vec_Ptr_t* vFuncs);

enum Type {
   LATTICE = 1, // A < A+x
   ORDER, // w1 < w2 < w3 ...
   IMPLY, // A < B => A+x < B+x

   R_IMPLY, // implied by required edge, including A < B => A-x < B-x
   REQUIRED, // function-dependent
   CHECK
};

class Edge; class Lattice; class Function; class Graph;
struct classcomp;

class Node
{
public:
   Node(unsigned long truth)
   {
      _bits = bitset<VAR_NUM> (truth);
      _inSize = 0;
      _outSize = 0;
   }

   unsigned long truth_position() { return _bits.to_ulong(); }
   int lvl() { return _bits.count(); }
   bool at(int i) { return _bits[i]; }
   unsigned long set(int i);
   unsigned long flip(int i);

   bool cover(Node* n) //check if there's line between two nodes
   { return (_bits^n->_bits).count()<=1; }

   void add_edge(Edge* e, bool in);
   void add_temp_edge(Edge* e, bool in);
   void remove_temp_edges();

   void topoStep(vector< vector<Node*> >& topo, int& s, int& j);

   string str();

   const vector<Edge*>& in() { return _in; }
   const vector<Edge*>& out() { return _out; }

private:
   bitset<VAR_NUM> _bits;
   vector<Edge*> _in;
   vector<Edge*> _out;
   int _inSize;
   int _outSize; // # of natural edges that should not be removed
};

class Edge
{
public:
   Edge(Lattice* lat, Node* from, Node* to) // used only for checking an edge exist or not
   {
      //_lat = lat;
      _from = from;
      _to = to;
      _type = CHECK;
   }

   Edge(Lattice* lat, Node* from, Node* to, Type type);

   bool operator==(const Edge& e1)
   {
      return (e1.from()==from() && e1.to()==to());
   }


   Node* from() const { return _from; }
   Node* to() const { return _to; }
   string str() { return _from->str()+"-->"+_to->str(); }

private:
   Node* _from;
   Node* _to;
   Type _type;
   //set<int> _func; // functions using it, only for _type = REQUIRED
   //Lattice* _lat;
};

struct classcomp // for hash table of Edge
{
   bool operator() (const Edge* lhs, const Edge* rhs) const
   {
      return ((unsigned long)(lhs->from()->truth_position())+(unsigned long)(lhs->to()->truth_position()<<VAR_NUM)) < ((unsigned long)(rhs->from()->truth_position())+(unsigned long)(rhs->to()->truth_position()<<VAR_NUM));
   }
};
class Lattice
{
public:
   Lattice();

   Node* at(int i) { return _nodeList[i]; }
   set<Edge*, classcomp> edges() const { return _edges; }

   const vector< vector< Node* > >& topo() { return _topological; }
   int topoSize() { return _topological.size(); }
   vector<Node*>& topoStep(int s) { return _topological[s]; }

   void add_implied_edge(Node* from, Node* to);

   //vector<Node*>& nodes(int lvl) { return _nodes.at(lvl); }
   
private:
   //vector< vector< Node* > > _nodes; // easy to retrieve nodes of N bits
   vector<Node*> _nodeList; // _nodeList[i] = Node* with value i
   vector< vector< Node* > > _topological; // topological order
   set<Edge*, classcomp> _edges;
};

class Function
{
public:
   Function(Thre_S* g);

   bool eval(unsigned long t);
   Thre_S* gate() { return _gate; }
   vector<Edge*>& required() { return _required; }

   void findRequiredEdges(Lattice* lat);
   void findRImply(Lattice* lat, Node* from, Node* to);

   void checkAllEdges(vector< vector<Node*> >& topo, Function* f, Graph* g, bool secondCheck, bool cleanUp);
   void checkAndAddAllEdges(vector< vector<Node*> >& topo, Function* f, Graph* g, bool secondCheck);
   bool findPath(Node* A, Node* B);
   //bool shouldCheck() { return _shouldCheck; }
   //void canSkip() { _shouldCheck = false; }
   void setAdj(Function* f) { _adj.insert(f); }
   void setImply(Function* f) { _imply.push_back(f); }
   bool compatible(Function* f) { return (_adj.find(f)!=_adj.end()); }

   void addAllEdges(vector< vector<Node*> >& topo);//, set<int>& knobs);
   void removeAllEdges(vector< vector<Node*> >& topo);
   void reorder(vector< vector<Node*> >& topo, int toStep, int toJ, int targetStep);

   void initD() { _d = _adj.size(); }
   void Dminus() { --_d; }
   void updateNeighborD(vector< vector<Function*> >& D, set<Function*>& _V);
   int d() { return _d; }

private:
   Thre_S* _gate;
   vector<int> _weights;
   int _threshold;
   vector<Edge*> _required;
   vector<Edge*> _R_imply;
   //bool _shouldCheck;
   set<Function*> _adj;
   vector<Function*> _imply;
   int _d; // for degeneracy ordering
};

class Graph
{
public:
   Graph(vector<Function*> functions)
   {
      for (int i=0; i<functions.size(); ++i)
         _V.insert(functions[i]);
   }

   void setConflict(Function* f1, Function* f2);
   void setImply(Function* f1, Function* f2);
   void setCompatible(Function* f1, Function* f2);

   void findCliques();
   bool findCliques_rec(vector<Function*>& R, vector<Function*>& P, vector<Function*>& X);
   void FindSubClique(vector<Function*>& clique, int& i);
   void reportCliques();
   void printV();

   //int minCoverLp();
   void minCoverHeuristic();
   int solveLp(vector<Function*>& funcs, vector<int>& ans);
   int solveLpDd(vector<Function*>& funcs, vector<int>& ans);
   void checkCompatibilityLP(Function* f1, Function* f2);

private:
   set<Function*> _V;
   set<Function*> _Vcopy;
   vector< vector<Function*> > _cliques;
   vector< vector<int> > _implementations;
};

#endif
