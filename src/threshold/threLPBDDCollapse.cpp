#ifndef THRELPBDD_CPP
#define THRELPBDD_CPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <climits>
#include "threshold.h"
#include "base/abc/abc.h"
#include "misc/vec/vec.h"
#include "bdd/extrab/extraBdd.h"
#include "bool/kit/kit.h"

#define nBddSizeMax INT_MAX/2
//#define symDetact
#define canonical

using namespace std;

bool lpbdd_Collapse_cpp(Thre_S* f_ori, Thre_S* l_ori, Thre_S* t);
DdNode * BDDmux_tree_traverse(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin);
DdNode * BDDmux_tree_traverse_front(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin, int ori_cur_weight);
DdNode * BDDdetermine(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin);
DdNode * BDDdetermineF(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin, int ori_cur_weight);
DdNode * BDDleaf0(DdManager * dd);
DdNode * BDDleaf1(DdManager * dd);
extern int max(Vec_Int_t* weights);
extern int min(Vec_Int_t* weights);
extern void threSort(Thre_S* t);
extern int solveLp(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, int n);
extern int solveLpCanonical(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, map<int,int>& symm);
extern int solveLpBlock(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, map<int,int>& symm);
bool bdd2lp(DdManager * dd, DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& minus);
bool bdd2lpDC(DdManager * dd, DdNode * bFunc, DdNode * bDC, vector< vector<int> >& A, vector<int>& B, vector<int>& minus);
void bdd_traverse(DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& row, int fComp, vector<int>& minus);
void bdd_traverse_onSet(DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& row, int fComp, vector<int>& minus);
void bdd_traverse_offSet(DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& row, int fComp, vector<int>& minus);
bool getDontCare(Thre_S* f, Thre_S* l, int k, unsigned * puTruth);
vector< vector<int> > Cudd_GetMinterm(DdManager * manager,DdNode * node, vector<int>& minus);
static void ddGetMintermAux(DdManager * dd,DdNode * node,vector<int>& list,vector< vector<int> >& minterms, vector<int>& minus);

extern void Vec_IntPopFront( Vec_Int_t* p );
extern Thre_S* copyThre(Thre_S* t); //copy constructer
extern void printGate(Thre_S* tObj);

extern int did;


extern "C" bool Th_Collapse_LPBDD(Thre_S* f_ori, Thre_S* l_ori, Thre_S* t)
{
   return lpbdd_Collapse_cpp(f_ori, l_ori, t);
}

extern "C" bool Th_Bdd2th_LP(DdManager * dd, DdNode * bFunc, Thre_S* t)
{
   vector< vector<int> > A;
   vector<int> vec,B,minus;
   vec.resize(dd->size+1);
   A.push_back(vec); 
   minus.resize(dd->size,0); 
   map<int,int> symm;
   for (int i=0; i<dd->size+1; ++i) //includes T
      symm[i]=i;

   if (!bdd2lp(dd,bFunc,A,B,minus))
   {
      //cout<<"not unate"<<endl;
      return 0;
   }
   //Cudd_PrintMinterm(dd,bFunc);
   
#ifdef symDetact
   Cudd_AutodynDisable( dd );
   Cudd_zddVarsFromBddVars( dd, 2 );
   Extra_SymmInfo_t * pSymms;
   if ( !Cudd_IsConstant(bFunc) )
   {
      pSymms = Extra_SymmPairsCompute( dd, bFunc );
      
      int fStart = 1;
      int nSize = pSymms->nVars;
      int i,k;
      int nVars = dd->size;
      int * pVarTaken;
      pVarTaken = ABC_ALLOC( int, nVars );
      memset( pVarTaken, 0, sizeof(int) * nVars );
      for ( i = 0; i < nSize; i++ )
      {
        // skip the variable already considered
        if ( pVarTaken[i] )
            continue;
        // find all the vars symmetric with this one
        for ( k = 0; k < nSize; k++ )
        {
            if ( k == i )
                continue;
            if ( pSymms->pSymms[i][k] == 0 )
                continue;
            // vars i and k are symmetric
            assert( pVarTaken[k] == 0 );
            // there is a new symmetry pair 
            if ( fStart == 1 )
            {  // start a new symmetry class
                fStart = 0;
                //printf( "  { %d ",i);
                // mark the var as taken
                pVarTaken[i] = 1;
            }
            symm[k]=i;
            //printf( "%d ",k );
            // mark the var as taken
            pVarTaken[k] = 1;
        }
        if ( fStart == 0 )
        {
            //printf( " }" );
            fStart = 1; 
        }   
      }   
   }
   Extra_SymmPairsDissolve( pSymms );
#endif

   //for (int i=0; i<dd->size; ++i)
      //cout<<i<<":"<<symm[i]<<endl;
   vector<int> ans;
   #ifdef canonical
   int res = solveLpCanonical(A,B,ans,symm);
   #else
   int res = solveLp(A,B,ans,0);
   #endif
   
      /*//////
      int res = solveLpBlock(A,B,ans,symm);
      vector<int> ans2;
      int res2 = solveLpCanonical(A,B,ans2,symm);
      if (res!=res2){
         cout<<"counter-example found"<<endl;
         for (int i=0; i<ans.size(); i++)
            cout<<ans[i]<<" ";
         cout<<endl;
         for (int i=0; i<ans.size(); i++)
            cout<<ans2[i]<<" ";
         cout<<endl;
         cout<<res<<" "<<res2<<endl;
         exit(1);
      }
      return 1;
      *///////
   
   if (res!=0) // solve failed
   {
      if (res==5)
      {
         //cout<<"not optimal"<<endl;
         /*res=solveLp(A,B,ans,0);
         if (res==0){
            cout<<"oh oh"<<endl;
            Cudd_PrintMinterm(dd,bFunc);
            for (int i=0; i<ans.size(); ++i)cout<<ans[i]<<" "; cout<<endl;
            cout<<endl;
         }*/
      }
      else if (res==6){
         cout<<"counter-example found"<<endl;
         exit(1);
      }
      else 
         cout<<"Abnormal error in solving lp..."<<endl;
      return 0;
   }
   else
   {
      if (!t)
      {
         return 1;
      }
      
      assert(ans.size()-1==minus.size());

      t->thre=ans.back();
      for (int i=0; i<ans.size()-1; ++i)
      {
         Vec_IntPush(t->weights,ans[i]);
         if (minus[i]==1)
         {
            Vec_IntWriteEntry(t->weights,i,-1*ans[i]);
            t->thre -= ans[i];
         }
      }
      for (int i=0; i<ans.size()-1; ++i)
         Vec_IntPush(t->Fanins,i);
      
      t->Fanouts=Vec_IntAlloc(1);
      Vec_IntPush(t->Fanouts,ans.size()-1);
      
#if 0
      if (1){
      printGate(t);
      //cout<<endl;
      }
#endif
      return 1;
   }
}


bool lpbdd_Collapse_cpp(Thre_S* f_ori, Thre_S* l_ori, Thre_S* t)
{
   if (f_ori==NULL || l_ori==NULL){
      cout<<"gate not exist"<<endl;
      return 0;
   }

   Thre_S* f=copyThre(f_ori);
   Thre_S* l=copyThre(l_ori);
   threSort(f); threSort(l);
   int k=-1;
   for (int i=0; i<Vec_IntSize(l->Fanins); ++i)
   {
      if (Vec_IntEntry(l->Fanins,i)==f->Id)
      {
         k=i;
         break;
      }
   }
   if (k==-1){
      cout<<"the gates are not connected"<<endl;
      delete f;
      delete l;
      return 0;
   }

   vector< pair<int,int> > sharedFanin;
   for (int i=0; i<Vec_IntSize(f->Fanins); ++i)
   {
      for (int j=0; j<Vec_IntSize(l->Fanins); ++j)
      {
         if (Vec_IntEntry(f->Fanins,i)==Vec_IntEntry(l->Fanins,j))
         {
            int a=j<k?j:j-1+Vec_IntSize(f->Fanins);
            int b=k+i;
            if (a<b) swap(a,b);
            sharedFanin.push_back(pair<int,int>(a, b));
            break;
         }
      }
   }
   
#if 0
      cout<<"original front:"<<endl;
      printGate(f_ori);
      cout<<"original later:"<<endl;
      printGate(l_ori);
      cout<<endl;
#endif
   DdManager * dd = Cudd_Init( Vec_IntSize(f->weights) + Vec_IntSize(l->weights)-1 , 0 , CUDD_UNIQUE_SLOTS , CUDD_CACHE_SLOTS , 0 );
//#define DC
#define nVM 15
#ifdef DC
   int nVarsMax=nVM; //should be the same as in getDontCare!
   unsigned * puTruth;
   int nBits = 1<<nVarsMax;
   int nWords = (nBits<=32)? 1: (nBits/32);
   puTruth = ABC_ALLOC(unsigned,nWords);
   DdNode* bDC;
   bool doDC=false;
   if (getDontCare(f, l, k, puTruth)) 
   {
      doDC=true;
      //did++; if (did==254) {did++; doDC=false;} else{
      bDC = Kit_TruthToBdd( dd, puTruth, Vec_IntSize(f->weights) + Vec_IntSize(l->weights)-1, 0 );  Cudd_Ref( bDC );
      //}
      //bDC = Cudd_Not(bDC); //don't care set
   }
#endif
   //Cudd_AutodynDisable( dd );
   vector< vector<int> > A,A2;
   vector<int> vec,B,minus,B2,minus2;
   vec.resize(Vec_IntSize(f->weights) + Vec_IntSize(l->weights));
   A.push_back(vec); A2.push_back(vec);
   minus.resize(dd->size,0); minus2.resize(dd->size,0);
   DdNode * bFunc = BDDmux_tree_traverse(dd,f,l,k,0,sharedFanin);

#ifdef DC
   bool oriLP=false;
   if (doDC) {

      DdNode*  bFuncOri = bFunc;
        if (bdd2lp(dd,bFuncOri,A2,B2,minus2)) oriLP=true;
#if 0
      if (did==254){
      DdNode*  bFuncOri = bFunc;
      cout<<"original front:"<<endl;
      printGate(f_ori);
      cout<<"original later:"<<endl;
      printGate(l_ori);
      Cudd_PrintMinterm(dd,bDC);
        if (bdd2lp(dd,bFuncOri,A2,B2,minus2)){
           cout<<"original lp:"<<endl;
      for (int i=0; i<minus2.size(); ++i)
         cout<<minus2[i]<<" "; cout<<endl;
      for (int i=0; i<k; ++i)
         cout<<Vec_IntEntry(l->Fanins,i)<<" ";
      for (int i=0; i<Vec_IntSize(f->Fanins); ++i)
         cout<<Vec_IntEntry(f->Fanins,i)<<" ";
      for (int i=k+1; i<Vec_IntSize(l->Fanins); ++i)
         cout<<Vec_IntEntry(l->Fanins,i)<<" ";
      cout<<endl;
      for (int i=0; i<A2.size(); ++i)
      {
         for (int j=0; j<A2[i].size(); ++j)
            cout<<setw(2)<<A2[i][j]<<" ";
         cout<<"| "<<B2[i]<<endl;
      }
      cout<<endl;
        } else cout<<"original not unate"<<endl;}
#endif
      
      if (!bdd2lpDC(dd,bFunc,bDC,A,B,minus))
      {
         //cout<<"dont care not unate"<<endl;
         Cudd_Quit(dd);
         delete f;
         delete l;
         return 0;
      }
   }
   else
   {
      if (!bdd2lp(dd,bFunc,A,B,minus))
      {
         Cudd_Quit(dd);
         delete f;
         delete l;
         return 0;
      }
   }
   ABC_FREE(puTruth);
#endif

#ifndef DC
   if (!bdd2lp(dd,bFunc,A,B,minus))
   {
      Cudd_Quit(dd);
      delete f;
      delete l;
      return 0;
   }
#endif

   Cudd_Quit(dd);
   for (int i=0; i<sharedFanin.size(); ++i)
   {
      for (int j=0; j<A.size(); ++j)
         A[j][sharedFanin[i].first]=0;
   }
#if 0
      cout<<"original front:"<<endl;
      printGate(f_ori);
      cout<<"original later:"<<endl;
      printGate(l_ori);
      for (int i=0; i<minus.size(); ++i)
         cout<<minus[i]<<" "; cout<<endl;
      for (int i=0; i<k; ++i)
         cout<<Vec_IntEntry(l->Fanins,i)<<" ";
      for (int i=0; i<Vec_IntSize(f->Fanins); ++i)
         cout<<Vec_IntEntry(f->Fanins,i)<<" ";
      for (int i=k+1; i<Vec_IntSize(l->Fanins); ++i)
         cout<<Vec_IntEntry(l->Fanins,i)<<" ";
      cout<<endl;
#endif
#if 0
      if (doDC && did==254){
      for (int i=0; i<A.size(); ++i)
      {
         for (int j=0; j<A[i].size(); ++j)
            cout<<setw(2)<<A[i][j]<<" ";
         cout<<"| "<<B[i]<<endl;
      }
      cout<<endl;}
#endif

   vector<int> ans,ans2;
   int res = solveLp(A,B,ans,0);
   //int res = solveLpCanonical(A,B,ans);
   
   if (res!=0) // solve failed
   {
      if (res!=5)
         cout<<"Abnormal error in solving lp..."<<endl;
#ifdef DC
      if (oriLP)
      {
         res = solveLp(A2,B2,ans2,0);
         if (res==0)
         {
            cout<<"DC failed but original solved"<<endl;
      for (int i=0; i<A2.size(); ++i)
      {
         for (int j=0; j<A2[i].size(); ++j)
            cout<<setw(2)<<A2[i][j]<<" ";
         cout<<"| "<<B2[i]<<endl;
      }
      cout<<endl;
            for(int i=0; i<ans2.size(); ++i) cout<<i<<" "; cout<<endl;
         }
      }
#endif
      delete f;
      delete l;
      return 0;
   }
   else
   {
      if (!t)
      {
         delete f;
         delete l;
         return 1;
      }
      /*for (int i=0; i<ans.size(); ++i)
         cout<<ans[i]<<" ";
      cout<<endl;*/
      
      //cout<<ans.size()<<" "<<minus.size()<<" "<<Vec_IntSize(f->Fanins)+Vec_IntSize(l->Fanins)<<endl; 
      assert(ans.size()-1==minus.size());

      t->thre=ans.back();
      for (int i=0; i<ans.size()-1; ++i)
      {
         Vec_IntPush(t->weights,ans[i]);
         if (minus[i]==1)
         {
            Vec_IntWriteEntry(t->weights,i,-1*ans[i]);
            t->thre -= ans[i];
         }
      }
      for (int i=0; i<k; ++i)
         Vec_IntPush(t->Fanins,Vec_IntEntry(l->Fanins,i));
      for (int i=0; i<Vec_IntSize(f->Fanins); ++i)
         Vec_IntPush(t->Fanins,Vec_IntEntry(f->Fanins,i));
      for (int i=k+1; i<Vec_IntSize(l->Fanins); ++i)
         Vec_IntPush(t->Fanins,Vec_IntEntry(l->Fanins,i));
      
      t->Fanouts=Vec_IntDup(l->Fanouts);
      
      sort(sharedFanin.begin(),sharedFanin.end());
      for (int i=sharedFanin.size()-1; i>=0; --i)
      {
         Vec_IntDrop(t->Fanins,sharedFanin[i].first);
         Vec_IntDrop(t->weights,sharedFanin[i].first);
      }
#if 0
      if (1){
      cout<<"original front:"<<endl;
      printGate(f_ori);
      cout<<"original later:"<<endl;
      printGate(l_ori);
      cout<<"merged gate:"<<endl;
      printGate(t);
      cout<<endl;}
#endif
      delete f;
      delete l;
      return 1;
   }
}

DdNode * BDDmux_tree_traverse(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin)
{
   int cur_weight=Vec_IntEntry(l->weights,0);
   Vec_IntPopFront(l->weights);
   
   if (k==0)
   {
      Thre_S* f2=copyThre(f);
      DdNode * bFunc = BDDmux_tree_traverse_front(dd,f2,l,k-1,lvl,sharedFanin,cur_weight);
      delete f2;
      return bFunc;
   }
   else
   {
      DdNode * bFunc, * bFunc0, * bFunc1, * bFuncC;
      
      int skip = -1;
      for (int i=0; i<sharedFanin.size(); ++i)
      {
         if (lvl==sharedFanin[i].first)
         {
            skip=sharedFanin[i].second;
            break;
         }
      }
      
      if (skip==-1)
      {
         bFuncC = dd->vars[lvl];
         Cudd_Ref(bFuncC);
      }
      else
      {
         bFuncC = dd->vars[skip];
         Cudd_Ref(bFuncC);
      }
      
      // 0 branch
      Thre_S* l2=copyThre(l);
      bFunc0 = BDDdetermine(dd,f,l2,k-1,lvl+1,sharedFanin);
      Cudd_Ref(bFunc0);
      delete l2;
      
      // 1 branch
      l->thre-=cur_weight;
      bFunc1 = BDDdetermine(dd,f,l,k-1,lvl+1,sharedFanin);
      Cudd_Ref(bFunc1);
      
      bFunc = Cudd_bddIte( dd, bFuncC, bFunc1, bFunc0 );   Cudd_Ref( bFunc );
      Cudd_RecursiveDeref( dd, bFunc0 );
      Cudd_RecursiveDeref( dd, bFunc1 );
      Cudd_RecursiveDeref( dd, bFuncC );
      
      return bFunc;
   }
}

DdNode * BDDmux_tree_traverse_front(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin, int ori_cur_weight)
{
   int cur_weight=Vec_IntEntry(f->weights,0);
   Vec_IntPopFront(f->weights);
   
   DdNode * bFunc, * bFunc0, * bFunc1, * bFuncC;
   int skip = -1;
   for (int i=0; i<sharedFanin.size(); ++i)
   {
      if (lvl==sharedFanin[i].first)
      {
         skip=sharedFanin[i].second;
         break;
      }
   }
   if (skip==-1)
   {
      bFuncC = dd->vars[lvl];
      Cudd_Ref(bFuncC);
   }
   else
   {
      bFuncC = dd->vars[skip];
      Cudd_Ref(bFuncC);
   }

   // 0 branch
   Thre_S* f2=copyThre(f);
   Thre_S* l2=copyThre(l);
   bFunc0 = BDDdetermineF(dd,f2,l2,k,lvl+1,sharedFanin,ori_cur_weight);
   delete f2;
   delete l2;
   
   // 1 branch
   f->thre-=cur_weight;
   bFunc1 = BDDdetermineF(dd,f,l,k,lvl+1,sharedFanin,ori_cur_weight);

   bFunc = Cudd_bddIte( dd, bFuncC, bFunc1, bFunc0 );   Cudd_Ref( bFunc );
   Cudd_RecursiveDeref( dd, bFunc0 );
   Cudd_RecursiveDeref( dd, bFunc1 );
   Cudd_RecursiveDeref( dd, bFuncC );
   
   return bFunc;
}

DdNode * BDDdetermine(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin)
{
   if (max(l->weights) < l->thre)
      return BDDleaf0(dd);
   else if (min(l->weights) >= l->thre)
      return BDDleaf1(dd);
   else
      return BDDmux_tree_traverse(dd,f,l,k,lvl,sharedFanin);
}

DdNode * BDDdetermineF(DdManager * dd, Thre_S* f, Thre_S* l, int k, int lvl, vector< pair<int,int> >& sharedFanin, int ori_cur_weight)
{
   if (max(f->weights) < f->thre) // 0 leaf of front --> 0 branch of late
   {
      while (Vec_IntSize(f->weights)>0)
      {
         lvl++;
         Vec_IntPopFront(f->weights);
      }
      return BDDdetermine(dd,f,l,k,lvl,sharedFanin);
   }
   else if (min(f->weights) >= f->thre)
   {
      l->thre-=ori_cur_weight;
      while (Vec_IntSize(f->weights)>0)
      {
         lvl++;
         Vec_IntPopFront(f->weights);
      }
      return BDDdetermine(dd,f,l,k,lvl,sharedFanin);
   }
   else
      return BDDmux_tree_traverse_front(dd,f,l,k,lvl,sharedFanin,ori_cur_weight);
}

DdNode * BDDleaf0(DdManager * dd)
{
   return Cudd_Not(dd->one);//dd->zero;
}

DdNode * BDDleaf1(DdManager * dd)
{
   return dd->one;
}

bool bdd2lp(DdManager * dd, DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& minus)
{
   assert(A[0].size()==minus.size()+1);
   // make positive unate
   for (int i=0; i<dd->size; ++i)
   {
      int unate = Extra_bddCheckUnateNaive(dd,bFunc,i);
      if (unate==0)
         return false;
      else if (unate==-1)
      {
         minus[i]=1;
      }
   }
   
   // go through each path and translate
   vector<int> row;
   row.resize(A[0].size(),2);
   bdd_traverse(bFunc,A,B,row,0,minus);
   
   A.erase(A.begin());
   return true;
}

bool bdd2lpDC(DdManager * dd, DdNode * bFunc, DdNode * bDC, vector< vector<int> >& A, vector<int>& B, vector<int>& minus)
{
   assert(A[0].size()==minus.size()+1);
   // make positive unate
   for (int i=0; i<dd->size; ++i)
   {
      int unate = Extra_bddCheckUnateNaive(dd,bFunc,i);
      if (unate==0)
         return false;
      else if (unate==-1)
      {
         minus[i]=1;
      }
   }
   
   vector<int> row;
   // care on-set
   DdNode * onSet = Cudd_bddAnd(dd,bFunc,bDC); Cudd_Ref(onSet);
   if (did==254) {cout<<"on-set:"<<endl; Cudd_PrintMinterm(dd,onSet);}
   row.resize(A[0].size(),2);
   bdd_traverse_onSet(onSet,A,B,row,0,minus);
   
   DdNode * offSet = Cudd_bddAnd(dd,Cudd_Not(bFunc),bDC); Cudd_Ref(offSet);
   if (did==254) {cout<<"off-set:"<<endl; Cudd_PrintMinterm(dd,offSet);}
   row.resize(0);
   row.resize(A[0].size(),2);
   bdd_traverse_offSet(offSet,A,B,row,0,minus);
   
   A.erase(A.begin());
   return true;
}

void bdd_traverse(DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& row, int fComp, vector<int>& minus)
{
   int numPi = Cudd_Regular( bFunc )->index;
   if (Cudd_IsComplement( bFunc ))
      fComp += 1;
   
   if ( Cudd_IsConstant( bFunc ) )
   {
#if 0            
      for (int i=0; i<row.size(); ++i)
         cout<<row[i]<<" ";
      cout<<endl;
#endif
      if (fComp%2 == 1) // 0 leaf
      {
         vector<int> nrow;
         nrow.resize(A[0].size(),-1);
         for (int i=0; i<row.size(); ++i)
            if ((row[i]==0&&!minus[i]) || (row[i]==1&&minus[i])) nrow[i]=0;
         nrow.back()=1;
         A.push_back(nrow);
         B.push_back(1);
      }
      else // 1 leaf
      {
         vector<int> nrow;
         nrow.resize(A[0].size(),0);
         for (int i=0; i<row.size(); ++i)
            if ((row[i]==1&&!minus[i]) || (row[i]==0&&minus[i])) nrow[i]=1;
         nrow.back()=-1;
         A.push_back(nrow);
         B.push_back(0);
      }
      return;
   }
   
   row[numPi]=1;
   bdd_traverse(Cudd_T(bFunc),A,B,row,fComp,minus);
   row[numPi]=0;
   bdd_traverse(Cudd_E(bFunc),A,B,row,fComp,minus);
   row[numPi]=2;
}

void bdd_traverse_onSet(DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& row, int fComp, vector<int>& minus)
{
   int numPi = Cudd_Regular( bFunc )->index;
   if (Cudd_IsComplement( bFunc ))
      fComp += 1;
   
   if ( Cudd_IsConstant( bFunc ) )
   {
#if 0
      for (int i=0; i<row.size(); ++i)
         cout<<row[i]<<" ";
      cout<<endl;
#endif
      if (fComp%2 == 0) // 1 leaf
      {
         vector<int> nrow;
         nrow.resize(A[0].size(),0);
         for (int i=0; i<row.size(); ++i)
            if ((row[i]==1&&!minus[i]) || (row[i]==0&&minus[i])) nrow[i]=1;
         nrow.back()=-1;
         A.push_back(nrow);
         B.push_back(0);
      }
      return;
   }
   
   row[numPi]=1;
   bdd_traverse_onSet(Cudd_T(bFunc),A,B,row,fComp,minus);
   row[numPi]=0;
   bdd_traverse_onSet(Cudd_E(bFunc),A,B,row,fComp,minus);
   row[numPi]=2;
}

void bdd_traverse_offSet(DdNode * bFunc, vector< vector<int> >& A, vector<int>& B, vector<int>& row, int fComp, vector<int>& minus)
{
   int numPi = Cudd_Regular( bFunc )->index;
   if (Cudd_IsComplement( bFunc ))
      fComp += 1;
   
   if ( Cudd_IsConstant( bFunc ) )
   {
#if 0
      cout<<"off set leaf "<<fComp<<": ";
      for (int i=0; i<row.size(); ++i)
         cout<<row[i]<<" ";
      cout<<endl;
#endif
      if (fComp%2 == 0) // 1 leaf
      {
         vector<int> nrow;
         nrow.resize(A[0].size(),-1);
         for (int i=0; i<row.size(); ++i)
            if ((row[i]==0&&!minus[i]) || (row[i]==1&&minus[i])) nrow[i]=0;
         nrow.back()=1;
         A.push_back(nrow);
         B.push_back(1);
      }
      return;
   }
   
   row[numPi]=1;
   bdd_traverse_offSet(Cudd_T(bFunc),A,B,row,fComp,minus);
   row[numPi]=0;
   bdd_traverse_offSet(Cudd_E(bFunc),A,B,row,fComp,minus);
   row[numPi]=2;
}

bool getDontCare(Thre_S* f, Thre_S* l, int k, unsigned* puTruth)
{
   assert(!(Vec_PtrSize(Th_aigMap)-1 < l->Id));

   int nVarsMax=nVM, nLevels=9; //5~15,1~9, verbose=0, veryVerbose=0
   Odc_Man_t * pManOdc;
   Abc_Obj_t * pNode = (Abc_Obj_t*)Vec_PtrEntry(Th_aigMap,l->Id);
   Vec_Ptr_t * vLeaves = Vec_PtrAlloc(15);
   int Entry,i,j;
   Vec_IntForEachEntry(l->Fanins,Entry,i){
      if (i==k) {
         Vec_IntForEachEntry(f->Fanins,Entry,j) 
            Vec_PtrPush(vLeaves,Vec_PtrEntry(Th_aigMap,Entry)); 
         continue;
      }
      // FIXME:
      //Vec_PtrPushUnique(vLeaves,Vec_PtrEntry(Th_aigMap,Entry));
      Vec_PtrPush(vLeaves,Vec_PtrEntry(Th_aigMap,Entry));
   }
   if (Vec_PtrSize(vLeaves)>nVarsMax) {
      Vec_PtrFree(vLeaves);
      return false;
   }

   //pManOdc = Abc_NtkDontCareAlloc( Vec_PtrSize(vLeaves), nLevels, 0,0 );
   pManOdc = Abc_NtkDontCareAlloc( nVarsMax, nLevels, 0,0 );
   Abc_NtkDontCareClear( pManOdc );
   int ret = Abc_NtkDontCareCompute( pManOdc, pNode, vLeaves, puTruth );
   //cout<<"ret"<<ret<<endl;
   //Abc_AigPrintNode(pNode);
   //if (ret!=0)
   //cout<<f->Id<<" "<<l->Id<<" "<<ret<<endl;
   Abc_NtkDontCareFree( pManOdc );
   Vec_PtrFree(vLeaves);
   return (ret!=0);
}

vector< vector<int> > Cudd_GetMinterm(DdManager * manager,DdNode * node, vector<int>& minus)
{
    vector<int> list;
    vector< vector<int> > minterms;
   
    list.resize(manager->size,2);

    ddGetMintermAux(manager,node,list,minterms,minus);
    return minterms;
} /* end of Cudd_PrintMinterm */

static void ddGetMintermAux(DdManager * dd,DdNode * node,vector<int>& list,vector< vector<int> >& minterms, vector<int>& minus)
{
    DdNode      *N,*Nv,*Nnv;
    int         index;
    N = Cudd_Regular(node);

    if (cuddIsConstant(N)) {
        /* Terminal case: Print one cube based on the current recursion
        ** path, unless we have reached the background value (ADDs) or
        ** the logical zero (BDDs).
        */
        if (node != dd->background && node != Cudd_Not(dd->one)) {
           assert(cuddV(node)==1);
           minterms.push_back(list);
        }
    } else {
        Nv  = cuddT(N);
        Nnv = cuddE(N);
        //cout<<N->index<<" "<<Cudd_IsComplement(node)<<" "<<minus[N->index]<<endl;
        //if ((Cudd_IsComplement(node)&&!minus[N->index]) || (!Cudd_IsComplement(node)&&minus[N->index])) {
        if (Cudd_IsComplement(node)) {
            Nv  = Cudd_Not(Nv);
            Nnv = Cudd_Not(Nnv);
        }
        index = N->index;
        list[index] = minus[index]? 1: 0;
        ddGetMintermAux(dd,Nnv,list,minterms,minus);
        list[index] = minus[index]? 0: 1;
        ddGetMintermAux(dd,Nv,list,minterms,minus);
        list[index] = 2;
    }
    return;
}
DdNode* th2bdd(Thre_S* f_ori, DdManager* dd);
extern "C" DdNode * Th2Bdd(Thre_S* f, DdManager * dd)
{
   return th2bdd(f, dd);
}

DdNode* th2bdd(Thre_S* f_ori, DdManager* dd)
{
   //printGate(f_ori);
   if (f_ori==NULL){
      cout<<"gate not exist"<<endl;
      return 0;
   }

   Thre_S* f=copyThre(f_ori);
   threSort(f);    

   vector< pair<int,int> > sharedFanin;
   DdNode* bFunc = BDDmux_tree_traverse(dd,0,f,-1,0,sharedFanin);

   //Cudd_PrintMinterm(dd,bFunc);
   return bFunc;
}
#endif

