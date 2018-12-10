#ifndef THRELP_CPP
#define THRELP_CPP

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
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

using namespace std;
#define M 1000
#define completeSolve 1
#define recursiveSolve 1 // if completeSolve==1, use recursive 2^n (1) or encode to constraints(0)

bool lp_Collapse_cpp(Thre_S* f_ori, Thre_S* l_ori, Thre_S* t);
void mux_tree_traverse(Thre_S* f, Thre_S* l, vector<int> selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, vector< pair<int,int> >& sharedFanin);
void mux_tree_traverse_front(Thre_S* f, Thre_S* l, vector<int> selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, int ori_cur_weight, vector< pair<int,int> >& sharedFanin);
void determine(Thre_S* f, Thre_S* l, vector<int>& selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, vector< pair<int,int> >& sharedFanin);
void determineF(Thre_S* f, Thre_S* l, vector<int> selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, int ori_cur_weight, vector< pair<int,int> >& sharedFanin);
void leaf0(vector<int>& selection, vector< vector<int> >& A, vector<int>& B);
void leaf1(vector<int>& selection, vector< vector<int> >& A, vector<int>& B);
int max(Vec_Int_t* weights);
int min(Vec_Int_t* weights);
void threSort(Thre_S* t);
int solveLpRescursive(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, vector< pair<int,int> >& sharedFanin, vector<int>& minus, int i);
int solveLpMultiConstraint(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, vector< pair<int,int> >& sharedFanin, vector<int>& minus);
extern int solveLp(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, int n);
extern int solveLpBlock(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, map<int,int>& symm);
extern int solveLpCanonical(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, map<int,int>& symm);
void Vec_IntPopFront( Vec_Int_t* p )
{
   assert(p->nSize>0);
   for (int i=0; i<p->nSize-1; ++i)
      p->pArray[i]=p->pArray[i+1];
   --p->nSize;
}

Thre_S* copyThre(Thre_S* t) //copy constructer
{
   Thre_S* cp = new Thre_S;
   cp->thre=t->thre;
   cp->Type=t->Type;
   cp->Ischoose=t->Ischoose;
   cp->Id=t->Id;
   cp->oId=t->oId;
   cp->nId=t->nId;
   cp->cost=t->cost;
   cp->level=t->level;
   cp->pName=t->pName;
   cp->weights=Vec_IntDup(t->weights);
   cp->Fanins=Vec_IntDup(t->Fanins);
   cp->Fanouts=Vec_IntDup(t->Fanouts);
   cp->pCopy=t->pCopy;
   return cp;
}

extern "C" bool Th_Collapse_LP(Thre_S* f_ori, Thre_S* l_ori, Thre_S* t)
{
   return lp_Collapse_cpp(f_ori, l_ori, t);
}

bool lp_Collapse_cpp(Thre_S* f_ori, Thre_S* l_ori, Thre_S* t)
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

   vector<int> minus;
   for (int i=0; i<k; ++i)
   {
      if (Vec_IntEntry(l->weights,i)<0) minus.push_back(1);
      else minus.push_back(0);
   }
   for (int i=0; i<Vec_IntSize(f->weights); ++i)
   {
      if (Vec_IntEntry(l->weights,k)<0)
      {
         if (Vec_IntEntry(f->weights,i)<0) minus.push_back(0);
         else minus.push_back(1);
      }
      else
      {
         if (Vec_IntEntry(f->weights,i)<0) minus.push_back(1);
         else minus.push_back(0);
      }
   }
   for (int i=k+1; i<Vec_IntSize(l->weights); ++i)
   {
      if (Vec_IntEntry(l->weights,i)<0) minus.push_back(1);
      else minus.push_back(0);
   }

   vector< pair<int,int> > sharedFanin;
   for (int i=0; i<Vec_IntSize(f->Fanins); ++i)
   {
      for (int j=0; j<Vec_IntSize(l->Fanins); ++j)
      {
         if (Vec_IntEntry(f->Fanins,i)==Vec_IntEntry(l->Fanins,j))
         {
            //cout<<"shared fanin"<<endl;
            int a=j<k?j:j-1+Vec_IntSize(f->Fanins);
            int b=k+i;
            if (a<b) swap(a,b);
            sharedFanin.push_back(pair<int,int>(a, b));
            break;
         }
      }
   }

   vector<int> selection;
   vector< vector<int> > A;
   vector<int> vec;
   vec.resize(Vec_IntSize(f->weights) + Vec_IntSize(l->weights));
   A.push_back(vec);
   vector<int> B;
   mux_tree_traverse(f,l,selection,k,A,B,minus,sharedFanin);
   A.erase(A.begin());
   for (int i=0; i<sharedFanin.size(); ++i)
   {
      for (int j=0; j<A.size(); ++j)
         A[j][sharedFanin[i].first]=0;
   }

#if 0
      if (sharedFanin.size()!=0){
         cout<<sharedFanin[0].first<<" "<<sharedFanin[0].second<<" "<<k<<endl;
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
   for (int i=0; i<A.size(); ++i)
   {
      for (int j=0; j<A[i].size(); ++j)
         cout<<setw(2)<<A[i][j]<<" ";
      cout<<"| "<<B[i]<<endl;
   }
      cout<<endl;}
#endif

   for(int j=0; j<A[0].size(); ++j)
   {
      bool skip=false;
      for (int m=0; m<sharedFanin.size(); ++m)
      {
         if (j==sharedFanin[m].second && minus[j]!=minus[sharedFanin[m].first])
         {
           skip=true;
           break;
         }
      }
      if (!completeSolve) skip=false;
      if (!skip)
      {
         for (int i=0; i<A.size(); ++i)
         {
            if (A[i][j]==2) A[i][j]=0;
            else if (A[i][j]==-2) A[i][j]=-1;
         }
      }
   } 
#if 0
   cout<<endl;
   for (int i=0; i<sharedFanin.size(); ++i) cout<<sharedFanin[i].second<<" "; cout<<endl;
   cout<<"original A"<<endl;
   for (int i=0; i<minus.size(); ++i) cout<<minus[i]<<" "; cout<<endl;
   for (int i=0; i<A.size(); ++i)
   {
      for (int j=0; j<A[i].size(); ++j)
         cout<<setw(2)<<A[i][j]<<" ";
      cout<<"| "<<B[i]<<endl;
   }
#endif
   vector<int> ans;
   int res;
   if (completeSolve && recursiveSolve)
      res = solveLpRescursive(A,B,ans,sharedFanin,minus,0);
   else if (completeSolve && !recursiveSolve)
      res = solveLpMultiConstraint(A,B,ans,sharedFanin,minus);
   else
      res = solveLp(A,B,ans,0);
   
   if (res!=0) // solve failed
   {
#ifdef checkCononical
      if (res==-1)
      {
         cout<<"original front:"<<endl;
         printGate(f_ori);
         cout<<"original later:"<<endl;
         printGate(l_ori);
      }
#endif
      if (res!=5)
         cout<<"Abnormal error in solving lp..."<<endl;
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
      //cout<<t->Id<<" built"<<endl;
      return 1;
   }
}

void mux_tree_traverse(Thre_S* f, Thre_S* l, vector<int> selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, vector< pair<int,int> >& sharedFanin)
{
   int cur_weight=Vec_IntEntry(l->weights,0);
   Vec_IntPopFront(l->weights);
   
   if (k==0)
   {
      Thre_S* f2=copyThre(f);
      mux_tree_traverse_front(f2,l,selection,k-1,A,B,minus,cur_weight,sharedFanin);
      delete f2;
   }
   else
   {
      int skip = -1;
      for (int i=0; i<sharedFanin.size(); ++i)
      {
         if (selection.size()==sharedFanin[i].first)
         {
            skip=sharedFanin[i].second;
            break;
         }
      }
      
      if (skip==-1)
      {
         // 0 branch
         if (cur_weight>0)
            selection.push_back(0);
         else
            selection.push_back(1);
         Thre_S* l2=copyThre(l);
         determine(f,l2,selection,k-1,A,B,minus,sharedFanin);
         delete l2;
         
         // 1 branch
         if (cur_weight>0)
            selection.back()=1;
         else
            selection.back()=0;
         l->thre-=cur_weight;
         determine(f,l,selection,k-1,A,B,minus,sharedFanin);
      }
      else
      {
         if (selection[skip]==2) // not determined
         {
            // 0 branch
            if (minus[skip]==0)
               selection[skip]=0;
            else
               selection[skip]=1;
            selection.push_back(0);
            //for (int i=0; i<selection.size();++i) cout<<selection[i]<<" "; cout<<endl;
            Thre_S* l2=copyThre(l);
            determine(f,l2,selection,k-1,A,B,minus,sharedFanin);
            delete l2;
         
            // 1 branch
            if (selection[skip]==0)
               selection[skip]=1;
            else
               selection[skip]=0;
            //for (int i=0; i<selection.size();++i) cout<<selection[i]<<" "; cout<<endl;
            l->thre-=cur_weight;
            determine(f,l,selection,k-1,A,B,minus,sharedFanin);
         }
         else
         {
            selection.push_back(0);
            //for (int i=0; i<selection.size();++i) cout<<selection[i]<<" "; cout<<endl;
            if (selection[skip]+minus[skip]==1) //chose 1 and not minus or chose 0 and minus
               l->thre-=cur_weight;
            determine(f,l,selection,k-1,A,B,minus,sharedFanin);
         }
      }
   }
}

void mux_tree_traverse_front(Thre_S* f, Thre_S* l, vector<int> selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, int ori_cur_weight, vector< pair<int,int> >& sharedFanin)
{
   int cur_weight=Vec_IntEntry(f->weights,0);
   Vec_IntPopFront(f->weights);
   
   int skip = -1;
   for (int i=0; i<sharedFanin.size(); ++i)
   {
      if (selection.size()==sharedFanin[i].first)
      {
         skip=sharedFanin[i].second;
         break;
      }
   }
   if (skip==-1)
   {
      // 0 branch
      if (ori_cur_weight<0)
      {
         if (cur_weight>0)
            selection.push_back(1);
         else
            selection.push_back(0);
      }
      else
      {
         if (cur_weight>0)
            selection.push_back(0);
         else
            selection.push_back(1);
      }
      Thre_S* f2=copyThre(f);
      Thre_S* l2=copyThre(l);
      determineF(f2,l2,selection,k,A,B,minus,ori_cur_weight,sharedFanin);
      delete f2;
      delete l2;
      
      // 1 branch
      if (selection.back()==1) selection.back()=0;
      else selection.back()=1;
      f->thre-=cur_weight;
      determineF(f,l,selection,k,A,B,minus,ori_cur_weight,sharedFanin);
   }
   else
   {
      selection.push_back(0);
      if (selection[skip]+minus[skip]==1) //chose 1 and not minus or chose 0 and minus
         f->thre-=cur_weight;
      determineF(f,l,selection,k,A,B,minus,ori_cur_weight,sharedFanin);
   }
}

void determine(Thre_S* f, Thre_S* l, vector<int>& selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, vector< pair<int,int> >& sharedFanin)
{
   if (max(l->weights) < l->thre)
      leaf0(selection,A,B);
   else if (min(l->weights) >= l->thre)
      leaf1(selection,A,B);
   else
      mux_tree_traverse(f,l,selection,k,A,B,minus,sharedFanin);
}

void determineF(Thre_S* f, Thre_S* l, vector<int> selection, int k, vector< vector<int> >& A, vector<int>& B, vector<int>& minus, int ori_cur_weight, vector< pair<int,int> >& sharedFanin)
{
   if (max(f->weights) < f->thre) // 0 leaf of front --> 0 branch of late
   { 
      while (Vec_IntSize(f->weights)>0)
      {
         selection.push_back(2);
         Vec_IntPopFront(f->weights);
      }
      determine(f,l,selection,k,A,B,minus,sharedFanin);
   }
   else if (min(f->weights) >= f->thre)
   {
      l->thre-=ori_cur_weight;
      while (Vec_IntSize(f->weights)>0)
      {
         selection.push_back(2);
         Vec_IntPopFront(f->weights);
      }
      determine(f,l,selection,k,A,B,minus,sharedFanin);
   }
   else
      mux_tree_traverse_front(f,l,selection,k,A,B,minus,ori_cur_weight,sharedFanin);
}

void leaf0(vector<int>& selection, vector< vector<int> >& A, vector<int>& B)
{
   vector<int> row;
   for (int i=0; i<selection.size(); ++i)
   {
      row.push_back(selection[i]*-1);
      //if (row.back()==-2)
         //row.back()=-1;
   }
   while (row.size()<A[0].size())
      row.push_back(-2);//-1);
   row.back()=1;
   A.push_back(row);
   B.push_back(1);
   
   /*cout<<"selection:";
   for (int i=0; i<selection.size(); ++i)
      cout<<selection[i]<<" ";
   cout<<"leaf0"<<endl<<"row: ";
   for (int i=0; i<row.size(); ++i)
      cout<<row[i]<<" ";
   cout<<"| 1"<<endl;*/
}

void leaf1(vector<int>& selection, vector< vector<int> >& A, vector<int>& B)
{
   vector<int> row;
   for (int i=0; i<selection.size(); ++i)
   {
      row.push_back(selection[i]);
      //if (row.back()==2)
         //row.back()=0;
   }
   while (row.size()<A[0].size())
      row.push_back(2);//0);
   row.back()=-1;
   A.push_back(row);
   B.push_back(0);
   
   /*cout<<"selection:";
   for (int i=0; i<selection.size(); ++i)
      cout<<selection[i]<<" ";
   cout<<"leaf1"<<endl<<"row: ";
   for (int i=0; i<row.size(); ++i)
      cout<<row[i]<<" ";
   cout<<"| 0"<<endl;*/
}

int max(Vec_Int_t * weights)
{
   int ans=0;
   
   for (int i=0; i<Vec_IntSize(weights); ++i)
   {
      if (Vec_IntEntry(weights,i)>0)
         ans+=Vec_IntEntry(weights,i);
   }
   return ans;
}

int min(Vec_Int_t* weights)
{
   int ans=0;
   for (int i=Vec_IntSize(weights)-1; i>=0; --i)
   {
      if (Vec_IntEntry(weights,i)<0)
         ans+=Vec_IntEntry(weights,i);
   }
   return ans;
}

const bool absSort(pair<int,int> a, pair<int,int> b)
{
   return abs(a.first)<abs(b.first);
}

void threSort(Thre_S* t)
{
   assert(Vec_IntSize(t->weights)==Vec_IntSize(t->Fanins));
   vector< pair<int,int> > Wpair;
   for (int i=0; i<Vec_IntSize(t->weights); ++i)
      Wpair.push_back(pair<int,int>(Vec_IntEntry(t->weights,i),Vec_IntEntry(t->Fanins,i)));
   sort(Wpair.begin(),Wpair.end(),absSort);
   for (int i=0; i<Vec_IntSize(t->weights); ++i)
   {
      Vec_IntWriteEntry(t->weights,i,Wpair[Vec_IntSize(t->weights)-i-1].first);
      Vec_IntWriteEntry(t->Fanins,i,Wpair[Vec_IntSize(t->weights)-i-1].second);
   }
}

int solveLpRescursive(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, vector< pair<int,int> >& sharedFanin, vector<int>& minus, int s)
{
   if (s==sharedFanin.size()) 
   {
#if 0
   for (int jj=0; jj<minus.size(); ++jj) cout<<minus[jj]<<" "; cout<<endl;
   for (int ii=0; ii<A.size(); ++ii)
   {
      for (int jj=0; jj<A[ii].size(); ++jj)
         cout<<setw(2)<<A[ii][jj]<<" ";
      cout<<"| "<<B[ii]<<endl;
   }
#endif
      return solveLp(A,B,ans,0);
   }
   
   int j=sharedFanin[s].second;
   int res; 
   vector< vector<int> > A2=A;
   for (int i=0; i<A.size(); ++i)
   {
      if (A[i][j]==2) A2[i][j]=0;
      else if (A[i][j]==-2) A2[i][j]=-1;
   }
   res = solveLpRescursive(A2,B,ans,sharedFanin,minus,s+1);
   if (res==0) return 0;
   if (minus[sharedFanin[s].first]==minus[sharedFanin[s].second]) return res;
   
   for (int i=0; i<A.size(); ++i)
   {
      if (A[i][j]==2) A2[i][j]=0;
      else if (A[i][j]==-2) A2[i][j]=-1;
      else if (A[i][j]==1) A2[i][j]=0;
      else if (A[i][j]==-1) A2[i][j]=0;
      else if (A[i][j]==0 && B[i]==0) A2[i][j]=1; //l-leaf
      else if (A[i][j]==0 && B[i]==1) A2[i][j]=-1; //0-leaf
      else assert(0);
   }
   if (minus[sharedFanin[s].second]==0) minus[sharedFanin[s].second]=1;
   else minus[sharedFanin[s].second]=0;
   res = solveLpRescursive(A2,B,ans,sharedFanin,minus,s+1);
   if (res==0) return 0;
   else
   {
      if (minus[sharedFanin[s].second]==0) minus[sharedFanin[s].second]=1;
      else minus[sharedFanin[s].second]=0;
      return res;
   }
}

int solveLpMultiConstraint(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, vector< pair<int,int> >& sharedFanin, vector<int>& minus)
{
   vector<int> sfId;
   for (int s=0; s<sharedFanin.size(); ++s)
   {
#if 0
   for (int i=0; i<A.size(); ++i)
   {
      for (int j=0; j<A[i].size(); ++j)
         cout<<setw(2)<<A[i][j]<<" ";
      cout<<"| "<<B[i]<<endl;
   }cout<<endl;
#endif
      if (minus[sharedFanin[s].first]==minus[sharedFanin[s].second]) continue;
      sfId.push_back(s);
      A.insert(A.end(),A.begin(),A.end());
      B.insert(B.end(),B.begin(),B.end());
      int j=sharedFanin[s].second;
      for (int i=0; i<A.size()/2; ++i)
      {
         if (A[i][j]==2) A[i][j]=0;
         else if (A[i][j]==-2) A[i][j]=-1;
         A[i].insert(A[i].end()-1,M);
      }
      for (int i=A.size()/2; i<A.size(); ++i)
      {
         if (A[i][j]==2) A[i][j]=0;
         else if (A[i][j]==-2) A[i][j]=-1;
         else if (A[i][j]==1) A[i][j]=0;
         else if (A[i][j]==-1) A[i][j]=0;
         else if (A[i][j]==0 && A[i].back()==-1) A[i][j]=1; //l-leaf
         else if (A[i][j]==0 && A[i].back()==1) A[i][j]=-1; //0-leaf
         else assert(0);
         A[i].insert(A[i].end()-1,-1*M);
         B[i]-=M;
      }
   }
#if 0
   for (int i=0; i<A.size(); ++i)
   {
      for (int j=0; j<A[i].size(); ++j)
         cout<<setw(2)<<A[i][j]<<" ";
      cout<<"| "<<B[i]<<endl;
   }cout<<endl;
#endif
   int res = solveLp(A,B,ans,sfId.size());
   if (res==0)
   {
      for (int i=0; i<sfId.size(); ++i)
      {
         if (ans[A[0].size()-sfId.size()-1]==1)
         {
            int s=sfId[i];
            if (minus[sharedFanin[s].second]==0) minus[sharedFanin[s].second]=1;
            else minus[sharedFanin[s].second]=0;
         }
         ans.erase(ans.begin()+(A[0].size()-sfId.size()-1));
      }
   }
   return res;
}


#endif
