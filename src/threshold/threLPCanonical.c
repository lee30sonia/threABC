#ifndef THRELPCAN_C
#define THRELPCAN_C

#include "base/abc/abc.h"
#include "threshold.h"
#include "misc/vec/vec.h"
#include "bdd/extrab/extraBdd.h"
#include "bool/kit/kit.h"

#define nBddSizeMax INT_MAX/2

typedef struct Abc_TtStore_t_  Abc_TtStore_t;
struct Abc_TtStore_t_
{
   int               nVars;
   int               nWords;
   int               nFuncs;
   word **           pFuncs;
};

extern Abc_TtStore_t * Abc_TtStoreLoad( char * pFileName, int nVarNum );

void Th_LibraryConstruct(char* FileName, int nVarNum)
{
   Abc_TtStore_t* p = Abc_TtStoreLoad(FileName,-1);
   if ( p==NULL )
   {
      printf( "File load failed.\n" );
      return;
   }
   if (current_TList==0)
      current_TList = Vec_PtrAlloc( p->nFuncs *2);

   p->nVars=nVarNum;
   DdManager * dd = Cudd_Init( p->nVars , 0 , CUDD_UNIQUE_SLOTS , CUDD_CACHE_SLOTS , 0 );
   DdNode* bFunc,*bFunc2;
   int offset,i,j,Entry;
   int n=0;
   for (i=0; i<p->nFuncs; ++i)
   {
      offset = i*p->nWords;
      bFunc = Kit_TruthToBdd(dd, p->pFuncs[offset],p->nVars,0); Cudd_Ref(bFunc);
      //printf("i=%d\n",i);
      //Cudd_PrintMinterm(dd,bFunc);
      Thre_S * tObj = Th_CreateObj( current_TList, Th_Node );
      if (Th_Bdd2th_LP(dd,bFunc,tObj))
      {
          ++n;
          Vec_IntSort(tObj->weights,0);
          int len = Vec_IntSize(tObj->weights);
          printf("[");
          Vec_IntForEachEntry(tObj->weights, Entry, j) {
             if (j<len-1)
               printf("%d, ", Entry);
             else
               printf("%d", Entry);
          }
          printf("; %d]\n", tObj->thre);
      }
      //printf("%d\n",Th_Bdd2th_LP(dd,bFunc,tObj));
/*
      bFunc2 = Cudd_Not(bFunc); Cudd_Ref(bFunc2);
      //printf("%d\n",i);
      //Cudd_PrintMinterm(dd,bFunc);

      tObj = Th_CreateObj( current_TList, Th_Node );
      if (Th_Bdd2th_LP(dd,bFunc2,tObj))
      {
          ++n;
          Vec_IntSort(tObj->weights,0);
          int len = Vec_IntSize(tObj->weights);
          printf("[");
          Vec_IntForEachEntry(tObj->weights, Entry, j) {
             if (j<len-1)
               printf("%d, ", Entry);
             else
               printf("%d", Entry);
          }
          printf("; %d]\n", tObj->thre);
      }
      //printf("%d\n",Th_Bdd2th_LP(dd,bFunc,tObj));
      Cudd_RecursiveDeref(dd,bFunc2);
*/
      Cudd_RecursiveDeref(dd,bFunc);
   }
   //printf("n=%d\n",n);
}
abctime accum;
// for quick canonical checking
void Th_Canonical_Check()
{
   //Thre_S* pEntry; int i;
   //Vec_PtrForEachEntry(Thre_S*, current_TList, pEntry, i)
      //printGate(pEntry);
      accum=0;
   DdManager * dd=Cudd_Init( Vec_IntSize(((Thre_S*)Vec_PtrEntry(current_TList,Vec_PtrSize(current_TList)-1))->weights)  , 0 , CUDD_UNIQUE_SLOTS , CUDD_CACHE_SLOTS , 0 );
   int i;
   int n=Vec_PtrSize(current_TList);
   for (i=1; i<=n; ++i){
      //printf("%d\n",i);
      Thre_S* tObj = Vec_PtrEntry(current_TList,n-i);
      if (tObj->Type != Th_Node && tObj->Type != Th_Po) continue;
      DdNode* bFunc=Th2Bdd(tObj, dd);
      //Thre_S* tObj2 = Th_CreateObj( current_TList, Th_Node );
      if (!Th_Bdd2th_LP(dd,bFunc,0)) printf("solve failed\n");
      //printGate(tObj2);
   }
   Abc_PrintTime( 1 , "solve time" ,  accum );
   Cudd_Quit(dd); 
}

#endif
