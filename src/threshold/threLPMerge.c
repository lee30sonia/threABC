#include <stdio.h>
#include <math.h>
#include "base/abc/abc.h"
#include "threshold.h"
#include "misc/vec/vec.h"
#include "test.c"



void
Th_CollapseNtk_LP( Vec_Ptr_t * TList , int fIterative , int fOutBound, int nMax )
{
   Thre_S * tObj;
   int i , j , FinId , sizeBeforeIter , sizeBeforeCollapse;
   //printf( "Th_CollapseNtk() : fIterative = %d\n" , fIterative );
   //printf( "Th_CollapseNtk() : fOutBound  = %d\n" , fOutBound  );
   do {
      Th_UnmarkAllNode();
      sizeBeforeIter = Vec_PtrSize( TList );
      while ( 1 ) {
         sizeBeforeCollapse = Vec_PtrSize( TList );
         Vec_PtrForEachEntry( Thre_S* , TList , tObj , i )
         {
            // Following nodes are skipped:
            if ( !tObj )                    continue; // NULL  node
            if ( tObj->Type != Th_Node )    continue; // PI/PO/CONST
            if ( tObj->nId == globalRef )   continue; // black node : those who have nId = 1
            if ( !Th_ObjNormalCheck(tObj) ) continue; // Abnormal node : const or 0-weight, collect and clean
            
            Vec_IntForEachEntry( tObj->Fanins , FinId , j )
            {
               if ( Th_CollapseNodes_LP( tObj , j , fOutBound ) ) {
                  // delete tObj`s j-fanin and all its fanouts
                  Th_DeleteClpObj( tObj , j );
                  break;
               }
               //printf("(%d) cannot be merged.\n", tObj->Id);
               // non-mergable node-> color = black
               if ( j == Vec_IntSize( tObj->Fanins ) - 1 ) tObj->nId = globalRef;
            }
         if (nMax!=0 && Vec_PtrSize(current_TList) >nMax/*175100*/) break;
            //printf("%d / %d\n", i, sizeBeforeMerge);
         }
         if (nMax!=0 && Vec_PtrSize(current_TList) >nMax) break;
         if ( sizeBeforeCollapse == Vec_PtrSize(TList) ) break;
      }
   } while ( fIterative && Vec_PtrSize( TList ) > sizeBeforeIter );
   
   //printf("merging process completed...\n");
}

int
Th_CollapseNodes_LP( const Thre_S * tObj2 , int nFanin , int fOutBound )
{
   assert( tObj2 );
   assert( nFanin >= 0 && nFanin < Vec_IntSize(tObj2->Fanins) );
   
   Thre_S * tObj1, *tObj3, *tObjMerge;
   Vec_Ptr_t* pVec;
   int i, Entry;
   
        // printGate(tObj2);
         //printf("%d %d\n",tObj2->Id,nFanin);
   tObj1 = Th_GetObjById( current_TList , Vec_IntEntry( tObj2->Fanins, nFanin ) );
   assert( tObj1 );

   if ( tObj1->Type != Th_Node || !Th_CheckMultiFoutCollapse_LP( tObj1 , fOutBound ) )
		return 0;
   //Vec_Ptr_t* vMerged; //**
   //if ( tObj1->Type != Th_Node ) //**
     // return 0; //**
   //else vMerged = Th_CheckMultiFoutCollapse_LP( tObj1 , fOutBound );  //**
   //if (!vMerged) //**
     // return 0; //**
	//else return Th_CalLPCollapse( tObj1, vMerged ); //**
	else return Th_CalLPCollapse( tObj1 );
}

int
//Vec_Ptr_t* //**
Th_CheckMultiFoutCollapse_LP( const Thre_S * tObj1 , int fOutBound )
{
	// controlling multi-fanout number
	//int foutBound = 30;
	if ( fOutBound == -1 ); // -1 --> no limit
	else if ( Vec_IntSize( tObj1->Fanouts ) > fOutBound ) return 0;

	Thre_S * tObj2, *pEntry;
	int RetValue , Entry , i, j;
	RetValue = 1;
   //Vec_Ptr_t* vMerged = Vec_PtrAlloc(8); //**
   Vec_IntForEachEntry( tObj1->Fanouts , Entry , i )
	{
		tObj2 = Th_GetObjById( current_TList , Entry );
		assert(tObj2);
      //Thre_S* tObjMerge = Th_CreateObjNoInsert(Th_Node); //**
      //Vec_PtrPush(vMerged, tObjMerge); //**
		if ( tObj2->nId == globalRef || tObj2->Type != Th_Node || !Th_Collapse_LPBDD(tObj1 , tObj2,/* tObjMerge */0) ) {
         //Vec_PtrForEachEntry(Thre_S*, vMerged, pEntry, j) //**
            //Th_DeleteObjNoInsert(pEntry); //**
		   RetValue = 0;
			break;
         //return 0; //**
		}
      /*else { 
         printGate(tObj1);
         printGate(tObj2);
         printGate(tObjMerge);
         if (Vec_IntSize(tObjMerge->Fanouts)>0) printGate(Vec_PtrEntry(current_TList,Vec_IntEntry(tObjMerge->Fanouts,0))); }*/
	}
	return /*vMerged;*/RetValue;
}

int 
Th_CalLPCollapse( const Thre_S * tObj1 )
//Th_CalLPCollapse( const Thre_S * tObj1, Vec_Ptr_t* vMerged ) //**
{
   //Vec_Int_t* vId = Vec_IntAlloc(8);
   //Vec_Int_t* vId2 = Vec_IntAlloc(8);
	Thre_S * tObj2, * tObjMerge;
	int Entry , i;
   Vec_IntForEachEntry( tObj1->Fanouts , Entry , i )
	{
      tObj2  = Th_GetObjById( current_TList , Entry );
		assert( tObj2 && tObj2->Type == Th_Node );
      tObjMerge = Th_CreateObj( current_TList, Th_Node );
		assert(Th_Collapse_LPBDD(tObj1 , tObj2, tObjMerge));
      Vec_PtrPush(Th_aigMap,Vec_PtrEntry(Th_aigMap,tObj2->Id));
      /*tObjMerge = Vec_PtrEntry(vMerged, i);
      Vec_IntPush(vId, Entry);
      tObjMerge->Id = Vec_PtrSize(current_TList);
      Vec_IntPush(vId2, tObjMerge->Id);
      Vec_PtrPush(current_TList,tObjMerge);
      tObjMerge->Fanouts = Vec_IntDup(tObj2->Fanouts);
      Th_LPPatchFanio(tObj1,tObj2,tObjMerge,vId,vId2);*/
      Th_KLPatchFanio(tObj1,tObj2,tObjMerge);
	}
	return 1;
}

void
Th_LPPatchFanio( const Thre_S * tObj1 , const Thre_S * tObj2 , const Thre_S * tObjMerge, Vec_Int_t* vId, Vec_Int_t* vId2 )
{
	Thre_S * tObjFanin , * tObjFanout;
	int nFanin , Entry , i, oriId;
	// connect fanins , tObj2 fanout part
	Vec_IntForEachEntry( tObj2->Fanouts , Entry , i )
	{
      tObjFanout = Th_GetObjById( current_TList , Entry );
	   assert(tObjFanout);
      if (Vec_IntFind(tObj1->Fanouts, Entry)!= -1)
      {
         oriId = Vec_IntFind(vId2,Entry);
         if (oriId != -1)
         {
            oriId = Vec_IntEntry(vId,oriId);
		      Vec_IntRemove( tObjMerge->Fanouts , oriId );
		      Vec_IntPush  ( tObjMerge->Fanouts , Entry );
         }
      }
		// Unmark a node if some of its fanins are merged
		if ( tObjFanout->nId == globalRef ) --(tObjFanout->nId);
	   nFanin = Th_ObjFanoutFaninNum( tObj2 , tObjFanout );
	   Vec_IntWriteEntry( tObjFanout->Fanins , nFanin , tObjMerge->Id );
	   assert( Vec_IntSize(tObjFanout->Fanins) == Vec_IntSize(tObjFanout->weights) );
	}
	// connect fanouts , tObj1 fanin part 
	Vec_IntForEachEntry( tObj1->Fanins , Entry , i )
	{
      tObjFanin = Th_GetObjById( current_TList , Entry );
	   assert(tObjFanin);
		Vec_IntRemove( tObjFanin->Fanouts , tObj1->Id );
		Vec_IntPush  ( tObjFanin->Fanouts , tObjMerge->Id );
	}
	// connect fanouts , tObj2 fanin part 
	Vec_IntForEachEntry( tObj2->Fanins , Entry , i )
	{
	   if ( Entry == tObj1->Id ) continue;
      tObjFanin = Th_GetObjById( current_TList , Entry );
	   assert(tObjFanin);
      if (Vec_IntFind(tObj1->Fanouts, Entry)!= -1)
	      Vec_IntWriteEntry( tObjMerge->Fanins , i , Entry );
		Vec_IntRemove( tObjFanin->Fanouts , tObj2->Id );
		Vec_IntPushUnique  ( tObjFanin->Fanouts , tObjMerge->Id );
	}
	assert( Vec_IntSize(tObjMerge->Fanins) == Vec_IntSize(tObjMerge->weights) );
}
