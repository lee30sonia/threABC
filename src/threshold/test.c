#include <stdio.h>
#include "base/abc/abc.h"
#include "threshold.h"


void printGate(Thre_S* tObj)
{
   //printf("%d\n",tObj->thre);
   int i=0;
   printf("id:%d\n",tObj->Id);
   printf("fanins: ");
   for (i=0; i<Vec_IntSize(tObj->Fanins); ++i){
   printf("%d/",Vec_IntEntry(tObj->Fanins,i));
   printf("%d ",Vec_IntEntry(tObj->weights,i));}
   printf("\nfanouts: ");
   for (i=0; i<Vec_IntSize(tObj->Fanouts); ++i)
   printf("%d ",Vec_IntEntry(tObj->Fanouts,i));
   printf("\nthreshold: %d\n",tObj->thre);
   //for (i=0; i<Vec_IntSize(tObj->weights); ++i)
   //printf("%d ",Vec_IntEntry(tObj->weights,i));
   //printf("|%d\n",tObj->thre);
}
