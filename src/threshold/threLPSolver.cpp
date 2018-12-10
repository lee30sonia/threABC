#ifndef THRELPSOLVER
#define THRELPSOLVER

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "../../lib/lp_solve_5.5.2.5_dev/lp_lib.h"
#include "base/abc/abc.h"

extern abctime accum;
using namespace std;
#define M 1000
#define wRestrict 255

int solveLp(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, int n)
{
   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0, ite=0;
   REAL *row = NULL;
   
   Ncol = A[0].size();
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=n; i<Ncol; ++i)
         set_int(lp,i+1,true);
      for (int i=0; i<n; ++i)
         set_binary(lp,i+1,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
   }
   
   int k=0;
   while ((ret == 0)&&(k<A.size())) {
      j = 0;
      for (int i=0; i<n; ++i)
      {
         if (A[k][Ncol-n+i]==0) continue;
         colno[j]=i+1;
         row[j++]=A[k][Ncol-n+i];
      }
      for (int i=n; i<A[k].size(); ++i)
      {
         if (A[k][i-n]==0) continue;
         colno[j]=i+1;
         row[j++]=A[k][i-n];
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, GE, B[k]))
         ret = 3;
      k++;
   }
#ifdef wRestrict
   k=0;
   while ((ret == 0)&&(k<A[0].size()-1)) {
      j = 0;
      for (int i=0; i<A[0].size()-1; ++i)
      {
         colno[j]=i+1;
         row[j++]=(i==k)?1:0;
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, LE, wRestrict))
         ret = 3;
      k++;
   }
#endif
   k=0;

   if(ret == 0) {
      set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
      
      /* set the objective function */
      j = 0;
      
      for (int i=0; i<A[0].size(); ++i)
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
      ite++;
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
      for(j = n; j < Ncol; j++)
      {
         //printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
         ans.push_back(row[j]);
         //cout<<row[j]<<" ";
         //assert(ans.back()>=0);
      }
      //cout<<endl;
      for(j = 0; j < n; j++)
      {
         ans.push_back(row[j]);
      }
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
   //cout<<ite<<endl; 
   return(ret);
}

int solveLpBlock(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, map<int,int>& symm)
{
   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0;
   REAL *row = NULL;
   REAL *row2 = NULL;
   int opti;

   //Ncol = A[0].size()*3-2; 
   Ncol = A[0].size()*5-4;
   //1~A[0].size: weights & T, 
   //A[0].size+1 ~ 2*A[0].size-1: yi (correspond to i-A[0].size())
   //2*A[0].size ~ 3*A[0].size-2: zi (correspond to i-2*A[0].size+1)
   //3*A[0].size-1 ~ 4*A[0].size-3: yi (correspond to i-3*A[0].size()+2)
   //4*A[0].size-2 ~ 5*A[0].size-4: zi (correspond to i-4*A[0].size+3)
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=1; i<=A[0].size(); ++i)
         set_int(lp,i,true);
      for (int i=A[0].size()+1; i<=Ncol; ++i)
         set_binary(lp,i,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      row2 = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
   }
   
   int k=0;
   while ((ret == 0)&&(k<A.size())) {
      j = 0;
      for (int i=0; i<A[k].size(); ++i)
      {
         if (A[k][i]==0) continue;
         colno[j]=i+1;
         row[j++]=A[k][i];
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, GE, B[k]))
         ret = 3;
      k++;
   }
#ifdef wRestrict
   k=0;
   while ((ret == 0)&&(k<A[0].size()-1)) {
      j = 0;
      for (int i=0; i<A[0].size()-1; ++i)
      {
         colno[j]=i+1;
         row[j++]=(i==k)?1:0;
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, LE, wRestrict))
         ret = 3;
      k++;
   }
#endif
   k=0;
   while (ret==0 && k<A[0].size())
   {
      if (symm[k]!=k)
      {
         //add k>symm[k]
         colno[0]=k+1; colno[1]=symm[k]+1;
         row[0]=1; row[1]=-1;
         if(!add_constraintex(lp, 2, row, colno, GE, 0))
            { ret = 3; break; }
         if (k<A[0].size()-1)
         {
            for (int i=k+1;i<A[0].size();++i)
            {
               if (symm[i]==symm[k])
               {
                  //add i>k
                  colno[0]=i+1; colno[1]=k+1;
                  row[0]=1; row[1]=-1;
                  if(!add_constraintex(lp, 2, row, colno, GE, 0))
                  { ret = 3; break; }
               }
            }
         }
      }
      ++k;
   }

   if(ret == 0) {
      set_add_rowmode(lp, FALSE); 
      
      /* set the objective function */
      j = 0;
      for (int i=0; i<A[0].size(); ++i)
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
      set_verbose(lp,NEUTRAL);// IMPORTANT);
      
      /* Now let lpsolve calculate a solution */
      ret = solve(lp);
      get_variables(lp, row2);
      for (int i=0; i<A[0].size()-1; ++i)
         cout<<row2[i]<<" ";
      cout<<"; "<<row2[A[0].size()-1]<<endl;
      if(ret == OPTIMAL)
         ret = 0;
      else
         ret = 5;
   }

   // add blocking constraint
   if (ret==0) {
      opti = get_objective(lp);
      get_variables(lp, row2);
      set_add_rowmode(lp, TRUE);
      for (int i=1; i<=A[0].size()-1; ++i)
      {
         colno[0]=i; colno[1]=i+A[0].size();
         row[0]=-1; row[1]=M;
         if(!add_constraintex(lp, 2, row, colno, GE, 1-row2[i-1]))
            { ret = 3; break; }
         colno[1]=i-1+2*A[0].size();
         row[0]=1;
         if(!add_constraintex(lp, 2, row, colno, GE, 1+row2[i-1]))
            { ret = 3; break; }
      }
   }
   if (ret==0)
   {
      for (int i=0; i<2*A[0].size()-2; ++i)
      {
         colno[i]=i+A[0].size()+1;
         row[i]=1;
      }
      if(!add_constraintex(lp, 2*A[0].size()-2, row, colno, LE, 2*(A[0].size()-1)-1))
         ret = 3;
   }
   
   if(ret == 0) {
      set_add_rowmode(lp, FALSE);
      /* set the object direction to minimize */
      set_minim(lp);
      ret = solve(lp);
      if(ret == OPTIMAL) {
         if (get_objective(lp)<=opti)
         {
            // report counter-example!!
            get_variables(lp, row2);
            for (int i=0; i<A[0].size()-1; ++i)
               cout<<row2[i]<<" ";
            cout<<"; "<<row2[A[0].size()-1]<<endl;
            ret = 6;
         }
         else ret = 0;
      }
      else
         ret = 5;
   }

   //3*A[0].size-1 ~ 4*A[0].size-3: yi (correspond to i-3*A[0].size()+2)
   //4*A[0].size-2 ~ 5*A[0].size-4: zi (correspond to i-4*A[0].size+3)
   if (ret==6) {
      get_variables(lp, row2);
      set_add_rowmode(lp, TRUE);
      for (int i=1; i<=A[0].size()-1; ++i)
      {
         colno[0]=i; colno[1]=i+3*A[0].size()-2;
         row[0]=-1; row[1]=M;
         if(!add_constraintex(lp, 2, row, colno, GE, 1-row2[i-1]))
            { ret = 3; break; }
         colno[1]=i+4*A[0].size()-3;
         row[0]=1;
         if(!add_constraintex(lp, 2, row, colno, GE, 1+row2[i-1]))
            { ret = 3; break; }
      }
   }
   if (ret==6){
      for (int i=0; i<2*A[0].size()-2; ++i)
      {
         colno[i]=i+3*A[0].size()-1;
         row[i]=1;
      }
      if(!add_constraintex(lp, 2*A[0].size()-2, row, colno, LE, 2*(A[0].size()-1)-1))
         ret = 3;
   }
   if (ret==6){
      set_add_rowmode(lp, FALSE);
      /* set the object direction to minimize */
      set_minim(lp);
      ret = solve(lp);
      if(ret == OPTIMAL) {
         cout<<get_objective(lp)<<endl;
         if (get_objective(lp)<=opti)
         {
            // report counter-example!!
            get_variables(lp, row2);
            for (int i=0; i<A[0].size()-1; ++i)
               cout<<row2[i]<<" ";
            cout<<"; "<<row2[A[0].size()-1]<<endl;
            ret = 6;
         }
         else ret = 0;
      }
      else
         ret = 5;
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

//#define debug
#define sym 0
#define ILP 1
int solveLpCanonical(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, map<int,int>& symm)
{
   lprec *lp2;
   REAL *row2 = NULL;
  /* 
   row2 = (REAL *) malloc(A[0].size()*2 * sizeof(*row2));
   lp2 = read_LP("solve.lp",NORMAL,"test");
      write_LP(lp2,stdout);
         int ret2 = solve(lp2);
   cout<<ret2<<endl;
   get_variables(lp2, row2);
   for(int i=0; i<A[0].size(); i++)
   cout<<row2[i]<<" ";cout<<endl;
   return 0;*/

   abctime clk;
   clk = Abc_Clock();

   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0, fixed_constraint = 0;
   float max_weight;
   REAL *row = NULL;
   int k=0, ite=0;
   set<int> free_weights;
   for (int i=1; i<=A[0].size()-1; ++i)
      free_weights.insert(i);

   Ncol = A[0].size()*2; 
   //1~A[0].size: weights & T, 
   //A[0].size+1: current max weight, 
   //A[0].size+2 ~ 2*A[0].size: yi for weights (correspond to i-A[0].size()-1)
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=1; i<=A[0].size()+1; ++i)
      {
         #if ILP
            set_int(lp,i,true);
         #else
            set_int(lp,i,false);
         #endif
      }
      for (int i=A[0].size()+2; i<=2*A[0].size(); ++i)
         set_binary(lp,i,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      row2 = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
      set_verbose(lp,NEUTRAL);// IMPORTANT);
   }
  
   // set given constraints 
   k=0;
   while ((ret == 0)&&(k<A.size())) {
      j = 0;
      for (int i=0; i<A[k].size(); ++i)
      {
         if (A[k][i]==0) continue;
         #if sym
            if (symm[i]!=i)
            {
               int a;
               for (a=0; a<j; a++)
               {
                  if (colno[a]==symm[i]+1)
                  {
                     row[a] += A[k][i];
                     break;
                  }
               }
               if (a==j)
               {
                  colno[j]=symm[i]+1;
                  row[j++]=A[k][i];
               }
            }
            else
            {
               colno[j]=i+1;
               row[j++]=A[k][i];
            }
         #else
            colno[j]=i+1;
            row[j++]=A[k][i];
         #endif
         
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, GE, B[k]))
         ret = 3;
      k++;
   }

   while ((ret==0)&&(free_weights.size()>0))
   {
      fixed_constraint = get_Nrows(lp);
      //1. minimize max weight
#ifdef debug 
      cout<<"minimize max weight"<<endl;
#endif
      set<int>::iterator it=free_weights.begin();
      while ((ret==0) && (it!=free_weights.end()))
      {
         #if sym
            if (symm[*it-1]!=*it-1) { ++it; continue; }
            colno[0]=symm[*it-1]+1; 
         #else
            colno[0]=*it;
         #endif
         colno[1]=A[0].size()+1;
         row[0]=-1; row[1]=1;
      
         if(!add_constraintex(lp, 2, row, colno, GE, 0))
            ret = 3;
         ++it;
      }
      if (ret==0)
      {
         colno[0]=A[0].size()+1;
         row[0]=1;
         if(!set_obj_fnex(lp, 1, row, colno))
            ret = 4;
      }
      if (ret==0)
      {
         set_add_rowmode(lp, FALSE);
         set_minim(lp);
#ifdef debug 
      write_LP(lp,stdout);
#endif
      
         //clk = Abc_Clock();
         ret = solve(lp);
         //accum += Abc_Clock()-clk;
         //Abc_PrintTime( 3 , "solve time" , Abc_Clock()-clk );
         ite++;
         //cout<<ret<<endl;
         if(ret == OPTIMAL)
            ret = 0;
         else
            ret = 5;
      }

      //2. maximize # of max weight
      if (ret==0)
      {
         get_variables(lp, row2);
         //for (j=0; j<A[0].size()+1; ++j) cout<<row2[j]<<" ";
         //cout<<endl; cout<<get_objective(lp)<<" "; cout<<endl;
#ifdef debug 
      cout<<"maximize # max weight"<<endl;
#endif
         //max_weight=get_objective(lp);
         max_weight=row2[A[0].size()];
         //cout<<"max weight: "<< get_objective(lp)<<endl;

      /*
         set_add_rowmode(lp, TRUE);
         row[0]=1; colno[0]=A[0].size()+1;
         if(!add_constraintex(lp,1,row,colno,EQ,max_weight)) ret=3;
         else 
         {
            set<int>::iterator it=free_weights.begin();
            while ((ret==0) && (it!=free_weights.end()))
            {
               colno[0]=(*it); colno[1]=(*it)+A[0].size()+1;
               row[0]=-1; row[1]=1;
               if(!add_constraintex(lp, 2, row, colno, GE, 1-max_weight))
                  ret = 3;
               else
               {
                  row[1]=max_weight;
                  if(!add_constraintex(lp, 2, row, colno, LE, 0))
                     ret = 3;
               }
               it++;
            }
         }
      }
      if (ret==0)
      {
         j=0;
         for (set<int>::iterator it=free_weights.begin(); it!=free_weights.end(); ++it)
         {
            colno[j]=(*it)+A[0].size()+1;
            row[j++]=1;
         }
      
         if(!set_obj_fnex(lp, j, row, colno))
            ret = 4;
      }
      if (ret==0)
      {
         set_add_rowmode(lp, FALSE);
         set_maxim(lp);
#ifdef debug 
      write_LP(lp,stdout);
#endif
         ret = solve(lp);
         if(ret == OPTIMAL)
            ret = 0;
         else
            ret = 5;
      */
      }


      //3. fix the max weights
      if (ret==0)
      {
#ifdef debug 
      cout<<"fix max weight"<<endl;
#endif
         j=get_Nrows(lp);
         //cout<<j<<" "<<fixed_constraint<<endl;
         while (j>fixed_constraint) del_constraint(lp,j--);
         set_add_rowmode(lp, TRUE);
         get_variables(lp, row2);
         for (j=1; j<=A[0].size()-1; ++j)
         {
         #if sym
            #if ILP
               if (row2[symm[j-1]]==max_weight)
            #else
               if ((row2[symm[j-1]]>=max_weight && row2[symm[j-1]]-max_weight<0.001) || (row2[symm[j-1]]<=max_weight && max_weight-row2[symm[j-1]]<0.001))
            #endif
         #else
            #if ILP
               if (row2[j-1]==max_weight)
            #else
               if ((row2[j-1]>=max_weight && row2[j-1]-max_weight<0.001) || (row2[j-1]<=max_weight && max_weight-row2[j-1]<0.001))
            #endif
         #endif
            {
               set<int>::iterator it = free_weights.find(j);
               if (it!=free_weights.end())
               {
                  free_weights.erase(it);
               #if sym
                  if (symm[j-1]!=j-1) continue;
               #endif
                  colno[0]=j; row[0]=1; 
                  if(!add_constraintex(lp, 1, row, colno, EQ, max_weight))
                  { ret = 3; break; }
               #if sym
                  ;
               #else
                  #if ILP
                     break;
                  #endif
               #endif
               }
            }
         }
#ifdef debug 
      write_LP(lp,stdout);
#endif
      //cout<<"# free weights: "<<free_weights.size()<<endl;
      }
      if (ret==0)
      {
         lp2=copy_lp(lp);
         delete_lp(lp);
         lp=lp2;
      }
   }
   //4. minimize T
   if (ret==0)
   {
#ifdef debug 
      cout<<"minimize T"<<endl;
#endif
      colno[0]=A[0].size();
      row[0]=1;
      if (!set_obj_fnex(lp, 1, row, colno))
         ret=4;
      
   }
   if (ret==0)
   {
      set_add_rowmode(lp, FALSE);
      set_minim(lp);
#ifdef debug 
      write_LP(lp,stdout);
#endif
      //clk = Abc_Clock();
      ret = solve(lp);
      //accum += Abc_Clock()-clk;
      //Abc_PrintTime( 1 , "solve time" , Abc_Clock()-clk );
      ite++;
      if(ret == OPTIMAL)
         ret = 0;
      else
         ret = 5;
   }
   
   if (ret==0)
   {
      get_variables(lp, row);
      for(j = 0; j < A[0].size(); ++j){
      #if sym
         //cout<<row[symm[j]]<<" ";
         ans.push_back(row[symm[j]]);
      #else
         //cout<<row[j]<<" ";
         ans.push_back(row[j]);
      #endif
      }//cout<<endl;
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
   //cout<<"ite: "<<ite<<endl; 
   //accum += Abc_Clock()-clk;
   return(ret);
}

#endif
