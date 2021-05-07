#include<string>
#include<algorithm>
#include<iostream>
#include<string.h>
#include<pthread.h>
using namespace std;

//Global variables --shared by all threads
int matchPenalty;
int misMatchPenalty;   
int gapPenalty;

struct gene_struct {
    string gene1;
    string gene2;
    int rowSize;
    int colSize;
    int **mat;
};

void* set_matrix_row(void *args)
{
   struct gene_struct *shared_block = (struct gene_struct *) args;
   for (int i = 0; i < shared_block->rowSize; i++)
   {
	shared_block->mat[i][0] = i * gapPenalty;
   }
   
   
   return NULL;
}

void *set_matrix_column(void *args)
{
   struct gene_struct *shared_block = (struct gene_struct *) args;
   for (int i = 0; i < shared_block->colSize; i++)
   {
	shared_block->mat[0][i] = i * gapPenalty;
   }
   
   
   return NULL;
}

void set_matrix_remaining(struct gene_struct shared_block)
{
    for (int i = 1; i < shared_block.rowSize; i++)
   {
	for (int j = 1; j < shared_block.colSize; j++)
	{
	   // Diagnol values will be filled serially, all other values will be computed by parallel threads (Long method --Plan B)
	  
	       if (shared_block.gene1[j - 1] == shared_block.gene2[i - 1])    //Similar gene values (A==A ||T==T)
	       {
		 shared_block.mat[i][j] = shared_block.mat[i - 1][j - 1]+matchPenalty;
	       }
	       else
		{
		   shared_block.mat[i][j] = max({ shared_block.mat[i - 1][j - 1] + misMatchPenalty , shared_block.mat[i - 1][j] + gapPenalty, shared_block.mat[i][j - 1] + gapPenalty });
		}
	   }   
   } 
}

void Print_Matrix(struct gene_struct shared_block)
{
	      cout<<"\n";
		for (int i = 0; i < shared_block.rowSize; i++)
		{
			for (int j = 0; j < shared_block.colSize; j++)
			{
				cout << shared_block.mat[i][j] << "\t";
			}
			cout << "\n";
		}
}

// For matrix traceback, split the matrix in two halves. Combine the results obtained by both threads. 
// Compare the output each time with series.cpp by executing the file like this: time ./obj1

int main()
{
        pthread_t ptid1, ptid2;
        string gene1, gene2;
        struct gene_struct shared_block;   
	cout << "Enter First Gene : ";
	cin >> gene1;
	cout << "Enter Second Gene : ";
	cin >> gene2;
	cout<<"Enter the value of Match Penalty" << endl;
	cin >>  matchPenalty;	
	cout<<"Enter the value of misMatchPenalty Penalty" << endl;
	cin >>  misMatchPenalty;	
	cout<<"Enter the value of gapPenalty Penalty" << endl;
	cin >>  gapPenalty;
	
        shared_block.gene1 = gene1;
        shared_block.gene2 = gene2;
        shared_block.rowSize = gene2.length()+1;
        shared_block.colSize = gene1.length()+1; 
        shared_block.mat = new int* [shared_block.rowSize];  //col = gene2+1
        for (int i = 0; i < shared_block.rowSize; i++)
        {
		shared_block.mat[i] = new int[shared_block.colSize];
             
	 }
		
        pthread_create(&ptid1, NULL, &set_matrix_row, (void *)&shared_block);
        pthread_create(&ptid2, NULL ,&set_matrix_column, (void *)&shared_block);
        pthread_join(ptid1,NULL); 
        pthread_join(ptid2,NULL);
        set_matrix_remaining(shared_block);
        Print_Matrix(shared_block);
}

