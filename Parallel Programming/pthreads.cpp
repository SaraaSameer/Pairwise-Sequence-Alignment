#include<string>
#include<algorithm>
#include<iostream>
#include<string.h>
#include<pthread.h>
using namespace std;

//Global variables --shared by all threads
int match_penalty;
int mismatch_penalty;   
int gap_penalty;
int minimum_penalty;

struct gene_struct {
    string gene1;
    string gene2;
    int rowSize;
	int *gene1Result,*gene2Result;
    int colSize;
    int **mat;
};

void* set_matrix_row(void *args)
{
   struct gene_struct *shared_block = (struct gene_struct *) args;
   for (int i = 0; i < shared_block->rowSize; i++)
   {
	shared_block->mat[i][0] = i * gap_penalty;
   }
   
   
   return NULL;
}

void *set_matrix_column(void *args)
{
   struct gene_struct *shared_block = (struct gene_struct *) args;
   for (int i = 0; i < shared_block->colSize; i++)
   {
	shared_block->mat[0][i] = i * gap_penalty;
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
		 		shared_block.mat[i][j] = shared_block.mat[i - 1][j - 1]+match_penalty;
	       }
	       else
			{
		  	     shared_block.mat[i][j] = max({ shared_block.mat[i - 1][j - 1] + mismatch_penalty , shared_block.mat[i - 1][j] + gap_penalty, shared_block.mat[i][j - 1] + gap_penalty });
			}
	   }   
   } 
}
void  Final_Resultant_Strings(struct gene_struct shared_block){
		int lenGene2 = shared_block.rowSize - 1;		//Length of gene2
		int lenGene1 = shared_block.colSize - 1;		//length of gene 1

		int maxLength = lenGene2 + lenGene1;

		int i = lenGene1;								//backtracking from the last value of matrix
		int j = lenGene2;
	

		int xpos = maxLength;		
		int ypos = maxLength;

		//RULES:
		//1) IF A MATCH THEN GO DIAGNOL
		//2) IF A MISMATCH GO TO NEAREST BLOCK VALUE.
		shared_block.gene1Result = new int[maxLength + 1];
		shared_block.gene2Result = new int[maxLength + 1];

		while (!(i == 0 || j == 0)) {
			if (shared_block.gene1[i - 1] == shared_block.gene2[j - 1]) {//If it is match then go diagnol.

				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
				
				i--;
				j--;
				
                minimum_penalty= minimum_penalty + match_penalty;
                
			}
			//IF LETTERS ARE NOT SAME THEN  WE NEED TO FIND THE HIGHEST VALUE OF NEIGHBOURS ( ISKO KARNE KELIYE JO UPAR SETMARIX() MAI MINUS KARA THA WOH PLUS KARDO

			else if (shared_block.mat[i - 1][j - 1] + mismatch_penalty == shared_block.mat[i][j]) {	//Diagnol case.

				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
			
				
				i--; j--;
				
				minimum_penalty= minimum_penalty + mismatch_penalty;
			}
			else if (shared_block.mat[i - 1][j] + gap_penalty == shared_block.mat[i][j]) {	//FOR LEFT CASE
			
				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];

				shared_block.gene2Result[ypos--] = (int)'_';
				i--;
				
				minimum_penalty= minimum_penalty + gap_penalty;
				
			}
			else if (shared_block.mat[i][j - 1] + gap_penalty == shared_block.mat[i][j]) {	//FOR UP CASE.

				shared_block.gene1Result[xpos--] = (int)'_';

				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];

				j--;
			
				minimum_penalty= minimum_penalty + gap_penalty;
			
			}
			else if (shared_block.mat[i - 1][j - 1] + match_penalty == shared_block.mat[i][j]){
				
			
				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
			
				
				i--; j--;
				minimum_penalty= minimum_penalty + mismatch_penalty;
				
			}

		}
		while (xpos > 0)
		{		

			if (i > 0) {
				
				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[--i];
				minimum_penalty= minimum_penalty + gap_penalty;
				
			}
	
			else shared_block.gene1Result[xpos--] = (int)'_';	//Filling the starting gaps 
		}
		while (ypos > 0)
		{
			if (j > 0){
				
				 shared_block.gene2Result[ypos--] = (int)shared_block.gene2[--j];
				 minimum_penalty= minimum_penalty + gap_penalty;
			}
			else {
			
			shared_block.gene2Result[ypos--] = (int)'_';	//Filling the starting gaps.
			
			}
		}
		int gapsEncountered = 1;
		for (int i = maxLength; i >= 1; i--)
		{
			if ( (char)shared_block.gene1Result[i] == '_' && (char)shared_block.gene2Result[i] == '_')		//Yeh gapsEncountered ki value waha tak laraha hai jahan tk dashes lagey wai hain in the starting.
			{
				gapsEncountered = i + 1;
				
				break;
			}
		}
		    cout <<endl<< "Minimum Penalty in aligning the genes = ";
	    	cout << minimum_penalty << "\n";
		
		    cout << "The Aligned Genes Are :" << std::endl;


	    for (i = gapsEncountered; i <= maxLength; i++)
	    {
	        cout<<(char)shared_block.gene1Result[i];
	    }
	    cout << "\n";
	    for (i = gapsEncountered; i <= maxLength; i++)
	    {
	        cout << (char)shared_block.gene2Result[i];
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
		cin >>  match_penalty;	
		cout<<"Enter the value of misMatchPenalty Penalty" << endl;
		cin >>  mismatch_penalty;	
		cout<<"Enter the value of gapPenalty Penalty" << endl;
		cin >>  gap_penalty;
	
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
		Final_Resultant_Strings(shared_block);
        Print_Matrix(shared_block);
}

