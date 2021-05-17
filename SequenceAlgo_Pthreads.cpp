#include<string>
#include<algorithm>
#include<iostream>
#include<string.h>
#include<pthread.h>
#include<fstream>
#include<iomanip>
using namespace std;

//Global variables --shared by all threads
int match_penalty;
int mismatch_penalty;   
int gap_penalty;
int minimum_penalty;

//Multiple arguments passing in threads via struct
struct gene_struct {
    string gene1;
    string gene2;
    int rowSize;
    int *gene1Result,*gene2Result;
    int colSize;
    int **mat;
    int currValue;
};

void *set_matrix_row(void *args)
{
   struct gene_struct *shared_block = (struct gene_struct *) args;
   for (int i = 0; i < shared_block->rowSize; i++)
   {
	shared_block->mat[i][0] = i * gap_penalty;
   }
   
   pthread_exit(0);
}

void *set_matrix_column(void *args)
{
   struct gene_struct *shared_block = (struct gene_struct *) args;
   for (int i = 0; i < shared_block->colSize; i++)
   {
	shared_block->mat[0][i] = i * gap_penalty;
   }
   
   pthread_exit(0);
}

void *set_matrix_remaining_column(void* args)
{
    struct gene_struct *shared_block = (struct gene_struct *) args;
    for(int i = shared_block->currValue+1; i < shared_block->colSize; i++)
    {
      if (shared_block->gene1[i- 1] == shared_block->gene2[shared_block->currValue  - 1])    //Similar gene values (A==A ||T==T)
	{
	   shared_block->mat[shared_block->currValue][i] = shared_block->mat[shared_block->currValue - 1][i - 1]+match_penalty;
	}
	
	else
	{
	   shared_block->mat[shared_block->currValue][i] = max({ shared_block->mat[shared_block->currValue - 1][i - 1] + mismatch_penalty , shared_block->mat[shared_block->currValue -1][i] + gap_penalty, shared_block->mat[shared_block->currValue][i - 1] + gap_penalty });
        }
    }
    
    pthread_exit(0);
}


void *set_matrix_remaining_row(void* args)
{
    struct gene_struct *shared_block = (struct gene_struct *) args;
    for(int i = shared_block->currValue+1; i < shared_block->rowSize; i++)
    {
       if(shared_block->gene1[shared_block->currValue- 1] == shared_block->gene2[i - 1])   
	{
	   shared_block->mat[i][shared_block->currValue] = shared_block->mat[i-1][shared_block->currValue - 1]+match_penalty;
	}
	
	else
	{
	   shared_block->mat[i][shared_block->currValue] = max({ shared_block->mat[i-1][shared_block->currValue - 1] + mismatch_penalty , shared_block->mat[i][shared_block->currValue - 1] + gap_penalty, shared_block->mat[i-1][shared_block->currValue] + gap_penalty });
        }
    }
    pthread_exit(0);
}

void set_matrix_diagnol(struct gene_struct shared_block)
{
    pthread_t th1, th2;
    int minimum = min(shared_block.rowSize , shared_block.colSize);
    for (int i = 1; i < minimum; i++)
    {
	if (shared_block.gene1[i - 1] == shared_block.gene2[i - 1])    //Similar gene values (A==A ||T==T)
	{
	    shared_block.mat[i][i] = shared_block.mat[i - 1][i - 1]+match_penalty;
	}
	else
	{
	    shared_block.mat[i][i] = max({ shared_block.mat[i - 1][i - 1] + mismatch_penalty , shared_block.mat[i - 1][i] + gap_penalty, shared_block.mat[i][i - 1] + gap_penalty });
        }
     
     shared_block.currValue = i;
     pthread_create(&th1, NULL, &set_matrix_remaining_row, (void *)&shared_block);
     pthread_create(&th2, NULL ,&set_matrix_remaining_column, (void *)&shared_block);
     pthread_join(th1,NULL); 
     pthread_join(th2,NULL); 
     
     }  
}
void  Final_Resultant_Strings(struct gene_struct shared_block){
		int lenGene2 = shared_block.rowSize - 1;		//Length of gene2
		int lenGene1 = shared_block.colSize - 1;		//length of gene 1
		int maxLength = lenGene2 + lenGene1;
		int i = lenGene1;							
		int j = lenGene2;
		int xpos = maxLength;		
		int ypos = maxLength;
		
		shared_block.gene1Result = new int[maxLength + 1];
		shared_block.gene2Result = new int[maxLength + 1];

		while (!(i == 0 || j == 0)) 
		{
			if (shared_block.gene1[i - 1] == shared_block.gene2[j - 1]) 
			       {

				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
				i--;
				j--;
                               minimum_penalty= minimum_penalty + match_penalty;
			        }

			 else if (shared_block.mat[i - 1][j - 1] + mismatch_penalty == shared_block.mat[i][j]) 
			       {	

				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
				i--; j--;
				minimum_penalty= minimum_penalty + mismatch_penalty;
				
			       }
			       
			 else if (shared_block.mat[i - 1][j] + gap_penalty == shared_block.mat[i][j]) 
			      {	//FOR LEFT CASE
			
				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)'_';
				i--;
				minimum_penalty= minimum_penalty + gap_penalty;
				
			      }
			      
			 else if (shared_block.mat[i][j - 1] + gap_penalty == shared_block.mat[i][j]) 
			      {	//FOR UP CASE.

				shared_block.gene1Result[xpos--] = (int)'_';
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
				j--;
				minimum_penalty= minimum_penalty + gap_penalty;
			
			     }
			     
			 else if (shared_block.mat[i - 1][j - 1] + match_penalty == shared_block.mat[i][j])
			     {
							
				shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
				shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];				
				i--; j--;
				minimum_penalty= minimum_penalty + mismatch_penalty;
			    
			    }
                    }
                    
		 while (xpos > 0)
		    {		

			if (i > 0) 
			{
				
			  shared_block.gene1Result[xpos--] = (int)shared_block.gene1[--i];
			  minimum_penalty= minimum_penalty + gap_penalty;
				
			}
	
			else shared_block.gene1Result[xpos--] = (int)'_';	//Filling the starting gaps 
			
		    }
		    
		while (ypos > 0)
		{
			
			if (j > 0)
			{
			
			  shared_block.gene2Result[ypos--] = (int)shared_block.gene2[--j];
			  minimum_penalty= minimum_penalty + gap_penalty;
			}
			else 
			{
			
			  shared_block.gene2Result[ypos--] = (int)'_';	//Filling the starting gaps.
			
			}
		}
		
		int gapsEncountered = 1;
		for (int i = maxLength; i >= 1; i--)
		{
			if ( (char)shared_block.gene1Result[i] == '_' && (char)shared_block.gene2Result[i] == '_')		
			{
				gapsEncountered = i + 1;
				break;
			}
		}
		
	    	cout<<endl<<"Step02: Deducing the alignment by tracing back the scoring matrix"<<endl<<endl;
		cout << "The Aligned Genes Are :" << endl;
	        for (i = gapsEncountered; i <= maxLength; i++)
	        {
	          cout<<(char)shared_block.gene1Result[i];
	        }
	        cout << "\n";
	        for (i = gapsEncountered; i <= maxLength; i++)
	        {
	          cout << (char)shared_block.gene2Result[i];
	        }
                
               cout <<endl<< "Minimum Penalty in aligning the genes = ";
	    	cout << minimum_penalty << "\n";
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

// Compare the output each time with series.cpp by executing the file like this: time ./obj1

int main()
{
        pthread_t ptid1, ptid2;
        string gene1, gene2;
        struct gene_struct shared_block; 
        cout<<"Reading Input from file..."<<endl<<endl;  
	fstream fin("Input.txt");
	fin>>gene1;
	fin>>gene2;
	fin>>match_penalty;
	fin>>mismatch_penalty;
	fin>>gap_penalty;
	
	
	cout<<"Gene01: "<<gene1<<endl;
	cout<<"Gene02: "<<gene2<<endl;
	cout<<"Match Penalty: "<<match_penalty<<endl;
	cout<<"Mismatch Penalty: "<<mismatch_penalty<<endl;
	cout<<"Gap Penalty: "<<gap_penalty<<endl;
	
        shared_block.gene1 = gene1;
        shared_block.gene2 = gene2;
        shared_block.rowSize = gene2.length()+1;
        shared_block.colSize = gene1.length()+1; 
        shared_block.mat = new int* [shared_block.rowSize];  //col = gene2+1
        for (int i = 0; i < shared_block.rowSize; i++)
        {
		shared_block.mat[i] = new int[shared_block.colSize];
             
	 }
	 
	clock_t start = clock();	
        pthread_create(&ptid1, NULL, &set_matrix_row, (void *)&shared_block);
        pthread_create(&ptid2, NULL ,&set_matrix_column, (void *)&shared_block);
        pthread_join(ptid1,NULL); 
        pthread_join(ptid2,NULL); 
        set_matrix_diagnol(shared_block);
	cout<<endl<<endl;
	cout<<"STEP 01: Designing scoring matrix by calculating penalties"<<endl<<endl;
        Print_Matrix(shared_block);
        Final_Resultant_Strings(shared_block);
        clock_t end = clock();
        cout << endl << "Program Execution Time: " << setprecision(4)<< double((end-start)/double(CLOCKS_PER_SEC))  << " seconds" << endl;
        fstream fout("Output_Time_Pthread.txt",ios::out);
        fout << double((end-start)/double(CLOCKS_PER_SEC));
        fout.close();
}

