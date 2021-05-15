#include<iostream>
#include<algorithm>
#include<string.h>
#include<string>
#include<fstream>
#include<omp.h>

using namespace std;

int matchPenalty;
int misMatchPenalty;   
int gapPenalty;
int minimum_penalty;

struct gene_struct {
    string gene1;
    string gene2;
    int rowSize;
    int colSize;
    int **mat;
    int *gene1Result;
    int *gene2Result;
};

void Set_Matrix(struct gene_struct shared_block){

	#pragma omp parallel num_threads(shared_block.rowSize)
	#pragma omp for  
	for (int i = 0 ; i < shared_block.rowSize ; i++)
	{
		shared_block.mat[i][0] = i * gapPenalty;
	}
	
	#pragma omp parallel num_threads(shared_block.colSize)
	#pragma omp for 
	for (int i = 0 ; i < shared_block.colSize ; i++)
	{
		shared_block.mat[0][i] = i * gapPenalty;
	}

	for (int i = 1 ; i < shared_block.rowSize ; i++)
	{
		#pragma omp parallel num_threads(shared_block.colSize)
		#pragma omp for 
		for (int j = 1; j < shared_block.colSize ; j++)
		{
			//cout << endl << "Thread: " << omp_get_thread_num() << endl;
			if (shared_block.gene1[j - 1] == shared_block.gene2[i - 1])    //Similar gene values (A==A ||T==T)
			{
				shared_block.mat[i][j] = shared_block.mat[i - 1][j - 1] + matchPenalty;
			}
			else
			{
				shared_block.mat[i][j] = max({ shared_block.mat[i - 1][j - 1] + misMatchPenalty , shared_block.mat[i - 1][j] + gapPenalty, shared_block.mat[i][j - 1] + gapPenalty });
			}
		}
	}
}
void Final_Resultant_Strings(struct gene_struct shared_block) 
{
	int lenGene2 = shared_block.rowSize - 1;		
	int lenGene1 = shared_block.colSize - 1;	

	int maxLength = lenGene2 + lenGene1;

	int i = lenGene1;								
	int j = lenGene2 ;

	int ypos = maxLength ;
	int xpos = maxLength;
				
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
        		minimum_penalty= minimum_penalty + matchPenalty;
		}
		else if (shared_block.mat[i - 1][j - 1] + misMatchPenalty == shared_block.mat[i][j]) 
		{
			shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
			shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
			i--; 
			j--;
			minimum_penalty= minimum_penalty + misMatchPenalty;
		}
		else if (shared_block.mat[i - 1][j] + gapPenalty == shared_block.mat[i][j]) 
		{	
			shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
			shared_block.gene2Result[ypos--] = (int)'_';
			i--;
			minimum_penalty= minimum_penalty + gapPenalty;
		}
		else if (shared_block.mat[i][j - 1] + gapPenalty == shared_block.mat[i][j]) 
		{	
			shared_block.gene1Result[xpos--] = (int)'_';
			shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
			j--;
			minimum_penalty= minimum_penalty + gapPenalty;
		}
		else if (shared_block.mat[i - 1][j - 1] + matchPenalty == shared_block.mat[i][j])
		{
			shared_block.gene1Result[xpos--] = (int)shared_block.gene1[i - 1];
			shared_block.gene2Result[ypos--] = (int)shared_block.gene2[j - 1];
			i--; 
			j--;
			minimum_penalty= minimum_penalty + misMatchPenalty;
		}
	}
	
	#pragma omp parallel num_threads(maxLength)
	{
		#pragma omp for 	
		
		for (int x = xpos ; x > 0 ; x--)
		{		
			if (i > 0) 
			{
				shared_block.gene1Result[x] = (int)shared_block.gene1[--i];
				minimum_penalty= minimum_penalty + gapPenalty;
			}

			else 
			{
				shared_block.gene1Result[x] = (int)'_';	//Filling the starting gaps 
			}
		}
	}
	#pragma omp parallel num_threads(maxLength)
	{
		#pragma omp for 
	
		for (int y = ypos ; y > 0 ; y--)
		{
			if (j > 0)
			{	
				 shared_block.gene2Result[y] = (int)shared_block.gene2[--j];
				 minimum_penalty= minimum_penalty + gapPenalty;
			}
			else 
			{
				shared_block.gene2Result[y] = (int)'_';
			
			}
		}
	}
	int gapsEncountered = 1;
	 
	for (int i = maxLength; i >= 1; i--)
	{
		if ((char)shared_block.gene1Result[i] == '_' && (char)shared_block.gene2Result[i] == '_')
		{
			gapsEncountered = i + 1;
			
			break;
		}
	}

	cout << endl << "Minimum Penalty in aligning the genes = " << minimum_penalty << endl;
	cout << "The Aligned Genes Are : \n";
	
	for (i = gapsEncountered; i <= maxLength; i++)
	{
		cout << (char)shared_block.gene1Result[i];
	}
	
	cout << endl;
	
	for (i = gapsEncountered; i <= maxLength; i++)
	{
		cout << (char)shared_block.gene2Result[i];
	}
	
	cout << endl;		
}	
void Print_Matrix(struct gene_struct shared_block)
{
	cout << endl;
	
	for (int i = 0; i < shared_block.rowSize; i++)
	{
		for (int j = 0; j < shared_block.colSize; j++)
		{
			cout << shared_block.mat[i][j] << "\t";
		}
		cout << "\n";
	}
}

int main(){
	
	double start; 
	double end; 
	start = omp_get_wtime();
	
	string gene1, gene2;
	struct gene_struct shared_block;
	
	fstream fin("Input.txt");
	fin >> gene1;
	fin >> gene2;
	fin >> matchPenalty;
	fin >> misMatchPenalty;
	fin >> gapPenalty;
	
	shared_block.gene1 = gene1;
	shared_block.gene2 = gene2;
	shared_block.rowSize = gene2.length()+1;
	shared_block.colSize = gene1.length()+1;

	shared_block.mat = new int* [shared_block.rowSize];
	
	#pragma omp parallel num_threads(shared_block.rowSize)
	{
		//printf("Number of Threads: %d\n" , omp_get_thread_num());
		#pragma omp for 
		for (int i = 0; i < shared_block.rowSize; i++)
		{
			shared_block.mat[i] = new int[shared_block.colSize];   
		}	
 	}
	Set_Matrix(shared_block);
	Final_Resultant_Strings(shared_block);
	Print_Matrix(shared_block);
	
	end = omp_get_wtime(); 
	cout << endl << "Work took " << (end-start) << "seconds" << endl;	
	return 0;
}

