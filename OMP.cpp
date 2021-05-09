#include<iostream>
#include<algorithm>
#include<string.h>
#include<string>
#include<omp.h>

using namespace std;

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

void Set_Matrix(struct gene_struct shared_block){

	#pragma omp parallel
	//1st row and 1st Column values
	#pragma omp for  
	for (int i = 0 ; i < shared_block.rowSize ; i++)
	{
		shared_block.mat[i][0] = i * gapPenalty;
	}
	#pragma omp for 
	for (int i = 0 ; i < shared_block.colSize ; i++)
	{
		shared_block.mat[0][i] = i * gapPenalty;
	}

	// Other matrix values

	for (int i = 1 ; i < shared_block.rowSize ; i++)
	{
		#pragma omp for 
		for (int j = 1; j < shared_block.colSize ; j++)
		{
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
void Print_Matrix(struct gene_struct shared_block)
{
	cout<<"\n";
	#pragma omp for 
	for (int i = 0; i < shared_block.rowSize; i++)
	{
		for (int j = 0; j < shared_block.colSize; j++)
		{
			cout << shared_block.mat[i][j] << "\t";
		}
		cout << "\n";
	}
}

int main(int argc, char *argv[]){
	
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
	 
	#pragma omp parallel
	
	printf("Number of Threads: %d\n" , omp_get_thread_num());
	
	shared_block.mat = new int* [shared_block.rowSize];  
	
	#pragma omp for 
	for (int i = 0; i < shared_block.rowSize; i++)
	{
		shared_block.mat[i] = new int[shared_block.colSize];   
	}
	
	
	Set_Matrix(shared_block);
	Print_Matrix(shared_block);
	return 0;
}

