#include<iostream>
#include<string.h>
#include<string>
#include<fstream>
#include<algorithm>
#include<iomanip>
#include<unistd.h>
using namespace std;

class Matrix
{
private:
	int x;
	int y;
	int** mat;
	string gene1;
	string gene2;
	int match_penalty;
	int mismatch_penalty;
	int gap_penalty;
	int minimum_penalty;
	
public:
	Matrix(int gene1Len, int gene2Len)
	{
		x = gene2Len + 1;			//gene2 length
		y = gene1Len + 1; 	//gene 1 length;
		mat = new int* [x];
		for (int i = 0; i < x; i++)
		{
			mat[i] = new int[y];
		}
		for (int i = 0; i < x; ++i) 
		{
			for (int j = 0; j < y; ++j) 
			{
				mat[i][j] = 0;
			}
		}
		match_penalty = 1;
		mismatch_penalty = 3;
		gap_penalty = 2;
		minimum_penalty=0;
	}

	void Print_Matrix()
	{
		for (int i = 0; i < x; i++)
		{
			for (int j = 0; j < y; j++)
			{
				cout << mat[i][j] << "\t";
			}
			cout << "\n";
		}
	}

	void setGenes(string gene1, string gene2)
	{
		this->gene1 = gene1;
		this->gene2 = gene2;
	}

	void setPenalty(int mismatch, int gp,int match)
	{
		mismatch_penalty = mismatch;
		gap_penalty = gp;
		match_penalty=match;
	}

	void setMatrix()
	{
		for (int i = 0; i < x; i++)
		{
			mat[i][0] = i * gap_penalty;
		}

		for (int i = 0; i < y; i++)
		{
			mat[0][i] = i * gap_penalty;
		}
		for (int i = 1; i < x; i++)
		{
			for (int j = 1; j < y; j++)
			{
				if (gene1[j - 1] == gene2[i - 1])    //Similar gene values (A==A ||T==T)
				{
					mat[i][j] = mat[i - 1][j - 1]+match_penalty;
				}
				else
				{
					mat[i][j] = max({ mat[i - 1][j - 1] + mismatch_penalty , mat[i - 1][j] + gap_penalty, mat[i][j - 1] + gap_penalty });
				}
			}
		}
	}
	
	void FinalResStrings() 
	{
		int lenGene2 = x - 1;		
		int lenGene1 = y - 1;	
		int maxLength = lenGene2 + lenGene1;
		int i = lenGene1;								
		int j = lenGene2;
		int xpos = maxLength;
		int ypos = maxLength;
						
		int* gene1Ans = new int[maxLength + 1];
		int* gene2Ans = new int[maxLength + 1];

		while (!(i == 0 || j == 0)) 
		{						
			if (this->gene1[i - 1] == this->gene2[j - 1]) 
			{		
				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
				i--;
				j--;
                		minimum_penalty= minimum_penalty + match_penalty;
			}
			else if (mat[i - 1][j - 1] + mismatch_penalty == mat[i][j]) 
			{	
				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
				i--; 
				j--;
				minimum_penalty= minimum_penalty + mismatch_penalty;
			}
			else if (mat[i - 1][j] + gap_penalty == mat[i][j]) 
			{
				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)'_';
				i--;
				minimum_penalty= minimum_penalty + gap_penalty;
			}
			else if (mat[i][j - 1] + gap_penalty == mat[i][j]) 
			{
				gene1Ans[xpos--] = (int)'_';
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
				j--;
				minimum_penalty= minimum_penalty + gap_penalty;
			}
			else if (mat[i - 1][j - 1] + match_penalty == mat[i][j])
			{
				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
				i--; 
				j--;
				minimum_penalty= minimum_penalty + mismatch_penalty;
			}
		}
			
		while (xpos > 0)
		{		
			if (i > 0) 
			{	
				gene1Ans[xpos--] = (int)gene1[--i];
				minimum_penalty= minimum_penalty + gap_penalty;	
			}
	
			else gene1Ans[xpos--] = (int)'_';	//Filling the starting gaps 
		}
		
		while (ypos > 0)
		{
			if (j > 0)
			{	
				gene2Ans[ypos--] = (int)gene2[--j];
				minimum_penalty= minimum_penalty + gap_penalty;
			}
			else 
			{
				gene2Ans[ypos--] = (int)'_';		
			}
		}
	
		int gapsEncountered = 1;
		for (int i = maxLength; i >= 1; i--)
		{
			if ( (char)gene1Ans[i] == '_' && (char)gene2Ans[i] == '_')		
			{
				gapsEncountered = i + 1;
			
				break;
			}
		}
	
		cout<<endl<<"Step02: Deducing the alignment by tracing back the scoring matrix"<<endl<<endl;
		cout << "The Aligned Genes Are :" << endl;
		for (i = gapsEncountered; i <= maxLength; i++)
		{
			cout<<(char)gene1Ans[i];
		}
		cout << endl;
		for (i = gapsEncountered; i <= maxLength; i++)
		{
			cout << (char)gene2Ans[i];
		}
		
		cout <<endl<< "Minimum Penalty in aligning the genes = ";
		cout << minimum_penalty << "\n";		
	}	
};

int main()
{
	string gene1;    
	string gene2;
	int matchPenalty;
	int misMatchPenalty;
	int gapPenalty;

	clock_t start = clock();
	cout<<"Reading Input from file..."<<endl<<endl;
	fstream fin("Input.txt");
	fin >> gene1;
	fin >> gene2;
	fin >> matchPenalty;
	fin >> misMatchPenalty;
	fin >> gapPenalty;
	
	cout<<"Gene01: "<<gene1<<endl;
	cout<<"Gene02: "<<gene2<<endl;
	cout<<"Match Penalty: "<<matchPenalty<<endl;
	cout<<"Mismatch Penalty: "<<misMatchPenalty<<endl;
	cout<<"Gap Penalty: "<<gapPenalty<<endl;
	
	Matrix dp(gene1.length(), gene2.length());
	dp.setGenes(gene1, gene2);
	dp.setPenalty(misMatchPenalty, gapPenalty,matchPenalty);
	dp.setMatrix();
	cout<<endl<<endl;
	cout<<"STEP 01: Designing scoring matrix by calculating penalties"<<endl<<endl;
	dp.Print_Matrix();
	dp.FinalResStrings();
	clock_t end = clock();
	cout << endl << "Program Execution Time: " << setprecision(4)<< (end-start)  << " seconds" << endl;
}
