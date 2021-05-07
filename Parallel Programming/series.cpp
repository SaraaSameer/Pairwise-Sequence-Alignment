#include<iostream>
#include<string.h>
#include<string>
#include<algorithm>

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
			mat[i] = new int[y];

		for (int i = 0; i < x; ++i) {
			for (int j = 0; j < y; ++j) {
				mat[i][j] = 0;
			}
		}

		//Default Penalties
		match_penalty = 1;
		mismatch_penalty = 3;
		gap_penalty = 2;
		minimum_penalty=0;
	}

	void Print_Matrix()
	{
	      cout<<"\n";
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
		//1st row and 1st Column values
		for (int i = 0; i < x; i++)
		{
			mat[i][0] = i * gap_penalty;
		}

		for (int i = 0; i < y; i++)
		{
			mat[0][i] = i * gap_penalty;
		}

		// Other matrix values

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
	
};

int main()
{
	
	 //string gene1 = "ACCA";
	 //string gene2 =  "CCA";      //Test case from https://www.youtube.com/watch?v=LhpGz5--isw  //Working BUT GIVING CORRECT MP
	 
	string gene1, gene2;

	// intialsing penalties of different types
	int matchPenalty;
	int misMatchPenalty ;
	int gapPenalty ;
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


        
	Matrix dp(gene1.length(), gene2.length());
	dp.setGenes(gene1, gene2);
	dp.setPenalty(misMatchPenalty, gapPenalty,matchPenalty);
	dp.setMatrix();
	dp.Print_Matrix();

}
