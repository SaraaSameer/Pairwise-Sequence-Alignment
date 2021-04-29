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

public:
	Matrix(int x_axis, int y_axis)
	{
		x = y_axis + 1;
		y = x_axis + 1;
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
		mismatch_penalty = -1;
		gap_penalty = -2;
	}

	void Print_Matrix()
	{
		for (int i = 0; i < x; i++)
		{
			for (int j = 0; j < y; j++)
			{
				cout << mat[i][j]<<"\t";
			}
			cout << "\n";
		}
	}

	 void setGenes(string gene1, string gene2)
	 {
		 this->gene1 = gene1;
		 this->gene2 = gene2;
	 }

	 void setPenalty(int mismatch, int gp)
	 {
		 mismatch_penalty = mismatch;
		 gap_penalty = gp;
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
				 if (gene1[i - 1] == gene2[j - 1])    //Similar gene values (A==A ||T==T)
				 {
					 mat[i][j] = mat[i - 1][j - 1] + match_penalty;
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
	string gene1 = "ATCGT";      //Remove case sensitivity
	string gene2 = "TGGTG";

	// intialsing penalties of different types
	int misMatchPenalty = -1;
	int gapPenalty = -2;

	Matrix dp(gene1.length(), gene2.length());
	dp.setGenes(gene1, gene2);
	dp.setPenalty(misMatchPenalty, gapPenalty);
	dp.setMatrix();
	dp.Print_Matrix();

}