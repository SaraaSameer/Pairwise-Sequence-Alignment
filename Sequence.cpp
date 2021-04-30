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
	 void FinalResStrings(){
	 	int lenGene2=x-1;		//doing minus 1 because upar length mai +1 kiya wa tha.
	 	int lenGene1=y-1;	//length of gene 1
	 	
	 	int maxLength=lenGene2+lenGene1;
	 	
	 	int i=lenGene1;								//backtracking from the last value of matrix
	 	int j=lenGene2;
//		cout<<"LENGTH OF GENE 1 : " <<lenGene1;

		int xpos=maxLength;	//maxlength =7
		int ypos=maxLength;
		
		//RULES:
		//1) IF A MATCH THEN GO DIAGNOL
		//2) IF A MISMATCH GO TO NEAREST BLOCK VALUE.
		
		
		int gene1Ans[maxLength+1],gene2Ans[maxLength+1];								//Final Answers!
		
//		cout<<"WOW"<<

		while( !(i == 0 || j == 0) ){								//Stopping condition; stop the search when it reaches matrix[0][0]
		
			if( this->gene1[i-1] == this->gene2[j-1] ){								//If it is match then go diagnol.
				
				gene1Ans[xpos--]=(int)this->gene1[i-1];
				gene2Ans[ypos--]=(int)this->gene2[i-1];
	
				i--; 
				j--;
			}
			//IF LETTERS ARE NOT SAME THEN  WE NEED TO FIND THE HIGHEST VALUE OF NEIGHBOURS ( ISKO KARNE KELIYE JO UPAR SETMARIX() MAI MINUS KARA THA WOH PLUS KARDO
			
			else if ( mat[i-1][j-1] + mismatch_penalty == mat[i][j]  ){	//Diagnol case.
			
				gene1Ans[xpos--] = (int)this->gene1[i - 1];
           	 	gene2Ans[ypos--] = (int)this->gene2[j - 1];
            
				i-- ; j-- ;
				
			}
			else if(mat[i-1][j] + gap_penalty == mat[i][j]){	//FOR LEFT CASE
			
				gene1Ans[xpos--]=(int)this->gene1[i-1];
				
				gene2Ans[ypos--]=(int)'_';
				i--;
			}
			else if(mat[i][j-1] + gap_penalty == mat[i][j]){	//FOR UP CASE.
			
				gene1Ans[xpos--]=(int)'_';
				
				gene2Ans[ypos--]=(int)this->gene2[i-1];
				
				j--;
			}
		}
		
		
	 }

};

int main()
{
	string gene1 = "TCG";      //Remove case sensitivity
	string gene2 = "ATCG";

	// intialsing penalties of different types
	int misMatchPenalty = -1;
	int gapPenalty = -2;

	Matrix dp(gene1.length(), gene2.length());
	dp.setGenes(gene1, gene2);
	dp.setPenalty(misMatchPenalty, gapPenalty);
	dp.setMatrix();
	dp.Print_Matrix();
	dp.FinalResStrings();

}
