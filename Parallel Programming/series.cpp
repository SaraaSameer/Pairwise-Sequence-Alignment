#include<iostream>
#include<string.h>
#include<string>
#include<algorithm>

using namespace std;


//https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-3/#nwRR       (Reference)

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
	void FinalResStrings() {
		int lenGene2 = x - 1;		//doing minus 1 because upar length mai +1 kiya wa tha.
		int lenGene1 = y - 1;	//length of gene 1

		int maxLength = lenGene2 + lenGene1;

		int i = lenGene1;								//backtracking from the last value of matrix
		int j = lenGene2;
	

		int xpos = maxLength;	//maxlength =7
		int ypos = maxLength;

		//RULES:
		//1) IF A MATCH THEN GO DIAGNOL
		//2) IF A MISMATCH GO TO NEAREST BLOCK VALUE.


							//Final Answers!
		int* gene1Ans = new int[maxLength + 1];
		int* gene2Ans = new int[maxLength + 1];


		while (!(i == 0 || j == 0)) {								//Stopping condition; stop the search when it reaches matrix[0][0]
			if (this->gene1[i - 1] == this->gene2[j - 1]) {								//If it is match then go diagnol.

				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
				
				i--;
				j--;
				
                minimum_penalty= minimum_penalty + match_penalty;
               
			}
			
			//IF LETTERS ARE NOT SAME THEN  WE NEED TO FIND THE HIGHEST VALUE OF NEIGHBOURS ( ISKO KARNE KELIYE JO UPAR SETMARIX() MAI MINUS KARA THA WOH PLUS KARDO

			else if (mat[i - 1][j - 1] + mismatch_penalty == mat[i][j]) {	//Diagnol case.

				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
			
				
				i--; j--;
				//cout << "Going diagnolly to mat [ " << i << " ]" << "[ " << j << " ]" << endl;
				minimum_penalty= minimum_penalty + mismatch_penalty;
			}
			else if (mat[i - 1][j] + gap_penalty == mat[i][j]) {	//FOR LEFT CASE
			
				gene1Ans[xpos--] = (int)this->gene1[i - 1];

				gene2Ans[ypos--] = (int)'_';
				i--;
				//cout << "Going diagnolly to mat [ " << i << " ]" << "[ " << j << " ]" << endl;
				minimum_penalty= minimum_penalty + gap_penalty;
				
			}
			else if (mat[i][j - 1] + gap_penalty == mat[i][j]) {	//FOR UP CASE.

				gene1Ans[xpos--] = (int)'_';

				gene2Ans[ypos--] = (int)this->gene2[j - 1];

				j--;
				//cout << "Going diagnolly to mat [ " << i << " ]" << "[ " << j << " ]" << endl;
				minimum_penalty= minimum_penalty + gap_penalty;
			
			}
			else if (mat[i - 1][j - 1] + match_penalty == mat[i][j]){
				
			
				gene1Ans[xpos--] = (int)this->gene1[i - 1];
				gene2Ans[ypos--] = (int)this->gene2[j - 1];
			
				
				i--; j--;
				minimum_penalty= minimum_penalty + mismatch_penalty;
				//cout << "Going diagnolly to mat [ " << i << " ]" << "[ " << j << " ]" << endl;
			}
		}
			
		while (xpos > 0)
		{		

			if (i > 0) {
				
				gene1Ans[xpos--] = (int)gene1[--i];
				minimum_penalty= minimum_penalty + gap_penalty;
				
			}
	
			else gene1Ans[xpos--] = (int)'_';	//Filling the starting gaps 
		}
		while (ypos > 0)
		{
			if (j > 0){
				
				 gene2Ans[ypos--] = (int)gene2[--j];
				 minimum_penalty= minimum_penalty + gap_penalty;
			}
			else {
			
			gene2Ans[ypos--] = (int)'_';	//Filling the starting gaps.
			
			}
		}
	
		int gapsEncountered = 1;
		for (int i = maxLength; i >= 1; i--)
		{
			if ( (char)gene1Ans[i] == '_' && (char)gene2Ans[i] == '_')		//Yeh gapsEncountered ki value waha tak laraha hai jahan tk dashes lagey wai hain in the starting.
			{
				gapsEncountered = i + 1;
				//minimum_penalty= minimum_penalty - gap_penalty;
				break;
			}
		}
	
			cout <<endl<< "Minimum Penalty in aligning the genes = ";
	    	cout << minimum_penalty << "\n";
		
		cout << "The Aligned Genes Are :" << std::endl;


	    for (i = gapsEncountered; i <= maxLength; i++)
	    {
	        cout<<(char)gene1Ans[i];
	    }
	    cout << "\n";
	    for (i = gapsEncountered; i <= maxLength; i++)
	    {
	        cout << (char)gene2Ans[i];
		}
				
	}	
};

int main()
{
	
	string gene1 = "ATCG";      //Working with MP=1 GP=-2 AND MSIMATCH=-1
	string gene2 = "TCG";

    //string gene1 = "ATGCT";     
    //string gene2 = "AGCT";        //Working with MP=2
	
    //string gene1 = "CA";	         //Working
    //string gene2 = "CG";
	
	//string gene1="ATCGT";            //Working
	//string gene2="TGGTG";
	
	//string gene1 = "GAGC";      
	//string gene2 = "GATC";       //TEST CASE FROM THAT ROUGH VIDEO      //Working

	//  string gene1 = "ACGGCTC";
	//  string gene2 =  "ATGGCCTC";      //Test case from https://www.youtube.com/watch?v=LhpGz5--isw  //Working BUT GIVING CORRECT MP
	 
	  //string gene1 = "CCATACGA";
	  //string gene2 = "CAGCTAGCG";       //CORNER CASE ...NOT WORKING
	  
      //string gene1 = "AGGGCT";
      //string gene2 = "AGGCA";         //FAILED TEST CASE

	// intialsing penalties of different types
	int matchPenalty;
	int misMatchPenalty ;
	int gapPenalty ;
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
	dp.FinalResStrings();


}
