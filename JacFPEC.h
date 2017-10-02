/*

 -------------------------------------------------------------- 
   Written by Parisa Babaheidarian, Basis Decomposition project: 
           Functions and Jacobians for the Speck transform
           Change number 10 in the code to the number of bases
           for other basis transforms. 
    -------------------------------------------------------------- 


*/

//#include<iostream>
//#include<string>
#include<vector>
#include<math.h>
#include <stdio.h> 

using namespace std;

class Essentials {
public:
	 vector<vector<double> > alpha;
	 vector<vector<double> > Wbrute;
	 vector<vector<double> >WeightMatrix;
	 double I0;
	Essentials(double x[], double B[2][100], double Weight [10][100]) {
		I0 = 2 * pow(10, 5);
		vector<double> tmp(2, 0);
		alpha.resize(1, tmp);
		for (int i = 0; i < 2; i++) {
			alpha[0][i]=x[i];
		}
		vector<double> temp(100, 0);
		Wbrute.resize(2,temp);
		for (int i = 0; i < 2; i++) {
		//	Wbrute.push_back({ B[i][0] });
			for (int j=0;j < 100; j++) {
				Wbrute[i][j]=B[i][j];
			
			}
		}
		WeightMatrix.resize(10, temp);
		for (int i = 0; i < 10; i++) {
			
			for (int j=0; j < 100; j++) {
				WeightMatrix[i][j]=Weight[i][j];

			}
		}
		/**/

	}

	// Computes C=A*B
	
	vector<vector<double>> Multi(vector<vector<double> > A, vector<vector<double> > B) {
	int M = A.size(); int N = B[0].size(); int P = B.size();
	vector<vector<double> > result(M, vector<double>(N, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			result[i][j] = 0;
			for (int k = 0; k < P; k++) {
				result[i][j] = result[i][j] + A[i][k] * B[k][j];
			}
		}


	}
	return result;
}

// Computes C=A'*B
vector<vector<double> > TMulti(vector<vector<double> > A, vector<vector<double> > B) {
	int M = A[0].size(); int N = B[0].size(); int P = B.size();
	vector<vector<double> > result(M, vector<double>(N, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			result[i][j] = 0;
			for (int k = 0; k < P; k++) {
				result[i][j] = result[i][j] + A[k][i] * B[k][j];
			}
		}


	}
	return result;
}

// Computes C=A*B';
vector<vector<double>> MultiT(vector<vector<double> > A, vector<vector<double> > B) {
	int M = A.size(); int N = B.size(); int P = B[0].size();
	vector<vector<double> > result(M, vector<double>(N, 0));
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			result[i][j] = 0;
			for (int k = 0; k < P; k++) {
				result[i][j] = result[i][j] + A[i][k] * B[j][k];
			}
		}


	}
	
	return result;
}


vector<vector<double>> X(vector<vector<double> > alpha, vector<vector<double> >  Wbrute, double I0, vector<vector<double> > WeightMatrix) {
	vector<vector<double>> temp1 = Multi(alpha, Wbrute);

	vector<vector<double> >  temp2(1, vector<double>(100, 0));
	for (int i = 0; i<100;i++) {
		temp2[0][i] = exp(-temp1[0][i]);
	}
	vector<vector<double>> temp3 = MultiT(temp2, WeightMatrix);
	for (int i = 0; i<10;i++) {
		temp3[0][i] = temp3[0][i] / (I0*0.1);
	}
	return temp3;
}

vector<vector<double> > JacX(vector<vector<double> > alpha, vector<vector<double> > Wbrute, double I0, vector<vector<double> >WeightMatrix)
{
	vector<vector<double> > temp1(10, vector<double>(10, 0));
	vector<vector<double> > temp2 = X(alpha, Wbrute, I0, WeightMatrix);

	for (int i = 0; i<10; i++) {
		temp1[i][i] = 1.0 / (temp2[0][i]);
	//	cout << temp2[0][i] << endl;
	}
	return temp1;
}

vector<vector<double>> T(vector<vector<double> > alpha, vector<vector<double> > Wbrute) {
	return Multi(alpha, Wbrute);
}

vector<vector<double> > Jacd(vector<vector<double> >  alpha, vector<vector<double> > Wbrute) {
	vector<vector<double>> temp1 = T(alpha, Wbrute);
	vector<vector<double> > temp2(100, vector<double>(100, 0));
	for (int i = 0; i<100; i++) {
		temp2[i][i] = -1 * exp(-temp1[0][i]);
		//cout << temp2[i][i] << endl;
	}
	return temp2;
}

  double * Grad(
		double bin1, double bin2, double bin3, double bin4, double bin5, double bin6, double bin7, double bin8, double bin9, double bin10)
	{
		vector<vector<double> > Sino(1, vector<double>(10, 0));
		Sino[0][0] = bin1; Sino[0][1] = bin2; Sino[0][2] = bin3; Sino[0][3] = bin4; Sino[0][4] = bin5;
		Sino[0][5] = bin6; Sino[0][6] = bin7; Sino[0][7] = bin8; Sino[0][8] = bin9; Sino[0][9] = bin10;

		vector<vector<double> >  g(1, vector<double>(10, 0));
		vector<vector<double> > temp1(10, vector<double>(2, 0));
        vector<vector<double> > test(100, vector<double>(100, 0));
        test=Jacd(alpha, Wbrute);
		temp1 = Multi(WeightMatrix, MultiT(Jacd(alpha, Wbrute), Wbrute));
        
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 2; j++) {
			temp1[i][j] = temp1[i][j] / (0.1*I0);

		}
	}
	
	temp1 = Multi(JacX(alpha, Wbrute, I0, WeightMatrix),temp1);
   
	vector<vector<double> > temp2 = T(alpha, Wbrute);
	for (int i = 0; i < temp2[0].size(); i++) {
		temp2[0][i] = exp(-temp2[0][i]);

	}
	vector<vector<double> > temp3 = MultiT(temp2, WeightMatrix);
	for (int i = 0; i < 10; i++) {

		temp3[0][i] = log(temp3[0][i] / (0.1*I0));

	}
	for (int i = 0; i < 10; i++) {

		g[0][i] = Sino[0][i] + temp3[0][i];
        
		
	}
    double s=0;
    for (int i = 0; i < 10; i++) {
      s=s+1*temp1[i][0];
    }
	vector<vector<double> > result(1,vector<double> (2,0)); 
    result= Multi(g,temp1);
    static double output[2];
	for (int i = 0; i < 2; i++) {
		output[i] = result[0][i];
	}
	return output;
	}

	double FunMain(
		double bin1, double bin2, double bin3, double bin4, double bin5, double bin6, double bin7, double bin8, double bin9, double bin10)
	{
		vector<vector<double> > Sino(1, vector<double>(10, 0));
		Sino[0][0] = bin1; Sino[0][1] = bin2; Sino[0][2] = bin3; Sino[0][3] = bin4; Sino[0][4] = bin5;
		Sino[0][5] = bin6; Sino[0][6] = bin7; Sino[0][7] = bin8; Sino[0][8] = bin9; Sino[0][9] = bin10;

		vector<vector<double>> g(1, vector<double>(10, 0.));
	vector<vector<double>> temp2;
	temp2 = T(alpha, Wbrute);
	vector<vector<double>> aux(1,vector<double>(100,0.));
	for (int i = 0; i < 100; i++) {
		temp2[0][i] = exp(-temp2[0][i]);
		aux[0][i] = temp2[0][i];
        
	}
	vector<vector<double>> temp3;
	temp3 = MultiT(aux, WeightMatrix);
	for (int i = 0; i < 10; i++) {
    
		temp3[0][i] = log(temp3[0][i] / (0.1*I0));

	}
	for (int i = 0; i < 10; i++) {
		g[0][i]= temp3[0][i]+Sino[0][i];
       
	}
    
	vector<vector<double>> temp4;
	temp4 = MultiT(g, g);
   
	return 0.5*temp4[0][0];
    
	
	
	} 
	/**/

};

