#include "mex.h"
#include "armadillo"
#include <iostream>
#include <random>
using namespace std;
using namespace arma;


// compile this in Matlab using
// mex -v CXXFLAGS="\$CFLAGS -std=gnu++0x" ABCmex.cpp
// run this using
// load('DifferentLatticeSpacings.mat', 'M24')
// [A,FR,FC]=ABCmex(0.25,0.0025,M24);
void matlab2arma(mat& A, const mxArray *mxdata){
// delete [] A.mem; // don't do this!
access::rw(A.mem)=mxGetPr(mxdata);
access::rw(A.n_rows)=mxGetM(mxdata); // transposed!
access::rw(A.n_cols)=mxGetN(mxdata);
access::rw(A.n_elem)=A.n_rows*A.n_cols;
};

void freeVar(mat& A, const double *ptr){
access::rw(A.mem)=ptr;
access::rw(A.n_rows)=1; // transposed!
access::rw(A.n_cols)=1;
access::rw(A.n_elem)=1;
};



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

    int ROWS=24;
    int COLS=32;
    int iterations=288;
    //double Pm=0.8;
    //double Pp=0.01;
    
    
/*
	mat M=zeros<mat>(2,2);
	M(0,0)=2.1;
	M(1,0)=3.4;
	M(0,1)=2.3;
	M(1,1)=2.45;
	*/
//	int ROWS=10;
//	int COLS=10;
//	int iterations=100;
	
        //mat output(1,1);
        //const double* outputmem=access::rw(output.mem);
        //matlab2arma(output,plhs[0]);
    
    if (nrhs != 3)
    mexErrMsgTxt("Incorrect number of input arguments");
    //if (nlhs != 4)
        if (nlhs != 3)
    mexErrMsgTxt("Incorrect number of output arguments");
    //if(!mxIsStruct(prhs[2]))
    //mexErrMsgTxt("Input must be a structure.");
    /*
    mat RowPos(1,1);
    const double* RowPosmem=access::rw(RowPos.mem);
    matlab2arma(RowPos,prhs[0]); // First create the matrix, then change it to point to the matlab data.

    mat ColPos(1,1);
    const double* ColPosmem=access::rw(ColPos.mem);
    matlab2arma(ColPos,prhs[0]); // First create the matrix, then change it to point to the matlab data.
    */
    double *p_Pm,*p_Pp;
    double Pm,Pp;
    //size_t mrows,ncols;

    if( !(mxGetN(prhs[2])==COLS && mxGetM(prhs[2])==ROWS) ) {
    mexErrMsgIdAndTxt( "MATLAB:inputNotCorrectDimension",
            "Input must be a matrix with correct dimension.");
    }
    
    /* The input must be a noncomplex scalar double.*/
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(mxGetM(prhs[0])==1 && mxGetN(prhs[0])==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:inputNotRealScalarDouble",
            "Input must be a noncomplex scalar double.");
    }

    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      !(mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:inputNotRealScalarDouble",
            "Input must be a noncomplex scalar double.");
    }


    p_Pm = mxGetPr(prhs[0]);
    p_Pp = mxGetPr(prhs[1]);
    Pm=*p_Pm;
    Pp=*p_Pp;

    mat Minit(1,1);
    const double* Minitmem=access::rw(Minit.mem);
    matlab2arma(Minit,prhs[2]); // First create the matrix, then change it to point to the matlab data.

    
    plhs[0] = mxCreateDoubleMatrix(ROWS, COLS, mxREAL);
    mat M(1,1);
    const double* Mmem=access::rw(M.mem);
    matlab2arma(M,plhs[0]);
    
    
    
    /*plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mat SS(1,1);
    const double* SSmem=access::rw(SS.mem);
    matlab2arma(SS,plhs[1]);
    
    plhs[1] = mxCreateDoubleMatrix(1, 24, mxREAL);
    mat totalMoves(1,1);
    const double* totalMovesmem=access::rw(totalMoves.mem);
    matlab2arma(totalMoves,plhs[2]);
    
    plhs[2] = mxCreateDoubleMatrix(2, 36, mxREAL);
    mat traj(1,1);
    const double* trajmem=access::rw(traj.mem);
    matlab2arma(traj,plhs[3]);
    */
    
    
    M=zeros<mat>(ROWS,COLS);

    //totalMoves=zeros<mat>(1,24);    
    //traj=zeros<mat>(2,36);    
    
    // THIS IS FOR WHEN WE WANT TO RANDOMLY INITIALISE A MATRIX WITH N CELLS.
   /* int N=10;
   
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator

	srand(time(0));

	int maxcells=ROWS*COLS;

	//int Y0=10; //initialise cells at top y0 rows of domain
	ivec RowPos(maxcells);
	ivec ColPos(maxcells);
    std::uniform_int_distribution<> distr1(0, ROWS-1); // define the range
	std::uniform_int_distribution<> distr2(0, COLS-1); // define the range


	RowPos(0)=distr1(eng);
	ColPos(0)=distr2(eng);
	int j=0;
	M(RowPos(j),ColPos(j))=1;
	while (j<N-1)
	{
		j=j+1;


	    RowPos(j)=distr1(eng);
	    ColPos(j)=distr2(eng);

	    if (M(RowPos(j),ColPos(j))==0)
		{
	    	M(RowPos(j),ColPos(j))=1;
   		}
	    else
	    {
	        j=j-1;
		}
	}*/
	
    
    int maxcells=ROWS*COLS;
	ivec RowPos(maxcells);
	ivec ColPos(maxcells);
    int N=0;
    int subN=0;
    
    vec CellIndices = zeros<vec>(5);
    
    for (int i=0; i<ROWS; i++) {
        for (int j=0; j<COLS; j++) {
            if (Minit(i,j)==1) {
                RowPos(N)=i;
                ColPos(N)=j;
                N=N+1;
                M(i,j)=1;
                //if (i>=24 && i<48 && j>=32 && j<64)
                //{
//                    if (subN<10)
//                    {
 //                       CellIndices(subN)=N-1;
 //                   }
 //                   subN++;
                    
                //}
            }
        }
    }
    
    if( N==0 ) {
    mexErrMsgIdAndTxt( "MATLAB:inputMatrixEmpty",
            "Input must not be empty.");
    }
	CellIndices=randi<vec>(5,distr_param(0,N-1));
    subN=5;
    
    
    plhs[1] = mxCreateDoubleMatrix(subN, 36, mxREAL);
    mat RowTraj(1,1);
    const double* RowTrajmem=access::rw(RowTraj.mem);
    matlab2arma(RowTraj,plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(subN, 36, mxREAL);
    mat ColTraj(1,1);
    const double* ColTrajmem=access::rw(ColTraj.mem);
    matlab2arma(ColTraj,plhs[2]);
    
    
    

    //M.print("M:");
	int oldN;
//    int trajIndex=rand()%N;
	ivec Permutation;
	vec Rands;
	for (int k=1;k<=iterations;k++)
	{
        if ((k-1)%8==0)
        {

            for (int cells=0;cells<subN;cells++)
            {
                //if (RowPos(CellIndices(cells))>=24 && RowPos(CellIndices(cells))<48 && ColPos(CellIndices(cells))>=32 && ColPos(CellIndices(cells))<64)
                //{
                        RowTraj(cells,(k-1)/8)=(double) RowPos(CellIndices(cells));
                        ColTraj(cells,(k-1)/8)=(double) ColPos(CellIndices(cells));
                //}
            }
            
        }

		oldN=N;
		Permutation=randi<ivec>(oldN,distr_param(0,oldN-1));
		Rands=randu<vec>(oldN);
		
		// Movement
		for (int j=0;j<oldN;j++)
		{
            if (Rands(j)<Pm)
            {
			//move left
			if (Rands(j) < Pm/4 && ColPos(Permutation(j))>0 && M(RowPos(Permutation(j)),ColPos(Permutation(j))-1)==0)
			{
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))--;
				ColPos(Permutation(j))--;
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))++;
                if (Permutation(j)<24) 
                {
    //                totalMoves(0,Permutation(j))++;
                }
			}
			//move right
			else if (Rands(j) > Pm/4 && Rands(j) < Pm/2 && ColPos(Permutation(j))+1 < COLS && M(RowPos(Permutation(j)),ColPos(Permutation(j))+1)==0)
			{
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))--;
				ColPos(Permutation(j))++;
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))++;
                if (Permutation(j)<24) 
                {
      //              totalMoves(0,Permutation(j))++;
                }
			}
			//move down
			else if (Rands(j) > Pm/2 && Rands(j) < 3*Pm/4 && RowPos(Permutation(j))>0 && M(RowPos(Permutation(j))-1,ColPos(Permutation(j)))==0)
			{
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))--;
				RowPos(Permutation(j))--;
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))++;
                if (Permutation(j)<24) 
                {
        //            totalMoves(0,Permutation(j))++;
                }
			}
			//move up
			else if (Rands(j) > 3*Pm/4 && Rands(j) < Pm && RowPos(Permutation(j))+1<ROWS && M(RowPos(Permutation(j))+1,ColPos(Permutation(j)))==0)
			{
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))--;
				RowPos(Permutation(j))++;
				M(RowPos(Permutation(j)),ColPos(Permutation(j)))++;
                if (Permutation(j)<24) 
                {
          //          totalMoves(0,Permutation(j))++;
                }
			}
            }
            else if (Rands(j)<Pm+Pp)
            {
            

	//	}

		// Proliferation
	//	Permutation=randi<ivec>(oldN,distr_param(0,oldN-1));
	//	Rands=randu<vec>(oldN);
	//	for (int j=0;j<oldN;j++)
	//	{
			//move left
			if (Rands(j) > Pm && Rands(j) < Pm+Pp/4 && ColPos(Permutation(j))>0 && M(RowPos(Permutation(j)),ColPos(Permutation(j))-1)==0)
			{
				N++;
				ColPos(N-1)=ColPos(Permutation(j))-1;
				RowPos(N-1)=RowPos(Permutation(j));
				M(RowPos(Permutation(j)),ColPos(Permutation(j))-1)++;
			}
			//move right
			else if (Rands(j) > Pm+Pp/4 && Rands(j) < Pm+Pp/2 && ColPos(Permutation(j))+1 < COLS && M(RowPos(Permutation(j)),ColPos(Permutation(j))+1)==0)
			{
				N++;
				ColPos(N-1)=ColPos(Permutation(j))+1;
				RowPos(N-1)=RowPos(Permutation(j));
				M(RowPos(Permutation(j)),ColPos(Permutation(j))+1)++;
			}
			//move down
			else if (Rands(j) > Pm+Pp/2 && Rands(j) < Pm+3*Pp/4 && RowPos(Permutation(j))>0 && M(RowPos(Permutation(j))-1,ColPos(Permutation(j)))==0)
			{
				N++;
				ColPos(N-1)=ColPos(Permutation(j));
				RowPos(N-1)=RowPos(Permutation(j))-1;
				M(RowPos(Permutation(j))-1,ColPos(Permutation(j)))++;
			}
			//move up
			else if (Rands(j) > Pm+3*Pp/4 && Rands(j) < Pm+Pp && RowPos(Permutation(j))+1<ROWS && M(RowPos(Permutation(j))+1,ColPos(Permutation(j)))==0)
			{
				N++;
				ColPos(N-1)=ColPos(Permutation(j));
				RowPos(N-1)=RowPos(Permutation(j))+1;
				M(RowPos(Permutation(j))+1,ColPos(Permutation(j)))++;
			}
            }
			//printf("%d , %d\n",accu(M),N);
		}
		//M.print("M:");
	}    
    
    
    
    // Compute Summary Statistics:   
    //SS=sum(sum(M));
    
    /*X=27;
    Y=36;
    N=length(RowPos);
    D=dist(RowPos);
    for i=1:Y-1
    c(i)=length(find(triu(D)==i));
    q(i)=c(i)/(X*X*(Y-i)*(N/(X*Y))*((N-1)/(X*Y-1)));
    end*/
    
    freeVar(RowTraj,RowTrajmem);
    freeVar(ColTraj,ColTrajmem);
//    freeVar(totalMoves,totalMovesmem);
    freeVar(Minit,Minitmem);
    freeVar(M,Mmem);
    //freeVar(SS,SSmem);
    //mxFree(pointer);
    //mxFree(pointer2);
    return;
}
