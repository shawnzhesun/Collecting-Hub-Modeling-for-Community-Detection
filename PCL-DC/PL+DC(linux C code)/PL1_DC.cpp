/********************************************************************************************************
   this program implements the PCL link model + Discriminative Content Model on attributes 
   by EM algorithm. In DC Model, the Nesterov Method is used. If there is no attributes
   available, we only use PCL, otherwise the PCL on link + DC on attributes are 
   used. 
   Author: Tianbao Yang
*******************************************************************************************************/


#include<iostream>
#include<string>
#include<cmath>
#include<fstream>
#include<cstdlib>
#include<algorithm>
#include<iterator>
#include<numeric>
#include<vector>
#include <iomanip> 
using namespace std;
struct link 
	{
	   public:
      link(){};
      link(int x, int y, double weight, int link_ind){i=x; j=y; w=weight; link_index=link_ind;}
	    int i;
	    int j;
	    double w;
      int link_index;
	} ;
class target_point
  {
     public :
        target_point(){};
		    target_point(int y, double weight, int link_ind){ind=y; w=weight; link_index=link_ind;}
	      int ind;
		    double w;
		    int link_index;
        ~target_point(){};
  } ;

#define  MAXROW 50000



void usage(); // output the useage 
int Read_File(ifstream &, ifstream &,  double ** W, double ** X, int n, int D, int flag_data); // Read link, attributes file, save link in W, attributes in X, n is the number of examples, D is the number of attributes 
void Initialize(double ** G, double ** Y,double * B, double ** R, double **Model, int n,int K,int D, int flag_data, int seed);  // Initial G,Y, B, R, Model 
int Iteration(double **W,double **X, double **G, double **Y, double *B, double **R, double **Model,int n, int K, int D, int N,double ep,double lamda, string stop, int verbosity, int flag_data, string save_file);
double objective1(vector<struct link> S, double **Y, double **R, double *col_sum_R, int n, int K);
void LogisticModel(double **X, double **G, double **Y, double **Model, double lamda, double ep, int N, string stop, int verbosity, int n, int K, int D);
void newton_solve_equation(double p0, double delta, double epsilon, int maxN, double *Ax, double *Cx, double dx, int K, double *returned);
double objective2(double **X, double **G, double **model, double **tempY, double lamda, int n, int K, int D);
void gradient2(double **X, double **G, double **tempY, double **model, double **grad, double lamda,int n, int K, int D);
double norm(double **grad, int K, int D);



int main(int argc, char * argv[])
{
   
  int  n, K,D, N, verbosity,seed; //number of examples,number of communities, number of attributes, number of iterations, level of ouptu, initial state for rand
	double ep, lamda; // precision for iteration, lamda for logistic model
	string link_file, data_file, stop, save_file;
	int flag_data; // indicate whether attributes are available or not

  /*default value */
	n=0;// number of nodes
	K=0;// number of attributes
	D=0;// number of clusters

	N=10000;// number of iterations
	verbosity=0;// verbosity level
	ep=0.00000001;
	lamda=0;
	flag_data=0;
	seed=1;// seed for rand number 
  link_file="";
	data_file="";
	stop="fulliteration";
	string cmd_file="";

	if(argc<5)
	{
		usage();
    return -1;
    n=9;// number of nodes
	  D=2;// number of attributes
	  K=2;// number of clusters
		N=1000;
		verbosity=1;
		seed=1;
		flag_data=1;
		link_file=string("toy.link");
    data_file=string("toy.data");
    save_file=string("toy");
	}
	else
    {
	    ifstream cmd_line(argv[1]);
	    if(!cmd_line)
	    {
		  cerr<<"unable to open cmd file: "<< argv[1]<<"-- bailing out!"<<endl;
		  return -1;
	    }
	    char str[100];
      while(cmd_line.getline(str,100))
	    { 
        string line(str);
		    size_t len=line.length();
        size_t epos=line.find_first_of(" ",0);
	    	string para_name=line.substr(0, epos);
		    string para_value=line.substr(epos+1, len-epos-1);
		   //cout<<para_name<<" "<<para_value<<endl;
		   if( strcmp(para_name.c_str(), "Link")==0 )
		   	link_file=para_value;
		   else if( strcmp(para_name.c_str(), "Attribute")==0 )
               data_file=para_value;
		   else if( strcmp(para_name.c_str(), "n")==0 )
		   	n=atoi(para_value.c_str());
		   else if( strcmp(para_name.c_str(), "D")==0 )
              D=atoi(para_value.c_str());
		   else if( strcmp(para_name.c_str(), "K")==0 )
               K=atoi(para_value.c_str());
		   else if( strcmp(para_name.c_str(), "ep")==0 && strcmp(para_value.c_str(), "")!=0)
		   {
		   	ep=atof(para_value.c_str());
		   	stop="converge";
		   }
		   else if( strcmp(para_name.c_str(), "N")==0 && strcmp(para_value.c_str(), "")!=0)
		   {
		       N=atoi(para_value.c_str());
		   	stop="fulliteration";
		   }
		   else if( strcmp(para_name.c_str(), "lamda")==0 )
                       lamda=atof(para_value.c_str());
 		   else if( strcmp(para_name.c_str(), "verbosity")==0 )
		       verbosity=atoi(para_value.c_str());
		   else if( strcmp(para_name.c_str(), "seed")==0 )
		       seed=atoi(para_value.c_str());
	     }
	}
	if (strcmp(link_file.c_str(),"")==0 && strcmp(data_file.c_str(), "")==0)
		cerr<<"no link and data"<<endl;
	else if (strcmp(data_file.c_str(), "")==0)
		flag_data=0;
	else
		flag_data=1;

	cout<<"Link: "<<link_file<<endl;
        cout<<"Attribute: "<<data_file<<endl;
	cout<<"n: "<<n<<endl;
	cout<<"D: "<<D<<endl;
        cout<<"K: "<<K<<endl;
	cout<<"N: "<<N<<endl;
	cout<<"ep: "<<ep<<endl;
	cout<<"lamda: "<<lamda<<endl;
	cout<<"verbosity: "<<verbosity<<endl;
	cout<<"seed: "<<seed<<endl;
	cout<<"stop: "<<stop<<endl;
	cout<<"flag_data: "<<flag_data<<endl;
    

   /* read link and attriubtes */
   
    // open files
  ifstream input1,input2;
	input1.open(argv[2]);
	//input1.open(link_file.c_str(), ifstream::in);
	input2.open(argv[3]);
  //input2.open(data_file.c_str(), ifstream::in);
  if (!input1)
	{
		cerr<<"unable to open link file: "<< link_file<<endl;
		return -1;
	}
  if (!input2)
	{
		cerr<<"no data file "<< data_file<<endl;
	}
   // allocate memeory
	double **W,**X;
	W=new double*[n+1];
	for(int i=0;i<=n; i++)
	   W[i]=new double[n+1];
    X=new double*[n+1];
	for(int i=0; i<=n; i++)
		X[i]=new double[D+1];

  // read file
  int status=Read_File(input1,input2, W, X,n, D, flag_data);
	/*cout<<"W"<<endl;
	for (int i=1; i<=n; i++){
	   copy(&(W[i][1]), &(W[i][n+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }
	if(flag_data==1)
	{
	  cout<<endl;
   	  for (int i=1; i<=n; i++){
	   copy(&(X[i][1]), &(X[i][D+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }    
	}*/
  
	/* initialize G  Y B R Model*/

	double **G=new double*[n+1];
	for(int i=0;i<=n;i++)
	   G[i]=new double[K+1];

	double **Y=new double*[n+1];
	for(int i=0;i<=n;i++)
	   Y[i]=new double[K+1];

	double *B=new double[n+1];

	double **R=new double*[n+1];
	for(int i=0;i<=n;i++)
	   R[i]=new double[K+1];    

	double **Model=new double*[K+1];
	   for(int k=0;k<=K;k++)
	     Model[k]=new double[D+1];   


   Initialize(G,Y, B, R,Model,n,K,D,flag_data,seed);
   //befor iteration
   /*cout << fixed << showpoint; 
   cout << setprecision(10);
   cout<<endl<<"before iteration"<<endl;
   cout<<endl<<"G"<<endl;
	for (int i=1; i<=n; i++)
	{
	   copy(&(G[i][1]), &(G[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	 }
	cout<<endl<<"B"<<endl;
	copy(&(B[1]), &(B[n+1]), ostream_iterator<double>(cout, "\n"));

	cout<<endl<<"Y"<<endl;
	for (int i=1; i<=n; i++)
	{
	   copy(&(Y[i][1]), &(Y[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }


	cout<<endl<<"R"<<endl;
	for (int i=1; i<=n; i++)
	{
	   copy(&(R[i][1]), &(R[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	 }*/
    
   /*  Iteration EM algorithm   */
   save_file=argv[4];
   int converge=Iteration(W,X,G,Y,B,R,Model,n,K,D,N,ep,lamda,stop,verbosity,flag_data, save_file);
   if(converge)
	   cout<<endl<<"Converge!"<<endl;
   else
	   cout<<endl<<"Fulliteration!"<<endl;
	   
   
  //after iteration
  string result_file;
  ofstream out;  
  
  // output G   
  result_file=save_file+".G";
  out.open(result_file.c_str());
  out << fixed << showpoint; 
  out << setprecision(10);
	for (int i=1; i<=n; i++){
	   copy(&(G[i][1]), &(G[i][K+1]), ostream_iterator<double>(out, " "));
	   out<<endl;
	  }
	// output Y   
	result_file=save_file+".Y";
  out.open(result_file.c_str());  
  out << fixed << showpoint; 
  out << setprecision(10);
	for (int i=1; i<=n; i++)
	{
	   copy(&(Y[i][1]), &(Y[i][K+1]), ostream_iterator<double>(out, " "));
	   out<<endl;
	  }
	  
	//output B  
  result_file=save_file+".B";
  out.open(result_file.c_str());  
  out << fixed << showpoint; 
  out << setprecision(10);
	copy(&(B[1]), &(B[n+1]), ostream_iterator<double>(out, "\n"));


  //output R
  result_file=save_file+".R";
  out.open(result_file.c_str());  
  out << fixed << showpoint; 
  out << setprecision(10);
	for (int i=1; i<=n; i++)
	{
	   copy(&(R[i][1]), &(R[i][K+1]), ostream_iterator<double>(out, " "));
	   out<<endl;
	  
	}
	if(flag_data==1)
	{  	// output model
		  result_file=save_file+".model";
      out.open(result_file.c_str());  
      out << fixed << showpoint; 
      out << setprecision(10);
       for (int k=1; k<=K; k++)
	   {
	     copy(&(Model[k][1]), &(Model[k][D+1]), ostream_iterator<double>(out, " "));
	     out<<endl;
	   }
	}

   delete [] W;
   delete [] X;
   delete [] G;
   delete [] Y;
   delete [] B;
   delete [] R;
   delete [] Model;

}


void usage()
{   
	cout<<"*******************************************************************************"<<endl;
	cout<<"usage: IPHITS.exe cmd_line_file link_file data_file save_file"<<endl;
	cout<<"*******************************************************************************"<<endl;
}


int Read_File(ifstream &input1,ifstream &input2, double **W, double **X,int n, int D, int flag_data)
{
   
	int i=1, j=1;
	char str[MAXROW];
	/* read link */
	cout<<"\nloading link data"<<endl;
	while(input1.getline(str,MAXROW)){
      string line(str);
	  size_t pos=0;
	  size_t len=line.length();
	  while(pos<len){
		  size_t epos=line.find_first_of(" ",pos);
		  string wstr;
		  if(epos!=string::npos){
		    wstr=line.substr(pos, epos-pos);
			pos=epos+1;
		  }
		  else{
            wstr=line.substr(pos, len-pos);
			pos=len;
		  }
		  W[i][j]=atof(wstr.c_str());
		  if(i==j)
			  W[i][j]=0;
          j=j+1;
	  }
	  if(j<n+1){
		  cerr<<"Link: the format of the file is not correct: number of columns"<<endl;
          return -1;
	  }
	  else{
		  j=1;
		  i=i+1;
	  }
	} 

	if(i<n+1)
	{
		cerr<<"Link: the format of the file is not correct: number of rows"<<endl;
		return -1;
	}

	// normalize W
	double sum=0.0;
	for (i=1; i<=n; i++)
	{   
		double row_sum=accumulate(&W[i][1], &W[i][n+1], 0.0);
		sum=sum+row_sum;
	}
	for(i=1; i<=n; i++)
	{
		 for(j=1;j<=n; j++)
			W[i][j]=W[i][j]/sum;
	}


  if (flag_data==0)
	{
		return 1;
	}
	/* read data */
	cout<<"\nloading attributes data"<<endl;
	i=1;
	j=1;
	while(input2.getline(str,MAXROW)){
      string line(str);
	  size_t pos=0;
	  size_t len=line.length();
	  while(pos<len){
		  size_t epos=line.find_first_of(" \n",pos);
		  string wstr;
		  if(epos!=string::npos){
		    wstr=line.substr(pos, epos-pos);
			pos=epos+1;
		  }
		  else{
            wstr=line.substr(pos, len-pos);
			pos=len;
		  }
		  X[i][j]=atof(wstr.c_str());
          j=j+1;
	  }
	  if(j<D+1){
		  cerr<<"Data: the format of the file is not correct: number of columns"<<endl;
          return -1;
	  }
	  else{
		  j=1;
		  i=i+1;
	  }
	}

	if(i<n+1)
	{
		cerr<<"Data: the format of the file is not correct: number of rows"<<endl;
		return -1;
	}

	return 1;
}



void Initialize(double ** G, double ** Y,double * B, double ** R, double **Model, const int n,const int K, const int D, int flag_data, int seed)
{
	srand(seed);
    for(int i=1;i<=n;i++)
	{
	    for (int k=1; k<=K;k++)
		 {
		   G[i][k]=double(rand())/double(RAND_MAX);
		   Y[i][k]=G[i][k];
	     }
	 }
	for (int i=1; i<=n; i++)
	{   
		double row_sum=accumulate(&G[i][1], &G[i][K+1], 0.0);
		if(row_sum>0)
		{
			for(int k=1;k<=K; k++)
			{
			  G[i][k]=G[i][k]/row_sum;
			  Y[i][k]=G[i][k];
			}	    
		}
	}

   for(int i=1;i<=n;i++)
   {
      B[i]=double(rand())/double(RAND_MAX);
   }
   for(int i=1; i<=n; i++)
   {
	   double sumB=accumulate(&B[1], &B[n+1], 0.0);
	   B[i]=B[i]/sumB;
	   for(int k=1; k<=K; k++)
		 R[i][k]=Y[i][k]*B[i];
   }  

   if(flag_data==1)
   {
	   for(int k=1;k<=K; k++)
		   for(int j=1; j<=D; j++)
			   Model[k][j]=0.0;
   }
}


int Iteration(double **W,double **X, double **G, double **Y, double *B, double **R, double **Model,const int n, const int K, const int D, int N,double ep, double lamda, string stop, int verbosity, int flag_data, string save_file)
{ /* implement EM iteratoin */
 
  // S array link
  int status=-1;
  vector<struct link> S;
  typedef vector<target_point> linklist;


  linklist *outlink= new linklist[n+1];
  linklist *inlink=new linklist[n+1];


  int linkindex=1;
  for(int i=1; i<=n; i++)
	  for(int j=1; j<=n; j++)
	  if (W[i][j]>0){
          struct link link_obj(i,j, W[i][j], linkindex);
		  target_point dp(j, W[i][j], linkindex);
		  target_point sp(i, W[i][j], linkindex);
		  S.push_back(link_obj);
		  outlink[i].push_back(dp);
		  inlink[j].push_back(sp);
		  linkindex++;
	  }
  int linkN=int(S.size());

  
  double **q=new double*[linkN+1];
  for(int i=1; i<=linkN; i++)
	  q[i]=new double[K+3];

  int *nin=new int[n+1];
  int *nout=new int[n+1];
  for(int i=1; i<=n; i++)
  {
	  nin[i]=int(inlink[i].size());
	  nout[i]=int(outlink[i].size());
  }

  /*for(vector<link>::iterator  it=S.begin(); it!= S.end(); it++)
		cout<<"link: "<<(*it).i<<"==>"<<(*it).j<<" "<<(*it).w<<endl;
  for(int i=1; i<=n; i++)
  {
	cout<<"node"<<i<<":==>";
    for(vector<target_point>::iterator  it=outlink[i].begin(); it!= outlink[i].end(); it++)
		cout<<"("<<(*it).ind<<" "<<(*it).w<<")";
    cout<<endl; 
    cout<<"node"<<i<<":<==";
    for(vector<target_point>::iterator  it=inlink[i].begin(); it!= inlink[i].end(); it++)
		cout<<"("<<(*it).ind<<" "<<(*it).w<<")";
    cout<<endl; 
  }*/
  

  double **Ain=new double*[n+1];
	  for(int i=0;i<=n;i++)
	   Ain[i]=new double[K+1];

  double **Aout=new double*[n+1];
	  for(int i=0;i<=n;i++)
	   Aout[i]=new double[K+1];

  double **A=new double*[n+1];
	  for(int i=0;i<=n;i++)
	   A[i]=new double[K+1];

  double **M=new double*[n+1];
	 for(int i=0;i<=n;i++)
	   M[i]=new double[K+1];

  double *col_sum_R=new double[K+1];
  double *col_sum_M=new double[K+1];

  double *obj=new double[N+1];

  //compute the column sum of R
	  
  for(int k=1; k<=K; k++)
	 {
		 col_sum_R[k]=0;
		 for(int i=1; i<=n; i++)
			 col_sum_R[k]=col_sum_R[k]+R[i][k];
	  }



  double obj_curt=objective1(S, Y, R, col_sum_R, n, K);
  double obj_last=obj_curt-1;
  int iteration=0;
  while(true)
  {
	 if(obj_curt>=obj_last && (obj_curt-obj_last)<ep)
	  {
		  status=1;
		  break;
	  }
	 if( strcmp(stop.c_str(), "fulliteration")==0 && iteration>=N)
	  {
		  status=0;
		  break;
	  }



     // compute q
	  for(int l=1; l<=linkN; l++)
	  {   int i=S[l-1].i;
	      int j=S[l-1].j;
		  double w=S[l-1].w;
		  q[l][1]=i;
		  q[l][2]=j;
		  double sum=0.0;
		  for(int k=1; k<=K; k++)
		  {
			  q[l][k+2]=Y[i][k]*R[j][k]/(col_sum_R[k]);
			  sum+=q[l][k+2];
		  }
		  for(int k=1; k<=K; k++)
			  q[l][k+2]/=sum;
	  }

	  // compute Ain Aout M
     for(int i=1; i<=n; i++)
	 {
		 for(int k=1; k<=K; k++)
		 {    Ain[i][k]=0.0;
		      Aout[i][k]=0.0;
			  for(vector<target_point>::iterator  it=inlink[i].begin(); it!= inlink[i].end(); it++)
		          // j=(*it).ind inlink of i, qjik=Y[jk]*Y[ik]*B[i]/sum(Y[ik]B[i])
			  {   int j=(*it).ind;
			      int link_index=(*it).link_index;
				  double qjik=q[link_index][k+2];
				  Ain[i][k]=Ain[i][k]+(*it).w*qjik;
			  }
			  for(vector<target_point>::iterator  it=outlink[i].begin(); it!= outlink[i].end(); it++)
		          // j=(*it).ind outlink of i, qijk=Y[ik]*Y[jk]*B[j]/sum(Y[jk]B[j])
			  {   int j=(*it).ind;
			      int link_index=(*it).link_index;
				  double qijk=q[link_index][k+2];
				  Aout[i][k]=Aout[i][k]+(*it).w*qijk;
			  }
			  M[i][k]=Aout[i][k]/col_sum_R[k];
			  A[i][k]=Ain[i][k]+Aout[i][k];
		 }
	 }
	/*cout<<endl<<"Ain"<<endl;
	for (int i=1; i<=n; i++){
	   copy(&(Ain[i][1]), &(Ain[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }
    cout<<endl<<"Aout"<<endl;
	for (int i=1; i<=n; i++){
	   copy(&(Aout[i][1]), &(Aout[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }
    cout<<endl<<"A"<<endl;
	for (int i=1; i<=n; i++){
	   copy(&(A[i][1]), &(A[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }
    cout<<endl<<"M"<<endl;
	for (int i=1; i<=n; i++){
	   copy(&(M[i][1]), &(M[i][K+1]), ostream_iterator<double>(cout, " "));
	   cout<<endl;
	  }*/

     // update G and B
	 // compute column sum of M
     for(int k=1; k<=K; k++)
	  {
		  col_sum_M[k]=0;
		  for(int i=1; i<=n; i++)
			  col_sum_M[k]=col_sum_M[k]+M[i][k];
	  }   

	 for(int i=1; i<=n ; i++)
	 {
		 if(nin[i]==0 && nout[i]==0)
		 {
			 B[i]=0.0;
			 for(int k=1; k<=K; k++){
				 G[i][k]=double(1)/double(K);
				 Y[i][k]=G[i][k];
			 }
		 }
		 if(nin[i]==0 && nout[i]>0)
		 {
			 B[i]=0;
			 double sum=accumulate(&(Aout[i][1]), &(Aout[i][K+1]), 0.0);
			 for(int k=1; k<=K; k++){
				 G[i][k]=Aout[i][k]/sum;
				 Y[i][k]=G[i][k];
			 }
		 }
		 if(nin[i]>0 && nout[i]==0)
		 {
		    double sum=0.0;
			for(int k=1; k<=K; k++)
			{
				G[i][k]=Ain[i][k]/col_sum_M[k];
				sum+=G[i][k];
			}
			for(int k=1; k<=K; k++){
				G[i][k]=G[i][k]/sum;
				Y[i][k]=G[i][k];
			}
			double deno=0.0;
			double Ain_sum=accumulate(&(Ain[i][1]), &(Ain[i][K+1]), 0.0);
			for(int k=1; k<=K; k++)
			 deno+=(col_sum_M[k]*G[i][k]);
			B[i]=Ain_sum/deno;
		 }
		 if(nin[i]>0 && nout[i]>0)
		 {   /*first get the parameters of equations involving b, \sum_k A_k/(C_k+x) */
			 double p0=0.0;// the returned values of netwon iteration 
			 double delta=0.0000000001,epsilon=0.0000000001;// error tolerance
			 int maxN=1000;// max numberof iterations 
			 double *Ax=new double[K+1];
			 double *Cx=new double[K+1];
			 double *returned=new double[4];
			 double dx;
			 double Aout_sum=accumulate(&(Aout[i][1]), &(Aout[i][K+1]), 0.0);
			 double Ain_sum=accumulate(&(Ain[i][1]), &(Ain[i][K+1]), 0.0);
			 double A_sum=accumulate(&(A[i][1]), &(A[i][K+1]),0.0);
			 dx=A_sum-Ain_sum;
			 for(int k=1; k<=K; k++)
			 {
				 Ax[k]=Aout_sum*A[i][k]/col_sum_M[k];
				 Cx[k]=Aout_sum/col_sum_M[k];
			 }
			 newton_solve_equation(p0, delta, epsilon, maxN, Ax, Cx, dx, K, returned);// sovle the equation
			 B[i]=returned[0];
			 for(int k=1;k<=K; k++)
				 G[i][k]=A[i][k]/(B[i]*col_sum_M[k]+Aout_sum);
			 double sum=accumulate(&(G[i][1]), &(G[i][K+1]), 0.0);
			 for(int k=1; k<=K; k++){
				 G[i][k]/=sum;
				 Y[i][k]=G[i][k];
			 }
			 delete [] Ax;
			 delete [] Cx;
			 delete [] returned;
		 }
	 }

	 // normalize B
	 double sum=accumulate(&(B[1]), &(B[n+1]), 0.0);
	 for(int i=1; i<=n; i++)
		 B[i]/=sum;

	 if(flag_data==1)
	 {
		double itep=1e-5;
    int    itN=1000;
    string stopp="converge";
		if(iteration<N-1)
		{
      itN=100;
			stopp="fulliteration";
		}
     LogisticModel(X, G, Y, Model, lamda, itep, itN, stopp, verbosity, n, K, D);
	 }

     // compute R
     for(int i=1; i<=n ; i++)
		 for(int k=1; k<=K; k++)
			 R[i][k]=B[i]*Y[i][k];

    //compute the column sum of R	  
    for(int k=1; k<=K; k++)
	 {
		 col_sum_R[k]=0;
		 for(int i=1; i<=n; i++)
			 col_sum_R[k]=col_sum_R[k]+R[i][k];
	  }


	 obj_last=obj_curt;
     obj_curt=objective1(S, Y,  R,col_sum_R, n, K);
	 iteration++;
     obj[iteration]=obj_curt;
     if(verbosity>=2)
            printf("1.....................iteration=%d,  first Objective=%.10f\n", iteration, obj_curt);
     else if (verbosity==1 && fmod(double(iteration), double(100))==0)
            printf("1.....................iteration=%d, first Objective=%.10f\n", iteration, obj_curt);    
  }

  delete [] outlink;
  delete [] inlink;
  delete [] Ain;
  delete [] Aout;
  delete [] A;
  delete [] M;
  delete [] col_sum_R;
  delete [] col_sum_M;
  delete [] nin;
  delete [] nout;
  delete [] q;
  string str=save_file+".obj";
  ofstream output(str.c_str());
  output.precision(15);
  copy(&(obj[1]), &(obj[iteration+1] ), ostream_iterator<double>(output, " \n"));
  cout.precision(15);
  if(verbosity==0)
	  cout<<"\nfinal objective:"<<obj[iteration]<<endl;
  delete [] obj;
  return status;
}


double objective1(vector<struct link> S, double **Y, double **R, double *col_sum_R,  int n, int K)
{
	double obj=0.0;
	for(vector<struct link>::iterator it=S.begin(); it!=S.end(); it++)
	{
		double logterm=0;
		int i=(*it).i;
		int j=(*it).j;
		for(int k=1; k<=K; k++)
			if(R[j][k]>0 && Y[i][k]>0)
			{
				logterm+=Y[i][k]*R[j][k]/col_sum_R[k];
			}
		obj+=(*it).w*log(logterm);
	}
	return obj;
}

void newton_solve_equation(double p0, double delta, double epsilon, int maxN, double *Ax, double *Cx,double dx,int K, double *returned)
 {
   int it=1;
   double obj_p0;
   double gra_p0;
   double err,relerr,p1;
   for(it=1; it<=maxN; it++)
   {
       obj_p0=0.0;
	   gra_p0=0.0;
	   for(int k=1; k<=K; k++){
		   obj_p0+=Ax[k]/(Cx[k]+p0);
		   gra_p0+=-Ax[k]/((Cx[k]+p0)*(Cx[k]+p0));
	   }
	   obj_p0=obj_p0-dx;
       if(abs(obj_p0)<epsilon)
		   break;
	   p1=p0-obj_p0/gra_p0;
	   err=fabs(p1-p0);
	   relerr=2*err/(fabs(p1)+delta);
	   p0=p1;
	   if(err<delta || relerr<delta)
		   break;
   }
   returned[0]=p0;
   returned[1]=err;
   returned[2]=it;
   returned[3]=obj_p0;
 }


void LogisticModel(double **X, double **G, double **Y, double **Model, double lamda, double ep, int N, string stop, int verbosity, int n, int K, int D)
{ 
	

	double **model_Yt=new double*[K+1];
	double **model_Yt_1=new double*[K+1];
	double **model_Xt=new double*[K+1];
	double **model_Xt_L=new double*[K+1];
	double **grad_Xt=new double*[K+1];
	for(int k=0; k<=K; k++){
		model_Yt[k]=new double[D+1];
		model_Yt_1[k]=new double[D+1];
		model_Xt[k]=new double[D+1];
		model_Xt_L[k]=new double[D+1];
		grad_Xt[k]=new double[D+1];
	}
	double **tempY=new double*[n+1];
	for(int i=0; i<=n; i++)
		tempY[i]=new double[K+1];

    for(int k=1; k<=K; k++)
		for(int d=1; d<=D; d++)
		{
           model_Yt[k][d]=0.0;
		   model_Yt_1[k][d]=0.0;
		}
    int iteration=1;
	double t_1=1;
	double t_2=0;
	double Lt_1=10;
	double obj_Yt=objective2(X, G, model_Yt, tempY, lamda, n, K, D);
	double obj_Yt_1=obj_Yt-1;

	while( true)
	{
	  if(obj_Yt>=obj_Yt_1 && (obj_Yt-obj_Yt_1)<=ep)
		  break;
	  if(strcmp(stop.c_str(), "fulliteration")==0 && iteration>=N)
		  break;
	  double alphat=(t_2-1)/t_1;
	  for(int k=1; k<=K; k++)
		  for(int d=1; d<=D; d++)
 	         model_Xt[k][d]=(1+alphat)*model_Yt[k][d]-alphat*model_Yt_1[k][d];
	  double obj_Xt=objective2(X, G, model_Xt, tempY, lamda, n, K, D);
	  gradient2(X, G, tempY, model_Xt, grad_Xt, lamda, n, K, D);
	  double L=Lt_1;
	  for(int k=1; k<=K; k++)
		  for(int d=1; d<=D; d++)
	         model_Xt_L[k][d]=model_Xt[k][d]+(1/L)*grad_Xt[k][d];
	  double obj_Xt_L=objective2(X, G, model_Xt_L, tempY, lamda, n, K, D);
	  double norm_grad_Xt=norm(grad_Xt, K, D);
	  while( obj_Xt_L< obj_Xt+(1/(2*L))*norm_grad_Xt)
	  {
		  L=2*L;
		  for(int k=1; k<=K; k++)
		    for(int d=1; d<=D; d++)
	          model_Xt_L[k][d]=model_Xt[k][d]+(1/L)*grad_Xt[k][d];
		  obj_Xt_L=objective2(X,G, model_Xt_L, tempY, lamda, n, K, D);
	  }
	  Lt_1=L;

	  for(int k=1; k<=K; k++)
		 for(int d=1; d<=D; d++)
	        model_Yt_1[k][d]=model_Yt[k][d];

  	  for(int k=1; k<=K; k++)
		  for(int d=1; d<=D; d++)
	         model_Yt[k][d]=model_Xt_L[k][d];
	  obj_Yt_1=obj_Yt;
	  obj_Yt=obj_Xt_L;

	  t_2=t_1;
	  t_1=(0.5)*(1+sqrt(1+4*t_2*t_2));
	  iteration++;
	  if(verbosity>=3)
		  printf("2.......iteration=%d, obj=%.10f\n", iteration, obj_Yt);
	}
   
   for(int k=1; k<=K; k++)
	   for(int d=1; d<=D; d++)
		   Model[k][d]=model_Yt[k][d];
   // compute Y
   double *A_row=new double[n+1];
   for(int i=1; i<=n; i++)
   {   
	   double maxA_row=0.0;
	   for(int k=1; k<=K; k++)
	   {   A_row[k]=0.0;
		   for(int d=1; d<=D; d++)
		     A_row[k]+=X[i][d]*Model[k][d];
		   if(A_row[k]>maxA_row)
			   maxA_row=A_row[k];
	   }
	   for(int k=1; k<=K; k++){
		   A_row[k]-=maxA_row;
		   Y[i][k]=exp(A_row[k]);
	   }
   }
   for(int i=1; i<=n; i++)
   {   
	   double sum=accumulate(&(Y[i][1]), &(Y[i][K+1]), 0.0);
	   for(int k=1; k<=K; k++)
		   Y[i][k]/=sum;
   }

   delete [] A_row;
   delete [] model_Yt;
   delete [] model_Yt_1;
   delete [] model_Xt;
   delete [] model_Xt_L;
   delete [] grad_Xt;
   delete [] tempY;
}


double objective2(double **X, double **G, double **model, double **tempY, double lamda, int n, int K, int D)
{
   // compute tempY
   double *A_row=new double[n+1];
   for(int i=1; i<=n; i++)
   {   
	   double maxA_row=0.0;
	   for(int k=1; k<=K; k++)
	   {   A_row[k]=0.0;
		   for(int d=1; d<=D; d++)
		     A_row[k]+=X[i][d]*model[k][d];
		   if(A_row[k]>maxA_row)
			   maxA_row=A_row[k];
	   }
	   for(int k=1; k<=K; k++){
		   A_row[k]-=maxA_row;
		   tempY[i][k]=exp(A_row[k]);
	   }
   }
   for(int i=1; i<=n; i++)
   {   
	   double sum=accumulate(&(tempY[i][1]), &(tempY[i][K+1]), 0.0);
	   for(int k=1; k<=K; k++)
		   tempY[i][k]/=sum;
   }
   double obj=0;
   for(int i=1; i<=n; i++)
	   for(int k=1; k<=K; k++)
		   if(G[i][k]>0)
			   obj+=G[i][k]*log(tempY[i][k]);
   if(lamda>0)
	   {
		  double Frob_norm=0;
		  for(int k=1; k<=K; k++)
		    for(int d=1; d<=D; d++)
			   Frob_norm+=model[k][d]*model[k][d];
	      obj-=(lamda/2)*Frob_norm;
       }

   delete [] A_row;
   return obj;
}

void gradient2(double **X, double **G, double **Y, double **model, double **grad, double lamda,int n, int K, int D)	
{
	for(int k=1; k<=K; k++)
		for(int d=1; d<=D; d++)
		{   grad[k][d]=0.0;
			for(int i=1; i<=n; i++)
				grad[k][d]+=(G[i][k]-Y[i][k])*X[i][d];
			grad[k][d]-=lamda*model[k][d];
		}    
}

double norm(double **grad, int K, int D)
{
	double Frob_norm=0.0;
	for(int k=1; k<=K; k++)
		for(int d=1; d<=D; d++)
			Frob_norm+=grad[k][d]*grad[k][d];
	return Frob_norm;
}
