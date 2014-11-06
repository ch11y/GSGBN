#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>

#include "solve_GSGBN.h" 

using namespace std ; 

#define SZ(x) (int)x.size()

char ss[5000000+ 10] ; 
//./GSBN -i input.txt  -l1 1.0 -l2 1. -num 5 -thr 0.01 
int main (int argc , char **argv ){ 
    
    clock_t start = clock() ;
    vector< pair< string , string > > parameter  ; 
    for(int i=1; i < argc ; i+=2)
        parameter.push_back(make_pair(argv[i] , argv[i+1])) ; 
    sort( parameter.begin() , parameter.end())   ;  
    if( SZ(parameter) != 5 ){ 
        cout << "Grouped sparse Gaussian Bayesian network (GSGBN) structure learning." << endl;
        cout << "Usage: " << endl;
        cout << "    GSGBN -i input.txt -l1 1.0 -l2 0.1 -num 10 -thr 0.1" <<endl; 
        cout << "Options: "<< endl;
        cout << "    -i,   an m * n matrix with m samples and n variables;" <<  endl;
        cout << "    -l1,  the regularization parameter lambda1 for sparsity;" << endl ;
        cout << "    -l2,  the regularization parameter lambda2*n is used for grouping;" << endl ;
        cout << "    -num, number of rounds for enumerating DAGs;" << endl; 
        cout << "    -thr, threshold for filtering" << endl;  
        exit( 1 ); 
    }

    FILE *fp = fopen(parameter[0].second.c_str() , "r" ) ; 
    assert( NULL != fp );     
    vector< VD > input; 
      
    while( fgets( ss , 5000000, fp ) != NULL ){ 
        string s = string( ss ) ; 
        stringstream sin(s) ; 
        VD tmp ; 
        double x;
        while( sin >> x )
            tmp.push_back(x);
        input.push_back( tmp ) ; 
    } 
    fclose(fp) ;

    cout <<"Number of Samples: " <<  input.size() << " Number of Variabels: " << input[0].size() << endl; 
    srand( time(NULL) );
     
    for(int i = 0 ; i < SZ(input) ; ++i) 
       swap( input[i], input[ rand()% (i+1) ] ) ; 
   
    vector < VD > B = vector< VD > ( input[0].size(), VD(input[0].size(), 0.0)), 
                  W = vector< VD > ( input[0].size(), VD(input[0].size(), 0.0));
    double lambda1 = atof(parameter[1].second.c_str() ), lambda2 = atof( parameter[2].second.c_str() ); 
    int num = atoi(parameter[3].second.c_str());
    double thr = atof(parameter[4].second.c_str());

    swap( lambda1, lambda2); 
    solve_GSGBN(input, B , W, lambda1, lambda2,num, thr) ;
    fp = fopen( (parameter[0].second+"parameter").c_str() , "w") ;  
    assert( NULL != fp ) ;  
    for(int i = 0 ; i < SZ(B) ; ++i){
        for(int j = 0 ; j < SZ(B[0]) ; ++j){
            if(fabs(B[i][j]) < thr ) B[i][j] = 0. ; 
            fprintf( fp , "%lf " , B[i][j] ) ; 
        }
        fprintf( fp , "\n" ) ; 
    }
    fclose(fp) ; 
    fp = fopen( (parameter[0].second+"similarity").c_str(),"w") ; 
    assert(NULL!=fp); 
    for(int i=0; i < SZ(W); ++i){
        for(int j=0; j < SZ(W[0]); ++j){
            fprintf(fp, "%lf ", W[i][j]); 
        }
        fprintf(fp,"\n"); 
    }
    fclose(fp); 
    clock_t end = clock() ; 
    cout << ( end - start ) / CLOCKS_PER_SEC<<" Seconds" << endl ; 
    return 0 ;
}
