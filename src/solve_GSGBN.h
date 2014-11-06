#ifndef __SOLVE_GSGBN_H_
#define __SOLVE_GSGBN_H_ 

#include <cstdio>
#include <cstring>
#include <map>
#include <iostream>
#include <algorithm>
#include <queue>
#include <cmath>
#include <set>
#include <bitset>
#include <climits>
#include <utility>
using namespace std ;  


#define SZ(x) (int)x.size()

typedef vector<double> VD ;  
typedef long long int64 ; 

#define square(x) ((x)*(x))

bool isconverge(vector < VD > & preB , vector < VD > & B);

inline bool among( double x, double a, double b , double y);

double solve_nonternary( vector < pair< pair <double, double >, bool  > >& vc, double& val, int flag );

double nonternary( const vector < VD > & X ,const vector < VD > &B , const vector < VD > & W ,const  double & lambda1 ,
        const double& lambda2, const vector <VD> & hatX, int from, int to, bool none); 

void solve_coordinate_descent( const vector < VD > & X , vector < VD > &B,const vector< VD>  &W, const double &lambda1,
        const double& lambda2,vector<int >& order,const vector< vector< bool > > &skeleton,vector < VD > & hatX, double drop_zero, bool none);

void coordinate_descent( const vector < VD > & X, vector < VD > & B, const vector < VD>  & W, const double& lambda1,
        const double & lambda2, vector<int>& order,const vector< vector<bool> > & skeleton,bool none) ; 

bool isdecrease( vector < VD > & tmp , vector < VD > & W , vector < VD > & nextW,const double & lambda1);

void solve_projection( const vector < VD > & X, vector < VD > & B , vector < VD > & W, const double& lambda1, const double & lambda2);

double cal_score(const vector < VD > & X ,const vector < VD > & B ,const vector < VD > & W, const double&lambda1, const double& lambda2,bool print); 

void order_search( const vector< VD > &X , vector <VD> &B,const vector< VD > & W, 
        const double & lambda1, const double &lambda2, vector<int>& order, const vector< vector<bool> > & skeleton,bool none);

void estimate_initial_order(const vector< VD > & X, vector< VD > & B, vector< VD > & W,  const double & lambda1, const double & lambda2,
        vector< int > &order, vector< vector<bool> > & skeleton );

void solve_GSGBN( const vector < VD > & X , vector < VD > & B , vector < VD > & W, const double & lambda1, const double & lambda2,int num, double thr);

#endif
