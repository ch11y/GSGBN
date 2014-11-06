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

#include "solve_GSGBN.h"
using namespace std ;  


#define SZ(x) (int)x.size()
const double eps =  1e-4 ; 
static double threshold ; 
typedef vector<double> VD ;  
typedef long long int64 ; 

#define SZ(x) (int)x.size() 
#define square(x) ((x)*(x))

bool isconverge(vector < VD > & preB , vector < VD > & B){
    double res = 0 ; 
    for(int i=0 ; i < SZ(preB); ++i){
        for(int j=0; j < SZ(preB[0]); ++j){ 
            res = max( res, fabs(preB[i][j] - B[i][j])) ; 
        }
    }
    return res < eps ;  
}

inline bool among ( double x, double a, double b , double y){ 
    if( fabs(a) < eps || fabs(b-x) < eps || fabs( b- y) <eps ) return false;  
    return b >= x && b <= y ; 
} 

double solve_nonternary ( vector < pair< pair <double, double >, bool  > >& vc, double& val, int flag ){
    int n  = SZ(vc) ;  
    sort(vc.begin(), vc.end()) ; 
    n = SZ(vc); 
    vector< double > L = vector< double > ( n ), R = vector< double > ( n ) , AL = vector< double > ( n ) , AR = vector< double > ( n ), SUMBL = vector< double > ( n ), SUMBR = vector < double > ( n ); 
    int pos = 0, posq = 0  ; 
    for(int i = 0 ; i < n ; ++i){ 
        if( vc[ i ].second ) {
            pos = i ;
            posq = i; 
            break ; 
        }
    }
    AL[0]=L[ 0 ] = 0  ; 
    double sum = 0 ; 
    for(int i = 1; i < n ; ++i){ 
       AL[i] = 0 ; 
       sum += !vc[i-1].second?vc[ i - 1].first.second:0.; 
       L[ i ] = L[i-1] + sum  * (vc[i].first.first-  vc[ i -1 ].first.first) ; 
       AL[i] += (pos<i)? vc[pos].first.second*square( vc[ pos].first.first - vc[i].first.first):0. ;
    } 
    sum =  0 ; 
    AR[n-1]=R[ n - 1 ] = 0; 
    for(int i = n-2; i >= 0 ; --i){ 
        AR[ i ] = 0 ;
        sum += !vc[i+1].second?vc[i+1].first.second:0. ;  
        R[ i ] = R[ i+ 1] + sum * (vc[i + 1].first.first - vc[ i].first.first )  ;  
        AR[ i ] += (pos>i)?vc[pos].first.second*square( vc[ pos].first.first - vc[ i ].first.first ):0.; 
    } 
    double least_val = -1000000.;
    double least = 1000000.0;
    for(int i = 0 ; i < n ; ++i){ 
        if(flag<0&& vc[ i ].second ) continue;
        double y = L[i] + R[ i ] + AL[i] +AR[i];
        if( y < least ) { 
            least=y ; 
            least_val = vc[i].first.first; 
        }
    } 
    double suma =2* vc[posq].first.second ; 
    for(int i=0; i < n; ++i){
        SUMBL[i] = !vc[i].second?vc[i].first.second:0. ;  
        SUMBL[i] +=i?SUMBL[ i - 1]:0. ; 
    }
    for(int i= n - 1; i >= 0 ; --i){  
        SUMBR[ i ] = !vc[i].second?(-vc[i].first.second):0. ;  
        SUMBR[ i ] += (i!=n-1)?SUMBR[ i + 1]:0;  
    } 
    double cut_val =-100000.;  
    for(int i=1 + (flag<0); i < n ; ++i){ 
      if( among ( vc[ i -1 ].first.first ,  suma,  vc[posq].first.first - (SUMBL[i-1]+SUMBR[i])/suma , vc[i].first.first ) ){ 
           cut_val=vc[ posq].first.first-(SUMBL[i-1]+SUMBR[i])/suma ;
           break;
      }
    }
    val = (cut_val>-100.)? cut_val:least_val;
    double ans = 0 ; 
    for(int i= 0 ; i < n ; ++i)
        ans += !vc[i].second?vc[ i ].first.second * abs( val - vc[i].first.first ):vc[ i ].first.second * square( val - vc[ i ].first.first ); 
    return ans; 
} 

double nonternary( const vector < VD > & X ,const vector < VD > &B , const vector < VD > & W ,const  double & lambda1 , const double& lambda2, const vector <VD> & hatX 
        , int from, int to, bool none = false  ){
       int n = SZ(X) , p = SZ(X[0]);  
       double a = 0 , b = 0; 
       vector < double > hat(n);  
       for(int i = 0 ; i < n  ; ++i){
           hat[i] = hatX[i][to]-X[i][from]*B[from][to]-X[ i ][ to ] ; 
           a += square(X[i][from]);  
           b += X[i][from ] * hat[ i ] ; 
       }
       b*=2;
       if(none){
           double res = -b /(2*a) ;
           return (fabs(res)< threshold)?0:res;
       }

       int flag = 1; 
       if( b > 0 && fabs(b) > eps ) {
           b *= -1 ; 
           flag = -1;
       } 
       vector < pair< pair < double , double >, bool > > vc;  
       vc.push_back(make_pair ( make_pair(-b/(2*a) ,a ) ,  true )) ;
       vc.push_back(make_pair( make_pair( 0., lambda2 ) , false)) ;
       for(int i= 0 ; i < p ; ++i){
           if( i == to ) continue;
           if(fabs(B[from][i]) < eps) vc[1].first.second += 2*lambda1*W[to][i]; 
           else vc.push_back(make_pair( make_pair( abs(B[from][i]), 2*lambda1 * W[ to ][ i ] ) , false)) ; 
       }
       double val1, val2 ; 
       double res1= solve_nonternary( vc, val1, fabs(flag));
       val1 *= flag; 
       if( flag > 0 ) 
           return val1;
       for(int i = 0 ; i < SZ(vc) ; ++i){ 
           if( vc[i].second == true ){ 
               vc[i].first.first *= -1; 
               break;
           }
       }
       double res2 = solve_nonternary( vc, val2, flag );  
       if( res1 > res2 ) { 
           swap( res1, res2);
           swap( val1, val2); 
       }
       return val1; 
}

void solve_coordinate_descent( const vector < VD > & X , vector < VD > &B,const vector< VD>  &W, const double &lambda1, const double& lambda2,
        vector<int >& order,const vector< vector< bool > > &skeleton,vector < VD > & hatX, double drop_zero, bool none = false){
    int n = SZ(X) , p = SZ(X[0]);
    vector <VD> preB = B; 
    for(int i=0 ; i < p ; ++i){ 
        double diff=0;
        int tot = 0 ; 
        do{
            ++tot;
            diff = 0; 
            for(int j=0; j < p ; ++j){
                if(i == j || order[j] > order[i] ||!skeleton[j][i]) continue;
                double tmp_ji = B[j][i];
                if( fabs(tmp_ji) < eps && 1.*rand() / INT_MAX > drop_zero) 
                    continue;
                double ans_ji = nonternary( X , B, W , lambda1, lambda2, hatX, j , i ,none);
                if( fabs(ans_ji) < threshold ) ans_ji=0;
                if(fabs(ans_ji-tmp_ji)<eps) continue;
                diff = max(diff, fabs(ans_ji - tmp_ji)); 
                B[j][i] = ans_ji; 
                for(int k=0; k < n ; ++k) 
                    hatX[k][i] += X[k][j]*(B[j][i]-tmp_ji); 
            }
        }while(!(diff<eps)&&tot < p); 
    } 
} 

void coordinate_descent( const vector < VD > & X, vector < VD > & B, const vector < VD>  & W,
    const double& lambda1, const double & lambda2, vector<int>& order,const vector< vector<bool> > & skeleton,bool none=false){
     
    int n = SZ(X),  p = SZ(X[0]);
    vector <VD> preB;
    int tot = p*p;
    int number_interation = 0 ;
    double drop_zero = 0. ;
    for(int i = 0 ; i < p; ++i){
        for(int j=0; j < p; ++j){ 
            if(order[i] < order[j] ) 
                B[j][i]=0; 
        }
    }
    vector < VD > hatX = X ;
    for(int i = 0 ; i < p ; ++i){ 
       for(int j=0;  j < n; ++j){ 
            hatX[j][i] = 0 ;  
            for(int k= 0 ; k < p ; ++k){
                if(order[k]>order[i]) continue;
                hatX[j][i] +=  X[j][k]*B[k][i]; 
            }
       }
    }
    do{
        preB = B ;
        solve_coordinate_descent(X, B, W , lambda1, lambda2, order, skeleton,hatX,exp(-drop_zero),none); 
        ++number_interation; 
        drop_zero += 0.25; 
    }while(number_interation <=tot && !isconverge( preB,B)) ;
}

bool isdecrease( vector < VD > & tmp , vector < VD > & W , vector < VD > & nextW,const double & lambda1){ 
    double res1 = 0 , res2 =  0 ; 
    for(int i = 0 ; i < SZ(tmp) ; ++i){ 
        for(int j=0 ; j < SZ(tmp) ; ++j){  
            if( i == j ) continue ;  
            res1 += tmp[i][j] * W[ i ][ j ] * lambda1 ;  
            res1 -=  log( W[i][j] ) ;
            res2 += tmp[i][j] * nextW[ i ][ j ] * lambda1 ; 
            res2 -=  log( nextW[i][j]) ; 
        }
    }
    return res1 - res2 > eps ; 
}


void solve_projection( const vector < VD > & X, vector < VD > & B , vector < VD > & W, 
        const double& lambda1, const double & lambda2){ 
    int p = SZ(X[0]) ;
    for(int i = 0 ; i <  p; ++i){
        for(int j = 0 ; j < p ; ++j){
            if( i == j ) continue ; 
            W[i][j]=1./(p-1); 
        }
    }

    double check = 0 ;
    for(int i = 0 ; i < p ; ++i)
        for(int j=0; j < p ; ++j)
            check += W[i][j];
    double delta[ 100+10] ;
    delta[ 0 ] = 0.1; 
    for(int i=1; i < 110 ; ++i)
        delta[ i ] =delta[i-1]* 0.95; 
    vector < VD > tmp = B ; 
    for(int i = 0 ; i < p ; ++i){
        for(int j = 0 ; j < p ; ++j){ 
            tmp[i][j] = 0  ;
            for(int k = 0 ; k < p ; ++k){  
                tmp[i][j] += abs( abs(B[k][i]) - abs(B[k][j])) ; 
            }
        }
    } 
    vector < vector < bool > > used ( SZ(W));
    for(int step = 0 ; step < 110; ++step){
       bool legal = true; 
       vector < VD > nextW = W ;
       int cnt = 0 ;
       do{
           ++cnt;  
           W = nextW ;  
           double sum_W = 0  ;
           double sub_W = 0 ;
           int cnt_W = 0 ; 
           for(int i = 0 ; i < p ; ++i){ 
               vector < bool > tmp( p , false) ; 
               used[ i ] = tmp ;
           }
           double res=0; 
           for(int i = 0 ; i < p ; ++i){
               for(int j = 0 ; j < p; ++j){
                   if(i==j) continue; 
                   res += exp(-W[i][j]); 
               }
           }
           for(int i = 0 ; i < p ; ++i){
               for(int j = 0 ;j < p ; ++j){ 
                   if( i == j )continue;
                   if(W[i][j] > eps){
                       nextW[i][j] -= delta[step] *  lambda1 * tmp[i][j] - delta[step]/W[i][j];  
                       ++cnt_W; 
                       used[ i ][ j ] = 1; 
                   }
                   else { 
                       sub_W += W[i][j];   
                   }
                   sum_W += nextW[i][j] ; 
               }
           }
           sum_W -= (p-sub_W)  ;
           double check = 0 ; 
           for(int i = 0 ; i < p ; ++i){
               for(int j = 0 ; j < p ; ++j){
                   if( i != j ){
                   if( used[i][j]) 
                   nextW[i][j]  -= (sum_W-sub_W) /cnt_W ;
                   check += nextW[i][j]; 
                   }
               }  
           } 
           for(int i = 0 ; i < p ;++i){
               for(int j=i+1 ; j < p ; ++j){ 
                   if(nextW[i][j] <= 0 || nextW[i][j] >=p ){
                       legal=0;
                       break ;
                   }
               }
           }
           for(int i = 0 ; i < p ; ++i){
               double tmp = 0 ;
               for(int j = 0 ; j < p ; ++j)
                   tmp += nextW[ i ][ j ] ; 
               tmp *= 2; 
               if( tmp *lambda1> lambda2 ) {
                   legal = 0 ; 
                   break; 
               }
           }
       }while( cnt < 500 &&legal && isdecrease( tmp, W, nextW,lambda1)) ; 
    }
}

double cal_score(const vector < VD > & X ,const vector < VD > & B ,const vector < VD > & W, 
        const double&lambda1, const double& lambda2,bool print=false){
    double now = 0 ; 
    int p = SZ(B);
    double w = 0 ; 
    for(int i = 0; i < SZ(X) ; ++i){
        for(int j=0; j < p ; ++j){ 
            double tmp = 0 ; 
            for(int k= 0 ; k < p ; ++k)
                tmp += X[i][k]*B[k][j];
            now += square(tmp - X[i][j]); 
        }   
    }
    if(print) cout <<"Squared Error: " << now << endl;
    double now1 = 0 , now2 = 0 ;
    for(int i = 0 ; i < p ; ++i){
        for(int j=0 ; j < p ; ++j){
            if( i == j ) continue; 
            now += lambda2 * abs( B[j][i]) ; 
            now2 +=lambda2 * abs( B[j][i]) ; 
            for(int k=0 ; k < p ; ++k){ 
                now += lambda1*W[i][j]*abs(abs(B[k][i])- abs(B[k][j])) ;
                now1 += lambda1 * W[i][j] * abs(abs(B[k][i] - abs(B[k][j]))) ; 
            }
            w -= log(W[i][j]); 
        }
    }
    if(print)  cout <<"Penalties: "<<  now1 << "  " << now2<< " W: "<< w << endl; 
    return now;
} 
void order_search( const vector< VD > &X , vector <VD> &B,const vector< VD > & W, 
        const double & lambda1, const double &lambda2, vector<int>& order, const vector< vector<bool> > & skeleton,bool none=true){ 
    for(int i = 0; i < SZ(B); ++i)
        for(int j=0; j < SZ(B); ++j)
            if(order[i]>order[j])
                 B[i][j]=0.;
    coordinate_descent(X, B ,W, lambda1, lambda2, order,skeleton);
    double best_score = cal_score(X,B,W,lambda1, lambda2); 
    int number_interation = 0 ; 
    vector<int> state(B.size()) ;
    for(int i = 0; i < SZ(state); ++i)
        state[i] = i ;
    int cnt = 1 ; 
    int unchange=0 ;
    bool dec = false; 
    do{
        int indexi = -1,indexj=-1; 
        vector < VD > nowB = B ; 
        vector < VD > preB = B ;
        
        for(int i = 0 ; i < SZ(B); ++i){
            for(int j = i + 1; j < SZ(B); j++){ 
                swap( order[state[i]],order[state[j]]); 
                B = preB; 
                coordinate_descent (X,B,W,lambda1,lambda2, order,skeleton,none); 
                double now_score = cal_score( X, B , W ,lambda1, lambda2); 
                if( fabs(now_score - best_score) > eps &&  now_score < best_score ){ 
                    best_score = now_score ;
                    nowB = B ; 
                    indexi = i ;
                    indexj = j ;
                    cnt = 1;
                    dec = true; 
                }
                swap(order[state[i]],order[state[j]]);
                break;
            }
        }
        if(dec) unchange = 0 ; 
        else ++unchange;
        dec = false;
        B = nowB; 
        if( indexi != -1 ){ 
            swap( order[state[indexi]], order[state[indexj]]);  
            coordinate_descent( X , B , W , lambda1, lambda2, order, skeleton); 
        }else{ 
            for(int i = 0 ; i < SZ(state); ++i)
                swap(state[i], state[ rand() % ( i + 1 ) ] ); 
            ++number_interation;
        }
    }while( number_interation < 5&&unchange <= 20 );  
} 

void estimate_initial_order(const vector< VD > & X, vector< VD > & B, vector< VD > & W,  const double & lambda1, const double & lambda2, vector< int > &order, vector< vector<bool> > & skeleton ){
    vector< VD > & orgB = B;
    int p = B.size(); 
    order = vector< int >(p);
    for(int i = 0 ; i < SZ(order); ++i){ 
        order[i] = i;
        swap(order[i], order[ rand() % ( i + 1 )]);
    }
    vector< int >  best_order = order;
    double best_score = 10000000000000.; 
    for(int i = 0 ;  i < p*p ; ++i){
       orgB = B; 
       coordinate_descent( X, orgB , W, lambda1, lambda2, order, skeleton,true); 
       double cur_score = cal_score(X,orgB,W,lambda1,lambda2); 
       if(cur_score < best_score){
           best_score = cur_score; 
           best_order = order; 
       }
       for(int j = 0 ; j < p ; ++j)
           swap(order[j], order[ rand()%(j+1)]);
    }
    order = best_order; 
}

void solve_GSGBN( const vector < VD > & X , vector < VD > & B , vector < VD > & W, const double & lambda1, const double & lambda2,int num, double thr){
    int p = SZ(B) ;   
    for(int i= 0 ; i < SZ(W) ; ++i){
        for(int j= 0 ; j < SZ(W[0]);++j){
            W[i][j] = 0 ;
            if( i == j ) continue ; 
            W[ i ] [ j ] = 1./(p-1) ; 
        } 
    }
    threshold = thr; 
    vector <VD> preB;
    int number_interation = 0;
    vector< int > order(SZ(B),0); 
    vector< vector<bool> > skeleton( SZ(B)) ;
    for(int i=0 ;i < SZ(B); ++i)
        skeleton[i] = vector< bool > (SZ(B), true);   
    do{ 
         preB = B ;
         ++number_interation ; 
         vector<VD> preW = W ; 
         coordinate_descent( X , B , W , lambda1, lambda2, order,skeleton); 
         solve_projection( X , B ,W , lambda1, lambda2) ;
     }while( number_interation < 20 && !isconverge( B , preB) ); 
     for(int i = 0 ; i < SZ(B); ++i){
         for(int j = 0 ; j < SZ(B); ++j){ 
             if( i >= j ) continue; 
             double res = 0 ;
             for(int k = 0 ; k < SZ(B); ++k)
                 res +=  abs(abs(B[k][i])-abs(B[k][j])); 
         }
     }
     vector< vector< bool > > skeleton_no = skeleton ;  
     double best_score = 10000000000.;   
     vector< VD > bestB = B ; 
     vector< vector< bool > > skeleton_have = skeleton ;
     number_interation = 0 ;
     vector< int > bestorder ;
     srand( time(NULL) ); 
     vector< VD > oldB = B ;
     do{ 
         for(int i = 0 ; i < p ; ++i)
             for(int j = i + 1 ; j < p ; ++j)
                 skeleton_have[i][j] = skeleton_have[j][i] = ( fabs(bestB[i][j])+ fabs( bestB[j][i])) > threshold; 
         skeleton = skeleton_have; 
         estimate_initial_order(X,B, W, lambda1,lambda2,order, skeleton); 
         skeleton = skeleton_no ; 
         order_search(X,B,W,lambda1,lambda2,order,skeleton, false); 
         double now_score = cal_score(X,B,W,lambda1,lambda2);  
         if( now_score < best_score ){ 
             bestB = B ;  
             best_score = now_score ;
             coordinate_descent( X, B , W, lambda1, lambda2, order, skeleton); 
             cal_score(X,B,W,lambda1,lambda2);bestorder = order ; 
         }
         ++number_interation ;  
         cout <<"Finised Interations: " << number_interation << endl; 
     }while( number_interation < num) ;
     B = bestB;
     order = bestorder; 
     
     for(int i = 0 ; i < SZ(B); ++i){
         for(int j=0; j < SZ(B); ++j){
             if(order[i] > order[j]) skeleton[i][j] = false; 
             else skeleton[i][j] = fabs(B[i][j])>threshold;
         } 
     }
     coordinate_descent(X,B,W,lambda1,lambda2,bestorder, skeleton, true); 
}
