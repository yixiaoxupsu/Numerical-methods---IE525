//
//  Yixiao_Xu_hw10.cpp
//  Yixiao_Xu_hw10
//
//  Created by yixiao Xu on 11/30/17.
//  Copyright © 2017 yixiao Xu. All rights reserved.
//


//  Created by yixiao Xu on 11/14/17.
//  Copyright © 2017 yixiao Xu. All rights reserved.
//  Based on Professor's code in lesson 6

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include <vector>
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price, downtick_prob;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions, total_division;
//normal distribution cdf
void write_data(char* file_name, vector<double> &data)
{
    ofstream output_file(file_name);
    
    for (int i=0; i<data.size(); i++) {
        output_file<<data[i]<<endl;
    }
    // write code that writes "data" to file_name.
}
double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};

float max(float a, float b) {
    return (b < a )? a:b;
}
double option_price_put_black_scholes_closed_Form(const double& S,       // spot (underlying) price
                                                  const double& K,       // strike (exercise) price,
                                                  const double& r,       // interest rate
                                                  const double& sigma,   // volatility
                                                  const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    
    double vanilla_option_price = -S*N(-d1) + K*exp(-r*time)*N(-d2);
    return vanilla_option_price;
}
//void write_data(char* file_name, vector<double> &data)
//{
//    ofstream output_file(file_name);
//
//    for (int i=0; i<data.size(); i++) {
//        output_file<<data[i]<<endl;
//    }
//    // write code that writes "data" to file_name.
//}

double american_put_option_bbs(double s0,double k,double r,double t,double sigma,int N) {
    double stock_price[2*N+1];
    double option_value[2*N+1];
    stock_price[0]=s0*pow(1.0/up_factor, N);
    
    for (int i=1;i<2*N+1;i++){
        stock_price[i]=stock_price[i-1]*up_factor;
    }
    
    for (int i=1;i<2*N;i+=2){
        option_value[i]=max(option_price_put_black_scholes_closed_Form(stock_price[i], k, r, sigma, t/no_of_divisions),strike_price-stock_price[i]);
    }
    
    for (int i=2;i<N+1;i++){
        for (int j=i;j<2*N+1-i;j+=2){
            option_value[j]=max(strike_price-stock_price[j],(option_value[j-1]*downtick_prob+option_value[j+1]*uptick_prob)/R);
            
        }
        
    }
    
    return option_value[N];
}

double american_put_option(double s0,double k,double r,double t,double sigma,int N) {
    
    double stock_price[2*N+1];
    double option_value[2*N+1];
    stock_price[0]=s0*pow(1.0/up_factor, N);
    
    for (int i=1;i<2*N+1;i++){
        stock_price[i]=stock_price[i-1]*up_factor;
    }
    
    for (int i=0;i<2*N+1;i+=2){
        option_value[i]=max(strike_price-stock_price[i], 0);
    }
    
    for (int i=1;i<N+1;i++){
        for (int j=i;j<2*N+1-i;j+=2){
            option_value[j]=max(strike_price-stock_price[j],(option_value[j-1]*downtick_prob+option_value[j+1]*uptick_prob)/R);
            
        }
        
    }
    
    return option_value[N];
}

int main (int argc, char* argv[])
{
    
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &total_division);
    vector <double> put_prices;
    vector <double> deltas;
    vector  <double> put_bbs_price;
    vector <double> put_richard;
    for (int i=25; i<=total_division;i=i+500){
        
        no_of_divisions=i;
        up_factor= exp(volatility*sqrt(expiration_time / (double)no_of_divisions));
        R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
        uptick_prob = (R - 1/up_factor)/(up_factor - 1/up_factor);
        downtick_prob= 1-uptick_prob;
        //        cout<<no_of_divisions<<" "<<up_factor<<" "<<R<<" "<<uptick_prob<<endl;
        //
        
        double put_price = american_put_option(initial_stock_price,strike_price,risk_free_rate,expiration_time,volatility,no_of_divisions);
        double temp1=american_put_option(initial_stock_price*up_factor,strike_price,risk_free_rate,expiration_time,volatility,no_of_divisions-1);
        double temp2=american_put_option(initial_stock_price/up_factor,strike_price,risk_free_rate,expiration_time,volatility,no_of_divisions-1);
        double delta=(temp1-temp2)/((up_factor-1.0/up_factor)*initial_stock_price);
        double put_price_bbs=american_put_option_bbs(initial_stock_price,strike_price,risk_free_rate,expiration_time,volatility,no_of_divisions);
        
        no_of_divisions=i*2;
        up_factor= exp(volatility*sqrt(expiration_time / (double)no_of_divisions));
        R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
        uptick_prob = (R - 1/up_factor)/(up_factor - 1/up_factor);
        downtick_prob= 1-uptick_prob;
        double put_price_richard=2*american_put_option_bbs(initial_stock_price,strike_price,risk_free_rate,expiration_time,volatility,no_of_divisions)-put_price_bbs;
        put_prices.push_back(put_price);
        put_bbs_price.push_back(put_price_bbs);
        put_richard.push_back(put_price_richard);
        deltas.push_back(delta);
        
        
        
        
    }
    
    
    write_data("put_price", put_prices);
    write_data("put_price_bbs", put_bbs_price);
    write_data("put_price_richard", put_richard);
    write_data("delta", deltas);
    
    
    //
}

