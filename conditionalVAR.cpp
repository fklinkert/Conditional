#include <iostream>
#include <thread>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iomanip>



using namespace std;
using namespace std::chrono;

int N = 500;
int T = 15;
int MIN_VAR;



void selectionSort(int a[], int n) {
    int i, j, min, temp;
    for (i = 0; i < n - 1; i++) {
        min = i;
        for (j = i + 1; j < n; j++)
            if (a[j] < a[min])
                min = j;
        temp = a[i];
        a[i] = a[min];
        a[min] = temp;
    }
}

void simulation(double f[T+1], int simulationA[], int size){
    
    int i = 0;
    
    double random;
    srand ((unsigned)time(NULL));
    while(i < size){
        for(int j = 0; j < T; j++){
            
            if(j==0){
                simulationA[i] = N;
            }
            random = ((double) rand() / (RAND_MAX));
            if(random < (double)(1.0/3))
                simulationA[i] -= lambda1((int)(simulationA[i]-(double)(f[j]*simulationA[i])));
            else if(random < (double)(2.0/3)){
                simulationA[i] -= lambda2((int)(simulationA[i]-(double)(f[j]*simulationA[i])));
            }
            else{
                simulationA[i] -= lambda3((int)(simulationA[i]-(double)(f[j]*simulationA[i])));
            }
        }
        i++;
    }
    
    selectionSort(simulationA, size);
    f[T] = simulationA[(int)(size*0.99)-1];
}


void generate_f(double f[T+1]){
    
    double random;
    srand ((unsigned)time(NULL));
    for(int j = 0; j < T; j++){
        random = ((double) rand() / (RAND_MAX));
        f[j] = random;
    }
    
}



void worker(int n, int k, double range, int simulation_size, vector<vector<double>> *fvec) {
    
    double current_f[T+1];
    for(int i = 0; i < T+1;i++){
        current_f[i] = (*fvec)[n][i];
    }
    
    int simulation1[simulation_size];
    
    for(int i = 0; i < range; i++){
        //increment_f(current_f,step);
        generate_f(current_f);
        simulation(current_f, simulation1,simulation_size);
        
        if(current_f[T] < MIN_VAR){
            for( int j = 0; j < T+1; j++){
                (*fvec)[n][j] = current_f[j];
            }
            MIN_VAR = current_f[T];
            cout << "VaR value updated to " << MIN_VAR << endl;
        }
    
    }
    
    
}


void workers(int threads, double maxS, int simulation_size, int simulation1[], vector<vector<double> > *fvec){
    int range = maxS/threads;
    thread ths[threads];
    
    for(int k = 0; k < threads; k++){
        if(range <= maxS - k*range)
            ths[k] = thread(worker, k, range*k, range, simulation_size, fvec);
        else
            ths[k] = thread(worker, k, range*k, maxS-k*range, simulation_size, fvec);
    }
    for(int k = 0; k < threads; k++){
        ths[k].join();
    }
}


void VaR_Policy(int K, int T,int threads, double f[T+1]){
    
    int simulation_size = 1000;
    int simulation1[simulation_size];
    int step = 0.1;
    
    double current_f[T+1];
    vector<vector<double> > fvec(threads);
    for(int i = 0; i < T;i++){
        for(int j = 0; j < threads; j++)
            fvec[j].push_back(step);
        current_f[i] = step;
        f[i] = step;
    }
    current_f[T] = K;
    f[T] = K;
    MIN_VAR = K;
    
    double maxS = 10000;
    
    if(threads > 1){
        
        for(int j = 0; j < threads; j++)
            fvec[j].push_back(current_f[T]);
        
        f[T] = current_f[T];
        MIN_VAR = current_f[T];
        
        
        workers(threads, maxS, simulation_size, simulation1, &fvec);
        
        for(int k = 0; k < threads; k++){
            if(fvec[k][T] < f[T]){
                for(int i = 0; i < T+1; i++){
                    f[i]=fvec[k][i];
                }
            }
        }
    }
    else{
        for(int i = 0; i < maxS; i++){
            
            generate_f(current_f);
            simulation(current_f, simulation1,simulation_size);
            
            if(current_f[T] < MIN_VAR){
                for( int j = 0; j < T+1; j++){
                    f[j] = current_f[j];
                }
                MIN_VAR = f[T];
                cout << "VaR value updated to " << MIN_VAR << endl;
            }
        }
        
    }
}




int main(int argc, const char * argv[]){
    
    auto start = high_resolution_clock::now();
    
    int NumberThreads;
    if(argc != 2){
        NumberThreads = 2;
    }
    else{
        string arg = argv[1];
        std::istringstream iss (arg);
        iss >> NumberThreads;
    }
    
    double f[T+1];
    
    VaR_Policy(N, T, NumberThreads, f);
    
    
    cout << endl << "The optimal policy has the values of \n(";
    
    for(int i = 0; i < T; i++){
        cout << f[i];
        if(i != T-1)
            cout << ",";
    }
    cout << ") for a VaR value of "<< f[T] << endl;
    
    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
    
    cout << "Time taken by function: "
    << duration.count() /1000000.0 << " seconds" << endl;
    
    
    
    return 0;
}
