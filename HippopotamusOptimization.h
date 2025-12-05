#ifndef HIPPOPOTAMUS_OPTIMIZATION_H
#define HIPPOPOTAMUS_OPTIMIZATION_H

#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
using namespace std;

template<typename TYPE>
class HippopotamusOptimization {
private:
    int searchAgents;
    int maxIterations;
    int dimension;
    vector<TYPE> lowerBound;
    vector<TYPE> upperBound;
    function<TYPE(const vector<TYPE>&)> fitnessFunction;
    
    random_device rd;
    mutable mt19937 gen;
    mutable uniform_real_distribution<TYPE> uniformDist;
    mutable normal_distribution<TYPE> normalDist;
    
    vector<vector<TYPE>> population;
    vector<TYPE> fitness;
    vector<TYPE> bestPosition;
    TYPE bestScore;
    vector<TYPE> convergenceCurve;
    
    vector<TYPE> levy(int n, int m, TYPE beta) const;
    void initializePopulation();
    void evaluatePopulation();
    void updateBest(int iteration);
    void phase1Exploration();
    void phase2Defense();
    void phase3Exploitation(int iteration);
    void boundaryCheck(vector<TYPE>& solution) const;
    
public:
    HippopotamusOptimization(int searchAgents, int maxIterations, 
                           const vector<TYPE>& lowerBound,
                           const vector<TYPE>& upperBound,
                           function<TYPE(const vector<TYPE>&)> fitnessFunc);
    
    void optimize();
    
    TYPE getBestScore() const { return bestScore; }
    const vector<TYPE>& getBestPosition() const { return bestPosition; }
    const vector<TYPE>& getConvergenceCurve() const { return convergenceCurve; }
    
    void displayResults() const;
};

template<typename TYPE>
HippopotamusOptimization<TYPE>::HippopotamusOptimization(
    int searchAgents, int maxIterations,
    const vector<TYPE>& lowerBound,
    const vector<TYPE>& upperBound,
    function<TYPE(const vector<TYPE>&)> fitnessFunc)
    : searchAgents(searchAgents), maxIterations(maxIterations),
      lowerBound(lowerBound), upperBound(upperBound),
      fitnessFunction(fitnessFunc), gen(rd()),
      uniformDist(0.0, 1.0), normalDist(0.0, 1.0) {
    
    dimension = lowerBound.size();
    population.resize(searchAgents, vector<TYPE>(dimension));
    fitness.resize(searchAgents);
    bestPosition.resize(dimension);
    bestScore = numeric_limits<TYPE>::max();
    convergenceCurve.reserve(maxIterations);
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::initializePopulation() {
    for(int i = 0; i < searchAgents; ++i){
        for(int j = 0; j < dimension; ++j){
            population[i][j] = lowerBound[j] + uniformDist(gen) * (upperBound[j] - lowerBound[j]);
        }
    }
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::evaluatePopulation() {
    for(int i = 0; i < searchAgents; ++i)
        fitness[i] = fitnessFunction(population[i]);
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::updateBest(int iteration) {
    auto minIt = min_element(fitness.begin(), fitness.end());
    int bestIndex = distance(fitness.begin(), minIt);
    
    if(iteration == 0 || *minIt < bestScore){
        bestScore = *minIt;
        bestPosition = population[bestIndex];
    }
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::phase1Exploration() {
    uniform_int_distribution<int> intDist1(1, 2);
    uniform_int_distribution<int> intDist2(0, 1);
    uniform_int_distribution<int> randGroupDist(1, searchAgents);
    uniform_int_distribution<int> alfaDist(1, 5);
    
    for(int i = 0; i < searchAgents / 2; ++i){
        int I1 = intDist1(gen);
        int I2 = intDist1(gen);
        int Ip1_1 = intDist2(gen);
        int Ip1_2 = intDist2(gen);
        int randGroupNum = randGroupDist(gen);
        
        // Create random group
        vector<int> indices(searchAgents);
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), gen);
        vector<int> randGroup(indices.begin(), indices.begin() + randGroupNum);
        
        // Calculate mean of random group
        vector<TYPE> meanGroup(dimension, 0.0);
        for(int idx : randGroup){
            for(int j = 0; j < dimension; ++j)
                meanGroup[j] += population[idx][j];
        }
        for(int j = 0; j < dimension; ++j)
            meanGroup[j] /= randGroup.size();
        
        // Create alfa vectors
        vector<vector<TYPE>> alfa(5, vector<TYPE>(dimension));

        TYPE alfa4 = uniformDist(gen);
        for(int j = 0; j < dimension; ++j){
            alfa[0][j] = I2 * uniformDist(gen) + (1 - Ip1_1);
            alfa[1][j] = 2 * uniformDist(gen) - 1;
            alfa[2][j] = uniformDist(gen);
            alfa[3][j] = I1 * uniformDist(gen) + (1 - Ip1_2);
            alfa[4][j] = alfa4;
        }
        
        vector<TYPE> A = alfa[alfaDist(gen) - 1];
        vector<TYPE> B = alfa[alfaDist(gen) - 1];
        
        // Phase 1 updates
        vector<TYPE> X_P1(dimension), X_P2(dimension);
        
        TYPE r1 = uniformDist(gen);
        for(int j = 0; j < dimension; ++j)
            X_P1[j] = population[i][j] + r1 * (bestPosition[j] - I1 * population[i][j]);
        
    TYPE T = exp(-static_cast<TYPE>(convergenceCurve.size()) / maxIterations); 
	if(T>0.6)
            for(int j = 0; j < dimension; ++j)
                X_P2[j] = population[i][j] + A[j] * (bestPosition[j] - I2 * meanGroup[j]);
        else{
            if(uniformDist(gen) > 0.5)
                for(int j = 0; j < dimension; ++j)
                    X_P2[j] = population[i][j] + B[j] * (meanGroup[j] - bestPosition[j]);
            else
                for(int j = 0; j < dimension; ++j)
                    X_P2[j] = (upperBound[j] - lowerBound[j]) * uniformDist(gen) + lowerBound[j];
        }
        
        boundaryCheck(X_P2);
        boundaryCheck(X_P1);

        TYPE F_P1 = fitnessFunction(X_P1);
        if(F_P1 < fitness[i]){
            population[i] = X_P1;
            fitness[i] = F_P1;
        }
        
        TYPE F_P2 = fitnessFunction(X_P2);
        if(F_P2 < fitness[i]){
            population[i] = X_P2;
            fitness[i] = F_P2;
        }
    }
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::phase2Defense() {
    uniform_real_distribution<TYPE> bDist(2.0, 4.0);
    uniform_real_distribution<TYPE> cDist(1.0, 1.5);
    uniform_real_distribution<TYPE> dDist(2.0, 3.0);
    uniform_real_distribution<TYPE> lDist(-2.0 * M_PI, 2.0 * M_PI);
    
    for(int i = searchAgents / 2; i < searchAgents; ++i){
        vector<TYPE> predator(dimension);
        for(int j = 0; j < dimension; ++j)
            predator[j] = lowerBound[j] + uniformDist(gen) * (upperBound[j] - lowerBound[j]);
        
        TYPE F_HL = fitnessFunction(predator);
        
        vector<TYPE> distance2Leader(dimension);
        for(int j = 0; j < dimension; ++j)
            distance2Leader[j] = abs(predator[j] - population[i][j]);
        
        TYPE b = bDist(gen);
        TYPE c = cDist(gen);
        TYPE d = dDist(gen);
        TYPE l = lDist(gen);
        
        auto RL = levy(searchAgents, dimension, 1.5);
        for(auto &x : RL)  x *= 0.05;

        vector<TYPE> X_P3(dimension);

        if(fitness[i] > F_HL){
            for(int j = 0; j < dimension; ++j)
                X_P3[j] = RL[i * dimension + j] * predator[j] + 
                         (b / (c - d * cos(l))) * (1.0 / (distance2Leader[j]+1e-9));
        } 
        else {
            for(int j = 0; j < dimension; ++j)
                X_P3[j] = RL[i * dimension + j] * predator[j] + 
                         (b / (c - d * cos(l))) * (1.0 / (2.0 * distance2Leader[j] + uniformDist(gen)));
        }
        
        boundaryCheck(X_P3);
        
        TYPE F_P3 = fitnessFunction(X_P3);
        if(F_P3 < fitness[i])
            population[i] = X_P3, fitness[i] = F_P3;
    }
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::phase3Exploitation(int iteration) {
    TYPE t = static_cast<TYPE>(iteration + 1);
    
    for(int i = 0; i < searchAgents; ++i){
        vector<TYPE> LO_LOCAL(dimension), HI_LOCAL(dimension);
        for(int j = 0; j < dimension; ++j){
            LO_LOCAL[j] = lowerBound[j] / t;
            HI_LOCAL[j] = upperBound[j] / t;
        }
        
        uniform_int_distribution<int> alfaDist(1, 3);
        vector<vector<TYPE>> alfa(3, vector<TYPE>(dimension));
        
        TYPE alfa1_scalar = uniformDist(gen);
        TYPE alfa2_scalar = normalDist(gen);

        for(int j = 0; j < dimension; ++j){
            alfa[0][j] = 2.0 * uniformDist(gen) - 1.0; 
            alfa[1][j] = alfa1_scalar;          
            alfa[2][j] = alfa2_scalar;          
        }
        
        vector<TYPE> D = alfa[alfaDist(gen) - 1];
        vector<TYPE> X_P4(dimension);
        
        TYPE r3 = uniformDist(gen);
        for(int j = 0; j < dimension; ++j)
            X_P4[j] = population[i][j] + r3 * (LO_LOCAL[j] + D[j] * (HI_LOCAL[j] - LO_LOCAL[j]));
        
        boundaryCheck(X_P4);
        
        TYPE F_P4 = fitnessFunction(X_P4);
        if(F_P4 < fitness[i]){
            population[i] = X_P4;
            fitness[i] = F_P4;
        }
    }
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::boundaryCheck(vector<TYPE>& solution) const {
    for(int j = 0; j < dimension; ++j){
        solution[j] = min(max(solution[j], lowerBound[j]), upperBound[j]);
    }
}

template<typename TYPE>
vector<TYPE> HippopotamusOptimization<TYPE>::levy(int n, int m, TYPE beta) const {
    TYPE num = std::tgamma(1.0 + beta) * sin(M_PI * beta / 2.0);
    TYPE den = std::tgamma((1.0 + beta) / 2.0) * beta * pow(2.0, (beta - 1.0) / 2.0);
    TYPE sigma_u = pow(num / den, 1.0 / beta);

    vector<TYPE> z(n * m);
    const TYPE eps = static_cast<TYPE>(1e-10);
    TYPE max_step = static_cast<TYPE>(1e6); // exemplo de clamp

    for(int i = 0; i < n * m; ++i){
        TYPE u = normalDist(gen) * sigma_u;
        TYPE v = normalDist(gen);
        TYPE denom = pow(std::max(std::abs(v), eps), static_cast<TYPE>(1.0 / beta));
        z[i] = u / denom;
        z[i] = std::max(std::min(z[i], max_step), -max_step);
    }
    return z;
}


template<typename TYPE>
void HippopotamusOptimization<TYPE>::optimize() {
    initializePopulation();
    evaluatePopulation();
    
    for(int iteration = 0; iteration < maxIterations; ++iteration){
        updateBest(iteration);
        
        phase1Exploration();
        phase2Defense();
        phase3Exploitation(iteration);
        
        convergenceCurve.push_back(bestScore);
        
        cout << "Iteration " << iteration + 1 << ": Best Cost = " 
                  << bestScore << endl;
    }
}

template<typename TYPE>
void HippopotamusOptimization<TYPE>::displayResults() const {
    cout << "\n=== Results ===" << endl;
    cout << "Best Score: " << bestScore << endl;
    cout << "Best Position: [";
    cout << setprecision(6) << fixed;
    for(size_t i = 0; i < bestPosition.size(); ++i){
        cout << bestPosition[i];
        if(i < bestPosition.size() - 1) cout << ", ";
    }
    cout << "]" << endl;
}

#endif // HIPPOPOTAMUS_OPTIMIZATION_H
