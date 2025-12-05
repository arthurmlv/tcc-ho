#include "HippopotamusOptimization.h"
#include "TestFunctions.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <limits>
#include <bits/stdc++.h>
using namespace std;

// Custom optimization problem example
template<typename TYPE>
TYPE customFunction(const vector<TYPE>& x) {
    // Example: Minimize sum of squares with penalty
	TYPE sum = 0;
	TYPE m = 10;
	TYPE i = 1;
	for (const auto &xi : x) {
		sum += sin(xi) * pow(sin((i * xi * xi) / std::numbers::pi), 2 * m);
		i += 1;
	}
	return -sum;
}

void clearInputBuffer() {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
}

int getUserChoice(int min, int max) {
    int choice;
    while(true){
        cout << "Enter your choice (" << min << "-" << max << "): ";
        if(cin >> choice && choice >= min && choice <= max){
            clearInputBuffer();
            return choice;
        } else {
            cout << "Invalid input! Please enter a number between " << min << " and " << max << "." << endl;
            clearInputBuffer();
        }
    }
}

struct OptimizationParams {
    int searchAgents;
    int maxIterations;
};

OptimizationParams getOptimizationParams() {
    OptimizationParams params;
    
    cout << "\n--- Optimization Parameters ---" << endl;
    
    while(true){
        cout << "Enter number of search agents (population size) [10-100]: ";
        if(cin >> params.searchAgents && clamp(params.searchAgents,10,100)==params.searchAgents)
            break;
        else{
            cout << "Invalid input! Please enter a number between 10 and 100." << endl;
            clearInputBuffer();
        }
    }
    
    while(true){
        cout << "Enter maximum iterations [50-2000]: ";
        if(cin >> params.maxIterations && clamp(params.maxIterations,50,2000)==params.maxIterations)
            break;
        else{
            cout << "Invalid input! Please enter a number between 50 and 2000." << endl;
            clearInputBuffer();
        }
    }
    
    clearInputBuffer();
    return params;
}

void runTestFunction() {
    cout << "\n=== Available Test Functions ===" << endl;
    cout << "1.  F1  - Sphere Function" << endl;
    cout << "2.  F2  - Schwefel's Problem 2.22" << endl;
    cout << "3.  F3  - Schwefel's Problem 1.2" << endl;
    cout << "4.  F4  - Schwefel's Problem 2.21" << endl;
    cout << "5.  F5  - Rosenbrock Function" << endl;
    cout << "6.  F6  - Step Function" << endl;
    cout << "7.  F7  - Quartic Function with Noise" << endl;
    cout << "8.  F8  - Schwefel Function" << endl;
    cout << "9.  F9  - Rastrigin Function" << endl;
    cout << "10. F10 - Ackley Function" << endl;
    cout << "11. F11 - Griewank Function" << endl;
    cout << "12. F12 - Penalized Function 1" << endl;
    cout << "13. F13 - Penalized Function 2" << endl;
    cout << "14. F14 - Foxholes Function" << endl;
    cout << "15. F15 - Kowalik Function" << endl;
    cout << "16. F16 - Six-Hump Camel Function" << endl;
    cout << "17. F17 - Branin Function" << endl;
    cout << "18. F18 - Goldstein-Price Function" << endl;
    cout << "19. F19 - Hartman 3 Function" << endl;
    cout << "20. F20 - Hartman 6 Function" << endl;
    cout << "21. F21 - Shekel 5 Function" << endl;
    cout << "22. F22 - Shekel 7 Function" << endl;
    cout << "23. F23 - Shekel 10 Function" << endl;
    
    int functionChoice = getUserChoice(1, 23);
    string functionName = "F" + to_string(functionChoice);
    
    cout << "\n--- Data Type Selection ---" << endl;
    cout << "1. Double precision (recommended)" << endl;
    cout << "2. Float precision (faster)" << endl;
    
    int precisionChoice = getUserChoice(1, 2);
    
    OptimizationParams params = getOptimizationParams();
    
    cout << "\n--- Starting Optimization ---" << endl;
    cout << "Function: " << functionName << endl;
    cout << "Search Agents: " << params.searchAgents << endl;
    cout << "Max Iterations: " << params.maxIterations << endl;
    cout << "Precision: " << (precisionChoice == 1 ? "Double" : "Float") << endl;
    
    auto start = chrono::high_resolution_clock::now();
    
    if(precisionChoice == 1){
        auto funcInfo = TestFunctions<double>::getFunctionInfo(functionName);
        cout << "Function Name: " << funcInfo.name << endl;
        cout << "Dimension: " << funcInfo.dimension << endl;
        
        HippopotamusOptimization<double> optimizer(
            params.searchAgents,
            params.maxIterations,
            funcInfo.lowerBound,
            funcInfo.upperBound,
            funcInfo.function
        );
        
        optimizer.optimize();
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        
        optimizer.displayResults();
        cout << "Execution time: " << duration.count() << " ms" << endl;
        
        cout << "\nDo you want to save the convergence curve? (y/n): ";
        char saveChoice;
        cin >> saveChoice;
        if(saveChoice == 'y' || saveChoice == 'Y'){
            const auto& curve = optimizer.getConvergenceCurve();
            string filename = functionName + "_convergence.txt";
            ofstream file(filename);
            if(file.is_open()){
                file << "Iteration,Best_Score" << endl;
                for(size_t i = 0; i < curve.size(); ++i)
                    file << i + 1 << "," << curve[i] << endl;
                file.close();
                cout << "Convergence curve saved to '" << filename << "'" << endl;
            }
        }
        
    } else {
        auto funcInfo = TestFunctions<float>::getFunctionInfo(functionName);
        cout << "Function Name: " << funcInfo.name << endl;
        cout << "Dimension: " << funcInfo.dimension << endl;
        
        HippopotamusOptimization<float> optimizer(
            params.searchAgents,
            params.maxIterations,
            funcInfo.lowerBound,
            funcInfo.upperBound,
            funcInfo.function
        );
        
        optimizer.optimize();
        
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        
        optimizer.displayResults();
        cout << "Execution time: " << duration.count() << " ms" << endl;
    }
}

void runCustomFunction() {
    cout << "\n=== Custom Optimization Problem ===" << endl;
    
    OptimizationParams params = getOptimizationParams();
    
    vector<double> lowerBounds(5,0);
    vector<double> upperBounds(5,acos(-1));

    cout << "\n--- Starting Custom Optimization ---" << endl;
    cout << "Problem: Minimize sum((x[i] - i)^2)" << endl;
    cout << "Dimension: 10" << endl;
    cout << "Bounds: [0, PI] for each variable" << endl;
    cout << "Search Agents: " << params.searchAgents << endl;
    cout << "Max Iterations: " << params.maxIterations << endl;
    
    auto start = chrono::high_resolution_clock::now();
    
    HippopotamusOptimization<double> optimizer(
        params.searchAgents,
        params.maxIterations,
        lowerBounds,
        upperBounds,
        customFunction<double>
    );
    
    optimizer.optimize();
    
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    
    optimizer.displayResults();
    cout << "Execution time: " << duration.count() << " ms" << endl;
}

int main() { 
    while(true){
        cout << "\n=== Main Menu ===" << endl;
        cout << "1. Optimize Test Functions (F1-F23)" << endl;
        cout << "2. Run Custom Optimization Problem" << endl;
        cout << "3. Exit" << endl;
        
        int choice = getUserChoice(1, 3);
        
        switch(choice){
            case 1:
                runTestFunction();
                break;
            case 2:
                runCustomFunction();
                break;
            case 3:
                return 0;
        }
        
        cout << "\nPress Enter to continue...";
        cin.get();
    }
}
