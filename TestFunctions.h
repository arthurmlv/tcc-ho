#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <string>
#include <random>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

template<typename TYPE>
class TestFunctions {
public:

    static TYPE F1(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        for (const auto& val : x) {
            sum += val * val;
        }
        return sum;
    }

    static TYPE F2(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        TYPE prod = 1;
        for (const auto& val : x) {
            TYPE absVal = std::abs(val);
            sum += absVal;
            prod *= absVal;
        }
        return sum + prod;
    }

    static TYPE F3(const std::vector<TYPE>& x) {
        TYPE totalSum = 0;
        size_t dim = x.size();
        for (size_t i = 0; i < dim; ++i) {
            TYPE currentSum = 0;
            for (size_t j = 0; j <= i; ++j) {
                currentSum += x[j];
            }
            totalSum += currentSum * currentSum;
        }
        return totalSum;
    }

    static TYPE F4(const std::vector<TYPE>& x) {
        TYPE maxVal = 0;
        for (const auto& val : x) {
            TYPE absVal = std::abs(val);
            if (absVal > maxVal) maxVal = absVal;
        }
        return maxVal;
    }

    static TYPE F5(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        size_t dim = x.size();
        for (size_t i = 0; i < dim - 1; ++i) {
            TYPE term1 = x[i+1] - x[i] * x[i];
            TYPE term2 = x[i] - 1;
            sum += 100 * term1 * term1 + term2 * term2;
        }
        return sum;
    }

    static TYPE F6(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        for (const auto& val : x) {
            TYPE valPlus = std::floor(val + 0.5);
            sum += valPlus * valPlus;
        }
        return sum;
    }

    static TYPE F7(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        size_t dim = x.size();
        for (size_t i = 0; i < dim; ++i) {
            sum += (TYPE)(i + 1) * std::pow(x[i], 4);
        }
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<TYPE> dis(0.0, 1.0);
        
        return sum + dis(gen);
    }

    static TYPE F8(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        for (const auto& val : x) {
            sum += -val * std::sin(std::sqrt(std::abs(val)));
        }
        return sum;
    }

    static TYPE F9(const std::vector<TYPE>& x) {
        TYPE sum = 0;
        size_t dim = x.size();
        for (const auto& val : x) {
            sum += val * val - 10 * std::cos(2 * M_PI * val);
        }
        return sum + 10 * (TYPE)dim;
    }

    static TYPE F10(const std::vector<TYPE>& x) {
        TYPE sumSq = 0;
        TYPE sumCos = 0;
        size_t dim = x.size();
        
        for (const auto& val : x) {
            sumSq += val * val;
            sumCos += std::cos(2 * M_PI * val);
        }

        return -20 * std::exp(-0.2 * std::sqrt(sumSq / dim)) 
               - std::exp(sumCos / dim) 
               + 20 + M_E;
    }

    // F11
    static TYPE F11(const std::vector<TYPE>& x) {
        TYPE sumSq = 0;
        TYPE prodCos = 1;
        size_t dim = x.size();

        for (size_t i = 0; i < dim; ++i) {
            sumSq += x[i] * x[i];
            prodCos *= std::cos(x[i] / std::sqrt((TYPE)(i + 1)));
        }
    
        return sumSq / 4000.0 - prodCos + 1;
    }

    static TYPE Ufun(const std::vector<TYPE>& x, TYPE a, TYPE k, TYPE m) {
        TYPE result = 0;
        for (const auto& val : x) {
            if(val>a)
                result += k * std::pow(val - a, m);
            else if(val<-a)
                result += k * std::pow(-val - a, m);
        }
        return result;
    }

    static TYPE F12(const std::vector<TYPE>& x) {
        size_t dim = x.size();
        
        std::vector<TYPE> y(dim);
        for(size_t i=0; i<dim; ++i)
            y[i] = 1.0 + (x[i] + 1.0) / 4.0;
        

        TYPE term1 = 10 * std::pow(std::sin(M_PI * y[0]), 2);
        
        TYPE term2 = 0;
        for(size_t i=0; i<dim-1; ++i) {
            TYPE dy = y[i] - 1.0;
            term2 += (dy * dy) * (1.0 + 10.0 * std::pow(std::sin(M_PI * y[i+1]), 2));
        }

        TYPE term3 = std::pow(y[dim-1] - 1.0, 2);
        TYPE uSum = Ufun(x, 10, 100, 4);

        return (M_PI / dim) * (term1 + term2 + term3) + uSum;
    }

    static TYPE F13(const std::vector<TYPE>& x) {
        size_t dim = x.size();

        TYPE term1 = std::pow(std::sin(3 * M_PI * x[0]), 2);

        TYPE term2 = 0;
        for(size_t i=0; i < dim - 1; ++i) {
            TYPE dx = x[i] - 1.0;
            term2 += (dx * dx) * (1.0 + std::pow(std::sin(3 * M_PI * x[i+1]), 2));
        }

        TYPE dxEnd = x[dim-1] - 1.0;
        TYPE term3 = (dxEnd * dxEnd) * (1.0 + std::pow(std::sin(2 * M_PI * x[dim-1]), 2));

        TYPE uSum = Ufun(x, 5, 100, 4);

        return 0.1 * (term1 + term2 + term3) + uSum;
    }

    static TYPE F14(const std::vector<TYPE>& x) {
        static const TYPE aS[2][25] = {
            {-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32},
            {-32,-32,-32,-32,-32,-16,-16,-16,-16,-16,0,0,0,0,0,16,16,16,16,16,32,32,32,32,32}
        };

        TYPE sumOuter = 0;
        for (int j=0; j<25; ++j) {
            TYPE sumInner = 0;
            for (int i=0; i<2; ++i) 
                sumInner += std::pow(x[i] - aS[i][j], 6);
            
            sumOuter += 1.0 / ((double)(j + 1) + sumInner);
        }
        
        return std::pow(1.0/500.0 + sumOuter, -1);
    }

    // F15: Kowalik
    static TYPE F15(const std::vector<TYPE>& x) {
        static const TYPE aK[11] = {.1957, .1947, .1735, .16, .0844, .0627, .0456, .0342, .0323, .0235, .0246};
        static const TYPE bK_orig[11] = {.25, .5, 1, 2, 4, 6, 8, 10, 12, 14, 16};
        
        TYPE sum = 0;
        for(int i = 0; i < 11; ++i) {
            TYPE invBk = 1.0 / bK_orig[i];
            TYPE invBk2 = invBk * invBk;
            TYPE numerator = x[0] * (invBk2 + x[1] * invBk);
            TYPE denominator = invBk2 + x[2] * invBk + x[3];
            
            TYPE term = aK[i] - (numerator / denominator);
            sum += term * term;
        }
        return sum;
    }

    static TYPE F16(const std::vector<TYPE>& x) {
        TYPE x1 = x[0];
        TYPE x2 = x[1];
        return 4*std::pow(x1,2) - 2.1*std::pow(x1,4) + std::pow(x1,6)/3.0 + x1*x2 - 4*std::pow(x2,2) + 4*std::pow(x2,4);
    }

    static TYPE F17(const std::vector<TYPE>& x) {
        TYPE x1 = x[0];
        TYPE x2 = x[1];
        TYPE term1 = x2 - (std::pow(x1,2) * 5.1 / (4*std::pow(M_PI,2))) + (5/M_PI)*x1 - 6;
        TYPE term2 = 10 * (1 - 1/(8*M_PI)) * std::cos(x1) + 10;
        
        return std::pow(term1, 2) + term2;
    }

    static TYPE F18(const std::vector<TYPE>& x) {
        TYPE x1 = x[0];
        TYPE x2 = x[1];
        
        TYPE part1 = 1 + std::pow(x1 + x2 + 1, 2) * (19 - 14*x1 + 3*std::pow(x1,2) - 14*x2 + 6*x1*x2 + 3*std::pow(x2,2));
        TYPE part2 = 30 + std::pow(2*x1 - 3*x2, 2) * (18 - 32*x1 + 12*std::pow(x1,2) + 48*x2 - 36*x1*x2 + 27*std::pow(x2,2));
        
        return part1 * part2;
    }

    static TYPE F19(const std::vector<TYPE>& x) {
        static const TYPE aH[4][3] = {
            {3, 10, 30}, {.1, 10, 35}, {3, 10, 30}, {.1, 10, 35}
        };
        static const TYPE cH[4] = {1, 1.2, 3, 3.2};
        static const TYPE pH[4][3] = {
            {.3689, .117, .2673}, {.4699, .4387, .747}, {.1091, .8732, .5547}, {.03815, .5743, .8828}
        };

        TYPE R = 0;
        for(int i = 0; i < 4; ++i) {
            TYPE sumInner = 0;
            for(int j = 0; j < 3; ++j) {
                sumInner += aH[i][j] * std::pow(x[j] - pH[i][j], 2);
            }
            R -= cH[i] * std::exp(-sumInner);
        }
        return R;
    }

    static TYPE F20(const std::vector<TYPE>& x) {
        static const TYPE aH[4][6] = {
            {10, 3, 17, 3.5, 1.7, 8},
            {.05, 10, 17, .1, 8, 14},
            {3, 3.5, 1.7, 10, 17, 8},
            {17, 8, .05, 10, .1, 14}
        };
        static const TYPE cH[4] = {1, 1.2, 3, 3.2};
        static const TYPE pH[4][6] = {
            {.1312, .1696, .5569, .0124, .8283, .5886},
            {.2329, .4135, .8307, .3736, .1004, .9991},
            {.2348, .1415, .3522, .2883, .3047, .6650},
            {.4047, .8828, .8732, .5743, .1091, .0381}
        };

        TYPE R = 0;
        for(int i = 0; i < 4; ++i) {
            TYPE sumInner = 0;
            for(int j = 0; j < 6; ++j) {
                sumInner += aH[i][j] * std::pow(x[j] - pH[i][j], 2);
            }
            R -= cH[i] * std::exp(-sumInner);
        }
        return R;
    }

    static TYPE Shekel(const std::vector<TYPE>& x, int m) {
        static const TYPE aSH[10][4] = {
            {4, 4, 4, 4}, {1, 1, 1, 1}, {8, 8, 8, 8}, {6, 6, 6, 6},
            {3, 7, 3, 7}, {2, 9, 2, 9}, {5, 5, 3, 3}, {8, 1, 8, 1},
            {6, 2, 6, 2}, {7, 3.6, 7, 3.6}
        };
        static const TYPE cSH[10] = {.1, .2, .2, .4, .4, .6, .3, .7, .5, .5};

        TYPE R = 0;
        for(int i = 0; i < m; ++i) {
            TYPE sumInner = 0;
            for(int j = 0; j < 4; ++j) {
                sumInner += std::pow(x[j] - aSH[i][j], 2);
            }
            R -= 1.0 / (sumInner + cSH[i]);
        }
        return R;
    }

    static TYPE F21(const std::vector<TYPE>& x) { return Shekel(x, 5); }

    static TYPE F22(const std::vector<TYPE>& x) { return Shekel(x, 7); }

    static TYPE F23(const std::vector<TYPE>& x) { return Shekel(x, 10); }


    struct FunctionInfo {
        std::vector<TYPE> lowerBound;
        std::vector<TYPE> upperBound;
        int dimension;
        std::function<TYPE(const std::vector<TYPE>&)> function;
        std::string name;
    };
    
    static FunctionInfo getFunctionInfo(const std::string& functionName) {
        FunctionInfo info;
        
        auto setBounds = [&](int dim, TYPE lb, TYPE ub) {
            info.dimension = dim;
            info.lowerBound = std::vector<TYPE>(dim, lb);
            info.upperBound = std::vector<TYPE>(dim, ub);
        };

        if(functionName == "F1"){
            setBounds(30, -100, 100); info.function = F1; info.name = "Sphere Function";
        }
        else if(functionName == "F2"){
            setBounds(30, -10, 10); info.function = F2; info.name = "Schwefel's Problem 2.21";
        }
        else if(functionName == "F3"){
            setBounds(30, -100, 100); info.function = F3; info.name = "Schwefel's Problem 2.22";
        }
        else if(functionName == "F4"){
            setBounds(30, -100, 100); info.function = F4; info.name = "Schwefel's Problem 2.21 (Max)";
        }
        else if(functionName == "F5"){
            setBounds(30, -30, 30); info.function = F5; info.name = "Rosenbrock Function";
        }
        else if(functionName == "F6"){
            setBounds(30, -100, 100); info.function = F6; info.name = "Step Function";
        }
        else if(functionName == "F7"){
            setBounds(30, -1.28, 1.28); info.function = F7; info.name = "Quartic Function with Noise";
        }
        else if(functionName == "F8"){
            setBounds(30, -500, 500); info.function = F8; info.name = "Schwefel Function";
        }
        else if(functionName == "F9"){
            setBounds(30, -5.12, 5.12); info.function = F9; info.name = "Rastrigin Function";
        }
        else if(functionName == "F10"){
            setBounds(30, -32, 32); info.function = F10; info.name = "Ackley Function";
        }
        else if(functionName == "F11"){
            setBounds(30, -600, 600); info.function = F11; info.name = "Griewank Function";
        }
        else if(functionName == "F12"){
            setBounds(30, -50, 50); info.function = F12; info.name = "Penalized Function 1";
        }
        else if(functionName == "F13"){
            setBounds(30, -50, 50); info.function = F13; info.name = "Penalized Function 2";
        }
        else if(functionName == "F14"){
            setBounds(2, -65.536, 65.536); info.function = F14; info.name = "Foxholes Function";
        }
        else if(functionName == "F15"){
            setBounds(4, -5, 5); info.function = F15; info.name = "Kowalik Function";
        }
        else if(functionName == "F16"){
            setBounds(2, -5, 5); info.function = F16; info.name = "Six-Hump Camel Function";
        }
        else if(functionName == "F17"){
            info.dimension = 2;
            info.lowerBound = {(TYPE)-5, (TYPE)0};
            info.upperBound = {(TYPE)10, (TYPE)15};
            info.function = F17; info.name = "Branin Function";
        }
        else if(functionName == "F18"){
            setBounds(2, -2, 2); info.function = F18; info.name = "Goldstein-Price Function";
        }
        else if(functionName == "F19"){
            setBounds(3, 0, 1); info.function = F19; info.name = "Hartman 3 Function";
        }
        else if(functionName == "F20"){
            setBounds(6, 0, 1); info.function = F20; info.name = "Hartman 6 Function";
        }
        else if(functionName == "F21"){
            setBounds(4, 0, 10); info.function = F21; info.name = "Shekel 5 Function";
        }
        else if(functionName == "F22"){
            setBounds(4, 0, 10); info.function = F22; info.name = "Shekel 7 Function";
        }
        else if(functionName == "F23"){
            setBounds(4, 0, 10); info.function = F23; info.name = "Shekel 10 Function";
        }
        else {
            // Default
            setBounds(30, -100, 100); info.function = F1; info.name = "Sphere Function (Default)";
        }
        
        return info;
    }
};

#endif