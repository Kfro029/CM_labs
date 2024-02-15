// RunK.cpp: îïðåäåëÿåò òî÷êó âõîäà äëÿ ïðèëîæåíèÿ.
//


#include<iostream>
#include <Eigen/Dense>
#include<cmath>
#include <array>
#include <vector>
#include <fstream>

double T3(double t, double y) {
    return pow(t, 3);
}

/* Ýòî òàáëèöà Áóò÷åðà äëÿ ìåòîäà Ðóíãå-Êóòòû 4 ïîðÿäêà. ß åå çàïîëíèë */
struct RK4Table {
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = { { {0., 0., 0., 0.}, {0.5, 0., 0., 0.}, {0., 0.5, 0., 0.}, {0., 0., 1., 0.} } };
    static constexpr std::array<double, stages> cColumn = { 0, 0.5, 0.5, 1 };
    static constexpr std::array<double, stages> bString = { 1. / 6, 1. / 3, 1. / 3, 1. / 6 };
};

class RightPartT3 {


public:

    static constexpr unsigned int dim = 1;  // ðàçìåðíîñòü çàäà÷è

    using Argument = double;  // òèï àðãóìåíòà, òèï t

    using State = Eigen::Vector<double, dim>;  // ñîñòîÿíèå

    struct StateAndArg {
        State state; // y
        Argument arg; // t
    };

    /*** Âû÷èñëÿåò ïðàâóþ ÷àñòü ÄÓ - ôóíêöèþ f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        

        return Eigen::Vector<double, dim>{std::pow(stateAndArg.arg, 3)};

    }
};


//ïðàâàÿ ÷àñòü äëÿ y'' = -y ...
class RightPartY {


public:

    static constexpr unsigned int dim = 2;  // ðàçìåðíîñòü çàäà÷è

    using Argument = double;  // òèï àðãóìåíòà, òèï t

    using State = Eigen::Vector<double, dim>;  // ñîñòîÿíèå

    struct StateAndArg {
        State state; // y, x
        Argument arg; // t
    };

    /*** Âû÷èñëÿåò ïðàâóþ ÷àñòü ÄÓ - ôóíêöèþ f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {

        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};

    }
};

// Ñèãíàòóðà äëÿ ìåòîäà èíòåãðèðîâàíèÿ:

template<typename Table, typename RHS>  // òàáëèöà áóò÷åðà è êëàññ ïðàâîé ÷àñòè f
std::vector<typename RHS::StateAndArg> integrate(const typename RHS::StateAndArg& initialState,
    const typename RHS::Argument& endTime,
    double step,
    const RHS& rhs) {

    unsigned int num = std::ceil((endTime - initialState.arg) / step);

    step = (endTime - initialState.arg) / num;


    Eigen::Matrix<double, RHS::dim, Table::stages> K;
    std::vector<typename RHS::StateAndArg> result;

    result.push_back(initialState);


    Eigen::Vector<double, RHS::dim> y_fix;
    Eigen::Vector<double, RHS::dim> Ksum;


    for (std::size_t i = 1; i < num + 1; ++i) {

        for (std::size_t j = 0; j < Table::stages; ++j) {

            y_fix = Eigen::Vector<double, RHS::dim>::Zero();

            for (std::size_t n = 0; n < Table::stages; ++n) {

                y_fix += step * Table::table[j][n] * K.col(n);
            }
            K.col(j) = rhs.calc({ result.back().state + y_fix, result.back().arg + step * Table::cColumn[j] });
        }
        Ksum = Eigen::Vector<double, RHS::dim>::Zero();
        for (std::size_t j = 0; j < Table::stages; ++j) {

            Ksum += K.col(j) * Table::bString[j];
        }
        result.push_back({ result.back().state + step * Ksum, result.back().arg + step });
    }
    return result;

};

/*
Çàäàíèå:

Ðåàëèçîâàòü èíòåãðèðîâàíèå äëÿ çàäà÷ Êîøè: y' = t^3, y(0) = 0, y'' = - y, y(0) = 0, y'(0) = 1.
Ïðîâåñòè èíòåãðèðîâàíèå íà îòðåçêå [0, 5] êëàññè÷åñêèì ìåòîäîì Ðóíãå-Êóòòû ÷åòâåðòîãî ïîðÿäêà
Ïîñòðîèòü çàâèñèìîñòü îøèáêè (ìàêñèìàëüíàÿ ðàíîñòü ïî âñåì òî÷êàì) îò øàãà èíòåãðèðîâàíèÿ â ëîãàðèôìè÷åñêîì ìàñøòàáå.
Âû÷èñëèòü óãîë íàêëîíà ïðÿìîé.
*/



int main()
{
    std::ofstream fout("T3.txt");



    RightPartT3 M;
    RightPartT3::StateAndArg init;
    init.arg = 0.;
    init.state(0) = 0.;

    

    for (double i = 2; i > 0.01; i -= 0.01) {
        auto numbers = integrate<RK4Table, RightPartT3>(init, 5., i, M);
        double max = 0;
        for (int j = 0; j < size(numbers); j++) {
            if (abs(numbers[j].state(0) - pow(numbers[j].arg, 4) / 4) > max) {
                max = abs(numbers[j].state(0) - pow(numbers[j].arg, 4) / 4);
            }
        }
        fout << i << " " << max << std::endl;
    }

    fout.close();


    /*std::ofstream fout("sin.txt");

    RightPartY Z;
    RightPartY::StateAndArg init1;
    init1.arg = 0.;
    init1.state = { 0., 1. };


    for (double i = 2; i > 0.01; i -=0.01) {
        auto numbers1 = integrate<RK4Table, RightPartY>(init1, 5., i, Z);
        double max1 = 0;
        for (int j = 0; j < size(numbers1); j++) {
            if (abs(numbers1[j].state(0) - sin(numbers1[j].arg)) > max1){
                max1 = abs(numbers1[j].state(0) - sin(numbers1[j].arg));
            };
        };
        fout << i << " " << max1 << std::endl;
    }
    
    fout.close();*/




    return 0;
}

