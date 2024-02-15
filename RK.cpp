// RunK.cpp: определяет точку входа для приложения.
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

/* Это таблица Бутчера для метода Рунге-Кутты 4 порядка. Я ее заполнил */
struct RK4Table {
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = { { {0., 0., 0., 0.}, {0.5, 0., 0., 0.}, {0., 0.5, 0., 0.}, {0., 0., 1., 0.} } };
    static constexpr std::array<double, stages> cColumn = { 0, 0.5, 0.5, 1 };
    static constexpr std::array<double, stages> bString = { 1. / 6, 1. / 3, 1. / 3, 1. / 6 };
};

class RightPartT3 {


public:

    static constexpr unsigned int dim = 1;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg {
        State state; // y
        Argument arg; // t
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        

        return Eigen::Vector<double, dim>{std::pow(stateAndArg.arg, 3)};

    }
};


//правая часть для y'' = -y ...
class RightPartY {


public:

    static constexpr unsigned int dim = 2;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg {
        State state; // y, x
        Argument arg; // t
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {

        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};

    }
};

// Сигнатура для метода интегрирования:

template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
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
Задание:

Реализовать интегрирование для задач Коши: y' = t^3, y(0) = 0, y'' = - y, y(0) = 0, y'(0) = 1.
Провести интегрирование на отрезке [0, 5] классическим методом Рунге-Кутты четвертого порядка
Построить зависимость ошибки (максимальная раность по всем точкам) от шага интегрирования в логарифмическом масштабе.
Вычислить угол наклона прямой.
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

