#include <array>
#include <iostream>
#include <cmath>
#include <fstream>


template<typename xType, typename yType, unsigned int N>



class NewtonInterpolator {
    std::array<xType, N> X;
    std::array<yType, N> Y;
public:
    NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values) noexcept {

        X = points;
        Y = values;


        // вычисляем разделенные разности
        for (int i = 0; i < N-1; i++) {
            for (int j = N - 1; j > i; j--) {
                Y[j] = (Y[j] - Y[j - 1]) / (X[j] - X[j - 1 - i]);
            }

        }

    }
    //Вычисляем полином Ньютона по схеме Горнера
    yType interpolate(const xType& x) const noexcept {
        yType sum = 0;
        for (int i = N-1; i > 0; i--) {
            sum += Y[i];
            sum *= (x - X[i - 1]);
        }
        sum += Y[0];
        return sum;
    }
    
   



};

int main() {
   
    const unsigned int N = 5; // число узлов
    std::ofstream fout("5points.txt");
    

    for (int k = 1; k > (-5); k--) {

        std::array<double, N> X;
        std::array<double, N> Y;

        // заполнение массива координат x, y
        for (int i = 0; i < N; i++) {
            X[i] = (pow(2,k) * i) / (N - 1);
            Y[i] = exp(X[i]);
        }

        //интерполянт на равномерном распределении узлов
        NewtonInterpolator<double, double, N> M(X, Y);
        

        //макс. ошбка для равномерного
        double fail_evenly = 0;
        double qurrent = 0;


        //ошибка как максимальная разница в 1000 точках, равномерно распред. для "равномерного" интерполянта
        for (int i = 0; i < 1000; i++) {
            qurrent = abs(exp(i / 999. * pow(2, k)) - M.interpolate(i / 999. * pow(2, k)));
            if (fail_evenly < qurrent) fail_evenly = qurrent;
        }
        


        // Для нулей полинома Чебышева
         
         
        std::array<double, N> Xch;
        std::array<double, N> Ych;
         
         
        // заполнение массива координат x, y
        for (int i = 0; i < N; i++) {
            Xch[i] = pow(2, k - 1) + pow(2, k - 1) * cos((2 * i + 1) * 3.14159265358979323846 / (2 * N));
            Ych[i] = exp(Xch[i]);
        }
        




        NewtonInterpolator<double, double, N> Mch(Xch, Ych);

        //макс. ошбка для Чебышевского распределения
        double fail_cheb = 0;
        double qurrent_cheb = 0;



        //ошибка как максимальная разница в 1000 точках, равномерно распред. для "Чебышевского" интерполянта
        for (int i = 0; i < 1000; i++) {
            qurrent_cheb = abs(exp(i / 999. * pow(2, k)) - Mch.interpolate(i / 999. * pow(2, k)));
            if (fail_cheb < qurrent_cheb) fail_cheb = qurrent_cheb;
        }
        
        fout << fail_evenly << " " <<  fail_cheb << " " << std::endl; 
        

    }
    fout.close();
    
    return 0;
}
