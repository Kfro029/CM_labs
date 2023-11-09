#include "Code.h"
#include <array>
#include <cmath>
#include <fstream>

	template<typename RealType, unsigned int N>
struct DerivativeCoef {

	RealType centralCoef;
	std::array<RealType, N> otherCoefs;
};

	template<typename RealType, unsigned int N>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
	//тут N-число нецентральных точек
	std::array<std::array<RealType, N + 2>, N+1> A;

	//Уравнение вида Aa=x

	//Заполняем 1ую строку и 1ый и (N+2)ой столбец матрицы A|x
	for (unsigned int i = 0; i < N+1; i++) {
		A[i][N+1] = 0;
		A[i][0] = 0;
	};
	A[1][N+1] = 1;
	A[0][0] = 1;
	

	//Заполняем остальные столбцы A
	for (unsigned int i = 0; i < N+1; i++) {
		for (unsigned int j = 1; j < N+1; j++) {
			A[i][j] = pow(points[j-1], i);
			
			
		};
		
	};


	
	
	//Метод Гаусса с выбором ведущего элемента
	std::array<RealType, N + 2> help;
	for (unsigned int i = 1; i < N+1; i++) {

		//Выбор максимального элемента
		int ind_max = i;
		for (unsigned int m = i + 1; m < N+1; m++) {
			if (std::abs(A[m][i]) > std::abs(A[ind_max][i])) {
				ind_max = m;
			}
				
		}
		//Переводим строку с макс. элеметом на i-ое место
		if (ind_max != i) {
			help = A[i];
			A[i] = A[ind_max];
			A[ind_max] = help;
		}
		



		//делим на A[i][i] строку
		RealType del = A[i][i];
		for (unsigned int j = N+1; j > i-1; j--) {
			A[i][j] = A[i][j] / A[i][i];
		}

		//Вычетаем из "верхних" строк i-ую
		for (unsigned int j = 0; j < i; j++) {
			RealType rule = A[j][i];
			for (unsigned int k = N + 1; k > i-1; k--) {
				A[j][k] = A[j][k] - A[i][k] * A[j][i];
			}
		}



		//Вычетаем из "нижних" строк i-ую
		for (unsigned int j = i + 1; j < N+1; j++) {
			RealType rule = A[j][i];
			for (unsigned int k = N + 1; k > i - 1; k--) {
				A[j][k] = A[j][k] - A[i][k] * A[j][i];
			}
		}
		
		for (unsigned int l = 0; l < N + 1; l++) {
		for (unsigned int y = 0; y < N + 2; y++) {
			std::cout << A[l][y] << " ";
		}
		std::cout << std::endl;
	};
	std::cout << std::endl;

	};
	
	

	

	DerivativeCoef<RealType, N> M;
	M.centralCoef = A[0][N];
	for (unsigned int i = 0; i < N; i++) {
		M.otherCoefs[i] = A[i+1][N+1];
	}
	return M;
	};


int main()
{
	double x_0 = 1.; //исходная точка

	

	//для N узлов
	const unsigned int N = 3; // всего точек

	std::array<double, N-1> points3 = {-1., 1.}; // коэффициенты при h в доп. точках

	std::ofstream fout("3points_new.txt");



	//Вывод на экран коэффициентов
	DerivativeCoef<double, N-1> A3 = calcDerivativeCoef(points3);
	std::cout << "Central coef " << A3.centralCoef << std::endl;
	std::cout << "Other coefs" << " ";
	for (unsigned int i = 0; i < N - 1; i++) {
		std::cout << A3.otherCoefs[i] << " ";
	}



	//Запись в файл 
	for (int i = 0; i > -16; i--) {
		//sum - это D(f)*h (в обозначениях с семинара) 
		double sum = A3.centralCoef* exp(x_0);
		for (int j = 0; j < N-1; j++) {
			sum += A3.otherCoefs[j] * exp(x_0 + points3[j] * pow(10., i));
		}
		sum = sum / pow(10, i); // теперь sum это D(f)
		double fail = abs(sum - exp(x_0));
		fout << pow(10., i) << " " << fail << std::endl;

	}


	fout.close();


	

	return 0;
}
