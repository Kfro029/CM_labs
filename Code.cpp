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
	std::array<std::array<RealType, N + 1>, N> A;

	//Уравнение вида Aa=x

	//Заполняем 1ую строку и 1ый и (N+1)ый столбец матрицы A|x
	for (unsigned int i = 0; i < N; i++) {
		A[i][N] = 0.;
		A[i][0] = 0.;
	};
	A[1][N] = 1.;
	A[0][0] = 1.;
	

	//Заполняем остальные столбцы A
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 1; j < N; j++) {
			A[i][j] = pow(points[j-1], i);
			
			
		};
		
	};
	
	//Метод Гаусса с выбором ведущего элемента
	std::array<RealType, N + 1> help;
	for (unsigned int i = 1; i < N; i++) {

		//Выбор максимального элемента
		int ind_max = i;
		for (unsigned int m = i + 1; m < N; m++) {
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
		for (unsigned int j = 0; j < N + 1 - i; j++) {
			A[i][i+j] = A[i][i+j] / del;
		}

		
		for (unsigned int j = 0; j < i; j++) {
			RealType rule = A[j][i];
			for (unsigned int k = 0; k < N + 1 - i; k++) {
				A[j][i + k] = A[j][i + k] - A[i][i + k] * rule;
			}
		}



		//Вычетаем из "нижних" строк i-ую
		for (unsigned int j = i + 1; j < N; j++) {
			RealType rule = A[j][i];
			for (unsigned int k = 0; k < N + 1 - i; k++) {
				A[j][i + k] = A[j][i + k] - A[i][i + k] * rule;
			}
		}
		
	};
	

	

	DerivativeCoef<RealType, N> M;
	M.centralCoef = A[0][N];
	for (unsigned int i = 0; i < N-1; i++) {
		M.otherCoefs[i] = A[i+1][N];
	}
	return M;
	};


int main()
{
	double x_0 = 1.; //исходная точка

	

	//для N узлов
	const unsigned int N = 6; // всего точек
	std::array<double, N> points3 = {-2., -1., 1., 2., 3. }; // коэффициенты при h в доп. точках

	std::ofstream fout("6points.txt");

	DerivativeCoef<double, N> A3 = calcDerivativeCoef(points3);


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
