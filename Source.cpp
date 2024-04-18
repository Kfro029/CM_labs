#include <fstream>
#include <cmath>
#include <vector>

int k = 3; // параметр k в ГУ
float L = 1.0;
int step = 10000; // число шагов в задаче (шаг дискретизации)

float rho(int n) {
	return(sin(k * 2 * 3.14 * (n / step) / L));
}



int main() {
	// fi - неизв, rho = sin(2pikx/L)
	// eps0 = 1

	// Из ГУ fi_1 = 2(pi)kh/L
	// fi_(n-1) = 2(pi)kh*cos(2pik)/L-pow(2(pi)kh/L, 2)sin(2pik)/2

	//Надо найти alpha, beta

	std::vector<float> al;
	std::vector<float> bt;

	std::vector<double> fi(step + 1, 0.0);
	//std::vector<float> fi; //искомый

	//fi.resize(step + 1);
	fi[step] = 0;
	fi[0] = 0;
	fi[1] = (L / step) * (2 * 3.14 * k / L);
	fi[step - 1] = (L / step)*(2 * 3.14 * k / L) * cos(2 * 3.14 * k) - pow(2 * 3.14 * k * (L / step) / L, 2) * sin(2 * 3.14 * k) / 2;

	al.push_back(-0.5);
	bt.push_back(rho(0));


	//прямой ход прогонки
	float y = 2.0;
	for (int i = 2; i < step - 1; i++) {
		y = 2 + al[al.size() - 1];
		al.push_back(-1 / y);
		bt.push_back((rho(i) - bt[bt.size() - 1]) / y);
	}
	y = 2 + al[al.size() - 1];
	bt.push_back((rho(step) - bt[bt.size() - 1]) / y);


	//обратная прогонка
	fi.reserve(step + 1);
	for (int i = step - 2; i > 1; i--) {
		fi[i] = al[i -1] * fi[i + 1] + bt[i - 1];
	}


	std::ofstream fout("3.txt");
	for (int i = 0; i < step + 1; i++) {
		fout << i << " " << fi[i] << std::endl;
	}

	fout.close();
	return 0;
}