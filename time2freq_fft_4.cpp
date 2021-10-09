#include <complex>
#include <iostream>
#include <valarray>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <ctime>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;

int n = 4140;
const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

vector<float> x;
vector<float> y;
vector<float> z;

float freq[4140];

Complex x_fft[4140];
Complex y_fft[4140];
Complex z_fft[4140];


void fft(CArray& x) // Cooley–Tukey FFT 
{
	const size_t N = x.size();
	if (N <= 1) return;

	CArray even = x[std::slice(0, N/2, 2)];
	CArray  odd = x[std::slice(1, N/2, 2)];

	fft(even);
	fft(odd);

	for (size_t k = 0; k < N/2; ++k)
	{
		Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
		x[k    ] = even[k] + t;
		x[k+N/2] = even[k] - t;
	}
}


void ifft(CArray& x) 
{
	x = x.apply(std::conj); 
	fft( x ); 
	x = x.apply(std::conj); 
	x /= x.size(); 
}

void read_file(){
	ifstream file;
	file.open("1602245833-2715-NAO7856.txt");
		
	char line[256];
	char size_data[256];

	vector<float> data_file;

	for(int i=0; i<n; i++){
		char read_data[256];
		file.getline (read_data,256);
			
    	int tamanho = strlen(read_data); 
    	char *token = strtok(read_data, ",");
    	int ind = 0;
		while(token != NULL) {
        	float float_data = atof(token);
			if(ind == 0) x.push_back(float_data);
			else if(ind == 1) y.push_back(float_data);
			else{
				z.push_back(float_data);
			}
		
        	token = strtok(NULL, ",");
			ind++;
		}
	}
	
}

void save_file(){

	ofstream outFile; 
	
	outFile.open("output.txt", ios::out); 
	if (! outFile){ 
		cout << "Arquivo output.txt nao pode ser aberto" << endl;
 		abort();
	} 
	
	for(int e=0; e<n; e++){
		//cout << real(x_fft[e]) << "," << real(y_fft[e]) << "," << real(z_fft[e]) << "," << freq[e] << endl;
		outFile << real(x_fft[e]) << "," << real(y_fft[e]) << "," << real(z_fft[e]) << "," << freq[e] << endl;
	}

}

void fft_data(vector<float> function_value,char eixo){
	
	Complex L = n;
	Complex vector_data[n];
	
	for(int p = 0; p<n; p++){
		vector_data[p] = function_value[p];
	}
	
	CArray data(vector_data, n);

	fft(data);

	if(eixo == 'x'){
		for (int i = 0; i < n; ++i)
		{

			x_fft[i] = (data[i]*(conj(data[i])))/L;		
			freq[i] = (1/(0.000655*n))*i; // 0.000655 valor de 2.715 segundo dividido por 4140 amostras 
			if(i == 0) x_fft[i] = 0;
		}
	}
	else if(eixo == 'y'){
		for (int i = 0; i < n; ++i)
		{
			y_fft[i] = (data[i]*(conj(data[i])))/L;
			freq[i] = (1/(0.000655*n))*i; // 0.000655 valor de 2.715 segundo dividido por 4140 amostras 
			if(i == 0) y_fft[i] = 0;
		}
	}
	else{
		for (int i = 0; i < n; ++i)
		{
			z_fft[i] = (data[i]*(conj(data[i])))/L;
			freq[i] = (1/(0.000655*n))*i; // 0.000655 valor de 2.715 segundo dividido por 4140 amostras 
			if(i == 0) z_fft[i] = 0;
		}
	}

}


int main()
{
	
	// O arquivo de saída output.txt esta com uma amplitude muito alta em zero, desta forma os valores de 0 Hz foram forçados a 0.
	// O código para o calculo da fft foi retirado e basado no disponibilizado neste link: https://tfetimes.com/c-fast-fourier-transform/
	// Para o caculo e definição dos vetores de espectro e frequecia utilizou-se o databook. http://databookuw.com/databook.pdf

	read_file(); 
	fft_data(x,'x');
	fft_data(y,'y');
	fft_data(z,'z');

	save_file();

}