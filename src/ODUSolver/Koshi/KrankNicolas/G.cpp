#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> // для использования функции exit()

void umn_matr(double x[5][5], double y[5] ,double z[5], int n, int m) {
    double s;
	for (int i = 0; i < n; ++i) {
        s = 0.0;
        for (int j = 0; j < m; ++j) {
            s = s + x[i][j] * y[j];
        }
        z[i] = s;
	}
}

void obr_matr(double aa[5][5], double bb[5], double dd[5][5], double xx[5] , int n) {
	double Mstu, s, s1;
	double a[5][11];
	int cona[5];
    int  m, m1; 
	int kn, km;
	m = n;
    m1 = m * 2 + 1;
	for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m1; ++j) {
            if (j < m)
                a[i][j] = aa[i][j];
            else
                a[i][j] = 0.0;
            a[i][m+i] = 1.0;
    // cout << a[i][j] << " ";
        }
	//cout << endl;// Строка завершается символом перехода на новую строку
	}

    for (int i = 0; i < n; ++i) {
        a[i][m1-1] =bb[i];
    }

    // Вспомогательный массив
	for (int i = 0; i < n; ++i) {
        cona[i] = m1;
    }
      

    for (int li = 0; li < n; ++li) {
        Mstu= 0.0;
        kn=0;
        km=0;
        for (int i = 0; i < n; ++i) {
            if (cona[i] == m1  ) {
                // максимальный элемент в матрице
                for (int j = 0; j < m; ++j) {
                    if (abs(a[i][j]) > Mstu) {
                        Mstu =abs( a[i][j]); kn = i, km = j;
                    }
                }
            }
        }
        // меняем строки
        for (int l1 = 0; l1 < m1; ++l1) {
            s = a[kn][l1];
            a[kn][l1] = a[km][l1]; 
            a[km][l1] = s;
        }
        cona[km] = km;
        s1 = a[km][km];
        // делим на мах строки
        for (int l2 = 0; l2 < m1; ++l2) {
            a[km][l2] = a[km][l2]/s1;
        }
                
        for (int i = 0; i < n; ++i) {
            if (cona[i] == m1){
                s1 = a[i][km];
                for (int j = 0; j < m1; ++j) {   
                    a[i][j] = a[i][j]- a[km][j]* s1;
                }
            }
            cona[i] = m1;
        }


        for (int i = 0; i < n; ++i) { 
            for (int j = 0; j < m; ++j) {
                dd[i][j] = a[i][j+5];
            }                 
        }

        for (int i = 0; i < n; ++i) {
            xx[i]=a[i][m1-1];
        }
    }
}

int main() {
	using namespace std;// ifstream используется для чтения содержимого файла.
    setlocale(LC_ALL, "rus");
	ifstream inf("file.txt");// Попытаемся прочитать содержимое файла SomeText.txt
    ofstream outf("fileout.txt");
	{
        	
		double a[5][5], bb[5], dd[5][5], x[5], y[5];	
        int n, m ;
		
   inf >> n;
   m=n; 

    for (int i = 0; i < n; ++i)
	{
		 for (int j = 0; j < n; ++j)// Выводим на экран строку i
   {
		inf >> a[i][j];
    cout << a[i][j] << " ";
    }
	cout << endl;// Строка завершается символом перехода на новую строку
	}

for (int i = 0; i < n; ++i)
{
	inf >> bb[i];
}

    inf.close();
   /////////

obr_matr(a, bb,  dd,  x ,  n);

 	for (int i = 0; i < n; ++i)
{   
    for (int j = 0; j < m; ++j)// Выводим на экран строку i
    {
    outf << dd[i][j] << " dd ";
    }
	outf << endl;// Строка завершается символом перехода на новую строку
}
for (int i = 0; i < n; ++i)
{
	outf << x[i] << " x ";
}
 outf << endl;// Строка завершается символом перехода на новую строку

 umn_matr(a, bb, y, n, m);

 for (int i = 0; i < n; ++i)
{
	outf << y[i] << " y ";
}
 outf << endl;// Строка завершается символом перехода на новую строку


 //system("pause");
	}
 
	return 0;// Когда inf выйдет из области видимости, то деструктор класса ifstream автоматически закроет наш файл
 
}
