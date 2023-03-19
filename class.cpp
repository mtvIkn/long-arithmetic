#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <cmath>

#define size 1000
#define 
using namespace std;

class math{
private:
    int *mas;
public:
    math() {
        mas = (int*)malloc(size*sizeof(int));
        memset(mas, 0, size*sizeof(mas));       
    }
    
    math(int a) {
        mas = (int*)malloc(size * sizeof(int));
        int count = 0;

        for (int i = 0; a; i++) {
            mas[i] = a % 10;
            a /= 10;
            count++;
        }
        mas[size - 1] = count; //в последнем элементе массива указал длину массива
    }                                                          
                             

    ~math() {
        delete [] mas;
        mas = NULL;
    }

    string show() {
        int p;
        for (p = size - 2; p > 0 && !mas[p]; p--);
        string r = "";
        for (int i = 0; i <= p; i++) {
            r += '0' + mas[p - i];
        }
        cout << r;
        return r;
    }

    const math operator - (const math& ptr){
        int length;
        int k=3; //дополнительная переменная для удобства, если цифра 3- числа равны
        int sum[size];//вспомогательный массив

        if (mas[size - 1] > ptr.mas[size - 1]) {
            length = mas[size - 1];
            k = 1;
        }
        else {
            if (ptr.mas[size - 1] > mas[size - 1]) {
                length = ptr.mas[size - 1];
                k = 2;
            }
            else {
                length = mas[size - 1];//тут числа равны
                for (int i = length - 1; i != 0; i--) {
                    if (mas[i] > ptr.mas[i]) {
                        k = 1;
                        break;
                    }
                    if (ptr.mas[i] > mas[i]) {
                        k = 2;
                        break;
                    }
                }

            }
        }

        if (k == 1) {
            for (int i = 0; i < length - 1; i++){
                if (i < (length - 1)){
                    mas[i + 1]--; 
                    sum[i] += 10 + mas[i]; 
                }
                else{
                    sum[i] += mas[i];
                }
                sum[i] -= ptr.mas[i]; 
                if (sum[i] / 10 > 0){
                    sum[i + 1]++; 
                    sum[i] %= 10; 
                }
            }
            for (int i = 0; i < length; i++) {
                mas[i] = sum[i];
            }
            return *mas;
        }
        else if (k == 2) {
            for (int i = 0; i < length - 1; i++) {
                if (i < (length - 1)) {
                    ptr.mas[i + 1]--;
                    sum[i] += 10 + ptr.mas[i];
                }
                else {
                    sum[i] += ptr.mas[i];
                }
                sum[i] -= mas[i];
                if (sum[i] / 10 > 0) {
                    sum[i + 1]++;
                    sum[i] %= 10;
                }
            }
            for (int i = 0; i < length; i++) {
                mas[i] = sum[i];
            }
            return *mas;
        }
        else {
            for (int i = 0; i < length; i++) {
                mas[i] = 0;
            }
            return *mas;
        }
    }

    const math operator*(const math& ptr) {
        int sum[size];//вспомогательный массив
        int length;

        length = mas[size - 1] + ptr.mas[size - 1] + 1;

        for (int i = 0; i < mas[size - 1]; i++)
            for (int j = 0; j < ptr.mas[size - 1]; j++)
                sum[i + j - 1] += mas[i] * ptr.mas[j];

        for (int i = 0; i < length; i++)
        {
            sum[i + 1] += sum[i] / 10;
            sum[i] %= 10;
        }

        while (sum[length] == 0)
            length--;

        for (int i = 0; i < length; i++) {
            mas[i] = sum[i];
        }
        return *mas;
    }
};

class fourie {
private:
    typedef complex <double> dft;  //для умножения фурье необходимы комплексные чилса
    dft *mas_fourie;
    dft *grade; //массив для предподсчёта степеней числа (на просторах интернета узнал, что такой метод ускоряет алгоритм)
    int count;
    int rev[size];
    bool invert;
    int count;
    dft res[size]; // результирующий массив
public:
    fourie() {
        mas_fourie = (dft*)malloc(size * sizeof(dft));
        memset(mas_fourie, 0, size * sizeof(mas_fourie));
    }

    fourie(dft a) {
        mas_fourie = (dft*)malloc(size * sizeof(dft));
        int count = 0;

        for (int i = 0; a; i++) {
            mas_fourie[i] = a % 10;
            a /= 10;
            count++;
        }
    }

    void fft(const fourie& F, bool invert) {
        revers();// вызываем функцию реверса 
        for (int i = 0; i < count; ++i)
            if (i < rev[i])
                swap(mas_fourie[i], mas_fourie[rev[i]]);

        for (int len = 2; len <= count; len <<= 1) {
            double ang = 2 * 3.14159 / len * (invert ? -1 : +1);
            int len2 = len >> 1;

            dft wlen(cos(ang), sin(ang));
            grade[0] = dft(1, 0);
            for (int i = 1; i < len2; ++i)
                grade[i] = grade[i - 1] * wlen;

            for (int i = 0; i < count; i += len) {
                dft t;
                dft *pu = mas_fourie + i;          //указатели на текущие элементы массивов
                dft *pv = mas_fourie + i + len2;     
                dft *pu_end = mas_fourie + i + len2;
                dft *pw = grade;
                for (; pu != pu_end; ++pu, ++pv, ++pw) {
                    t = *pv * *pw;
                    *pv = *pu - t;
                    *pu += t;
                }
            }
        }

        if (invert) {
            for (int i = 0; i < count; ++i) {
                mas_fourie[i] /= count;
            }              
        }           
    }

    void revers() {                         //предподсчет реверса битов(также на просторах интернета узнал, что такой метод ускоряет алгоритм)
        for (int i = 0; i < count; ++i) {
            rev[i] = 0;
            for (int j = 0; j < log(count); ++j)
                if (i & (1 << j))
                    rev[i] |= 1 << (log(count) - 1 - j);
        }
    }

    void multiply(const fourie& F) {        
        fft(*mas_fourie, false), fft(F.mas_fourie, false);
        for (int i = 0; i < count; ++i)
            *mas_fourie[i] *= F.mas_fourie[i];
        fft(*mas_fourie, true);

        for (int i = 0; i < size; ++i)
            res[i] = int(*mas_fourie.real() + 0.5);

        int N = 0;
        for (int i = 0; i < size; ++i) {
            res[i] += N;
            N = res[i] / 10;
            res[i] %= 10;
        }
    }

    string show() {
        int p;
        for (p = size - 2; p > 0 && !res[p]; p--);
        string r = "";
        for (int i = 0; i <= p; i++) {
            r += '0' + res[p - i];
        }
        cout << r;
        return r;
    }

    ~fourie() {
        delete[] mas_fourie;
        mas_fourie = NULL;
    }
};

int main()
{
    setlocale(LC_ALL, "Russian");

}