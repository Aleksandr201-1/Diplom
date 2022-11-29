//reverse rang addCol addRow det
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstdint>
#include <vector>
#include <fstream>

struct MatrixSize {
    uint64_t n, m;
};

template <class T>
class Matrix {
    public:
        //конструктор пустой матрицы
        Matrix ();
        //конструктор единичной матрицы n x n
        Matrix (uint64_t n);
        //конструктор нулевой матрицы n строк на m столбцов
        Matrix (uint64_t n, uint64_t m);
        //конструктор копирования
        Matrix (const Matrix<T> &m);
        //конструктор преобразования вектора в матрицу n x m
        Matrix (uint64_t n, uint64_t m, const std::vector<T> &vec);
        //стандартный деструктор
        ~Matrix () = default;

        //текущие размеры матрицы
        MatrixSize size () const;
        //изменение размеров матрицы
        void resize (uint64_t newN, uint64_t newM, T el = 0);
        //преобразовать всю матрицу в одномерный массив
        const std::vector<T> &toVector () const;

        //унарные + и -
        const Matrix<T> operator+ () const;
        const Matrix<T> operator- () const;

        //арифметические операции
        template <class N>
        friend const Matrix<N> operator+ (const Matrix<N> &m1, const Matrix<N> &m2);
        template <class N>
        friend const Matrix<N> operator- (const Matrix<N> &m1, const Matrix<N> &m2);
        template <class N>
        friend const Matrix<N> operator* (const Matrix<N> &m1, const Matrix<N> &m2);
        //template <class N>
        //friend const N operator* (const Matrix<N> &m1, const Matrix<N> &m2);
        template <class N, class V>
        friend const Matrix<N> operator* (const Matrix<N> &m1, const V &m2);
        template <class N, class V>
        friend const Matrix<N> operator* (const V &m1, const Matrix<N> &m2);
        template <class N>
        friend const Matrix<N> operator/ (const Matrix<N> &m1, const Matrix<N> &m2);
        template <class N, class V>
        friend const Matrix<N> operator/ (const Matrix<N> &m1, const V &m2);

        //операторы сравнения
        template <class N>
        friend bool operator== (const Matrix<N> &m1, const Matrix<N> &m2);
        template <class N>
        friend bool operator!= (const Matrix<N> &m1, const Matrix<N> &m2);

        //оператор присваивания
        Matrix<T> &operator= (const Matrix<T> &mat);

        //проверка матрицы на квадратность
        bool isSquare () const;
        //проверка матрицы на диагональность
        bool isDiagonal () const;
        //проверка матрицы на симметричность
        bool isSymmetric () const;
    
        //получение элемента по индексу
        T at (uint64_t i, uint64_t j) const;
        //получение элемента по индексу
        T &operator() (uint64_t i, uint64_t j);
        //получение элемента по индексу
        T operator() (uint64_t i, uint64_t j) const;
        //добавить столбец
        void addCol (uint64_t i, const std::vector<T> &vec);
        //добавить строку
        void addRow (uint64_t i, const std::vector<T> &vec);
        //поменять местами колонки i и j
        void swapCols (uint64_t i, uint64_t j);
        //поменять местами строки i и j
        void swapRows (uint64_t i, uint64_t j);
        //транспонирование матрицы
        Matrix<T> transp () const;
        //обратная матрица
        Matrix<T> reverse () const;
        //определитель матрицы
        T det () const;
        //след матрицы
        T trace () const;
        //возведение в степень
        Matrix<T> pow (int64_t n) const;
        //ранг матрицы
        uint64_t rang () const;

        //ввод-вывод
        template <class N>
        friend std::istream &operator>> (std::istream &input, Matrix<N> &m);
        template <class N>
        friend std::ostream &operator<< (std::ostream &output, const Matrix<N> &m);

        //чтение-запись для файлов
        template <class N>
        friend std::ifstream &operator>> (std::ifstream &file, Matrix<N> &m);
        template <class N>
        friend std::ofstream &operator<< (std::ofstream &file, const Matrix<N> &m);
    private:
        std::vector<T> buff;
        uint64_t n, m;
};

//конструкторы

template <class T>
Matrix<T>::Matrix () : n(0), m(0) {}

template <class T>
Matrix<T>::Matrix (uint64_t n) : buff(n * n, T(0)), n(n), m(n) {
    for (uint64_t i = 0; i < n; ++i) {
        buff[i * n + i] = T(1);
    }
}

template <class T>
Matrix<T>::Matrix (uint64_t n, uint64_t m) : buff(n * m, T(0)), n(n), m(m) {}

template <class T>
Matrix<T>::Matrix (const Matrix<T> &matrix) {
    buff = matrix.buff;
    n = matrix.n;
    m = matrix.m;
}

template <class T>
Matrix<T>::Matrix (uint64_t n, uint64_t m, const std::vector<T> &vec) : buff(vec), n(n), m(m) {}

//текущие размеры матрицы
//n - количество строк
//m - количество столбцов
template <class T>
MatrixSize Matrix<T>::size () const {
    return {n, m};
}

//изменение размеров матрицы
template <class T>
void Matrix<T>::resize (uint64_t newN, uint64_t newM, T el) {
    if (newN == n && newM == m) {
        return;
    }
    std::vector<T> newBuff(newN * newM, el);
    uint64_t lowestN = std::min(n, newN), lowestM = std::min(m, newM);
    for (uint64_t i = 0; i < lowestN; ++i) {
        for (uint64_t j = 0; j < lowestM; ++j) {
            newBuff[i * newM + j] = buff[i * m + j];
        }
    }
    buff = newBuff;
    n = newN;
    m = newM;
}

//преобразовать всю матрицу в одномерный массив
template <class T>
const std::vector<T> &Matrix<T>::toVector () const {
    return buff;
}

//унарные + и -
template <class T>
const Matrix<T> Matrix<T>::operator+ () const {
    Matrix<T> matrix(*this);
    return matrix;
}

template <class T>
const Matrix<T> Matrix<T>::operator- () const {
    Matrix<T> matrix(*this);
    for (uint64_t i = 0; i < matrix.buff.size(); ++i) {
        matrix.buff[i] *= T(-1);
    }
    return matrix;
}

//арифметические операции

//сложение матриц
template <class T>
const Matrix<T> operator+ (const Matrix<T> &m1, const Matrix<T> &m2) {
    if (m1.n != m2.n || m1.m != m2.m) {
        throw std::logic_error("Matrices have different sizes");
    }
    Matrix<T> ans(m1);
    for (uint64_t i = 0; i < ans.buff.size(); ++i) {
        ans.buff[i] += m2.buff[i];
    }
    return ans;
}

//вычитание матриц
template <class T>
const Matrix<T> operator- (const Matrix<T> &m1, const Matrix<T> &m2) {
    if (m1.n != m2.n || m1.m != m2.m) {
        throw std::logic_error("Matrices have different sizes");
    }
    Matrix<T> ans(m1);
    for (uint64_t i = 0; i < ans.buff.size(); ++i) {
        ans.buff[i] -= m2.buff[i];
    }
    return ans;
}

//умножение матриц
template <class T>
const Matrix<T> operator* (const Matrix<T> &m1, const Matrix<T> &m2) {
    if (m1.m != m2.n) {
        throw std::logic_error("Matrices have different sizes");
    }
    Matrix<T> ans(m1.n, m2.m);
    for (uint64_t i = 0; i < ans.n; ++i) {
        for (uint64_t j = 0; j < ans.m; ++j) {
            T el = 0;
            for(uint64_t k = 0; k < m1.m; ++k) {
                //el += m1(i, k) * m2(k, j);
                el += m1.buff[i * m1.m + k] * m2.buff[k * m2.m + j];
            }
            ans.buff[i * ans.m + j] = el;
            //ans(i, j) = el;
        }
    }
    return ans;
}

// //умножение вектора-строки на вектор-столбец
// template <class T>
// const T operator* (const Matrix<T> &m1, const Matrix<T> &m2) {
//     if (m1.m != m2.n) {
//         throw std::logic_error("Matrices have different sizes");
//     }
//     if (m1.n != m2.m || m1.n != 1) {
//         throw std::logic_error("Matrices are not vectors");
//     }
//     T ans = 0;
//     for(uint64_t i = 0; i < m1.m; ++i) {
//         ans += m1.buff[i] * m2.buff[i];
//     }
//     return ans;
// }

//умножение матрицы на число
template <class T, class V>
const Matrix<T> operator* (const Matrix<T> &m1, const V &m2) {
    Matrix<T> ans(m1);
    T el = m2;
    for (uint64_t i = 0; i < ans.buff.size(); ++i) {
        ans.buff[i] *= el;
    }
    return ans;
}

//умножение матрицы на число
template <class T, class V>
const Matrix<T> operator* (const V &m1, const Matrix<T> &m2) {
    return m2 * m1;
}

//деление матриц
template <class T>
const Matrix<T> operator/ (const Matrix<T> &m1, const Matrix<T> &m2) {
    return m1 * m2.reverse();
}

//деление матрицы на число
template <class T, class V>
const Matrix<T> operator/ (const Matrix<T> &m1, const V &m2) {
    Matrix<T> ans(m1);
    T el = m2;
    for (uint64_t i = 0; i < ans.buff.size(); ++i) {
        ans.buff[i] /= el;
    }
    return ans;
}

//проверка на равенство
template <class T>
bool operator== (const Matrix<T> &m1, const Matrix<T> &m2) {
    if (m1.n != m2.n || m1.m != m2.m) {
        return false;
    }
    return m1.buff == m2.buff;
}

//проверка на неравенство
template <class T>
bool operator!= (const Matrix<T> &m1, const Matrix<T> &m2) {
    return !(m1 == m2);
}

//оператор присваивания
template <class T>
Matrix<T> &Matrix<T>::operator= (const Matrix<T> &mat) {
    buff = mat.buff;
    n = mat.n;
    m = mat.m;
    return *this;
}

//проверка матрицы на квадратность
template <class T>
bool Matrix<T>::isSquare () const {
    return n == m;
}

//проверка матрицы на диагональность
template <class T>
bool Matrix<T>::isDiagonal () const {
    if (n != m) {
        return false;
    }
    T null(0);
    for (uint64_t i = 0; i < n; ++i) {
        if (buff[i * n + i] == null) {
            return false;
        }
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = i + 1; j < n; ++j) {
            if (buff[i * n + j] != null || buff[j * n + i] != null) {
                return false;
            } 
        }
    }
    return true;
}

//проверка матрицы на симметричность
template <class T>
bool Matrix<T>::isSymmetric () const {
    if (n != m) {
        return false;
    }
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 1; j < n - i; ++j) {
            if (buff[i * n + j] != buff[j * n + i]) {
                return false;
            }
        }
    }
    return true;
}

//получение элемента по индексу
template <class T>
T Matrix<T>::at (uint64_t i, uint64_t j) const {
    if (i > n || j > m) {
        std::logic_error("Matrix: out of range");
    }
    return buff[i * m + j];
}

//получение элемента по индексу
template <class T>
T &Matrix<T>::operator() (uint64_t i, uint64_t j) {
    if (i >= n || j >= m) {
        std::logic_error("Matrix: out of range");
    }
    return buff[i * m + j];
}

//получение элемента по индексу
template <class T>
T Matrix<T>::operator() (uint64_t i, uint64_t j) const {
    if (i >= n || j >= m) {
        std::logic_error("Matrix: out of range");
    }
    return buff[i * m + j];
}

//добавить столбец
template <class T>
void Matrix<T>::addCol (uint64_t i, const std::vector<T> &vec) {

}

//добавить строку
template <class T>
void Matrix<T>::addRow (uint64_t i, const std::vector<T> &vec) {
    if (vec.size() > m || i > n) {
        throw std::logic_error("");
    }
    buff.resize(buff.size() + m);
}

//поменять местами колонки i и j
template <class T>
void Matrix<T>::swapCols (uint64_t i, uint64_t j) {
    if (i >= m || j >= m) {
        throw std::out_of_range("Operation \"swapCols\" out of range.");
    }
    for (uint64_t k = 0; k < n; ++k) {
        std::swap(buff[k * m + i], buff[k * m + j]);
    }
}

//поменять местами строки i и j
template <class T>
void Matrix<T>::swapRows (uint64_t i, uint64_t j) {
    if (i >= n || j >= n) {
        throw std::out_of_range("Operation \"swapRows\" out of range.");
    }
    for (uint64_t k = 0; k < m; ++k) {
        std::swap(buff[i * m + k], buff[j * m + k]);
    }
}

//транспонирование матрицы
template <class T>
Matrix<T> Matrix<T>::transp () const {
    Matrix<T> tr(m, n);
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            tr.buff[j * n + i] = buff[i * m + j];
        }
    }
    return tr;
}

//обратная матрица
template <class T>
Matrix<T> Matrix<T>::reverse () const {
    T determinator = det();
    if (determinator == 0) {
        throw std::runtime_error("Can\'t create reversed matrix");
    }
    Matrix<T> tmp(*this);
    return tmp;
}

//определитель матрицы
// 1 2 3 4  1   0
// 2 3 4 5  4:  0 0 0 0   
// 3 4 5 6  6:  -1 -2 -3 -1 -2 -1 ;; 1 2 4   5 3 1   6 3 2   6 5 4
// 4 5 6 7  4:  4 5 6 7
//not working
template <class T>
T Matrix<T>::det () const {
    if (n != m) {
        throw std::logic_error("Operation \"det\" available only for square matrix");
    }
    T ans = 1, null = 0;
    Matrix<T> tmp(*this);
    //создаём ступеньки
    for (uint64_t i = 0; i < n; ++i) {
        if (tmp.buff[i * n + i] == null) {
            uint64_t l = i;
            while (l < n) {
                if (tmp.buff[i * n + l] != null) {
                    break;
                }
                ++l;
            }
            if (l == n) {
                return null;
            } else {
                tmp.swapCols(i, l);
                ans *= T(-1);
            }
        }
        for (uint64_t j = i + 1; j < n; ++j) {
            if (tmp.buff[j * n + i] == null) {
                continue;
            }
            T c = -tmp.buff[j * n + i] / tmp.buff[i * n + i];
            for (uint64_t k = i; k < n; ++k) {
                tmp.buff[j * n + k] += c * tmp.buff[i * n + k];
            }
        }
    }
    for (uint64_t i = 0; i < n; ++i) {
        ans *= tmp.buff[i * n + i];
    }
    return ans;
}

//след матрицы
template <class T>
T Matrix<T>::trace () const {
    if (n != m) {
        throw std::logic_error("Operation \"trace\" available only for square matrix");
    }
    T ans = 0;
    for (uint64_t i = 0; i < n; ++i) {
        ans += buff[i * n + i];
    }
    return ans;
}

//возведение в степень
template <class T>
Matrix<T> Matrix<T>::pow (int64_t n) const {
    if (this->n != this->m) {
        throw std::logic_error("Operation \"pow\" available only for square matrix");
    }
    if (n < 0) {
        return this->reverse().pow(-n);
    }
    if (n == 1) {
        return *this;
    }
    if (n == 0) {
        return Matrix<T>(this->n);
    }
    Matrix<T> ans = Matrix<T>(this->n), A = *this;
    while (n != 0) {
        if (n & 1) {
            ans = ans * A;
        }
        A = A * A;
        n /= 2;
    }
    return ans;
}

//ранг матрицы
template <class T>
uint64_t Matrix<T>::rang () const {
    uint64_t ans = 0;
    T null = 0;
    Matrix<T> tmp(*this);
    //создаём ступеньки
    for (uint64_t i = 0; i < m; ++i) {
        if (tmp.buff[i * n + i] == null) {
            uint64_t l = i;
            while (l < n) {
                if (tmp.buff[i * n + l] != null) {
                    break;
                }
                ++l;
            }
            if (l == n) {
                return null;
            } else {
                tmp.swapCols(i + 1, l + 1);
                ans *= T(-1);
            }
        }
        for (uint64_t j = i + 1; j < n; ++j) {
            if (tmp.buff[j * n + i] == null) {
                continue;
            }
            T c = -tmp.buff[j * n + i]/tmp.buff[i * n + i];
            for (uint64_t k = i; k < n; ++k) {
                tmp.buff[j * n + k] += c * tmp.buff[i * n + k];
            }
        }
    }
    for (uint64_t i = 0; i < n; ++i) {

    }
    return ans;
}

//ввод-вывод
template <class T>
std::istream &operator>> (std::istream &input, Matrix<T> &m) {
    for (uint64_t i = 0; i < m.n; ++i) {
        for (uint64_t j = 0; j < m.m; ++j) {
            input >> m.buff[i * m.m + j];
        }
    }
    return input;
}

template <class T>
std::ostream &operator<< (std::ostream &output, const Matrix<T> &m) {
    for (uint64_t i = 0; i < m.n; ++i) {
        for (uint64_t j = 0; j < m.m - 1; ++j) {
            output << m.buff[i * m.m + j] << ' ';
        }
        output << m.buff[i * m.m + m.m - 1] << '\n';
    }
    return output;
}

//чтение-запись для файлов
template <class T>
std::ifstream &operator>> (std::ifstream &input, Matrix<T> &m) {
    input.read(reinterpret_cast<char*>(&m.n), sizeof(m.n));
    input.read(reinterpret_cast<char*>(&m.m), sizeof(m.m));
    m.buff.resize(m.m * m.n);
    for (uint64_t i = 0; i < m.m * m.n; ++i) {
        input.read(reinterpret_cast<char*>(&m.buff[i]), sizeof(T));
    }
    return input;
}

template <class T>
std::ofstream &operator<< (std::ofstream &output, const Matrix<T> &m) {
    output.write(reinterpret_cast<const char*>(&m.n), sizeof(m.n));
    output.write(reinterpret_cast<const char*>(&m.m), sizeof(m.m));
    for (uint64_t i = 0; i < m.buff.size(); ++i) {
        output.write(reinterpret_cast<const char*>(&m.buff[i]), sizeof(T));
    }
    return output;
}

#endif