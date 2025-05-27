#ifndef LONG_DOUBLE_HPP
#define LONG_DOUBLE_HPP

#include <iomanip>
#include <vector>
#include <iostream>

class LongDouble {
    private:
        static const uint64_t BASE = 1000000000;
        static const uint64_t RADIX = 9;
        static const uint64_t MAX_LEN = 10; //количество чисел после запятой
        void DeleteLeadingZeros (); //удаление лишних нулей спереди и сзади
        void Simplify ();
        void Set (std::string str); //строку в число
    public:
        //конструкторы
        LongDouble ();
        LongDouble (std::string str); //строку в число
        LongDouble (const LongDouble& num); //копирование
        ~LongDouble ();

        //унарные + и -
        const LongDouble operator+ () const;
        const LongDouble operator- () const;

        //арифметические операции
        friend const LongDouble operator+ (const LongDouble& num1, const LongDouble& num2);
        friend const LongDouble operator- (const LongDouble& num1, const LongDouble& num2);
        friend const LongDouble operator* (const LongDouble& num1, const LongDouble& num2);
        friend const LongDouble operator/ (const LongDouble& num1, const LongDouble& num2);

        //особые функции
        LongDouble Inverse () const;
        static LongDouble Pow (const LongDouble& num, int64_t n);
        //static LongDouble Pow  (const LongDouble& num, const LongDouble& n);
        //static LongDouble Sqrt (const LongDouble& num);
        //static LongDouble Sin  (const LongDouble& num1, const LongDouble& num2);
        //static LongDouble Cos  (const LongDouble& num1, const LongDouble& num2);
        //static LongDouble Tg   (const LongDouble& num1, const LongDouble& num2);
        //static LongDouble Ctg  (const LongDouble& num1, const LongDouble& num2);
        //static LongDouble Log  (const LongDouble& num1, const LongDouble& num2);
        //static LongDouble Ln   (const LongDouble& num);

        //операторы сравнения
        friend bool operator== (const LongDouble& num1, const LongDouble& num2);
        friend bool operator!= (const LongDouble& num1, const LongDouble& num2);
        friend bool operator>  (const LongDouble& num1, const LongDouble& num2);
        friend bool operator<= (const LongDouble& num1, const LongDouble& num2);
        friend bool operator<  (const LongDouble& num1, const LongDouble& num2);
        friend bool operator>= (const LongDouble& num1, const LongDouble& num2);

        //операторы присваивания
        LongDouble& operator= (const LongDouble& num);
        LongDouble& operator= (std::string str);

        //ввод-вывод
        friend std::istream& operator>> (std::istream& input, LongDouble& num);
        friend std::ostream& operator<< (std::ostream& output, const LongDouble& num);
    private:
        std::vector<uint64_t> data; //массив чисел
        int64_t exp; //экспонента
        bool isNegative; //знак
};

const LongDouble Pi("3.14");
const LongDouble E("2.71");

#endif