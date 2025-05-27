#ifndef LONG_INT_HPP
#define LONG_INT_HPP

#include <iostream>
#include <vector>
#include <iomanip>

namespace Math {
    class LongInt {
        private:
            static const uint64_t BASE = 1e9;
            static const uint64_t RADIX = 9;
            void DeleteLeadingZeros ();
        public:
            LongInt ();
            LongInt (uint64_t num);
            LongInt (const std::string& str);
            LongInt (const LongInt& num);
            ~LongInt ();

            void Set (const std::string& str);

            friend const LongInt operator+ (const LongInt& num1, const LongInt& num2);
            friend const LongInt operator- (const LongInt& num1, const LongInt& num2);
            friend const LongInt operator* (const LongInt& num1, const LongInt& num2);
            friend const LongInt operator/ (const LongInt& num1, const LongInt& num2);

            friend const LongInt operator+ (const LongInt& num1, uint64_t num2);
            friend const LongInt operator- (const LongInt& num1, uint64_t num2);
            friend const LongInt operator* (const LongInt& num1, uint64_t num2);
            friend const LongInt operator/ (const LongInt& num1, uint64_t num2);

            static LongInt Pow (const LongInt& num1, const LongInt& num2);

            friend bool operator== (const LongInt& num1, const LongInt& num2);
            friend bool operator!= (const LongInt& num1, const LongInt& num2);
            friend bool operator>  (const LongInt& num1, const LongInt& num2);
            friend bool operator<= (const LongInt& num1, const LongInt& num2);
            friend bool operator<  (const LongInt& num1, const LongInt& num2);
            friend bool operator>= (const LongInt& num1, const LongInt& num2);

            friend bool operator== (const LongInt& num1, uint64_t num2);
            friend bool operator!= (const LongInt& num1, uint64_t num2);

            LongInt& operator= (const LongInt& num);

            friend std::istream& operator>> (std::istream& input, LongInt& num);
            friend std::ostream& operator<< (std::ostream& output, const LongInt& num);
        private:
            std::vector<uint64_t> data;
            bool isNegative;
    };
}

#endif