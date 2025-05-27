#include <Math/LongInt.hpp>

namespace Math {

    void LongInt::DeleteLeadingZeros () {
        while (data.size() > 1 && data.back() == 0) {
            data.pop_back();
        }
        if (data.size() == 1) {
            isNegative = false;
        } 
    }

    LongInt::LongInt () {
        data.push_back(0);
        isNegative = false;
    }

    LongInt::LongInt (uint64_t num) {
        while (num > 0) {
            data.push_back(num % BASE);
            num = num / BASE;
        }
        isNegative = false;
    }

    LongInt::LongInt (const std::string& str) {
        Set(str);
    }

    LongInt::LongInt (const LongInt& num): data(num.data), isNegative(num.isNegative) {}

    LongInt::~LongInt () {}

    void LongInt::Set (const std::string& str) {
        data.clear();
        isNegative = false;
        if (str.empty()) {
            data.push_back(0);
        } else {
            if (str[0] == '-') {
                isNegative = true;
            }
            for (int64_t i = str.length(); i > 0; i -= RADIX) {
                if (i < RADIX) {
                    data.push_back(atoi(str.substr(0, i).c_str()));
                } else {
                    data.push_back(atoi(str.substr(i - RADIX, RADIX).c_str()));
                }
            }
            DeleteLeadingZeros();
        }
    }

    const LongInt operator+ (const LongInt& num1, const LongInt& num2) {
        uint64_t carry = 0;
        uint64_t size = std::max(num1.data.size(), num2.data.size());
        LongInt res;
        res.data.resize(size);
        for (uint64_t i = 0; i < size; ++i) {
            uint64_t sum = carry;
            if (i < num2.data.size()) {
                sum += num2.data[i];
            }
            if (i < num1.data.size()) {
                sum += num1.data[i];
            }
            carry = sum / LongInt::BASE;
            res.data[i] = sum % LongInt::BASE;
        }
        if (carry) {
            res.data.push_back(carry);
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator- (const LongInt& num1, const LongInt& num2) {
        if (num1 < num2) {
            throw std::logic_error("Error: trying to subtract bigger number from smaller");
        }
        uint64_t size = std::max(num2.data.size(), num1.data.size());
        uint64_t carry = 0;
        LongInt res;
        res.data.resize(size);
        for (uint64_t i = 0; i < size; ++i) {
            int64_t diff = num1.data[i] - carry;
            if (i < num2.data.size()) {
                diff -= num2.data[i];
            }
            if (diff < 0) {
                carry = 1;
                diff += LongInt::BASE;
            } else {
                carry = 0;
            }
            res.data[i] = diff % LongInt::BASE;
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator* (const LongInt& num1, const LongInt& num2) {
        LongInt res;
        res.data.resize(num1.data.size() + num2.data.size(), 0);
        for (uint64_t i = 0; i < num1.data.size(); ++i) {
            for (uint64_t j = 0, carry = 0; j < num2.data.size() || carry; ++j) {
                uint64_t current = res.data[i + j] + num1.data[i] * (j < num2.data.size() ? num2.data[j] : 0) + carry;
                res.data[i + j] = current % LongInt::BASE;
                carry = current / LongInt::BASE;
            }
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator/ (const LongInt& num1, const LongInt& num2) {
        if (num2 == 0) {
            throw std::logic_error("Error: trying to divide by zero");
        }
        LongInt res, curr;
        res.data.resize(num1.data.size());
        for (uint64_t i = num1.data.size() - 1; i < num1.data.size(); --i) {
            curr.data.insert(curr.data.begin(), num1.data[i]);
            curr.DeleteLeadingZeros();
            int64_t x = 0, l = 0, r = LongInt::BASE;
            while (l <= r) {
                int64_t m = (l + r) / 2;
                LongInt tmp = num2 * LongInt(m);
                if (tmp <= curr) {
                    x = m;
                    l = m + 1;
                } else {
                    r = m - 1;
                }
            }
            res.data[i] = x;
            curr = curr - num2 * LongInt(x);
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator+ (const LongInt& num1, uint64_t num2) {
        LongInt res(num1);
        uint64_t i = 0;
        while (num2 > 0) {
            if (i + 1 == res.data.size()) {
                res.data.push_back(0);
            }
            res.data[i] += num2 % LongInt::BASE;
            if (res.data[i] > LongInt::BASE) {
                res.data[i] -= LongInt::BASE;
                ++res.data[i + 1];
            }
            ++i;
            num2 = num2 / LongInt::BASE;
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator- (const LongInt& num1, uint64_t num2) {
        if (num1.data.size() == 1 && num1.data[0] < num2) {
            throw std::logic_error("Error: trying to subtract bigger number from smaller");
        }

        LongInt res(num1);
        uint64_t i = 0, curr = 0;
        while (num2 > 0) {
            curr = num2 % LongInt::BASE;
            if (res.data[i] >= curr) {
                res.data[i] -= curr;
            } else {
                res.data[i] = LongInt::BASE - (curr - res.data[i]);
                --res.data[i + 1];
            }
            ++i;
            num2 = num2 / LongInt::BASE;
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator* (const LongInt& num1, uint64_t num2) {
        LongInt res(num1);
        uint64_t carry = 0;
        for (uint64_t i = 0; i < num1.data.size() || carry > 0; ++i) {
            uint64_t currDigit = carry;
            if (i == num1.data.size()) {
                res.data.push_back(0);
            } else {
                currDigit += num1.data[i] * num2;
            }
            res.data[i] = currDigit % LongInt::BASE;
            carry = currDigit / LongInt::BASE;
        }
        res.DeleteLeadingZeros();
        return res;
    }

    const LongInt operator/ (const LongInt& num1, uint64_t num2) {
        if (num2 == 0) {
            throw std::logic_error("Error: trying to divide by zero");
        }
        LongInt res(num1);
        uint64_t carry = 0;
        for (uint64_t i = res.data.size() - 1; i < res.data.size(); --i) {
            uint64_t currDigit = carry * LongInt::BASE + res.data[i];
            res.data[i] = currDigit / num2;
            carry = currDigit % num2;
        }
        res.DeleteLeadingZeros();
        return res;
    }

    LongInt LongInt::Pow (const LongInt& num1, const LongInt& num2) {
        if (num1 == 0 && num2 == 0) {
            throw std::logic_error("Error: 0^0 is uncertain");
        }
        if (num1 == 1) {
            return num1;
        }
        LongInt res(1), x(num1), n(num2);
        while (n != 0) {
            if (n.data[0] & 1) {
                res = res * x;
            }
            x = x * x;
            n = n / 2;
        }
        return res;
    }

    bool operator== (const LongInt& num1, const LongInt& num2) {
        return num1.data == num2.data;
    }

    bool operator!= (const LongInt& num1, const LongInt& num2) {
        return !(num1 == num2);
    }

    bool operator> (const LongInt& num1, const LongInt& num2) {
        if (num1.data.size() != num2.data.size()) {
            return num1.data.size() > num2.data.size();
        }
        int64_t size = num1.data.size();
        for (int64_t i = size - 1; i >= 0; --i) {
            if (num1.data[i] != num2.data[i]) {
                return num1.data[i] > num2.data[i];
            }
        }
        return false;
    }

    bool operator<= (const LongInt& num1, const LongInt& num2) {
        return !(num1 > num2);
    }

    bool operator< (const LongInt& num1, const LongInt& num2) {
        return !(num1 > num2 || num1 == num2);
    }

    bool operator>= (const LongInt& num1, const LongInt& num2) {
        return !(num1 < num2);
    }

    bool operator== (const LongInt& num1, uint64_t num2) {
        uint64_t i = 0;
        if (num1.data.size() == 1 && num2 == num1.data[i]) {
            return true;
        }
        while (num2 > 0) {
            if (num1.data.size() == i || num1.data[i] != num2 % LongInt::BASE) {
                return false;
            }
            ++i;
            num2 = num2 / LongInt::BASE;
        }
        return num1.data.size() == i;
    }

    bool operator!= (const LongInt& num1, uint64_t num2) {
        return !(num1 == num2);
    }

    LongInt& LongInt::operator= (const LongInt& num) {
        if (this == &num) {
            return *this;
        }
        data.resize(num.data.size());
        for (uint64_t i = 0; i < this->data.size(); ++i) {
            data[i] = num.data[i];
        }
        return *this;
    }

    std::istream& operator>> (std::istream& input, LongInt& num) {
        std::string str;
        input >> str;
        num.Set(str);
        return input;
    }

    std::ostream& operator<< (std::ostream& output, const LongInt& num) {
        int64_t begin = static_cast<int64_t>(num.data.size()) - 1;
        output << num.data.back();
        for (int64_t i = begin - 1; i >= 0; --i) {
            output << std::setfill('0') << std::setw(LongInt::RADIX) << num.data[i];
        }
        return output;
    }
}