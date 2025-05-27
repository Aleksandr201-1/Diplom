#include <Math/LongDouble.hpp>

void LongDouble::DeleteLeadingZeros () {
    uint64_t n = (uint64_t) std::max((int64_t)1, exp);
    while (data.size() > n && data[data.size() - 1] == 0) {
        data.erase(data.end() - 1);
    }
    while (data.size() > 1 && data[0] == 0) {
        data.erase(data.begin());
        --exp;
    }
    while (data.size() > 1 && data[data.size() - 1] == 0) {
        data.erase(data.end() - 1);
    }
    if (data.size() == 1 && data[0] == 0) {
        exp = 1;
        isNegative = false;
    }// 343.7000 2
    //while (data[0] % 10 == 0 && data[0] != 0) {
        //--exp;
    //    data[0] /= 10;
    //}
    //data[0] *= 10;
    while (data[0] != 0 && data[0] % 10 == 0) {
        //--exp;
        data[0] /= 10;
    }
}

void LongDouble::Simplify () {
    uint64_t i = 0, rad = RADIX + 1;
    while (data[i] % rad == 0) {
        if (data[i] == 0) {
            data.erase(data.begin());
            exp += RADIX;
        } else {
            data[i] /= rad;
            exp += 1;
        }
    }
}

void LongDouble::Set (std::string str) {
    isNegative = false;
    data.clear();
    if (str.empty()) {
        data.push_back(0);
        exp = 0;
    } else {
        if (str[0] == '-') {
            isNegative = true;
            str.erase(0, 1);
        }
        exp = str.find_first_of('.', 0);
        if (exp != std::string::npos) {
            str.erase(exp, 1);
            //exp -= str.size();
            --exp;
        } else {
            exp = str.size() - 1;
        }
        uint64_t j = str.size() - 1;
        while (j > 0 && str[j] == '0') {
            --j;
        }
        str.erase(str.begin() + j + 1, str.end());
        if (exp == 0 && str[0] == '0') {
            uint64_t i = 1;
            while (i < str.size() && str[i] == '0') {
                --exp;
                ++i;
            }
            --exp;
            str.erase(str.begin(), str.begin() + i);
        }
        if (str.empty()) {
            str += '0';
        }
        std::cout << "exp: " << exp << " for " << str << "\n";
        for (int64_t i = str.size(); i > 0; i -= RADIX) {
            if (i <= RADIX) {
                data.push_back(atoi(str.substr(0, i).c_str()));
            } else {
                data.push_back(atoi(str.substr(i - RADIX, RADIX).c_str()));
            }
        }
        //DeleteLeadingZeros();
    }
}

LongDouble::LongDouble () {
    data.push_back(0);
    exp = 1;
    isNegative = false;
}

LongDouble::LongDouble (std::string str) {
    Set(str);
}

LongDouble::LongDouble (const LongDouble& num) {
    data = num.data;
    exp = num.exp;
    isNegative = num.isNegative;
}

LongDouble::~LongDouble () {}

const LongDouble LongDouble::operator+ () const {
    return LongDouble(*this);
}

const LongDouble LongDouble::operator- () const {
    LongDouble copy(*this);
    copy.isNegative = !copy.isNegative;
    return copy;
}

const LongDouble operator+ (const LongDouble& num1, const LongDouble& num2) {
    if (num1.isNegative == num2.isNegative) {
        int64_t exp1 = num1.exp;
        int64_t exp2 = num2.exp;
        int64_t exp = std::max(exp1, exp2);
        std::vector<uint64_t> d1(num1.data);
        std::vector<uint64_t> d2(num2.data);
        while (exp1 != exp) {
            if (exp1 + LongDouble::RADIX < exp) {
                d1.insert(d1.begin(), 0);
                exp1 += LongDouble::RADIX;
            } else {
                ++exp1;
            }
        }
        while (exp2 != exp) {
            if (exp2 + LongDouble::RADIX < exp) {
                d2.insert(d2.begin(), 0);
                exp2 += LongDouble::RADIX;
            } else {
                ++exp2;
            }
        }
        uint64_t size = std::max(d1.size(), d2.size());
        while (d1.size() != size) {
            d1.push_back(0);
        }
        while (d2.size() != size) {
            d2.push_back(0);
        }
        uint64_t len = 1 + size;
        LongDouble res;
        res.isNegative = num1.isNegative;
        res.data = std::vector<uint64_t>(len, 0);
        for (uint64_t i = 0; i < size; ++i) {
            res.data[i + 1] = d1[i] + d2[i];
        }
        for (uint64_t i = len - 1; i > 0; --i) {
            res.data[i - 1] += res.data[i] / LongDouble::BASE;
            res.data[i] %= LongDouble::BASE;
        }
        res.exp = exp + 1;
        res.DeleteLeadingZeros();
        return res;
    }
    if (num1.isNegative) {
        std::cout << num2 << " - " << (-num1) << "\n";
    }
    return num1.isNegative ? num2 - (-num1) : num1 - (-num2);
}

const LongDouble operator- (const LongDouble& num1, const LongDouble& num2) {
    //std::cout << "wee wee\n";
    if (num1.isNegative == num2.isNegative && num1.isNegative == false) { // если боа числа положительны
        bool cmp = num1 > num2; // получаем флаг того, больше ли первое число
        // if (cmp) {
        //     std::cout << "fizz\n";
        // } else {
        //     std::cout << "buzz\n";
        // }
        uint64_t exp1 = cmp ? num1.exp : num2.exp; // сохраняем экспоненту большего числа
        uint64_t exp2 = cmp ? num2.exp : num1.exp; // сохраняем экспоненту меньшего числа
        uint64_t exp = std::max(exp1, exp2); // определяем максимальную экспоненту, чтобы к ней привести числа

        std::vector<uint64_t> d1(cmp ? num1.data : num2.data); // запоминаем вектор цифр большего числа
        std::vector<uint64_t> d2(cmp ? num2.data : num1.data); // запоминаем вектор цифр меньшего числа

        // выравниваем экспоненты как при сложении (добавляя нули вначале числа)
        while (exp1 != exp) {
            if (exp1 + LongDouble::RADIX < exp) {
                d1.insert(d1.begin(), 0);
                exp1 += LongDouble::RADIX;
            } else {
                ++exp1;
            }
        }

        while (exp2 != exp) {
            if (exp2 + LongDouble::RADIX < exp) {
                d2.insert(d2.begin(), 0);
                exp2 += LongDouble::RADIX;
            } else {
                ++exp2;
            }
        }

        uint64_t size = std::max(d1.size(), d2.size()); // определяем максимальный размер

        // добавляем нули в конец векторов цифр
        while (d1.size() != size) {
            d1.push_back(0);
        }

        while (d2.size() != size) {
            d2.push_back(0);
        }

        uint64_t len = 1 + size;

        LongDouble res; // создаём число для результата

        res.isNegative = !cmp; // знак будет 1, если первое больше второго, и -1, если первое меньше второго
        res.data = std::vector<uint64_t>(len, 0); // создаём вектор из len нулей
    //std::cout << "wee wee\n";
        for (uint64_t i = 0; i < size; ++i) {
            if (d1[i] < d2[i]) {
                d1[i] += LongDouble::BASE;
                --d1[i + 1];
            }
            res.data[i + 1] = d1[i] - d2[i]; // вычитаем соответствующие разряды
        }

        // обрабатываем переполнения
        // for (uint64_t i = len - 1; i > 0; --i) {
        //     if (res.data[i] < 0) { // если текущий разряд стал меньше нуля
        //         res.data[i] += LongDouble::BASE; // занимаем у предыдущего, прибавляя 10 к текущему
        //         --res.data[i - 1]; // уменьшаем на 1 предыдущий разряд
        //     }
        // }

        res.exp = exp + 1; // восстанавливаем экспоненту
        res.DeleteLeadingZeros(); // удаляем лишнии нули

        return res; // возвращаем результат
    }

    if (num1.isNegative == num2.isNegative && num1.isNegative == true) {
        //std::cout << "fizz\n";
        return (-num2) - (-num1);
    }
    // if (sign == -1 && x.sign == -1)
    //     return (-x) - (-(*this)); // если оба числа отрицательны, то из второго с обратным знаком вычитаем первое с обратным знаком
    
    // return *this + (-x); // если знаки разные, прибавляем к первому отрицательное второе
    return num1 + (-num2);
}

const LongDouble operator* (const LongDouble& num1, const LongDouble& num2) {
    LongDouble ans;
    ans.data.resize(num1.data.size() + num2.data.size(), 0);
    ans.exp = num1.exp + num2.exp;
    ans.isNegative = num1.isNegative != num2.isNegative;
    for (uint64_t i = 0; i < num1.data.size(); ++i) {
        for (uint64_t j = 0, carry = 0; j < num2.data.size() || carry; ++j) {
            uint64_t current = ans.data[i + j] + num1.data[i] * (j < num2.data.size() ? num2.data[j] : 0) + carry;
            ans.data[i + j] = current % LongDouble::BASE;
            carry = current / LongDouble::BASE;
        }
    }
    ans.DeleteLeadingZeros();
    return ans;
}

const LongDouble operator/ (const LongDouble& num1, const LongDouble& num2) {
    return num1 * num2.Inverse();
}

LongDouble LongDouble::Inverse () const {
    if (data.size() == 1 && data[0] == 0)
        throw std::string("LongDouble LongDouble::inverse() - division by zero!"); // делить на ноль нельзя, поэтому бросим исключение

    LongDouble x(*this); // скопируем число,
    x.isNegative = false; // сделав его положительным

    LongDouble d("1"), zero("0"); // создадим то, что будем делить на x
    std::cout << "dd: " << d << "\n";
    LongDouble res; // создадит объект для результата
    res.isNegative = isNegative; // знак результата совпадёт со знаком числа
    res.exp = 1; // начнём с единичной экспоненты
    //res.data = std::vector<int>(); // создадим пустой вектор для цифр обратного элемента

    // пока число меньше 1
    while (x < d) {
        ++x.exp; // будем увеличивать его экспоненту (умножать на 10 фактически)
        ++res.exp; // и заодно экспоненту результата
    }
    
    // дальше сдлеаем число d большим x, также умножая его на 10, чтобы получить число 100...0
    while (d < x) {
        ++d.exp;
    }
    std::cout << "RRRR: " << x << " " << d << "\n";
    res.exp -= d.exp - 1; // подсчитаем реальное количество цифр в целой части

    uint64_t numbers = 0; // количество уже вычисленных цифр дробной части
    uint64_t totalNumbers = MAX_LEN + std::max((int64_t) 0, res.exp); // количество цифр с учётом целой части

    do {
        uint64_t div = 0; // будущая цифра
        
        // считаем, сколько раз нужно вычесть x из d
        while (d >= x) {
            ++div; 
            d = d - x;
            std::cout << "iter: " << div << " " << d << "\n";
        }

        // увеличиваем остаток в 10 раз
        d.exp++; 
        d.DeleteLeadingZeros();
        res.data.push_back(div); // записываем сформированную цифру
        ++numbers; // увеличиваем число вычисленных цифр
    } while (d != zero && numbers < totalNumbers); // считаем до тех пор, пока не дойдём до нулевого остатка или пока не превысим точность
    res.DeleteLeadingZeros();
    return res; // возвращаем результат
    //return LongDouble("1.0");
}

LongDouble LongDouble::Pow (const LongDouble& num, int64_t n) {
    LongDouble res("1.0"), x(num);
    int64_t exp = num.exp;
    if (n < 0) {
        n = -n;
        res = res.Inverse();
    }
    while (n != 0) {
        if (n & 1) {
            res = res * x;
            res.exp += exp;
        }
        exp *= 2;
        x = x * x;
        n = n >> 1;
    }
    return res;
}

bool operator== (const LongDouble& num1, const LongDouble& num2) {
    return num1.isNegative == num2.isNegative && num1.exp == num2.exp && num1.data == num2.data;
}

bool operator!= (const LongDouble& num1, const LongDouble& num2) {
    return !(num1.isNegative == num2.isNegative && num1.exp == num2.exp && num1.data == num2.data);
}

bool operator>  (const LongDouble& num1, const LongDouble& num2) {
    if (num1.isNegative != num2.isNegative) {
        return num2.isNegative;
    }
    if (num1.exp != num2.exp) {
        return num1.exp > num2.exp;
    }
    uint64_t size = std::max(num1.data.size(), num2.data.size());
    for (uint64_t i = size - 1; i > 0; --i) {
        if (num1.data[i] != num2.data[i]) {
            return num1.data[i] > num2.data[i] == !num1.isNegative;
        }
    }
    return false;
}

bool operator<= (const LongDouble& num1, const LongDouble& num2) {
    return !(num1 > num2);
}

bool operator<  (const LongDouble& num1, const LongDouble& num2) {
    if (num1.isNegative != num2.isNegative) {
        return num1.isNegative;
    }
    if (num1.exp != num2.exp) {
        return num1.exp < num2.exp;
    }
    uint64_t size = std::max(num1.data.size(), num2.data.size());
    for (uint64_t i = size - 1; i > 0; --i) {
        if (num1.data[i] != num2.data[i]) {
            return num1.data[i] < num2.data[i] == !num2.isNegative;
        }
    }
    return false;
}

bool operator>= (const LongDouble& num1, const LongDouble& num2) {
    return !(num1 < num2);
}

LongDouble& LongDouble::operator= (const LongDouble& num) {
    if (this == &num) {
        return *this;
    }
    data = num.data;
    exp = num.exp;
    isNegative = num.isNegative;
    return *this;
}

LongDouble& LongDouble::operator= (std::string str) {
    this->Set(str);
    return *this;
}

std::istream& operator>> (std::istream& input, LongDouble& num) {
    std::string str;
    input >> str;
    num.Set(str);
    return input;
}

std::ostream& operator<< (std::ostream& output, const LongDouble& num) {
    if (num.isNegative) {
        output << "-";
    }
    char oldFiller = output.fill('0');
    int64_t exp = num.exp, back = num.data.size() - 1;
    if (exp < 0) {
        output << "0.";
        for (int64_t i = -exp - 1; i > 0; --i) {
            output << "0";
        }
        output << num.data.back();
        --back;
    } else {
        if (exp > LongDouble::RADIX) {
            output << num.data[back];
            uint64_t n = num.data[back];
            --back;
            while (n != 0) {
                n /= 10;
                --exp;
            }
            //exp -= LongDouble::RADIX;
        }
        while (exp > LongDouble::RADIX) {
            output << std::setw(LongDouble::RADIX) << num.data[back];
            --back;
            exp -= LongDouble::RADIX;
        }
        uint64_t n = 1, n2 = LongDouble::BASE; 
        int64_t i = exp + 1;
        while (n < num.data[back]) {
            n *= 10;
            //--i;
        }
        for (; i > 0 && n > 1; --i) {
            n /= 10;
        }
        //std::cout << "first: " << num.data[back] << " n1: " << n << " n2: " << n2 << "\n";
        // if (n > num.data[back]) {
        //     output << num.data[back];
        //     --back;
        //     for (i = exp - i - 1; i > 0; --i) {
        //         n /= 10;
        //     }
        // }
        if (i != 0 && num.data.size() > 1) {
            output << num.data[back];
            --back;
            n = LongDouble::BASE;
            for (; i > 0; --i) {
                n /= 10;
            }
            //n2 /= n1;
        }
        //std::cout << "first: " << num.data[back] << "n: " << n << "\n";
        output << num.data[back] / n << "." << (n != 1 ? num.data[back] % n : 0);
        // if (n != 1) {
        //     output << num.data[back] % n;
        // } else {
        //     output << "0";
        // }
        --back;
    }
    for (int64_t i = back; i >= 0; --i) {
        output << std::setw(LongDouble::RADIX) << num.data[i];
    }
    output.fill(oldFiller);
    return output;
}