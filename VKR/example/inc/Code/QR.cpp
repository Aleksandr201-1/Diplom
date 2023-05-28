template <class T>
T getMainEl (const Matrix<T> &matrix, uint64_t i) {
    T ans = matrix(i, i), tmp = 0;
    T sign = ans >= T(0) ? T(1) : T(-1);
    for (uint64_t j = i; j < matrix.size().n; ++j) {
        tmp += matrix(j, i) * matrix(j, i);
    }
    return ans + sign * std::sqrt(tmp);
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>> QR (const Matrix<T> &matrix) {
    uint64_t n = matrix.size().n;
    Matrix<T> Q(n), R(matrix), v(n, 1);
    auto pred = [] (T el) -> bool {
        return el == T(0);
    };
    for (uint64_t i = 0; i < n - 1; ++i) {
        for (uint64_t j = 0; j < i; ++j)
            v(j, 0) = 0;
        v(i, 0) = getMainEl(R, i);
        for (uint64_t j = i + 1; j < n; ++j) {
            v(j, 0) = R(j, i);
        }
        if (std::count_if(v.toVector().cbegin(), v.toVector().cend(), pred) == n) {
            ++i;
            continue;
        }
        Matrix<T> H = Matrix<T>(n) - 2 * ((v * v.transp()) / (v.transp() * v)(0, 0));
        R = H * R;
        Q = Q * H;
    }
    for (uint64_t i = 0; i < n; ++i)
        for (uint64_t j = 0; j < i; ++j)
            R(i, j) = T(0);
    return std::make_tuple(Q, R);
}