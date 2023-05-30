#ifndef KOSHI_TASK_HPP
#define KOSHI_TASK_HPP

#include <vector>
#include <General/General.hpp>
#include <General/FuncMaker.hpp>
#include <ODUSolver/Task.hpp>

uint64_t getOrder (const std::string &task);

class KoshiTask : public Task {
    public:
        KoshiTask ();
        ~KoshiTask ();

        void setTaskInfo(const std::vector<std::string> &system, uint64_t order, float128_t X0, float128_t Xn);
        void setTaskInfo(const std::string &ode, uint64_t order, const std::vector<float128_t> &Y0, float128_t X0, float128_t Xn);
        void setSystemInfo(const std::vector<std::string> &system, uint64_t order, float128_t X0, float128_t Xn);

        std::tuple<float128_t, float128_t> getBorders () const;
    private:
        float128_t X0, Xn;
};

#endif