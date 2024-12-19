#ifndef KOSHI_TASK_HPP
#define KOSHI_TASK_HPP

#include <vector>
#include <Math/FuncMaker.hpp>
#include <ODUSolver/Task.hpp>

uint64_t getOrder (const std::string &task);

class KoshiTask : public Task {
    public:
        KoshiTask ();
        ~KoshiTask ();

        void setTaskInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn);
        void setTaskInfo(const std::string &ode, uint64_t order, const std::vector<double> &Y0, double X0, double Xn);
        void setSystemInfo(const std::vector<std::string> &system, uint64_t order, double X0, double Xn);

        std::tuple<double, double> getBorders () const;
    private:
        double X0, Xn;
};

#endif