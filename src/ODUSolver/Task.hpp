#ifndef TASK_HPP
#define TASK_HPP

#include <map>
#include <algorithm>
#include <General/General.hpp>
#include <General/Enum.hpp>

const uint64_t MIN_H_SIZE = 100;
const uint64_t MAX_H_SIZE = 10;

enum class TaskType {
    KOSHI,
    KOSHI_SYSTEM,
    CHEMICAL,
    ERROR
};

std::string taskTypeToString (TaskType type);

TaskType stringToTaskType (const std::string &str);

class Task {
    public:
        Task ();
        virtual ~Task ();
        const std::vector<std::function<float128_t(const std::vector<float128_t> &)>> &getODE () const;
        const std::vector<float128_t> &getY0 () const;
    protected:
        std::vector<std::function<float128_t(const std::vector<float128_t> &)>> ode_system;
        std::vector<float128_t> Y;
};

#endif