#ifndef TASK_HPP
#define TASK_HPP

#include <map>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include <functional>
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
        const std::vector<std::function<double(const std::vector<double> &)>> &getODE () const;
        const std::vector<double> &getY0 () const;
    protected:
        std::vector<std::function<double(const std::vector<double> &)>> ode_system;
        std::vector<double> Y;
};

#endif