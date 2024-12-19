#include "Task.hpp"

const std::map<TaskType, std::string> task_types = {
    {TaskType::KOSHI, "Koshi"},
    {TaskType::KOSHI_SYSTEM, "KoshiSystem"},
    {TaskType::CHEMICAL, "Chemical"}
};

std::string taskTypeToString (TaskType type) {
    return enumToString(type, task_types);
}

TaskType stringToTaskType (const std::string &str) {
    return stringToEnum(str, task_types);
}

Task::Task () {}

Task::~Task () {}

const std::vector<std::function<double(const std::vector<double> &)>> &Task::getODE () const {
    return ode_system;
}

const std::vector<double> &Task::getY0 () const {
    return Y;
}