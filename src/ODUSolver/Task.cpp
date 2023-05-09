#include "Task.hpp"

const std::map<std::string, TaskType> types = {
    {"Koshi", TaskType::KOSHI},
    {"KoshiSystem", TaskType::KOSHI_SYSTEM},
    {"Chemical", TaskType::CHEMICAL}
};

std::string taskTypeToString (TaskType type) {
    auto check = [type] (const auto &pair) -> bool {
        return pair.second == type;
    };
    auto result = std::find_if(types.begin(), types.end(), check);
    if (result != types.end()) {
        return result->first;
    }
    return "NotAType";
}

TaskType stringToTaskType (const std::string &type) {
    auto it = types.find(type);
    if (it != types.end()) {
        return it->second;
    }
    return TaskType::NOT_A_TYPE;
}

Task::Task () {}

Task::~Task () {}

const std::vector<std::function<float128_t(const std::vector<float128_t> &)>> &Task::getODE () const {
    return ode_system;
}

const std::vector<float128_t> &Task::getY0 () const {
    return Y;
}