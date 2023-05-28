#include <Input/InputHandle.hpp>

void help (const std::string &name) {
    std:: cout << "Usage: " << name << " [KEYS] [OPTIONS]\n"
                  "\t-m, --method\t[METHOD]\tРешать задачу методом [METHOD]\n"
                  "\t\tПоддерживаемые методы:\n\t\tRunge-Kutta, Cheskino, Dorman-Prince, Falberg, Gauss, Lobatto, L Stable Diagonal, Merson, Rado\n"
                  "\t-o, --order\t[ORDER]\t\tРешать задачу методом [METHOD] с порядком точности [ORDER]\n"
                  "\t-w, --way\t[WAY]\t\tРешать задачу методом [METHOD] с порядком точности [ORDER] и способом реализации [WAY]\n"
                  "\t-i, --iter\t[ALGORYTHM]\tИспользовать алгоритм [ALGORYTHM] для решения системы уравнений в жёстких методах\n"
                  "\t\tПоддерживаемые алгоритмы:\n\t\tSI, Zeidel, Newton\n"
                  "\t-a, --approx\t[NUM]\t\tПродолжать процесс итераций в жёстких методах пока не будет достигнута точность [NUM]\n"
                  "\t-f, --file\t[FILE]\t\tВзять ввод из файла [FILE]\n"
                  "\t-tt, --task_type\t[TYPE]\t\tРешать задачу заданного типа\n"
                  "\t\tПоддерживаемые типы задач:\n\t\tKoshi, KoshiSystem, Chemical\n";
}

std::string readLine (std::istream &input) {
    std::string str;
    while (str.empty()) {
        std::getline(input, str);
        if (!str.empty() && str[0] == '#') {
            str = "";
        }
    }
    return str;
}

void argsHandler (const std::vector<std::string> &args, ReportInfo &info) {
    bool noArgument = false;
    uint64_t idx = 0;
    try {
        for (uint64_t i = 1; i < args.size(); ++i) {
            if (args[i] == "--help" || args[i] == "-h") {
                help(args[0]);
                exit(0);
            } else if (args[i] == "--method" || args[i] == "-m") {
                if (args.size() > i + 1) {
                    info.method = stringToSolveMethod(args[i + 1]);
                    if (info.method == SolveMethod::ERROR) {
                        throw std::logic_error("argsHandler: solve method \"" + args[i + 1] + "\" not presented");
                    }
                    ++i;
                } else {
                    noArgument = true;
                    idx = i;
                }
            } else if (args[i] == "--order" || args[i] == "-o") {
                if (args.size() > i + 1) {
                    info.order = std::stoull(args[i + 1]);
                    ++i;
                } else {
                    noArgument = true;
                    idx = i;
                }
            } else if (args[i] == "--way" || args[i] == "-w") {
                if (args.size() > i + 1) {
                    info.way = std::stoull(args[i + 1]);
                    ++i;
                } else {
                    noArgument = true;
                    idx = i;
                }
            } else if (args[i] == "--iter" || args[i] == "-i") {
                if (args.size() > i + 1) {
                    if (args[i + 1] == "Newton") {
                        info.algo = IterationAlgo::NEWTON;
                    } else if (args[i + 1] == "Zeidel") {
                        info.algo = IterationAlgo::ZEIDEL;
                    } else if (args[i + 1] == "SI") {
                        info.algo = IterationAlgo::SIMPLE_ITERATION;
                    } else {
                        throw std::logic_error("argsHandler: iteration method \"" + args[i + 1] + "\" not presented");
                    }
                    ++i;
                } else {
                    noArgument = true;
                    idx = i;
                }
            } else if (args[i] == "--approx" || args[i] == "-a") {
                if (args.size() > i + 1) {
                    info.approx = std::stod(args[i + 1]);
                    ++i;
                } else {
                    noArgument = true;
                    idx = i;
                }
            } else if (args[i] == "--file" || args[i] == "-f") {
                if (args.size() > i + 1) {
                    info.fileInput = args[i + 1];
                    ++i;
                } else {
                    noArgument = true;
                    idx = i;
                }
            } else if (args[i] == "--task_type" || args[i] == "-tt") {
                if (args.size() > i + 1) {
                    info.type = stringToTaskType(args[i + 1]);
                    ++i;
                    if (info.type == TaskType::ERROR) {
                        throw std::logic_error("argsHandler: task type \"" + args[i + 1] + "\" not presented");
                    }
                } else {
                    noArgument = true;
                    idx = i;
                }
            }
            if (noArgument) {
                break;
            }
        }
        if (noArgument) {
            throw std::logic_error("argsHandler: no argument given after \"" + args[idx] + "\" key");
        }
    } catch (std::exception &exp) {
        std::cerr << exp.what() << "\n";
        exit(1);
    }
}