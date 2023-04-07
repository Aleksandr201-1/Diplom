#include "ReportGenerator.hpp"

std::vector<std::vector<double>> getAnaliticSolution (const std::vector<double> &X, const std::vector<FunctionalTree> &func) {
    std::vector<std::vector<double>> Ya(func.size(), std::vector<double>(X.size()));
    for (uint64_t i = 0; i < func.size(); ++ i) {
        for (uint64_t j = 0; j < X.size(); ++j) {
            Ya[i][j] = func[i](X[j]);
        }
    }
    return Ya;
}

std::tuple<double, double> getAnaliticCompare (const std::vector<double> &Yn, const std::vector<double> &Ya) {
    double max_miss = 0, average_miss = 0;
    for (uint64_t i = 0; i < Yn.size(); ++i) {
        double diff = std::abs(Yn[i] - Ya[i]);
        max_miss = std::max(max_miss, diff);
        average_miss += diff;
    }
    average_miss /= Yn.size();
    return std::make_tuple(max_miss, average_miss);
}

void printTable (const std::vector<std::vector<double>> &Yn, const std::vector<std::vector<double>> &Ya, std::ostream &out) {
    const uint64_t WIDTH = 12;
    const std::string space = "______________";
    const std::vector<double> &X = Yn[0];
    out << "Таблица:\n";
    for (uint64_t i = 0; i < Yn.size() + Ya.size(); ++i) {
        out << space;
    }
    for (uint64_t i = 0; i < Yn.size() + Ya.size() - 1; ++i) {
        out << "_";
    }
    //out << " ____________________________________________\n";
    out << "| " << std::setw(WIDTH) << "X";
    for (uint64_t i = 0; i < Ya.size(); ++i) {
        out << " | " << std::setw(WIDTH) << "Y nu" << i + 1;
        out << " | " << std::setw(WIDTH) << "Y an" << i + 1;
    }
    for (uint64_t i = Ya.size(); i < Yn.size() - 1; ++i) {
        out << " | " << std::setw(WIDTH) << "Y nu" << i + 1;
    }
    out << " |\n";
    for (uint64_t i = 0; i < Yn.size() + Ya.size(); ++i) {
        out << "|" << space;
    }
    out << "|\n";
    //out << "|______________|______________|______________|\n";
    for (uint64_t i = 0; i < X.size(); ++i) {
        out << "| " << std::setw(WIDTH) << X[i];
        for (uint64_t j = 0; j < Ya.size(); ++j) {
            out << " | " << std::setw(WIDTH) << Yn[j + 1][i];
            out << " | " << std::setw(WIDTH) << Ya[j][i];
        }
        for (uint64_t j = Ya.size(); j < Yn.size() - 1; ++j) {
            out << " | " << std::setw(WIDTH) << Yn[j + 1][i];
        }
        out << " |\n";
        for (uint64_t i = 0; i < Yn.size() + Ya.size(); ++i) {
            out << "|"  << space;
        }
        out << "|\n";
        //out << "|______________|______________|______________|\n";
    }
}

std::string toPgfplotsString (const std::string &str) {
    std::string new_str = str;
    //sinv(sinv()()) sinv() (sinv(()))
    //sin(deg(deg(agrg))) + sin()
    //0 - 2 * cos(x) ^ 2 + cos(x)
    //(sinv)(\(.+\)) //sinv(deg$2) 
    std::vector<std::string> to_rework = {"sin(", "cos(", "tan(", "cot("};
    for (auto el : to_rework) {
        //new_str = std::regex_replace(new_str, std::regex("(" + el + ")(\\(.+\\))"), el + "(deg$2)");
        uint64_t i = 0;
        while (i < new_str.size()) {
            uint64_t idx = new_str.find(el, i), brace = 0, start = idx + el.size(), end = start - 1;
            std::string sub_str;
            if (idx  == -1) {
                break;
            }
            std::cout << "idx = " << idx << "\n";
            do {
                //sub_str += new_str[end];
                if (new_str[end] == '(') {
                    ++brace;
                } else if (new_str[end] == ')') {
                    --brace;
                }
                ++end;
            } while (brace != 0);
            //sub_str = "deg(" + toPgfplotsString(new_str.substr(start, end - start)) + ")";
            sub_str = "deg(" + new_str.substr(start, end - start) + ")";
            std::cout << "sub str: " << el + sub_str << "\n";
            //new_str = std::regex_replace(new_str, std::regex(el + "\(.\)"), "deg(");
            i = start + el.size();
            new_str = new_str.substr(0, start) + sub_str + new_str.substr(end, new_str.size() - end);
        }
    }
    return new_str;
}

void generateTEXT (const ReportInfo &info, std::ostream &out) {
    std::string divider = "============================\n";
    double X0 = info.task.X0, Xn = info.task.Xn;
    out << "Отчёт по решению системы ОДУ\n";
    out << divider;
    out << "Задача:\n";
    for (uint64_t i = 0; i < info.input_task.size(); ++i) {
        out << info.input_task[i] << "\n";
    }
    out << "Порядок задачи: " << info.task.order << "\n"
           "h = " << info.h << "\n"
           "X in [" << X0 << ", " << Xn << "]\n";
    out << divider;
    out << "Метод решения: " << solveMethodToString(info.method) << "\n"
           "Порядок точности: " << info.order << "\n"
           "Способ: " << info.way << "\n"
           "Метод итерации: ";
    switch (info.algo) {
        case IterationAlgo::EXPLICIT_STEP:
            std::cout << "явный шаг\n";
            break;
        case IterationAlgo::SIMPLE_ITERATION:
            std::cout << "метод простых итераций\n";
            break;
        case IterationAlgo::ZEIDEL:
            std::cout << "метод Зейделя\n";
            break;
        case IterationAlgo::NEWTON:
            std::cout << "метод Ньютона\n";
            break;
    }
    out << divider;
    out << "Таблица Бутчера:\n" << info.butcher;
    out << divider;
    out << "Жёсткость задачи: " << info.tough_coeff << "\n"
           "Задача " << (info.tough_coeff > 100.0 ? "" : "не ") << "жёсткая\n";
    out << divider;
    out << "Решение задачи\n";
    //auto &X = info.solution.first, &Yn = info.solution.second;
    const auto &Yn = info.solution;
    auto Ya = getAnaliticSolution(info.solution[0], info.analitic);
    printTable(info.solution, Ya, out);
    for (uint64_t i = 0; i < Ya.size(); ++i) {
        double max_miss, average_miss;
        std::tie(max_miss, average_miss) = getAnaliticCompare(Yn[i + 1], Ya[i]);
        out << "\nСреднее отклонение от аналитического решения для функции " << i << ": " << average_miss << "\n"
               "Максимальное отклонение от аналитического решения для функции " << i << ": " << max_miss << "\n";
    }
    // if (!std::isnan(info.analitic(0))) {
    //     double max_miss, average_miss;
    //     std::tie(max_miss, average_miss) = getAnaliticCompare(Yn, Ya);
    //     out << "\nСреднее отклонение от аналитического решения: " << average_miss << "\n"
    //            "Максимальное отклонение от аналитического решения: " << max_miss << "\n";
    // }
    out << divider;
    for (uint64_t i = 1; i < Yn.size(); ++i) {
        auto func = LeastSquareMethod(Yn[0], Yn[i], 3);
        out << "Приближающий полином 3й степени для функции " << i << ": " << LSMToText(func) << "\n";
    }
}

void generateTEX (const ReportInfo &info, std::ostream &out) {
    auto time = std::chrono::system_clock::now();
    std::time_t date = std::chrono::system_clock::to_time_t(time);
    std::string date_str(std::ctime(&date));
    date_str.pop_back();

    //подключение нужных пакетов
    out << "\\documentclass[a4paper,14pt]{extarticle}\n"
           "\\usepackage[left=1.5cm,right=1.5cm,top=2cm,bottom=2cm]{geometry}\n"
           "\\usepackage[T2A]{fontenc}\n"
           "\\usepackage[utf8x]{inputenc}\n"
           "\\usepackage[russian]{babel}\n"
           "\\usepackage{amsmath}\n"
           "\\usepackage{graphicx}\n"
           "\\usepackage[table]{xcolor}\n"
           "\\usepackage{titlesec}\n"
           "\\usepackage{longtable}\n"
           "\\usepackage{tocloft}\n"
           "\\usepackage{pgfplots}\n"
           "\\pgfplotsset{width=10cm, compat=1.16}\n\n";

    //изменение формата заголовков
    out << "\\titleformat{\\section}\n"
           "{\\normalfont\\bfseries}{}{0pt}{\\fontsize{16}{12}\\selectfont}\n"
           "\\titlespacing*{\\section}{\\parindent}{5ex}{1.0ex}\n"
           "\\titleformat{\\subsection}\n"
           "{\\normalfont\\bfseries}{}{0pt}{\\fontsize{14}{12}\\selectfont}\n"
           "\\titlespacing*{\\subsection}{\\parindent}{0.5ex}{1ex}\n\n";

    //многоточие в содержании
    out << "\\setcounter{secnumdepth}{0}\n"
           "\\setlength{\\cftbeforesecskip}{1mm}\n"
           "\\setlength{\\cftbeforesubsecskip}{1mm}\n"
           "\\setlength{\\parindent}{0pt}\n"
           "\\renewcommand\\cftsecdotsep{\\cftdot}\n"
           "\\renewcommand\\cftsubsecdotsep{\\cftdot}\n\n";

    //титульная страница
    out << "\\title{Отчёт по решению системы ОДУ}\n"
           "\\author{LaTeX-версия}\n"
           "\\date{";
    out << date_str << "}\n\n"
           "\\begin{document}\n\n"
           "\\maketitle\n\n"
           "\\tableofcontents\n"
           "\\pagebreak\n\n";

    //Задача
    std::vector<std::string> var_list;
    double X0 = info.task.X0, Xn = info.task.Xn;
    uint64_t idx = info.input_task[0].find('=');
    var_list.push_back("x");
    var_list.push_back("y");
    for (uint64_t i = 0; i < info.task.order; ++i) {
        var_list.push_back(var_list.back());
        var_list.back() += "'";
    }
    out << "\\section{Задача}\n\n"
           "$$\n"
           "\\begin{cases}\n";
    out << "\t" << FunctionalTree(info.input_task[0].substr(0, idx), var_list).toString(FunctionalTree::Style::LATEX) << " = 0\\\\\n";
    for (uint64_t i = 1; i < info.input_task.size(); ++i) {
        //auto func = FunctionalTree(info.input_task[i]);
        out << "\t" << info.input_task[i] << "\\\\\n";
    }
    out << "\tx \\in [" << X0 << ", " << Xn << "]\n"
           "\\end{cases}\n"
           "$$\n\n"
           "Порядок задачи: " << info.task.order << "\n\n"
           "Начальный размер шага: " << info.h << "\n\n";

    //Метод решения
    out << "\\section{Метод решения}\n\n"
           "Метод: " << solveMethodToString(info.method) << "\\\\\n"
           "Порядок точности: " << info.order << "\\\\\n"
           "Способ: " << info.way << "\\\\\n"
           "Метод итерации: ";
    switch (info.algo) {
        case IterationAlgo::EXPLICIT_STEP:
            out << "явный шаг\n\n";
            break;
        case IterationAlgo::SIMPLE_ITERATION:
            out << "метод простых итераций\n\n";
            break;
        case IterationAlgo::ZEIDEL:
            out << "метод Зейделя\n\n";
            break;
        case IterationAlgo::NEWTON:
            out << "метод Ньютона\n\n";
            break;
    }

    //Таблица Бутчера
    uint64_t n = info.butcher.size().n, m = info.butcher.size().m;
    out << "\\section{Таблица Бутчера}\n\n"
           "\\begin{table}[h]\n"
           "\\centering\n"
           "\\begin{tabular}{|c||";
    for (uint64_t i = 1; i < m; ++i) {
        out << "c|";
    }
    out << "}\n";
    for (uint64_t i = 0; i < n; ++i) {
        out << "\\hline\n" << info.butcher(i, 0);
        for (uint64_t j = 1; j < m; ++j) {
            out << (i >= m - 1 ? " & \\cellcolor{lightgray} " : " & ") << info.butcher(i, j);
        }
        out << "\\\\\n";
    }
    out << "\\hline\n"
           "\\end{tabular}\n"
           "\\end{table}\n\n";

    //Жёсткость
    out << "\\section{Жёсткость}\n\n"
           "Коэффициент жёсткости задачи: " << info.tough_coeff << "\\\\\n" << "Задача " << (info.tough_coeff >= 99.0 ? "" : "не ") << "жёсткая\n\n";

    //Решение задачи
    const auto &Yn = info.solution;
    auto Ya = getAnaliticSolution(info.solution[0], info.analitic);
    //auto &X = info.solution.first, &Yn = info.solution.second;
    //auto Ya = getAnaliticSolution(X, info.analitic);
    //bool use_analitic = !std::isnan(info.analitic(X0));
    std::vector<std::vector<double>> table;
    //table.push_back(Yn);
    //if (use_analitic) {
    //    table.push_back(Ya);
    //}
    out << "\\section{Решение задачи}\n\n"
           "\\begin{longtable}{||m{3cm}|";
    for (uint64_t i = 0; i < Ya.size(); ++i) {
        out << "|m{3cm}|" << "m{3cm}|";
    }
    out << "|";
    for (uint64_t i = Ya.size(); i < Yn.size() - 1; ++i) {
        out << "m{3cm}||";
    }
    out << "}\n"
           "\\hline\n"
           "\\cellcolor{lightgray} X";
    for (uint64_t i = 0; i < Ya.size(); ++i) {
        out << " & \\cellcolor{lightgray} $Y_{numeric" << i + 1 << "}$ & \\cellcolor{lightgray} $Y_{analitic" << i + 1 << "}$";
    }
    for (uint64_t i = Ya.size(); i < Yn.size() - 1; ++i) {
        out << " & \\cellcolor{lightgray} $Y_{numeric" << i + 1 << "}$";
    }
    out << "\\\\\n";
    for (uint64_t i = 0; i < Yn[0].size(); ++i) {
        out << "\\hline\n" << Yn[0][i];
        for (uint64_t j = 0; j < Ya.size(); ++j) {
            out << " & " << Yn[j + 1][i] << " & " << Ya[j][i];
        }
        for (uint64_t j = Ya.size(); j < Yn.size() - 1; ++j) {
            out << " & " << Yn[j + 1][i];
        }
        out << "\\\\\n";
    }
    out << "\\hline\n"
           "\\end{longtable}\n\n";
    // if (use_analitic) {
    //     double max_miss, average_miss;
    //     std::tie(max_miss, average_miss) = getAnaliticCompare(Yn, Ya);
    //     out << "Среднее отклонение от аналитического решения: " << average_miss << "\n\n";
    //     out << "Максимальное отклонение от аналитического решения: " << max_miss << "\n\n";
    // }
    for (uint64_t i = 0; i < Ya.size(); ++i) {
        double max_miss, average_miss;
        std::tie(max_miss, average_miss) = getAnaliticCompare(Yn[i + 1], Ya[i]);
        out << "Среднее отклонение от аналитического решения для функции " << i + 1 << ": " << average_miss << "\n\n"
               "Максимальное отклонение от аналитического решения для функции " << i + 1 << ": " << max_miss << "\n\n";
    }

    //Приближающий полином
    out << "\\section{Приближающий полином}\n\n";
    for (uint64_t i = 1; i < Yn.size(); ++i) {
        auto func = LeastSquareMethod(Yn[0], Yn[i], 3);
        out << "Приближающий полином 3й степени для функции " << i << ": $" << LSMToText(func) << "$\n\n";
    }
    //auto func = LeastSquareMethod(X, Yn, 3);
    //out << "Приближающий полином 3й степени: $" << LSMToText(func) << "$\n\n";

    //График
    std::vector<double> Xpoints;
    uint64_t step_count = 5;
    double h = (Xn - X0) / step_count;
    for (uint64_t i = 0; i < step_count; ++i) {
        Xpoints.push_back(X0 + i * h);
    }
    Xpoints.push_back(Xn);
    out << "\\section{График}\n\n"
           "\\begin{tikzpicture}\n"
           "\\begin{axis}[\n"
           "\txlabel={$x$},\n"
           "\tylabel={$y$},\n";
    out << "\txmin=" << X0 << ", xmax=" << Xn <<",\n";
    out << "\txtick={" << Xpoints[0];
    for (uint64_t i = 1; i < step_count + 1; ++i) {
        out << "," << Xpoints[i];
    }
    out << "},\n"
           "\tlegend pos=outer north east,\n"
           "\tymajorgrids=true,\n"
           "\tgrid style=dashed,\n"
           "]\n\n";
    //графики приближённого решения
    uint64_t count;
    if (info.multigraph) {
        count = Yn.size();
    } else {
        count = 2;
    }
    for (uint64_t i = 1; i < count; ++i) {
        out << "\\addplot[\n"
            "\tcolor=blue,\n"
            "\tmark=square,\n"
            "\tmark size=0.5pt\n"
            "]\n"
            "coordinates {\n";
        for (uint64_t j = 0; j < Yn[0].size(); ++j) {
            out << "(" << Yn[0][j] << "," << Yn[i][j] << ")";
        }
        out << "\n"
            "};\n"
            "\\addlegendentry{Приближённое решение " << i << "}\n\n";
    }
    //графики аналитического решения
    if (info.multigraph) {
        count = Ya.size();
    } else {
        count = std::max(Ya.size(), 0ul);
    }
    for (uint64_t i = 0; i < count; ++i) {
        out << "\\addplot[\n"
               "\tdomain=";
        out << X0 << ":" << Xn << ",\n"
               "\tcolor=red,\n"
               "\tsamples=100\n"
               "]{";
        out << toPgfplotsString(info.analitic[i].toString(FunctionalTree::Style::DEFAULT));
        out << "};\n"
               "\\addlegendentry{Аналитическое решение " << i + 1 << "}\n\n";
    }

    out <<  "\\end{axis}\n"
            "\\end{tikzpicture}\n\n";

    //конец документа
    out << "\\end{document}";
}

void generatePDF (const ReportInfo &info, std::ostream &out) {}

ReportInfo::ReportInfo () {
    // analitic = [] (double x) -> double {
    //     return x / 0;
    // };
    multigraph = false;
    //analitic.push_back(FunctionalTree("x / 0", std::vector<std::string>{"x"}));
}

ReportInfo::~ReportInfo () {}

void generateReport (const ReportInfo &info, ReportType type, std::ostream &out) {
    switch (type) {
        case ReportType::TXT:
            generateTEXT(info, out);
            break;
        case ReportType::TEX:
            generateTEX(info, out);
            break;
        case ReportType::PDF:
            generatePDF(info, out);
            break;
        default:
            break;
    }
}