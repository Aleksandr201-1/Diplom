#include <GUI/MainWindow.hpp>
#include <GUI/ui_MainWindow.h>

struct SolveStruct {
    SolveMethod method;
    uint64_t order, way;
};

const std::map<QString, SolveStruct> strToExplicitSolve = {
    {"Рунге-Кутта 3 порядка, 1й способ", {SolveMethod::RUNGE_KUTTA, 3, 1}},
    {"Рунге-Кутта 3 порядка, 2й способ", {SolveMethod::RUNGE_KUTTA, 3, 2}},
    {"Рунге-Кутта 3 порядка, 3й способ", {SolveMethod::RUNGE_KUTTA, 3, 3}},
    {"Рунге-Кутта 3 порядка, 4й способ", {SolveMethod::RUNGE_KUTTA, 3, 4}},
    {"Рунге-Кутта 4 порядка, 1й способ (классический)", {SolveMethod::RUNGE_KUTTA, 4, 1}},
    {"Рунге-Кутта 4 порядка, 2й способ", {SolveMethod::RUNGE_KUTTA, 4, 2}},
    {"Рунге-Кутта 4 порядка, 3й способ", {SolveMethod::RUNGE_KUTTA, 4, 3}},
    {"Рунге-Кутта 4 порядка, 4й способ", {SolveMethod::RUNGE_KUTTA, 4, 4}},
    {"Рунге-Кутта 6 порядка, 1й способ", {SolveMethod::RUNGE_KUTTA, 6, 1}},
    {"Рунге-Кутта 6 порядка, 2й способ", {SolveMethod::RUNGE_KUTTA, 6, 2}}
};
const std::map<QString, SolveStruct> strToEmbeededSolve = {
    {"Ческино 2(3) порядка, 1й способ", {SolveMethod::CHESKINO, 2, 1}},
    {"Фалберг 2(3) порядка, 1й способ", {SolveMethod::FALBERG, 2, 1}},
    {"Фалберг 2(3) порядка, 2й способ", {SolveMethod::FALBERG, 2, 2}},
    {"Мерсон 3(4) порядка, 1й способ", {SolveMethod::MERSON, 3, 1}},
    {"Дорман-Принс 4(5) порядка, 1й способ", {SolveMethod::DORMAN_PRINCE, 4, 1}}
};
const std::map<QString, SolveStruct> strToImplicitSolve = {
    {"Лобатто 2 порядка, 1й способ", {SolveMethod::LOBATTO, 2, 1}},
    {"Лобатто 2 порядка, 2й способ", {SolveMethod::LOBATTO, 2, 2}},
    {"Лобатто 2 порядка, 3й способ", {SolveMethod::LOBATTO, 2, 3}},
    {"Лобатто 2 порядка, 4й способ", {SolveMethod::LOBATTO, 2, 4}},
    {"Лобатто 4 порядка, 1й способ", {SolveMethod::LOBATTO, 4, 1}},
    {"Лобатто 4 порядка, 2й способ", {SolveMethod::LOBATTO, 4, 2}},
    {"Лобатто 4 порядка, 3й способ", {SolveMethod::LOBATTO, 4, 3}},
    {"Лобатто 4 порядка, 4й способ", {SolveMethod::LOBATTO, 4, 4}},
    {"Лобатто 6 порядка, 1й способ", {SolveMethod::LOBATTO, 6, 1}},
    {"Лобатто 6 порядка, 2й способ", {SolveMethod::LOBATTO, 6, 2}},
    {"Лобатто 6 порядка, 3й способ", {SolveMethod::LOBATTO, 6, 3}},
    {"Лобатто 6 порядка, 4й способ", {SolveMethod::LOBATTO, 6, 4}},
    {"Гаусс 2 порядка, 1й способ", {SolveMethod::GAUSS, 2, 1}},
    {"Гаусс 4 порядка, 1й способ", {SolveMethod::GAUSS, 4, 1}},
    {"Гаусс 6 порядка, 1й способ", {SolveMethod::GAUSS, 6, 1}},
    {"Л-стабильный 3 порядка, 1й способ", {SolveMethod::L_STABLE_DIAGONAL, 3, 1}},
    {"Л-стабильный 4 порядка, 1й способ", {SolveMethod::L_STABLE_DIAGONAL, 4, 1}},
    {"Радо 1 порядка, 1й способ", {SolveMethod::RADO, 1, 1}},
    {"Радо 3 порядка, 1й способ", {SolveMethod::RADO, 3, 1}},
    {"Радо 5 порядка, 1й способ", {SolveMethod::RADO, 5, 1}}
};

const std::map<QString, IterationAlgo> strToIteration = {
    {"Метод Ньютона", IterationAlgo::NEWTON},
    {"Метод Зейделя", IterationAlgo::ZEIDEL},
    {"Метод простой итерации", IterationAlgo::SIMPLE_ITERATION}
    //{"Явный шаг", IterationAlgo::EXPLICIT_STEP}
};

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
    chemSys.initFromFile("./test/ChemicTest/bufermm.txt");
    chemSys.addAdditive("M", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
    //Y0.resize(chemSys.getSubstanceList().size() , 0.0);
    Y0.resize(toUnderlying(TaskType::ERROR));
    Y0[toUnderlying(TaskType::KOSHI)].resize(2, 0.0);
    Y0[toUnderlying(TaskType::CHEMICAL)].resize(chemSys.getSubstanceList().size() , 0.0);
    
    ui->setupUi(this);
    ui->stackedWidget->setCurrentWidget(ui->taskPage);
    QGraphicsScene *scene = new QGraphicsScene(ui->graphicsView);

    QPen pen(Qt::green);//Просто выбираем цвет для карандашика
    //scene->addLine(0,90,180,90,pen);//x
    //scene->addLine(90,0,90,180,pen);//y
    ui->graphicsView->setScene(scene);
    ui->graphicsView->scale(5.0,5.0);

    ui->solutionMethodBox->clear();
    for (auto &el : strToExplicitSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    ui->iterMethodBox->addItem("Явный шаг");

    ui->tableWidget->setColumnCount(7); // Указываем число колонок
    ui->tableWidget->setShowGrid(true); // Включаем сетку
    // Разрешаем выделение только одного элемента
    ui->tableWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    // Разрешаем выделение построчно
    ui->tableWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    // Устанавливаем заголовки колонок
    QStringList horzHeaders;
    horzHeaders << "Метод" << "Итер. алг." << "Время (мс)" << "Шаги" << "Сред. погрешность" << "Макс. погрешность" << "Используемая точность";
    ui->tableWidget->setHorizontalHeaderLabels(horzHeaders);
    // Растягиваем последнюю колонку на всё доступное пространство
    ui->tableWidget->horizontalHeader()->setStretchLastSection(true);
    //header->setVisible(true)
    ui->tableWidget->horizontalHeader()->setVisible(true);

    ui->chemTableWidget->setColumnCount(2);
    ui->chemTableWidget->setShowGrid(true);
    ui->chemTableWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    ui->chemTableWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    horzHeaders.clear();
    horzHeaders << "Название" << "Атомная масса (кг/моль)";
    ui->chemTableWidget->setHorizontalHeaderLabels(horzHeaders);
    ui->chemTableWidget->horizontalHeader()->setStretchLastSection(true);
    ui->chemTableWidget->horizontalHeader()->setVisible(true);
    auto substancesList = chemSys.getSubstanceList();
    auto substanceMasses = chemSys.getSubstanceMasses();
    for (uint64_t i = 0; i < substancesList.size(); ++i) {
        ui->chemTableWidget->insertRow(i);
        ui->chemTableWidget->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(substancesList[i])));
        ui->chemTableWidget->setItem(i, 1, new QTableWidgetItem(QString::fromStdString(std::to_string(substanceMasses[substancesList[i]]))));
        ui->substanceConcBox->addItem(QString::fromStdString(substancesList[i]));
    }
    ui->reactionTypeBox->addItems({"Изотерма, rho = const", "Адиабата, rho = const", "Изотерма, P = const", "Адиабата, P = const"});

    ui->reactionListWidget->addItem("1) H2 + O2 <==> 2OH\nA = 17000000.000000\nn = 0.000000\nE = 24044.000000\n");
    //ui->reactionListWidget->addItem("2) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
    //ui->reactionListWidget->addItem("3) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
    //ui->reactionListWidget->addItem("4) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
    //ui->reactionListWidget->addItem("5) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
    //ui->reactionListWidget->addItem("6) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
    //ui->reactionListWidget->addItem("7) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
    //ui->reactionListWidget->addItem("8) H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");

    std::string reactionStr = ui->reactInputLineEdit->text().toStdString() + "<==>" + ui->reactOutputLineEdit->text().toStdString();
    float128_t A = ui->ASpinBox->value() , n = ui->nSpinBox->value(), E = ui->ESpinBox->value();
    chemSys.addReaction(reactionStr, A, n, E);

    //ui->reactionListWidget->item(0)->setText(QString::fromStdString("qefaef\n"));
    //ui->reactionListWidget->item(0)->setData(9);
    //ui->reactionListWidget->item(0)->setData(ui->reactionListWidget->item(0)->text().replace("1)", "2)"));
    //ui->reactionListWidget->item(0)->text().replace("1)", "2)");
}

MainWindow::~MainWindow() {
    delete ui;
}


void MainWindow::on_solutionChoiceButton_clicked() {
    ui->stackedWidget->setCurrentWidget(ui->solutionPage);
}

void MainWindow::on_backToTaskInputButton_clicked() {
    switch(taskType) {
        case TaskType::KOSHI:
            ui->stackedWidget->setCurrentWidget(ui->KoshiPage);
            break;
        case TaskType::KOSHI_SYSTEM:
            break;
        case TaskType::CHEMICAL:
            ui->stackedWidget->setCurrentWidget(ui->ChemicPage);
            break;
        default:
            break;
    }
}

void MainWindow::on_reportButton_clicked() {

}

//Рассчёт
void MainWindow::on_calcButton_clicked() {
    float128_t X0, Xn, h, approx;
    uint64_t order, way, time = 0;
    Matrix<float128_t> butcher;
    std::vector<std::string> inputLines;
    std::vector<float128_t> &Ykoshi = Y0[toUnderlying(TaskType::KOSHI)];
    SolveMethod method;
    IterationAlgo iteration;
    SolveStruct solveStruct;
    std::string analiticLine;
    std::chrono::time_point <std::chrono::system_clock> startT, endT;
    
    // std::cout << "Its calculating time\n";
    // std::cout << "task: " << ui->inputLineEdit->text().toStdString() << "\n";
    // std::cout << "analitic: " << ui->analiticLineEdit->text().toStdString() << "\n";
    // std::cout << "borders: " << ui->borderASpinBox->value() << " " << ui->borderBSpinBox->value() << "\n";
    // std::cout << "step: " << ui->hStepBox->value() << "\n";
    // std::cout << "first init: " << ui->stepBoxY0->value() << " " << ui->stepBoxY1->value() << "\n";
    // std::cout << "solve method: " << ui->solutionMethodBox->currentText().toStdString() << "\n";
    // std::cout << "iteration method: " << ui->iterMethodBox->currentText().toStdString() << "\n";

    iteration = strToIteration.find(ui->iterMethodBox->currentText())->second;
    X0 = ui->borderASpinBox->value();
    Xn = ui->borderBSpinBox->value();
    h = ui->hStepBox->value();
    approx = ui->approxSpinBox->value();
    analiticLine = ui->analiticLineEdit->text().toStdString();
    Ykoshi = {ui->stepBoxY0->value(), ui->stepBoxY1->value()};

    if (ui->explicitButton->isChecked()) {
        //std::cout << "explicit\n";
        solveStruct = strToExplicitSolve.find(ui->solutionMethodBox->currentText())->second;
    } else if (ui->embeededButton->isChecked()) {
        //std::cout << "embeeded\n";
        solveStruct = strToEmbeededSolve.find(ui->solutionMethodBox->currentText())->second;
    } else {
        //std::cout << "implicit\n";
        solveStruct = strToImplicitSolve.find(ui->solutionMethodBox->currentText())->second;
    }
    method = solveStruct.method;
    order = solveStruct.order;
    way = solveStruct.way;
    butcher = createButcherTable(method, order, way);

    std::cout << "getting task\n";
    koshi.setTaskInfo(ui->inputLineEdit->text().toStdString(), getOrder(ui->orderLineEdit->text().toStdString()), Ykoshi, X0, Xn);
    std::cout << "solving\n";
    startT = std::chrono::system_clock::now();
    auto solution = KoshiSolver(method, koshi, butcher, h, iteration, approx);
    endT = std::chrono::system_clock::now();
    chemSys.printInfo(std::cout);
    std::cout << "=====Конец расчёта=====\n";
    time += std::chrono::duration_cast<duration_t>(endT - startT).count();
    std::cout << "Time: " << time << "\n";

    QGraphicsScene *scene = ui->graphicsView->scene();
    scene->clear();
    for(uint64_t i = 0; i < solution[0].size(); ++i) {
        scene->addEllipse(solution[0][i], solution[1][i], 1, 1);
    }
    uint64_t idx = ui->tableWidget->rowCount();
    ui->tableWidget->insertRow(idx);
    //ui->tableWidget->setItem(idx, 0, new QTableWidgetItem(QString::number(idx + 1)));
    ui->tableWidget->setItem(idx, 0, new QTableWidgetItem(ui->solutionMethodBox->currentText()));
    ui->tableWidget->setItem(idx, 1, new QTableWidgetItem(ui->iterMethodBox->currentText()));
    ui->tableWidget->setItem(idx, 2, new QTableWidgetItem(QString::number(time)));
    ui->tableWidget->setItem(idx, 3, new QTableWidgetItem(QString::number(solution[0].size())));
    if (taskType == TaskType::KOSHI && analiticLine != "") {
        FunctionalTree func(analiticLine, std::vector<std::string>{"x"});
        auto analitic = getAnaliticSolution(solution[0], {func});
        float128_t max_miss, average_miss;
        std::tie(max_miss, average_miss) = getAnaliticCompare(solution[1], analitic[0]);
        ui->tableWidget->setItem(idx, 4, new QTableWidgetItem(QString::fromStdString(std::to_string(average_miss))));
        ui->tableWidget->setItem(idx, 5, new QTableWidgetItem(QString::fromStdString(std::to_string(max_miss))));
    } else {
        ui->tableWidget->setItem(idx, 4, new QTableWidgetItem("-"));
        ui->tableWidget->setItem(idx, 5, new QTableWidgetItem("-"));
    }
    ui->tableWidget->setItem(idx, 6, new QTableWidgetItem(QString::fromStdString(std::to_string(approx))));

    std::cout << solution.back().back() << "\n";
    //std::cout << koshi. << "\n";
    //таблица с аналитическим решением
    //метод, время (мс), количестов шагов, средняя погрешность, максимальная погрешность, используемая точность
    //таблица без аналитического решения
    //метод, время (мс), количестов шагов, используемая точность
}

void MainWindow::on_toKoshiTask_clicked() {
    taskType = TaskType::KOSHI;
    //ui->tableWidget->clear();
    ui->tableWidget->setRowCount(0);
    ui->graphicsView->scene()->clear();
    ui->stackedWidget->setCurrentWidget(ui->KoshiPage);
}

void MainWindow::on_toChemicalTask_clicked() {
    taskType = TaskType::CHEMICAL;
    //ui->tableWidget->clear();
    ui->tableWidget->setRowCount(0);
    ui->graphicsView->scene()->clear();
    ui->stackedWidget->setCurrentWidget(ui->ChemicPage);
}

void MainWindow::on_aboutProgramm_clicked() {

}

void MainWindow::on_explicitButton_clicked() {
    ui->solutionMethodBox->clear();
    ui->iterMethodBox->clear();
    for (auto &el : strToExplicitSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    ui->iterMethodBox->addItem("Явный шаг");
}

void MainWindow::on_embeededButton_clicked() {
    ui->solutionMethodBox->clear();
    ui->iterMethodBox->clear();
    for (auto &el : strToEmbeededSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    ui->iterMethodBox->addItem("Явный шаг");
}

void MainWindow::on_implicitButton_clicked() {
    ui->solutionMethodBox->clear();
    ui->iterMethodBox->clear();
    for (auto &el : strToImplicitSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    for (auto &el : strToIteration) {
        ui->iterMethodBox->addItem(el.first);
    }
}

void MainWindow::on_KoshiToMenuButton_clicked() {
    ui->stackedWidget->setCurrentWidget(ui->taskPage);
}

void MainWindow::on_solutionChoiceButton2_clicked() {
    ui->stackedWidget->setCurrentWidget(ui->solutionPage);
}

void MainWindow::on_ChemicToMenuButton_clicked() {
    ui->stackedWidget->setCurrentWidget(ui->taskPage);
}

void MainWindow::on_reactionNumSpinBox_valueChanged(int arg1) {
    if (arg1 > chemSys.getCount() + 1) {
        ui->reactionNumSpinBox->setValue(chemSys.getCount() + 1);
        QMessageBox messageBox;
        messageBox.critical(0,"Ошибка","Некорректный номер реакции!");
        messageBox.setFixedSize(500,200);
        return;
    }
}

void MainWindow::on_addReactionButton_clicked() {
    uint64_t idx = ui->reactionNumSpinBox->value();
    std::string reactionStr = ui->reactInputLineEdit->text().toStdString() + " <==> " + ui->reactOutputLineEdit->text().toStdString();
    float128_t A = ui->ASpinBox->value() , n = ui->nSpinBox->value(), E = ui->ESpinBox->value();
    QString str = QString::fromStdString(std::to_string(idx) + ") " + reactionStr + "\nA = " + std::to_string(A) + "\nn = " + std::to_string(n) + "\nE = " + std::to_string(E) + "\n");
    try {
        if (idx == chemSys.getCount() + 1) {
            //ui->reactionListWidget->addItem(") H2 + O2 <==> 2OH\nA = 1.700e+10; n = 0.0; E = 24044.0\n");
            chemSys.addReaction(reactionStr, A, n, E);
            ui->reactionListWidget->addItem(str);
            ui->reactionNumSpinBox->setValue(idx + 1);
        } else {
            chemSys.changeReaction(idx - 1, reactionStr, A, n, E);
            ui->reactionListWidget->item(idx - 1)->setText(str);
        }
    } catch (std::exception &exp) {
        QMessageBox messageBox;
        messageBox.critical(0,"Ошибка","Некорректный ввод реакции!");
        messageBox.setFixedSize(500,200);
        return;
    }
    //chemSys.addReaction(reactionStr, A, n, E);
    //ChemicalReaction reaction(reactionStr, A, n, E);
    // if (idx > reactions.size()) {
    //     reactions.resize(idx);
    // }
}

void MainWindow::on_deleteReactionButton_clicked() {
    uint64_t idx = ui->reactionNumSpinBox->value();
    if (idx == chemSys.getCount() + 1) {
        QMessageBox messageBox;
        messageBox.critical(0,"Ошибка","Некорректный номер реакции!");
        messageBox.setFixedSize(500,200);
        return;
    }
    chemSys.deleteReaction(idx - 1);
    ui->reactionListWidget->removeItemWidget(ui->reactionListWidget->takeItem(idx - 1));

    for (uint64_t i = idx - 1; i < ui->reactionListWidget->count(); ++i) {
        QString str = ui->reactionListWidget->item(i)->text();
        str.remove(0, str.indexOf(")"));
        str = QString::fromStdString(std::to_string(i + 1)) + str;
        ui->reactionListWidget->item(i)->setText(str);
    }
    if (idx != ui->reactionNumSpinBox->minimum() && idx == chemSys.getCount() + 1) {
        ui->reactionNumSpinBox->setValue(idx - 1);
    }
}

void MainWindow::on_substanceConcBox_currentIndexChanged(int index) {
    //std::cout << "index " << index << "\n" << "conc size " << Y0.size() << "\n";
    ui->concSpinBox->setValue(Y0[toUnderlying(TaskType::CHEMICAL)][index]);
    //Y0[index] = ui->concSpinBox->value();
}

void MainWindow::on_concSpinBox_valueChanged(double arg1) {
    uint64_t index = ui->substanceConcBox->currentIndex();
    Y0[toUnderlying(TaskType::CHEMICAL)][index] = arg1;
}

void MainWindow::on_createBucherTable_clicked() {
    ui->stackedWidget->setCurrentWidget(ui->ButcherTablePage);
}

void MainWindow::on_ButcherToMenuButton_clicked() {
    ui->stackedWidget->setCurrentWidget(ui->taskPage);
}