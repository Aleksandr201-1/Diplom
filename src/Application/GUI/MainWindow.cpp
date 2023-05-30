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
    , ui(new Ui::MainWindow)
{
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
    ui->tableWidget->setHorizontalHeaderLabels({"Метод", "Время (мс)", "Шаги", "Сред. погрешность",
                                                "Макс. погрешность", "Используемая точность"});
    // Растягиваем последнюю колонку на всё доступное пространство
    ui->tableWidget->horizontalHeader()->setStretchLastSection(true);
    // Скрываем колонку под номером 0
    ui->tableWidget->hideColumn(0);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_solutionChoiceButton_clicked()
{
    ui->stackedWidget->setCurrentWidget(ui->solutionPage);
}

void MainWindow::on_backToKoshiButton_clicked()
{
    ui->stackedWidget->setCurrentWidget(ui->KoshiPage);
}

void MainWindow::on_reportButton_clicked()
{

}

//Рассчёт
void MainWindow::on_calcButton_clicked()
{
    float128_t X0, Xn, h, approx;
    uint64_t order, way;
    Matrix<float128_t> butcher;
    std::vector<std::string> inputLines;
    std::vector<float128_t> Y0;
    KoshiTask koshi;
    SolveMethod method;
    IterationAlgo iteration;
    SolveStruct solveStruct;
    std::cout << "Its calculating time\n";

    std::cout << "task: " << ui->inputLineEdit->text().toStdString() << "\n";
    std::cout << "analitic: " << ui->analiticLineEdit->text().toStdString() << "\n";
    std::cout << "borders: " << ui->borderASpinBox->value() << " " << ui->borderBSpinBox->value() << "\n";
    std::cout << "step: " << ui->hStepBox->value() << "\n";
    std::cout << "first init: " << ui->stepBoxY0->value() << " " << ui->stepBoxY1->value() << "\n";
    std::cout << "solve method: " << ui->solutionMethodBox->currentText().toStdString() << "\n";
    std::cout << "iteration method: " << ui->iterMethodBox->currentText().toStdString() << "\n";

    iteration = strToIteration.find(ui->iterMethodBox->currentText())->second;
    X0 = ui->borderASpinBox->value();
    Xn = ui->borderBSpinBox->value();
    h = ui->hStepBox->value();
    approx = ui->approxSpinBox->value();
    Y0 = {ui->stepBoxY0->value(), ui->stepBoxY1->value()};
    if (ui->explicitButton->isChecked()) {
        std::cout << "explicit\n";
        solveStruct = strToExplicitSolve.find(ui->solutionMethodBox->currentText())->second;
        QMessageBox messageBox;
        messageBox.critical(0,"Ошибка","Некорректный ввод метода!");
        messageBox.setFixedSize(500,200);
        return;
    } else if (ui->embeededButton->isChecked()) {
        std::cout << "embeeded\n";
        solveStruct = strToEmbeededSolve.find(ui->solutionMethodBox->currentText())->second;
    } else {
        std::cout << "implicit\n";
        solveStruct = strToImplicitSolve.find(ui->solutionMethodBox->currentText())->second;
    }
    method = solveStruct.method;
    order = solveStruct.order;
    way = solveStruct.way;
    butcher = createButcherTable(method, order, way);

    std::cout << "getting task\n";
    koshi.setTaskInfo(ui->inputLineEdit->text().toStdString(), getOrder(ui->orderLineEdit->text().toStdString()), Y0, X0, Xn);
    std::cout << "solving\n";
    auto solution = KoshiSolver(method, koshi, butcher, h, iteration, approx);

    QGraphicsScene *scene = ui->graphicsView->scene();
    scene->clear();
    for(uint64_t i = 0; i < solution[0].size(); ++i) {
        scene->addEllipse(solution[0][i], solution[1][i], 1, 1);
    }

    std::cout << solution.back().back() << "\n";
    //std::cout << koshi. << "\n";
    //таблица с аналитическим решением
    //метод, время (мс), количестов шагов, средняя погрешность, максимальная погрешность, используемая точность
    //таблица без аналитического решения
    //метод, время (мс), количестов шагов, используемая точность
}

void MainWindow::on_toKoshiTask_clicked()
{
    ui->stackedWidget->setCurrentWidget(ui->KoshiPage);
    //ui->stackedWidget->setCurrentIndex(1);
}

void MainWindow::on_toChemicalTask_clicked()
{

}

void MainWindow::on_aboutProgramm_clicked()
{

}

void MainWindow::on_explicitButton_clicked()
{
    ui->solutionMethodBox->clear();
    ui->iterMethodBox->clear();
    for (auto &el : strToExplicitSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    ui->iterMethodBox->addItem("Явный шаг");
}

void MainWindow::on_embeededButton_clicked()
{
    ui->solutionMethodBox->clear();
    ui->iterMethodBox->clear();
    for (auto &el : strToEmbeededSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    ui->iterMethodBox->addItem("Явный шаг");
}

void MainWindow::on_implicitButton_clicked()
{
    ui->solutionMethodBox->clear();
    ui->iterMethodBox->clear();
    for (auto &el : strToImplicitSolve) {
        ui->solutionMethodBox->addItem(el.first);
    }
    for (auto &el : strToIteration) {
        ui->iterMethodBox->addItem(el.first);
    }
}