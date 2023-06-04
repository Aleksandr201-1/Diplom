#ifndef MAIN_WINDOW_HPP
#define MAIN_WINDOW_HPP

#include <QMainWindow>
#include <QMessageBox>
//#include <qcustomplot.h>
//#include <GUI/QCustomPlot.hpp>
//#include <qwt_plot_curve.h>
//#include <qwt_plot.h>
//#include <qwt_plot_grid.h>
#include <ODUSolver/Koshi/KoshiSolver.hpp>
#include <ODUSolver/Chemical/ChemicalSolver.hpp>
#include <PDFReporter/ReportGenerator.hpp>
#include <General/ButcherTable.hpp>
#include <iostream>
#include <map>
#include <chrono>

using duration_t = std::chrono::milliseconds;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_solutionChoiceButton_clicked();

    void on_backToTaskInputButton_clicked();

    void on_reportButton_clicked();

    void on_calcButton_clicked();

    void on_toKoshiTask_clicked();

    void on_toChemicalTask_clicked();

    void on_aboutProgramm_clicked();

    void on_explicitButton_clicked();

    void on_embeededButton_clicked();

    void on_implicitButton_clicked();

    void on_KoshiToMenuButton_clicked();

    void on_solutionChoiceButton2_clicked();

    void on_ChemicToMenuButton_clicked();

    void on_reactionNumSpinBox_valueChanged(int arg1);

    void on_addReactionButton_clicked();

    void on_deleteReactionButton_clicked();

    void on_substanceConcBox_currentIndexChanged(int index);

    void on_concSpinBox_valueChanged(double arg1);

    void on_createBucherTable_clicked();

    void on_ButcherToMenuButton_clicked();

private:
    TaskType taskType;
    ChemicalSystem chemSys;
    KoshiTask koshi;
    //std::vector<ChemicalReaction> reactions;
    std::vector<std::vector<float128_t>> Y0;
    QGraphicsEffect *scene;
    std::vector<QString> explicitMethods, implicitMethods, embeededMethods, iterationMethods;
    Ui::MainWindow *ui;
};

#endif