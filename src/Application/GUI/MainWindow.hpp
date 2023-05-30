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
#include <General/ButcherTable.hpp>
#include <iostream>
#include <map>

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

    void on_backToKoshiButton_clicked();

    void on_reportButton_clicked();

    void on_calcButton_clicked();

    void on_toKoshiTask_clicked();

    void on_toChemicalTask_clicked();

    void on_aboutProgramm_clicked();

    void on_explicitButton_clicked();

    void on_embeededButton_clicked();

    void on_implicitButton_clicked();

private:
    QGraphicsEffect *scene;
    std::vector<QString> explicitMethods, implicitMethods, embeededMethods, iterationMethods;
    Ui::MainWindow *ui;
};

#endif