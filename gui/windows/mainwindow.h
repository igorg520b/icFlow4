#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QFileDialog>
#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>
#include <QVBoxLayout>
#include <QTreeWidget>
#include <QProgressBar>
#include <QMenu>
#include <QList>
#include <QDebug>
#include <QComboBox>
#include <QMetaEnum>
#include <QDir>

#include <QtCharts>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>

#include <vtkOBJReader.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkVersion.h>
#include <vtkWindowToImageFilter.h>
//#include <vtkPointSource.h>
//#include <vtkLineSource.h>
//#include <vtkOBBTree.h>
#include <vtkPolyDataMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAxesActor.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include <vtkAreaPicker.h>
#include <vtkPointPicker.h>
#include <vtkProp3DCollection.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCallbackCommand.h>
#include <vtkInteractorStyleRubberBandPick.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdint>

#include "objectpropertybrowser.h"
//#include "preferences_gui.h"

//#include "modelcontroller.h"
#include "modelcontrollerinterface.h"
#include "backgroundworker.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::MainWindow *ui;

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void showEvent( QShowEvent* event ) override;
    void closeEvent( QCloseEvent* event ) override;

    static void PickCallbackFunction(vtkObject* caller,
                              long unsigned int vtkNotUsed(eventId),
                              void* vtkNotUsed(clientData),
                              void* vtkNotUsed(callData));


private slots:

    void on_action_quit_triggered();

    void comboboxIndexChanged_visualizations(int index);

    void progress_updated();
    void render_results();
    void background_worker_paused();


    void on_action_simulation_single_step_triggered();

    void on_action_camera_reset_triggered();

    void on_action_simulation_start_triggered(bool checked);

    void updateGUI();   // when simulation is started/stopped or when a step is advanced



private:
//    PreferencesGUI prefsGUI;
//    icy::ModelController controller;   // simulation algorithms
    ModelControllerTest testController;
    BackgroundWorker *worker;

    QString m_sSettingsFile = "ic4_config";
    QLabel *statusLabel;                    // statusbar
    QLabel *statusLabelAttempt;             // attempt # at advancing a time step
    QLabel *labelStepCount;
    QComboBox *comboBox_visualizations;

    // splitter and the right side window
    QSplitter *splitter_main;
    QHBoxLayout *right_side_layout;
    QWidget *right_side_container;

    // properties
    ObjectPropertyBrowser *pbrowser;

    // VTK
    QVTKOpenGLNativeWidget *qt_vtk_widget;

    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkPointPicker> pointPicker;

    const QString wtitle = "icFlow4: Finite Element Simulation";

};
#endif // MAINWINDOW_H
