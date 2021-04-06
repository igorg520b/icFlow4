#include <algorithm>
#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::~MainWindow() {delete ui;}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    worker = new BackgroundWorker(&modelController);

//    connect(&controller.model, SIGNAL(requestGeometryUpdate()), SLOT(render_results()));
    connect(&modelController, SIGNAL(stepAborted()),SLOT(updateGUI()));
    connect(&modelController, SIGNAL(stepCompleted()), SLOT(updateGUI()));
    connect(worker, SIGNAL(workerPaused()), SLOT(background_worker_paused()));

    // property browser
    pbrowser = new ObjectPropertyBrowser(this);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);


    renderer->SetBackground(colors->GetColor3d("White").GetData());
    renderer->AddActor(modelController.model.actor_selected_nodes);
    renderer->AddActor(modelController.model.actor_boundary);
    renderer->AddActor(modelController.model.actor_mesh);
//

    renderWindow->AddRenderer(renderer);
    renderWindow->GetInteractor()->SetInteractorStyle(specialSelector2D);
    specialSelector2D->mw = this;

    renderWindow->GetInteractor()->SetPicker(pointPicker);

    // right frame
    right_side_container = new QWidget;
    right_side_layout = new QHBoxLayout;
    right_side_container->setLayout(right_side_layout);
    right_side_layout->setSpacing(0);
    right_side_layout->setMargin(0);
    right_side_layout->addWidget(qt_vtk_widget);

    // splitter
    splitter_main = new QSplitter(Qt::Orientation::Horizontal);
    splitter_main->addWidget(pbrowser);
    splitter_main->addWidget(right_side_container);
    setCentralWidget(splitter_main);

    // toolbar - comboboxes
    comboBox_visualizations = new QComboBox();
    ui->toolBar->addWidget(comboBox_visualizations);

    // populate combobox
    QMetaEnum qme = QMetaEnum::fromType<icy::Model::VisOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_visualizations->addItem(qme.key(i));

    connect(comboBox_visualizations, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [=](int index){ comboboxIndexChanged_visualizations(index); });


    // slider
    labelStepCount = new QLabel();
//    ui->toolBar->addWidget(labelStepCount);

    // statusbar
    statusLabel = new QLabel("-");

    statusLabelAttempt = new QLabel(" --- ");

    QSizePolicy sp;
    sp.setHorizontalPolicy(QSizePolicy::Fixed);
    statusLabelAttempt->setSizePolicy(sp);
    statusLabelAttempt->setFixedWidth(100);

    labelStepCount->setSizePolicy(sp);
    labelStepCount->setFixedWidth(100);

    ui->statusbar->addWidget(statusLabel);
    ui->statusbar->addPermanentWidget(statusLabelAttempt);
    ui->statusbar->addPermanentWidget(labelStepCount);

    // read/restore saved settings
    QSettings settings(m_sSettingsFile);
    splitter_main->restoreState(settings.value("splitter_main").toByteArray());

    vtkCamera* camera = renderer->GetActiveCamera();
    camera->ParallelProjectionOn();

    QVariant var = settings.value("camData");
    if(!var.isNull()) {
        QByteArray arr = var.toByteArray();
        double *vec = (double*)arr.constData();
        camera->SetPosition(vec[0],vec[1],vec[2]);
        camera->SetFocalPoint(vec[3],vec[4],vec[5]);
        camera->SetViewUp(vec[6],vec[7],vec[8]);
        camera->SetViewAngle(vec[9]);
    }
    camera->Modified();

    renderer->AddActor(scalarBar);
    scalarBar->SetLookupTable(modelController.model.hueLut);

    scalarBar->SetMaximumWidthInPixels(130);
    scalarBar->SetBarRatio(0.07);
    scalarBar->SetMaximumHeightInPixels(400);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar->GetPositionCoordinate()->SetValue(0.90,0.015, 0.0);
    scalarBar->SetLabelFormat("%.1e");
    scalarBar->GetLabelTextProperty()->BoldOff();
    scalarBar->GetLabelTextProperty()->ItalicOff();
    scalarBar->GetLabelTextProperty()->ShadowOff();
    scalarBar->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);
//    scalarBar->GetLabelTextProperty()->SetFontFamilyToTimes();

    renderWindow->Render();

    comboBox_visualizations->setCurrentIndex(settings.value("vis_option").toInt());
}

void MainWindow::showEvent( QShowEvent*)
{
    pbrowser->setActiveObject(&modelController.prms);
    QSettings settings(m_sSettingsFile);
//    comboBox_visualizations->setCurrentIndex(settings.value("vis_option").toInt());
    updateGUI();
    renderWindow->Render();
}


void MainWindow::on_action_quit_triggered() { this->close(); }

void MainWindow::closeEvent( QCloseEvent* event )
{
    // save settings and stop simulation
    QSettings settings(m_sSettingsFile);
    qDebug() << "MainWindow: closing main window; " << settings.fileName();

    settings.setValue("splitter_main", splitter_main->saveState());

    double data[10];
    renderer->GetActiveCamera()->GetPosition(&data[0]);
    renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    renderer->GetActiveCamera()->GetViewUp(&data[6]);
    data[9]=renderer->GetActiveCamera()->GetViewAngle();

    QByteArray arr((char*)&data[0], sizeof(double)*10);
    settings.setValue("camData", arr);
    settings.setValue("vis_option", comboBox_visualizations->currentIndex());

    // kill backgroundworker
    worker->Finalize();
    event->accept();
}

void MainWindow::on_action_simulation_start_triggered(bool checked)
{
    if(!worker->running && checked){
        qDebug() << "start button - starting";
        statusLabel->setText("starting simulation");
        modelController.Prepare();

        worker->Resume();
    }
    else if(worker->running && !checked)
    {
        statusLabel->setText("pausing simulation");
        qDebug() << "start button - pausing";
        worker->Pause();
        ui->action_simulation_start->setEnabled(false);
    }

}

void MainWindow::background_worker_paused()
{
    // enable the "Start" button
    qDebug() << "MainWindow::background_worker_paused()";
    ui->action_simulation_start->blockSignals(true);
    ui->action_simulation_start->setEnabled(true);
    ui->action_simulation_start->setChecked(false);
    ui->action_simulation_start->blockSignals(false);
    statusLabel->setText("stopped");
}

void MainWindow::render_results()
{
    modelController.model.UnsafeUpdateGeometry();
    renderWindow->Render();
}


void MainWindow::on_action_simulation_single_step_triggered()
{
    qDebug() << "take one step";
    modelController.Prepare();
    modelController.Step();
}

void MainWindow::on_action_camera_reset_triggered()
{
    qDebug() << "MainWindow::on_action_camera_reset_triggered()";
    vtkCamera* camera = renderer->GetActiveCamera();
    camera->ParallelProjectionOn();
    camera->SetClippingRange(1e1,1e3);
    camera->SetFocalPoint(0.0, 0.0, 0.0);
    camera->SetPosition(0.0, 0.0, 50.0);
    camera->SetViewUp(0.0, 1.0, 0.0);
    camera->Modified();
    renderWindow->Render();
}



void MainWindow::updateGUI()
{
    bool r = worker->running;
    ui->action_simulation_single_step->setEnabled(!r);

    labelStepCount->setText(QString{"step: %1"}.arg(modelController.currentStep));

    statusLabelAttempt->setText(QString{"%1"}.arg(modelController.timeStepFactor, 6, 'f', 3, '0'));

    render_results();

}



void MainWindow::comboboxIndexChanged_visualizations(int index)
{
    qDebug() << "comboboxIndexChanged_visualizations " << index;
    modelController.model.ChangeVisualizationOption((icy::Model::VisOpt)index);
    scalarBar->SetVisibility(index != 0);
    renderWindow->Render();
}








