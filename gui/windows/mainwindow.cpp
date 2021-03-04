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
    renderer->AddActor(modelController.model.actor_mesh);

    renderWindow->AddRenderer(renderer);
//    renderWindow->GetInteractor()->SetInteractorStyle(interactorStyleRB2D);
    renderWindow->GetInteractor()->SetInteractorStyle(specialSelector2D);

    renderWindow->GetInteractor()->SetPicker(pointPicker);
/*    vtkSmartPointer<vtkCallbackCommand> pickCallback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    pickCallback->SetCallback(MainWindow::PickCallbackFunction);
    pointPicker->AddObserver(vtkCommand::EndPickEvent, pickCallback);
    pickCallback->SetClientData((void*)this);
*/

    vtkSmartPointer<vtkCallbackCommand> selectionChangedCallback =
      vtkSmartPointer<vtkCallbackCommand>::New();
    selectionChangedCallback->SetCallback (MainWindow::SelectionChangedCallbackFunction);
    interactorStyleRB2D->AddObserver ( vtkCommand::SelectionChangedEvent, selectionChangedCallback );


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
    statusLabelAttempt->setFixedWidth(70);

    labelStepCount->setSizePolicy(sp);
    labelStepCount->setFixedWidth(70);

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
    renderWindow->Render();

//    prefsGUI.LoadState(settings);
//    ui->action_Load_Recent_on_Startup->setChecked(prefsGUI.LoadLastScene);

//    comboBox_visualizations->setCurrentIndex(prefsGUI.VisualizationOption);
}

void MainWindow::showEvent( QShowEvent*)
{
    pbrowser->setActiveObject(&modelController.prms);
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

//    prefsGUI.LoadLastScene = ui->action_Load_Recent_on_Startup->isChecked();
//    prefsGUI.SaveState(settings);

    // kill backgroundworker
    worker->Finalize();
    event->accept();
}

void MainWindow::on_action_simulation_start_triggered(bool checked)
{
    if(!worker->running && checked){
        qDebug() << "start button - starting";
        statusLabel->setText("starting simulation");
        // controller.Prepare();
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
//    controller.model.UnsafeUpdateGeometry(controller.ts.SimulationTime, controller.prms);
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

//    labelStepCount->setText(QString::number(modelController.currentStep));
    labelStepCount->setText(QString{"step: %1"}.arg(modelController.currentStep));

    //    statusLabel->setText(r ? "running" : "paused");
//    if(!r) ui->action_simulation_start->setEnabled(true);

}



void MainWindow::comboboxIndexChanged_visualizations(int index)
{
    /*
    controller.model.floes_vtk.UnsafeUpdateValues(controller.model.floes.nodes.get(),
                                                  controller.model.floes.elems.get(),
                                                  controller.prms.temporal_attenuation,
                                                  (icy::FloeVisualization::VisOpt)index);
    prefsGUI.VisualizationOption = index;
    renderWindow->Render();
    scalarBar->SetVisibility(prefsGUI.ShowScalarBar && prefsGUI.VisualizationOption!=0);
*/
}


// callback functions

void MainWindow::PickCallbackFunction(vtkObject* caller,
                          long unsigned int vtkNotUsed(eventId),
                          void* clientData,
                          void* vtkNotUsed(callData))
{
    vtkPointPicker* pp = static_cast<vtkPointPicker*>(caller);
    MainWindow *mw = (MainWindow*)clientData;

    vtkIdType id = pp->GetPointId();
    qDebug() << "pointpicker vtkCommand::EndPickEvent " << id;
    std::cout << "pointpicker" <<std::endl;

    /*
    mw->controller.model.floes_vtk.UnsafeUpdateSelection(mw->controller.model.floes.nodes.get(),
                                                     pp->GetPointId());
    mw->renderWindow->Render();
    */
}


void MainWindow::SelectionChangedCallbackFunction ( vtkObject* vtkNotUsed(caller),
  long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* callData )
{
  unsigned int* rect = reinterpret_cast<unsigned int*> ( callData );
  unsigned int pos1X = rect[0];
  unsigned int pos1Y = rect[1];
  unsigned int pos2X = rect[2];
  unsigned int pos2Y = rect[3];

  std::cout << "vtkCommand::SelectionChangedEvent " << "Start x: " << pos1X << " Start y: " << pos1Y
            << " End x: " << pos2X << " End y: " << pos2Y << std::endl;
}






