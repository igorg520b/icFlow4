#ifndef IMPLICITMODEL5_H
#define IMPLICITMODEL5_H

#include <vector>
#include <chrono>
#include <QObject>
#include <QDebug>

//#include "model.h"

#include "modelcontrollerinterface.h"

namespace icy { class ModelController; }

class icy::ModelController : public QObject, public ModelControllerInterface
{
    Q_OBJECT

public:

    void Prepare(void) override;
    bool Step(void) override;
    void RequestAbort(void) override;

//    icy::Model model;
//    SimParams prms;

    int currentStep = 0;

private:
    bool abortRequested = false;
    void Aborting();       // called before exiting Step() if aborted

signals:
    void stepCompleted();
    void stepAborted();
};

#endif // IMPLICITMODEL5_H
