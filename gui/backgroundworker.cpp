#include "backgroundworker.h"

BackgroundWorker::BackgroundWorker()
{
    this->start();
}

// resume the worker thread
void BackgroundWorker::Resume()
{
//    controller->Preapare();
    condition.wakeOne();
}

// cancel current step and pause the worker thread
void BackgroundWorker::Pause()
{
    if(!running) return;
    timeToPause = true;
//    model->RequestAbort();
}

// exit the worker thread
void BackgroundWorker::Finalize()
{
    qDebug() << "BackgroundWorker::Finalize()";
//    controller->RequestAbort();
    kill=true;
    condition.wakeOne();
    bool result = wait();
    qDebug() << "BackgroundWorker::Finalize() terminated";
    qDebug() << "QThread wait() returns " << result;
    running = false;
}

void BackgroundWorker::run()
{
    while(!kill) {
        if (timeToPause) {
            mutex.lock();
            running = false;
            emit workerPaused();
            condition.wait(&mutex);
            timeToPause = false;
            running = true;
            mutex.unlock();
        }
        if(kill) break;

        qDebug() << "sleeping for 1 second";
        sleep(1);

//        bool result = controller->Step();
//        if(!result) timeToPause = true;
    }
    qDebug() << "BackgroundWorker::run() terminated";
}
