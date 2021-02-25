#include "modelcontrollerinterface.h"
#include <chrono>
#include <thread>
#include <iostream>

void ModelControllerTest::Prepare(void)
{
    abortRequested = false;
    counter = 0;
}

bool ModelControllerTest::Step(void)
{
    counter++;
    std::cout << "ModelControllerTest::Step " << counter << std::endl;
    for(int i=0;i<30;i++) {
        if(!abortRequested) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::cout << i << " substep of step " << counter << std::endl;
        } else {
            std::cout << "abort requested at " << i << " substep of step " << counter << std::endl;
            break;
        }
    }
    return(counter < 3);
}

void ModelControllerTest::RequestAbort(void)
{
    abortRequested = true;
}
