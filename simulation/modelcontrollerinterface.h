#ifndef MODELCONTROLLERINTERFACE_H
#define MODELCONTROLLERINTERFACE_H


class ModelControllerInterface
{
public:
    virtual void Prepare(void) = 0;
    virtual bool Step(void) = 0;
    virtual void RequestAbort(void) = 0;
};

class ModelControllerTest : public ModelControllerInterface
{
public:
    int counter = 0;
    bool abortRequested = false;

    void Prepare(void) override;
    bool Step(void) override;
    void RequestAbort(void) override;
};


#endif // MODELCONTROLLERINTERFACE_H
