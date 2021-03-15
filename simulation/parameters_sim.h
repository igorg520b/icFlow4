#ifndef P_SIM_H
#define P_SIM_H

#include <QObject>
#include <QDebug>

#include <Eigen/Core>
#include <iostream>

// variables related to the formulation of the model

namespace icy { class SimParams; }

class icy::SimParams : public QObject
{
    Q_OBJECT

    // general
    Q_PROPERTY(int s_MaxSteps MEMBER MaxSteps NOTIFY propertyChanged)

    // integration
    Q_PROPERTY(double in_InitialTimeStep MEMBER InitialTimeStep NOTIFY propertyChanged)
    Q_PROPERTY(double in_ConvergenceEpsilon MEMBER ConvergenceEpsilon NOTIFY propertyChanged)
    Q_PROPERTY(double in_ConvergenceCutoff MEMBER ConvergenceCutoff NOTIFY propertyChanged)
    Q_PROPERTY(int in_MinIter MEMBER MinIter NOTIFY propertyChanged)
    Q_PROPERTY(int in_MaxIter MEMBER MaxIter NOTIFY propertyChanged)

    // material parameters and physical constants
    Q_PROPERTY(double p_Gravity MEMBER Gravity NOTIFY propertyChanged)
    Q_PROPERTY(double p_Density MEMBER Density NOTIFY propertyChanged)
    Q_PROPERTY(double p_YoungsModulus MEMBER YoungsModulus NOTIFY propertyChanged)
    Q_PROPERTY(double p_PoissonsRatio MEMBER PoissonsRatio NOTIFY propertyChanged)

    Q_PROPERTY(double p_Thickness MEMBER Thickness NOTIFY propertyChanged)

    // meshing
    Q_PROPERTY(double s_ElemSize MEMBER CharacteristicLength NOTIFY propertyChanged)

public:
    int MaxSteps, MinIter, MaxIter;
    double InitialTimeStep;
    double Gravity, Density, PoissonsRatio, YoungsModulus, Thickness;
    double CharacteristicLength;

    double ConvergenceEpsilon, ConvergenceCutoff;

    SimParams() { Reset(); }

    void Reset() {
        MaxSteps = 20000;
        InitialTimeStep = 0.005;

        // material parameters and physical constants
        Gravity = 9.81;
        Density = 1;
        Thickness = 0.1;

        PoissonsRatio = 0.3;
        YoungsModulus = 50;

        CharacteristicLength = 0.05;
        ConvergenceEpsilon = 1e-2;
        ConvergenceCutoff = 1e-7;

        MinIter = 3;
        MaxIter = 10;

        emit propertyChanged();
    }


signals:
    void propertyChanged();
};


#endif
