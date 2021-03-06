#include "equationofmotionsolver.h"
#include <iostream>

EquationOfMotionSolver::EquationOfMotionSolver()
{
    r = MSK_makeenv(&env, NULL);
    if (r != MSK_RES_OK) throw std::runtime_error("makeenv");

    constexpr int numvar = 10000, numcon = 0;
    r = MSK_maketask(env, numcon, numvar, &task);
    if (r != MSK_RES_OK) throw std::runtime_error("maketask");

    r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
    if (r != MSK_RES_OK) throw std::runtime_error("linkfunctotaskstream");
}

EquationOfMotionSolver::~EquationOfMotionSolver()
{
    MSK_deletetask(&task);
    MSK_deleteenv(&env);
}


void MSKAPI EquationOfMotionSolver::printstr(void *, const char str[])
{
    std::cout << str << std::endl;
}

//MSK_putcfix


void EquationOfMotionSolver::ClearAndResize(std::size_t N)
{
    cfix = 0;
}


void EquationOfMotionSolver::AddElementToStructure(int row, int column)
{

}

void EquationOfMotionSolver::CreateStructure()
{

}

// creating the values array
void EquationOfMotionSolver::AddToQ(const int row, const int column, const Eigen::Matrix2d &mat)
{

}

void EquationOfMotionSolver::AddToC(const int idx, const Eigen::Vector2d &vec)
{

}

void EquationOfMotionSolver::AddToConstTerm(double c)
{

}

void EquationOfMotionSolver::Solve()
{

}

void EquationOfMotionSolver::TestSolve()
{

}

void EquationOfMotionSolver::AdjustCurrentGuess(int idx, Eigen::Vector2d &vec)
{

}


