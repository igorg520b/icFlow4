#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <QFileInfo>
#include <QObject>
#include <QMutex>

#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkSphereSource.h>

#include <vtkLabeledDataMapper.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkIdFilter.h>
#include <vtkCellCenters.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkVertexGlyphFilter.h>

#include "parameters_sim.h"
#include "mesh.h"
#include "equationofmotionsolver.h"

#include <Eigen/Core>

namespace icy { class Model; class Node; class Element;}

class icy::Model : public QObject
{
    Q_OBJECT

public:    

    // visualization options
    enum VisOpt { none, elem_area, energy_density };
    Q_ENUM(VisOpt)

    Model();
    void Reset(SimParams &prms);

    void InitialGuess(SimParams &prms, double timeStep);
    void AssembleAndSolve(SimParams &prms, double timeStep);
    void AcceptTentativeValues(double timeStep);
    void UnsafeUpdateGeometry();    // called from the main thread
    void UpdateValues();
    void ChangeVisualizationOption(VisOpt option);  // called from the main thread

    icy::Mesh mesh;

    vtkNew<vtkActor> actor_mesh;
    vtkNew<vtkActor> actor_selected_nodes;

private:
    EquationOfMotionSolver eqOfMotion;

    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry
    VisOpt VisualizingVariable = VisOpt::none;

    vtkNew<vtkLookupTable> hueLut;
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkPoints> points;

    // 2D mesh
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkCellArray> cellArray;
    vtkNew<vtkDataSetMapper> dataSetMapper;

    vtkNew<vtkVertexGlyphFilter> glyph_filter;
    vtkNew<vtkPolyDataMapper> glyph_mapper;
    vtkNew<vtkPolyData> poly_data;
    vtkNew<vtkIntArray> glyph_int_data;
    vtkNew<vtkLookupTable> glyph_hueLut;


    // visualizing variables
    vtkNew<vtkDoubleArray> visualized_values;
    vtkIdType selectedPointId = -1;

    void InitializeLUT(int table);

    std::size_t freeNodeCount=-1;



signals:
    void requestGeometryUpdate(); // this goes to the main thread, which calls UnsafeUpdateGeometry()
    void propertyChanged();

private:
    static constexpr float lutArrayBands[6][3] =
{{0.084154, 0.075059, 0.095523},
    {0.0882734, 0.623682,0.0375036},
    {0.801724, 0.0748252, 0.0492262},
    {0.745543, 0.423584, 0.52455},
    {0.756931, 0.53206, 0.359658},
    {0.384161, 0.395036, 0.599407} };

    static constexpr float lutArrayTerrain[51][3] =
    {{0.54938, 0.772213, 0.848103},
     {0.54349, 0.751689, 0.811536},
     {0.5376, 0.731165, 0.77497},
     {0.53171, 0.71064, 0.738403},
     {0.525821, 0.690116, 0.701836},
     {0.523788, 0.674231, 0.670908},
     {0.526385, 0.663912, 0.646744},
     {0.528982, 0.653593, 0.622581},
     {0.531579, 0.643275, 0.598418},
     {0.534176, 0.632956, 0.574255},
     {0.542585, 0.630953, 0.560901},
     {0.551575, 0.629781, 0.548628},
     {0.560566, 0.628609, 0.536356},
     {0.569556, 0.627438, 0.524083},
     {0.579925, 0.628775, 0.515402},
     {0.592706, 0.634504, 0.513005},
     {0.605487, 0.640233, 0.510609},
     {0.618268, 0.645962, 0.508213},
     {0.631049, 0.651691, 0.505817},
     {0.644928, 0.660899, 0.509458},
     {0.65905, 0.67088, 0.514441},
     {0.673172, 0.680861, 0.519424},
     {0.687294, 0.690842, 0.524407},
     {0.701252, 0.701293, 0.530766},
     {0.71477, 0.712998, 0.540793},
     {0.728289, 0.724704, 0.55082},
     {0.741807, 0.736409, 0.560848},
     {0.755325, 0.748114, 0.570875},
     {0.767446, 0.759552, 0.583215},
     {0.779043, 0.77089, 0.596421},
     {0.79064, 0.782228, 0.609627},
     {0.802237, 0.793566, 0.622834},
     {0.813354, 0.804566, 0.636368},
     {0.822309, 0.81404, 0.651377},
     {0.831264, 0.823514, 0.666387},
     {0.84022, 0.832989, 0.681397},
     {0.849175, 0.842463, 0.696406},
     {0.8563, 0.850209, 0.711857},
     {0.862381, 0.856968, 0.727562},
     {0.868461, 0.863726, 0.743266},
     {0.874541, 0.870485, 0.75897},
     {0.880373, 0.876977, 0.774624},
     {0.883714, 0.880806, 0.789783},
     {0.887055, 0.884635, 0.804943},
     {0.890396, 0.888464, 0.820102},
     {0.893737, 0.892293, 0.835261},
     {0.895825, 0.894749, 0.849094},
     {0.896869, 0.896062, 0.86182},
     {0.897913, 0.897375, 0.874547},
     {0.898956, 0.898687, 0.887273},
     {0.9, 0.9, 0.9}};
};

#endif // MESHCOLLECTION_H
