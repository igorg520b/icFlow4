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
    enum VisOpt { none, elem_area, energy_density, stress_xx, stress_yy, stress_hydrostatic, non_symm_measure,
                ps1, ps2, shear_stress, volume_change};
    Q_ENUM(VisOpt)

    Model();
    void Reset(SimParams &prms);

    void InitialGuess(SimParams &prms, double timeStep, double timeStepFactor);
    bool AssembleAndSolve(SimParams &prms, double timeStep);    // return what solver returns
    void AcceptTentativeValues(double timeStep);
    void UnsafeUpdateGeometry();    // called from the main thread
    void UpdateValues();
    void ChangeVisualizationOption(VisOpt option);  // called from the main thread

    icy::Mesh mesh;

    vtkNew<vtkActor> actor_mesh;
    vtkNew<vtkActor> actor_selected_nodes;
    vtkNew<vtkActor> actor_boundary;
    vtkNew<vtkActor> actor_indenter;
    vtkNew<vtkActor> actor_indenter_intended;

    EquationOfMotionSolver eqOfMotion;
    vtkNew<vtkLookupTable> hueLut;

private:

    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry
    VisOpt VisualizingVariable = VisOpt::none;

    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkPoints> points;
    vtkNew<vtkPoints> points_manip;

    // 2D mesh
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkCellArray> cellArray;
    vtkNew<vtkDataSetMapper> dataSetMapper;

    // nodes drawn separately
    vtkNew<vtkVertexGlyphFilter> glyph_filter;
    vtkNew<vtkPolyDataMapper> glyph_mapper;
    vtkNew<vtkPolyData> poly_data;
    vtkNew<vtkIntArray> glyph_int_data;
    vtkNew<vtkLookupTable> glyph_hueLut;

    // boundary
    vtkNew<vtkUnstructuredGrid> ugrid_boundary;
    vtkNew<vtkCellArray> cellArray_boundary;
//    vtkNew<vtkPolyDataMapper> boundary_mapper;
    vtkNew<vtkDataSetMapper> dataSetMapper_boundary;

    // indenter
    vtkNew<vtkPoints> points_indenter;
    vtkNew<vtkUnstructuredGrid> ugrid_indenter;
    vtkNew<vtkCellArray> cellArray_indenter;
    vtkNew<vtkDataSetMapper> dataSetMapper_indenter;

    // indenter-intended
    vtkNew<vtkPoints> points_indenter_intended;
    vtkNew<vtkUnstructuredGrid> ugrid_indenter_intended;
    vtkNew<vtkCellArray> cellArray_indenter_intended;
    vtkNew<vtkDataSetMapper> dataSetMapper_indenter_intended;

    // visualizing variables
    vtkNew<vtkDoubleArray> visualized_values;
    vtkIdType selectedPointId = -1;

    void InitializeLUT(int table=1);

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


    static constexpr float lutArrayTemperatureAdj[51][3] =
{{0.770938, 0.951263, 0.985716}, {0.788065, 0.959241,
  0.986878}, {0.805191, 0.96722, 0.98804}, {0.822318, 0.975199,
  0.989202}, {0.839445, 0.983178, 0.990364}, {0.856572, 0.991157,
  0.991526}, {0.872644, 0.995552, 0.98386}, {0.887397, 0.995466,
  0.965157}, {0.902149, 0.99538, 0.946454}, {0.916902, 0.995294,
  0.927751}, {0.931655, 0.995208, 0.909049}, {0.946408, 0.995123,
  0.890346}, {0.961161, 0.995037, 0.871643}, {0.975913, 0.994951,
  0.85294}, {0.990666, 0.994865, 0.834237}, {0.996257, 0.991758,
  0.815237}, {0.994518, 0.986234, 0.795999}, {0.992779, 0.98071,
  0.77676}, {0.99104, 0.975186, 0.757522}, {0.989301, 0.969662,
  0.738283}, {0.987562, 0.964138, 0.719045}, {0.985823, 0.958614,
  0.699807}, {0.984084, 0.953089, 0.680568}, {0.982345, 0.947565,
  0.66133}, {0.97888, 0.936201, 0.641773}, {0.974552, 0.921917,
  0.622058}, {0.970225, 0.907633, 0.602342}, {0.965897, 0.893348,
  0.582626}, {0.961569, 0.879064, 0.562911}, {0.957242, 0.86478,
  0.543195}, {0.952914, 0.850496, 0.52348}, {0.948586, 0.836212,
  0.503764}, {0.944259, 0.821927, 0.484048}, {0.939066, 0.801586,
  0.464871}, {0.933626, 0.779513, 0.445847}, {0.928186, 0.757441,
  0.426823}, {0.922746, 0.735368, 0.4078}, {0.917306, 0.713296,
  0.388776}, {0.911866, 0.691223, 0.369752}, {0.906426, 0.669151,
  0.350728}, {0.900986, 0.647078, 0.331704}, {0.895546, 0.625006,
  0.312681}, {0.889975, 0.597251, 0.298625}, {0.884388, 0.568785,
  0.285191}, {0.8788, 0.54032, 0.271756}, {0.873212, 0.511855,
  0.258322}, {0.867625, 0.483389, 0.244888}, {0.862037, 0.454924,
  0.231453}, {0.856449, 0.426459, 0.218019}, {0.850862, 0.397993,
  0.204584}, {0.845274, 0.369528, 0.19115}};
};

#endif // MESHCOLLECTION_H
