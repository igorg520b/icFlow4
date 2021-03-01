#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"


icy::Model::Model()
{
    ugrid->SetPoints(points);


    dataSetMapper->SetInputData(ugrid);
    dataSetMapper->UseLookupTableScalarRangeOn();
    dataSetMapper->SetLookupTable(hueLut);

    actor_mesh->SetMapper(dataSetMapper);
    actor_mesh->GetProperty()->VertexVisibilityOff();
    actor_mesh->GetProperty()->EdgeVisibilityOn();
    actor_mesh->PickableOn();
    actor_mesh->GetProperty()->SetColor(0.54938, 0.772213, 0.848103);
//    actor_mesh->GetProperty()->SetEdgeColor(91.0/255.0, 116.0/255.0, 145.0/255.0);
    actor_mesh->GetProperty()->SetLineWidth(1.5);
    actor_mesh->GetProperty()->LightingOff();
    actor_mesh->GetProperty()->ShadingOff();
    actor_mesh->GetProperty()->SetInterpolationToFlat();

    visualized_values->SetName("visualized_values");
}

void icy::Model::Reset(SimParams &prms)
{
    mesh.Reset(prms.CharacteristicLength);
    UnsafeUpdateGeometry();
}

void icy::Model::AssembleAndSolve(SimParams &prms, double timeStep) {}
void icy::Model::GetResultFromSolver(double timeStep) {}
void icy::Model::AcceptTentativeValues(SimParams &prms) {}


void icy::Model::UnsafeUpdateGeometry()
{
    points->SetNumberOfPoints(mesh.nodes.size());
    vtkIdType count=0;
    double x[3];
    for(icy::Node &nd : mesh.nodes)
    {
        x[0]=nd.xn.x();
        x[1]=nd.xn.y();
        x[2]=0;
        points->SetPoint(count++, x);
    }
    points->Modified();


    cellArray->Reset();
    // create ugrid
    vtkIdType pts2[3];
    for(icy::Element &tr : mesh.elems)
    {
        for(int j=0;j<3;j++) pts2[j] = tr.nds[j]->id;
        cellArray->InsertNextCell(3, pts2);
    }
    ugrid->SetCells(VTK_TRIANGLE, cellArray);

    // if(selectedPointId >= 0) UnsafeUpdateSelection(nodes, -1);
}


void icy::Model::ChangeVisualizationOption(VisOpt option)
{


}



void icy::Model::InitializeLUT(int table=1)
{
    const int n = 51;
    hueLut->SetNumberOfTableValues(n);

    if(table==1)
    for ( int i=0; i<n; i++)
            hueLut->SetTableValue(i, (double)lutArrayTerrain[i][0],
                    (double)lutArrayTerrain[i][1],
                    (double)lutArrayTerrain[i][2], 1.0);
}
