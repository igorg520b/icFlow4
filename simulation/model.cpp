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
//    actor_mesh->GetProperty()->SetVertexColor(0.1,0.9,0.1);
//    actor_mesh->GetProperty()->SetPointSize(5);

    actor_mesh->GetProperty()->EdgeVisibilityOn();
//    actor_mesh->GetProperty()->SetLineWidth(1.0);
    actor_mesh->GetProperty()->SetColor(0.84938, 0.872213, 0.848103);
    actor_mesh->GetProperty()->SetEdgeColor(91.0/255.0, 116.0/255.0, 145.0/255.0);
    actor_mesh->GetProperty()->LightingOff();
    actor_mesh->GetProperty()->ShadingOff();
    actor_mesh->GetProperty()->SetInterpolationToFlat();
    actor_mesh->PickableOn();

    visualized_values->SetName("visualized_values");


    // selected and pinned points
    glyph_hueLut->SetNumberOfTableValues(6);
    for ( int i=0; i<6; i++)
            glyph_hueLut->SetTableValue(i, (double)lutArrayBands[i][0], (double)lutArrayBands[i][1], (double)lutArrayBands[i][2], 1.0);


    actor_selected_nodes->GetProperty()->SetPointSize(5);
    poly_data->SetPoints(points);
    glyph_filter->SetInputData(poly_data);
    glyph_mapper->SetInputConnection(glyph_filter->GetOutputPort());
    glyph_mapper->UseLookupTableScalarRangeOn();
    glyph_mapper->SetLookupTable(glyph_hueLut);

    actor_selected_nodes->SetMapper(glyph_mapper);
    actor_selected_nodes->GetProperty()->SetColor(0.1, 0.1, 0.1);
    glyph_int_data->SetName("glyph_int_data");





}

void icy::Model::Reset(SimParams &prms)
{
    mesh.Reset(prms.CharacteristicLength);
    UnsafeUpdateGeometry();
}




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

//    glyph_mapper->Modified();
    glyph_int_data->SetNumberOfValues(mesh.nodes.size());

    //for testing
    for(int i=0;i<mesh.nodes.size();i++)
    {
        int value = 0;
        if(mesh.nodes[i].pinned)
            value = mesh.nodes[i].selected ? 2 : 1;
        glyph_int_data->SetValue(i,value);
    }
    glyph_hueLut->SetTableRange(-0.5,5.5);

    glyph_mapper->SetScalarModeToUsePointData();
    poly_data->GetPointData()->AddArray(glyph_int_data);
    poly_data->GetPointData()->SetActiveScalars("glyph_int_data");
//    ugrid->GetPointData()->AddArray(visualized_values);
//            ugrid->GetPointData()->SetActiveScalars("visualized_values");
//            dataSetMapper->SetScalarModeToUsePointData();

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


void icy::Model::AssembleAndSolve(SimParams &prms, double timeStep) {}

void icy::Model::GetResultFromSolver(double timeStep) {}

void icy::Model::AcceptTentativeValues(SimParams &prms) {}
