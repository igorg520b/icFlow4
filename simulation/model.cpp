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


    poly_data->SetPoints(points);
    glyph_filter->SetInputData(poly_data);
    glyph_mapper->SetInputConnection(glyph_filter->GetOutputPort());
    glyph_mapper->UseLookupTableScalarRangeOn();
    glyph_mapper->SetLookupTable(glyph_hueLut);

    actor_selected_nodes->SetMapper(glyph_mapper);
    actor_selected_nodes->GetProperty()->SetColor(0.1, 0.1, 0.1);
    actor_selected_nodes->GetProperty()->SetPointSize(5);
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

    glyph_int_data->SetNumberOfValues(mesh.nodes.size());

    for(std::size_t i=0;i<mesh.nodes.size();i++)
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
}


void icy::Model::ChangeVisualizationOption(VisOpt option)
{
    qDebug() << "icy::Model::ChangeVisualizationOption " << option;
    if(VisualizingVariable == option) return; // option did not change
    VisualizingVariable = option;


    if(VisualizingVariable == VisOpt::elem_area)
    {
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->AddArray(visualized_values);
        ugrid->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper->SetScalarModeToUseCellData();
        dataSetMapper->ScalarVisibilityOn();
    }
    else
    {
        dataSetMapper->ScalarVisibilityOff();
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->RemoveArray("visualized_values");
        return;
    }


    UpdateValues();

}

void icy::Model::UpdateValues()
{
    if(mesh.nodes.size()==0)
    {
        dataSetMapper->ScalarVisibilityOff();
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->RemoveArray("visualized_values");
        return;
    }

    vtk_update_mutex.lock();
    switch(VisualizingVariable)
    {
        case elem_area:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].area_initial);
        break;
    default:
        break;
    }
    vtk_update_mutex.unlock();

    visualized_values->Modified();

    double minmax[2];
    visualized_values->GetValueRange(minmax);
    hueLut->SetTableRange(minmax[0], minmax[1]);
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

void icy::Model::InitialGuess(double timeStep)
{
    std::size_t nNodes = mesh.nodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes[i];
        if(nd.pinned) continue;
        nd.xt = nd.xn + nd.vn*timeStep;
    }

    freeNodeCount = 0;
    for(icy::Node &nd : mesh.nodes)
    {
        if(nd.pinned) nd.eqId = -1;
        else nd.eqId = freeNodeCount++;
    }

}


void icy::Model::AssembleAndSolve(SimParams &prms, double timeStep)
{
    eqOfMotion.ClearAndResize(freeNodeCount);
}

void icy::Model::GetResultFromSolver(double timeStep)
{

}

void icy::Model::AcceptTentativeValues(double timeStep)
{
    vtk_update_mutex.lock();
    std::size_t nNodes = mesh.nodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes[i];
        if(nd.pinned) continue;
        Eigen::Vector2d dx = nd.xt-nd.xn;
        nd.xn = nd.xt;
        nd.vn = dx/timeStep;
        // nd.xn.x() += 0.001; // for testing

    }
    vtk_update_mutex.unlock();
}
