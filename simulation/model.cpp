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


    poly_data->SetPoints(points_manip);
    glyph_filter->SetInputData(poly_data);
    glyph_mapper->SetInputConnection(glyph_filter->GetOutputPort());
    glyph_mapper->UseLookupTableScalarRangeOn();
    glyph_mapper->SetLookupTable(glyph_hueLut);

    actor_selected_nodes->SetMapper(glyph_mapper);
    actor_selected_nodes->GetProperty()->SetColor(0.1, 0.1, 0.1);
    actor_selected_nodes->GetProperty()->SetPointSize(5);
    glyph_int_data->SetName("glyph_int_data");
    InitializeLUT(2);

    ugrid_boundary->SetPoints(points);
    dataSetMapper_boundary->SetInputData(ugrid_boundary);
    actor_boundary->SetMapper(dataSetMapper_boundary);
    actor_boundary->GetProperty()->EdgeVisibilityOn();
    actor_boundary->GetProperty()->VertexVisibilityOn();
    actor_boundary->GetProperty()->SetColor(0,0,0.4);
    actor_boundary->GetProperty()->SetEdgeColor(0,0,0.4);
    actor_boundary->GetProperty()->SetVertexColor(0.4,0,0);
    actor_boundary->GetProperty()->SetPointSize(5);
    actor_boundary->GetProperty()->SetLineWidth(3);


    // indenter

    ugrid_indenter->SetPoints(points_indenter);
    dataSetMapper_indenter->SetInputData(ugrid_indenter);

    actor_indenter->SetMapper(dataSetMapper_indenter);
    actor_indenter->GetProperty()->EdgeVisibilityOn();
    actor_indenter->GetProperty()->VertexVisibilityOn();
    actor_indenter->GetProperty()->SetColor(0.3,0,0.4);
    actor_indenter->GetProperty()->SetEdgeColor(0.3,0,0.4);
    actor_indenter->GetProperty()->SetVertexColor(0.4,0,0);
    actor_indenter->GetProperty()->SetPointSize(5);
    actor_indenter->GetProperty()->SetLineWidth(3);

    // indenter_intended

    ugrid_indenter_intended->SetPoints(points_indenter_intended);
    dataSetMapper_indenter_intended->SetInputData(ugrid_indenter_intended);

    actor_indenter_intended->SetMapper(dataSetMapper_indenter_intended);
    actor_indenter_intended->GetProperty()->EdgeVisibilityOn();
    actor_indenter_intended->GetProperty()->VertexVisibilityOn();
    actor_indenter_intended->GetProperty()->SetColor(0.3,0,0.4);
    actor_indenter_intended->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_indenter_intended->GetProperty()->SetVertexColor(0.04,0,0);
    actor_indenter_intended->GetProperty()->SetPointSize(2);
    actor_indenter_intended->GetProperty()->SetLineWidth(1);

}

void icy::Model::Reset(SimParams &prms)
{
    mesh.Reset(prms.CharacteristicLength);
    UnsafeUpdateGeometry();
}




void icy::Model::UnsafeUpdateGeometry()
{
    points->SetNumberOfPoints(mesh.nodes.size());
    points_manip->SetNumberOfPoints(mesh.nodes.size());
    vtkIdType count=0;
    double x[3];
    for(icy::Node &nd : mesh.nodes)
    {
        x[0]=nd.xn.x();
        x[1]=nd.xn.y();
        x[2]=0;
        points->SetPoint(count, x);
        if(nd.selected)
        {
            x[0]=nd.intended_position.x();
            x[1]=nd.intended_position.y();
        }
        points_manip->SetPoint(count, x);
        count++;
    }
    points->Modified();
    points_manip->Modified();


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

    // boundary
    cellArray_boundary->Reset();
    for(unsigned i=0;i<mesh.boundary.size();i++)
    {
        int idx1 = mesh.boundary[i].first;
        int idx2 = mesh.boundary[i].second;
        pts2[0]=idx1;
        pts2[1]=idx2;
        cellArray_boundary->InsertNextCell(2, pts2);
    }
    ugrid_boundary->SetCells(VTK_LINE, cellArray_boundary);

    // indenter
    points_indenter->SetNumberOfPoints(mesh.nodes_indenter.size());
    count=0;
    for(icy::Node &nd : mesh.nodes_indenter)
    {
        x[0]=nd.xn.x();
        x[1]=nd.xn.y();
        x[2]=0;
        points_indenter->SetPoint(count, x);
        count++;
    }
    points_indenter->Modified();

    cellArray_indenter->Reset();
    for(unsigned i=0;i<mesh.boundary_indenter.size();i++)
    {
        int idx1 = mesh.boundary_indenter[i].first;
        int idx2 = mesh.boundary_indenter[i].second;
        pts2[0]=idx1;
        pts2[1]=idx2;
        cellArray_indenter->InsertNextCell(2, pts2);
    }
    ugrid_indenter->SetCells(VTK_LINE, cellArray_indenter);

    // indenter_intended
    points_indenter_intended->SetNumberOfPoints(mesh.nodes_indenter.size());
    count=0;
    for(icy::Node &nd : mesh.nodes_indenter)
    {
        x[0]=nd.intended_position.x();
        x[1]=nd.intended_position.y();
        x[2]=0;
        points_indenter_intended->SetPoint(count, x);
        count++;
    }
    points_indenter_intended->Modified();

    cellArray_indenter_intended->Reset();
    for(unsigned i=0;i<mesh.boundary_indenter.size();i++)
    {
        int idx1 = mesh.boundary_indenter[i].first;
        int idx2 = mesh.boundary_indenter[i].second;
        pts2[0]=idx1;
        pts2[1]=idx2;
        cellArray_indenter_intended->InsertNextCell(2, pts2);
    }
    ugrid_indenter_intended->SetCells(VTK_LINE, cellArray_indenter_intended);


    if(VisualizingVariable != VisOpt::none) UpdateValues();
}


void icy::Model::ChangeVisualizationOption(VisOpt option)
{
    qDebug() << "icy::Model::ChangeVisualizationOption " << option;
    if(VisualizingVariable == option) return; // option did not change
    VisualizingVariable = option;


    if(VisualizingVariable == VisOpt::none)
    {
        dataSetMapper->ScalarVisibilityOff();
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->RemoveArray("visualized_values");
        return;
    }
    else
    {
        ugrid->GetPointData()->RemoveArray("visualized_values");
        ugrid->GetCellData()->AddArray(visualized_values);
        ugrid->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper->SetScalarModeToUseCellData();
        dataSetMapper->ScalarVisibilityOn();
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

    case energy_density:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].strain_energy_density);
        break;

    case stress_xx:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].CauchyStress(0,0));
        break;

    case stress_yy:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].CauchyStress(1,1));
        break;

    case stress_hydrostatic:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].CauchyStress.trace()/2);
        break;

    case non_symm_measure:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++)
        {
            Element &elem = mesh.elems[i];
            double value1 = elem.CauchyStress(0,1)-elem.CauchyStress(1,0);
            double value = value1*value1;
            visualized_values->SetValue(i, value);
        }
        break;

    case ps1:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].principal_stress1);
        break;

    case ps2:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].principal_stress2);
        break;

    case shear_stress:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].max_shear_stress);
        break;

    case volume_change:
        visualized_values->SetNumberOfValues(mesh.elems.size());
        for(size_t i=0;i<mesh.elems.size();i++) visualized_values->SetValue(i, mesh.elems[i].volume_change);
        break;

    default:
        break;
    }
    vtk_update_mutex.unlock();

    visualized_values->Modified();

    double minmax[2];
    visualized_values->GetValueRange(minmax);
    hueLut->SetTableRange(minmax[0], minmax[1]);
    if(VisualizingVariable == volume_change)
    {
        double range = minmax[1]-minmax[0];
        hueLut->SetTableRange(1-range*0.75,1+range*0.75);
    }

}


void icy::Model::InitializeLUT(int table)
{
    const int n = 51;
    hueLut->SetNumberOfTableValues(n);

    if(table==1)
    for ( int i=0; i<n; i++)
            hueLut->SetTableValue(i, (double)lutArrayTerrain[i][0],
                    (double)lutArrayTerrain[i][1],
                    (double)lutArrayTerrain[i][2], 1.0);
    else if(table==2)
        for ( int i=0; i<n; i++)
                hueLut->SetTableValue(i, (double)lutArrayTemperatureAdj[i][0],
                        (double)lutArrayTemperatureAdj[i][1],
                        (double)lutArrayTemperatureAdj[i][2], 1.0);
}

void icy::Model::InitialGuess(SimParams &prms, double timeStep, double timeStepFactor)
{
    std::size_t nNodes = mesh.nodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes[i];
        if(nd.pinned)
        {
            if(nd.selected) nd.xt = (1/timeStepFactor)*nd.intended_position + (1-1/timeStepFactor)*nd.xn;
            nd.vn=Eigen::Vector2d::Zero();
            nd.x_hat = nd.xn;
        }
        else
        {
            nd.x_hat = nd.xn + timeStep*nd.vn;
            nd.x_hat.y() -= prms.Gravity*timeStep*timeStep;
            nd.xt = nd.xn;
        }
    }

    freeNodeCount = 0;
    for(icy::Node &nd : mesh.nodes)
    {
        if(nd.pinned) nd.eqId = -1;
        else nd.eqId = freeNodeCount++;
    }

    nNodes = mesh.nodes_indenter.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes_indenter[i];
        nd.xt = (1/timeStepFactor)*nd.intended_position + (1-1/timeStepFactor)*nd.xn;
        nd.vn=Eigen::Vector2d::Zero();
        nd.x_hat = nd.xn;
    }
}


bool icy::Model::AssembleAndSolve(SimParams &prms, double timeStep)
{
    eqOfMotion.ClearAndResize(freeNodeCount);

    std::size_t nElems = mesh.elems.size();
    std::size_t nNodes = mesh.nodes.size();

#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++) mesh.elems[i].AddToSparsityStructure(eqOfMotion);

    eqOfMotion.CreateStructure();

    // assemble
    bool mesh_iversion_detected = false;
#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++)
    {
        bool result = mesh.elems[i].ComputeEquationEntries(eqOfMotion, prms, timeStep);
        if(!result) mesh_iversion_detected = true;
    }

    if(mesh_iversion_detected) return false; // mesh inversion

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++) mesh.nodes[i].ComputeEquationEntries(eqOfMotion, prms, timeStep);

    mesh.DetectContactPairs(prms.InteractionDistance);

    // solve
    bool result = eqOfMotion.Solve();

    // pull
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes[i];
        Eigen::Vector2d delta_x;
        if(!nd.pinned) {
            eqOfMotion.GetTentativeResult(nd.eqId, delta_x);
            nd.xt+=delta_x;
        }
    }

    return result;
}

void icy::Model::AcceptTentativeValues(double timeStep)
{
    vtk_update_mutex.lock();
    std::size_t nNodes = mesh.nodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes[i];
        if(!nd.pinned)
        {
            Eigen::Vector2d dx = nd.xt-nd.xn;
            nd.vn = dx/timeStep;
        }
        nd.xn = nd.xt;
    }

    nNodes = mesh.nodes_indenter.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node &nd = mesh.nodes_indenter[i];
        Eigen::Vector2d dx = nd.xt-nd.xn;
        nd.vn = dx/timeStep;
        nd.xn = nd.xt;
    }

    vtk_update_mutex.unlock();
}

