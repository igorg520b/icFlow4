#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>





icy::Mesh::Mesh()
{
    InitializeLUT(2);
//    visualized_values->SetName("visualized_values");

    ugrid_deformable->SetPoints(points_deformable);
    dataSetMapper_deformable->SetInputData(ugrid_deformable);
    dataSetMapper_deformable->UseLookupTableScalarRangeOn();
    dataSetMapper_deformable->SetLookupTable(hueLut);

    actor_mesh_deformable->SetMapper(dataSetMapper_deformable);
    actor_mesh_deformable->GetProperty()->VertexVisibilityOff();
    actor_mesh_deformable->GetProperty()->EdgeVisibilityOn();
    actor_mesh_deformable->GetProperty()->SetColor(0.84938, 0.872213, 0.848103);
    actor_mesh_deformable->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
    actor_mesh_deformable->GetProperty()->LightingOff();
    actor_mesh_deformable->GetProperty()->ShadingOff();
    actor_mesh_deformable->GetProperty()->SetInterpolationToFlat();
    actor_mesh_deformable->PickableOff();
    actor_mesh_deformable->GetProperty()->SetLineWidth(1);

    // boundary that encompasses all objects
    ugrid_boundary_all->SetPoints(points_deformable);
    dataSetMapper_boundary_all->SetInputData(ugrid_boundary_all);
    actor_boundary_all->SetMapper(dataSetMapper_boundary_all);
    actor_boundary_all->GetProperty()->EdgeVisibilityOn();
    actor_boundary_all->GetProperty()->VertexVisibilityOn();
    actor_boundary_all->GetProperty()->SetColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetEdgeColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetVertexColor(0.4,0,0);
    actor_boundary_all->GetProperty()->SetPointSize(5);
    actor_boundary_all->GetProperty()->SetLineWidth(3);
    actor_boundary_all->PickableOff();

    /*

    // indenter's intended position
    ugrid_indenter_intended->SetPoints(points_indenter_intended);
    dataSetMapper_indenter_intended->SetInputData(ugrid_indenter_intended);

    actor_boundary_intended_indenter->SetMapper(dataSetMapper_indenter_intended);
    actor_boundary_intended_indenter->GetProperty()->EdgeVisibilityOn();
    actor_boundary_intended_indenter->GetProperty()->VertexVisibilityOn();
    actor_boundary_intended_indenter->GetProperty()->SetColor(0.3,0,0.4);
    actor_boundary_intended_indenter->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_boundary_intended_indenter->GetProperty()->SetVertexColor(0.04,0,0);
    actor_boundary_intended_indenter->GetProperty()->SetPointSize(2);
    actor_boundary_intended_indenter->GetProperty()->SetLineWidth(1);

    // collisions
    ugrid_collisions->SetPoints(points_collisions);
    mapper_collisions->SetInputData(ugrid_collisions);
    actor_collisions->SetMapper(mapper_collisions);
    actor_collisions->GetProperty()->EdgeVisibilityOn();
    actor_collisions->GetProperty()->VertexVisibilityOn();
    actor_collisions->GetProperty()->SetColor(0.3,0,0.4);
    actor_collisions->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_collisions->GetProperty()->SetVertexColor(0.04,0,0);
    actor_collisions->GetProperty()->SetPointSize(2);
    actor_collisions->GetProperty()->SetLineWidth(1);
    */
}


void icy::Mesh::Reset(double CharacteristicLengthMax)
{
    brick.GenerateBrick(CharacteristicLengthMax);
    indenter.GenerateIndenter(CharacteristicLengthMax);
    allMeshes.clear();
    allMeshes.push_back(&brick);
    allMeshes.push_back(&indenter);
    RegenerateVisualizedGeometry();
}

void icy::Mesh::RegenerateVisualizedGeometry()
{
    allNodes.clear();
    allElems.clear();
    allBoundaryEdges.clear();
    unsigned count = 0;
    freeNodeCount = 0;

    for(MeshFragment *mf : allMeshes)
    {
        for(unsigned i=0;i<mf->nodes.size();i++)
        {
            Node *nd = &mf->nodes[i];
            nd->globId = count++;
            if(nd->pinned) nd->eqId=-1;
            else nd->eqId=freeNodeCount++;
            allNodes.push_back(nd);
        }

        for(unsigned i=0;i<mf->elems.size();i++) allElems.push_back(&mf->elems[i]);
        allBoundaryEdges.insert(allBoundaryEdges.end(), mf->boundary_edges.begin(), mf->boundary_edges.end());
    }

    points_deformable->SetNumberOfPoints(allNodes.size());
    cellArray_deformable->Reset();

    // create ugrid
    for(icy::Element *tr : allElems)
    {
        vtkIdType pts[3] = {tr->nds[0]->globId, tr->nds[1]->globId, tr->nds[2]->globId};
        cellArray_deformable->InsertNextCell(3, pts);
    }
    ugrid_deformable->SetCells(VTK_TRIANGLE, cellArray_deformable);


    cellArray_boundary_all->Reset();

    for(auto edge : allBoundaryEdges)
    {
        vtkIdType pts[2] = {edge.first->globId, edge.second->globId};
        cellArray_boundary_all->InsertNextCell(2, pts);
    }

    ugrid_boundary_all->SetCells(VTK_LINE, cellArray_boundary_all);

}





double icy::Mesh::SegmentPointDistance(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d P, Eigen::Vector2d &D, double &t)
{
    Eigen::Vector2d seg = B-A;
    Eigen::Vector2d v = P-A;
    t = v.dot(seg)/seg.squaredNorm();
    t = std::clamp(t, 0.0, 1.0);
    D = A+seg*t;
    double dist = (D-P).norm();
    return dist;
}

void icy::Mesh::ChangeVisualizationOption(VisOpt option)
{
    /*
    qDebug() << "icy::Model::ChangeVisualizationOption " << option;
    if(VisualizingVariable == option) return; // option did not change
    VisualizingVariable = option;


    if(VisualizingVariable == VisOpt::none)
    {
        dataSetMapper_deformable->ScalarVisibilityOff();
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        return;
    }
    else
    {
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->AddArray(visualized_values);
        ugrid_deformable->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUseCellData();
        dataSetMapper_deformable->ScalarVisibilityOn();
    }
    UpdateValues();
    */
}






void icy::Mesh::UnsafeUpdateGeometry()
{
    for(icy::Node *nd : allNodes) points_deformable->SetPoint((vtkIdType)nd->globId, nd->xn.data());
    points_deformable->Modified();



/*
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
*/
    // boundary

    /*
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

    // collisions
    ugrid_collisions->Reset();
    ugrid_collisions->SetPoints(points_collisions);
    unsigned nCollisions1 = mesh.indenter_boundary_vs_deformable_nodes.size();
    unsigned nCollisions2 = mesh.deformable_boundary_vs_indenter_nodes.size();
    points_collisions->SetNumberOfPoints((nCollisions1+nCollisions2)*2);

    for(unsigned i=0;i<nCollisions1;i++)
    {
        Interaction &c = mesh.indenter_boundary_vs_deformable_nodes[i];
        x[0]=c.P.x();
        x[1]=c.P.y();
        x[2]=0;
        points_collisions->SetPoint(i*2+0, x);
        x[0]=c.D.x();
        x[1]=c.D.y();
        x[2]=0;
        points_collisions->SetPoint(i*2+1, x);

        pts2[0]=i*2+0;
        pts2[1]=i*2+1;
        ugrid_collisions->InsertNextCell(VTK_LINE, 2, pts2);
    }

    for(unsigned i=0;i<nCollisions2;i++)
    {
        Interaction &c = mesh.deformable_boundary_vs_indenter_nodes[i];
        x[0]=c.P.x();
        x[1]=c.P.y();
        x[2]=0;
        points_collisions->SetPoint(nCollisions1*2+i*2+0, x);
        x[0]=c.D.x();
        x[1]=c.D.y();
        x[2]=0;
        points_collisions->SetPoint(nCollisions1*2+i*2+1, x);

        pts2[0]=nCollisions1*2+i*2+0;
        pts2[1]=nCollisions1*2+i*2+1;
        ugrid_collisions->InsertNextCell(VTK_LINE, 2, pts2);
    }
    points_collisions->Modified();
*/

}


void icy::Mesh::DetectContactPairs(double distance_threshold)
{
/*
    indenter_boundary_vs_deformable_nodes.clear();
    deformable_boundary_vs_indenter_nodes.clear();
    collision_interactions.clear();

    for(unsigned b_idx=0; b_idx<indenter.boundary_edges.size(); b_idx++)
        for(unsigned n_idx=0; n_idx<brick.boundary_nodes.size(); n_idx++)
        {
            Interaction i;
            i.ndA_idx = boundary_indenter[b_idx].first;
            i.ndB_idx = boundary_indenter[b_idx].second;
            i.ndP_idx = deformable_boundary_nodes[n_idx];

            i.ndA = &nodes_indenter[i.ndA_idx];
            i.ndB = &nodes_indenter[i.ndB_idx];
            i.ndP = &nodes[i.ndP_idx];

            i.A = i.ndA->xt;
            i.B = i.ndB->xt;
            i.P = i.ndP->xt;

            i.dist = SegmentPointDistance(i.A, i.B, i.P, i.D, i.t);
            if(i.dist <= distance_threshold) {
                indenter_boundary_vs_deformable_nodes.push_back(i);
                collision_interactions.push_back(i);
            }
        }

    for(unsigned b_idx=0; b_idx<boundary.size(); b_idx++)
        for(unsigned n_idx=0; n_idx<nodes_indenter.size(); n_idx++)
        {
            Interaction i;
            i.ndA_idx = boundary[b_idx].first;
            i.ndB_idx = boundary[b_idx].second;
            i.ndP_idx = n_idx;

            i.ndA = &nodes[i.ndA_idx];
            i.ndB = &nodes[i.ndB_idx];
            i.ndP = &nodes_indenter[i.ndP_idx];

            i.A = i.ndA->xt;
            i.B = i.ndB->xt;
            i.P = i.ndP->xt;

            i.dist = SegmentPointDistance(i.A, i.B, i.P, i.D, i.t);
            if(i.dist <= distance_threshold) {
                deformable_boundary_vs_indenter_nodes.push_back(i);
                collision_interactions.push_back(i);
            }
        }

//    qDebug() << "icy::Mesh::DetectContactPairs(): ib_dn " << indenter_boundary_vs_deformable_nodes.size();
//    qDebug() << "icy::Mesh::DetectContactPairs(): db_in " << deformable_boundary_vs_indenter_nodes.size();
*/
}

void icy::Mesh::UpdateValues()
{
    /*
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
*/
}

void icy::Mesh::InitializeLUT(int table)
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
