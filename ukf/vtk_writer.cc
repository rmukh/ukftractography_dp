/**
 * \file vtk_writer.cc
 * \brief implementation of vtk_writer.h
 * \todo The case differentiation in the beginning is very hackish..
*/

#include "vtkVersion.h"
#include "ukf_types.h"
#include "vtk_writer.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "ISignalData.h"
#include "utilities.h"
#include "ukffiber.h"
#include "itksys/SystemTools.hxx"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkLine.h"
#include "vtkCellArray.h"
#include "vtkStringArray.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtkDataObject.h"

#include "git_version.h"

VtkWriter::VtkWriter(const ISignalData *signal_data, bool write_tensors) : _signal_data(signal_data),
                                                                           _transform_position(true),
                                                                           _scale_glyphs(0.01),
                                                                           _write_tensors(write_tensors),
                                                                           _eigenScaleFactor(1),
                                                                           _writeBinary(true),
                                                                           _writeCompressed(true)
{
  _p_m1 = 0;
  _p_m2 = 1;
  _p_m3 = 2;

  _p_l1 = 3;
  _p_l2 = 4;
  _p_l3 = 4;

  _tensor_space = 7;

  // this also upon initialization of writer, its the same for all
  ukfMatrixType i2r = _signal_data->i2r();
  vec3_t voxel = _signal_data->voxel();
  // Factor out the effect of voxel size
  _sizeFreeI2R << i2r(0, 0) / voxel[2], i2r(0, 1) / voxel[1], i2r(0, 2) / voxel[0],
      i2r(1, 0) / voxel[2], i2r(1, 1) / voxel[1], i2r(1, 2) / voxel[0],
      i2r(2, 0) / voxel[2], i2r(2, 1) / voxel[1], i2r(2, 2) / voxel[0];
}

void VtkWriter::PopulateFibersAndTensors(vtkPolyData *polyData,
                                         const std::vector<UKFFiber> &fibers)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  size_t num_fibers = fibers.size();
  size_t num_points = 0;

  for (size_t i = 0; i < num_fibers; ++i)
    num_points += fibers[i].position.size();

  for (size_t i = 0; i < num_fibers; ++i)
  {
    size_t fiber_size = fibers[i].position.size();
    for (size_t j = 0; j < fiber_size; ++j)
    {
      vec3_t current = PointConvert(fibers[i].position[j]);
      points->InsertNextPoint(current[0], current[1], current[2]);
    }
  }
  polyData->SetPoints(points);

  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  vtkIdType counter = 0;
  for (size_t i = 0; i < num_fibers; ++i)
  {
    size_t fiber_size = fibers[i].position.size();
    vtkIdType *ids = new vtkIdType[fiber_size];
    for (size_t j = 0; j < fiber_size; ++j)
    {
      ids[j] = counter;
      counter++;
    }
    vtkSmartPointer<vtkLine> curLine = vtkSmartPointer<vtkLine>::New();
    curLine->Initialize(fiber_size, ids, points);
    lines->InsertNextCell(curLine);
    delete[] ids;
    ids = NULL;
  }
  polyData->SetLines(lines);

  /////
  // Dataset attribute part starts
  /////

  if (_write_tensors)
  {
    counter = 0;
    vtkPointData *pointData = polyData->GetPointData();

    mat33_t D;
    for (int local_tensorNumber = 1; local_tensorNumber <= 3; ++local_tensorNumber)
    {
      vtkSmartPointer<vtkFloatArray> curTensor = vtkSmartPointer<vtkFloatArray>::New();
      curTensor->SetNumberOfComponents(9);
      curTensor->Allocate(num_points * 9);
      {
        std::stringstream ss;
        ss << "tensor" << local_tensorNumber;
        curTensor->SetName(ss.str().c_str());
      }
      for (size_t i = 0; i < num_fibers; i++)
      {
        const size_t fiber_size = fibers[i].position.size();
        for (size_t j = 0; j < fiber_size; ++j)
        {
          const ukfStateVector &state = fibers[i].state[j];
          State2Tensor(state, D, local_tensorNumber);
          ukfPrecisionType tmp[9];
          for (unsigned ii = 0, v = 0; ii < 3; ++ii)
          {
            for (unsigned jj = 0; jj < 3; ++jj, ++v)
              tmp[v] = D(ii, jj);
          }
          curTensor->InsertNextTuple(tmp);
        }
      }
      const size_t idx = pointData->AddArray(curTensor);
      pointData->SetActiveAttribute(idx, vtkDataSetAttributes::TENSORS);
    }
  }
}

void VtkWriter::PopulateFibersDirs(vtkPolyData *polyData,
                                   const std::vector<UKFFiber> &fibers)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  size_t num_fibers = fibers.size();
  size_t num_points = 0;
  for (size_t i = 0; i < num_fibers; ++i)
    num_points += fibers[i].position.size();

  for (size_t i = 0; i < num_fibers; ++i)
  {
    size_t fiber_size = fibers[i].position.size();
    for (size_t j = 0; j < fiber_size; ++j)
    {
      vec3_t current = PointConvert(fibers[i].position[j]);
      points->InsertNextPoint(current[0], current[1], current[2]);
    }
  }
  polyData->SetPoints(points);

  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  // Do loop for all fibers
  vtkIdType counter = 0;
  for (size_t i = 0; i < num_fibers; ++i)
  {
    size_t fiber_size = fibers[i].position.size();

    for (size_t h = 0; h < (size_t)(fiber_size / 2); ++h)
    {
      vtkIdType *ids = new vtkIdType[2];
      for (size_t j = 0; j < 2; ++j)
      {
        // each id should be unique
        ids[j] = counter;
        counter++;
      }

      vtkSmartPointer<vtkLine> curLine = vtkSmartPointer<vtkLine>::New();
      curLine->Initialize(2, ids, points);
      lines->InsertNextCell(curLine);
      delete[] ids; // remove indecies for the next fiber
      ids = NULL;
    }
  }
  polyData->SetLines(lines);
}

void VtkWriter::WritePolyData(vtkSmartPointer<vtkPolyData> pd, const char *filename) const
{
  // Add version information to polydata as vtkStringArray
  std::stringstream header_string;
  header_string << "UKF_GIT_HASH:" << UKF_GIT_HASH;
  vtkSmartPointer<vtkStringArray> version_info = vtkSmartPointer<vtkStringArray>::New();
  version_info->SetNumberOfValues(1);
  version_info->SetValue(0, header_string.str());
  version_info->SetName("UKF_VERSION_INFO");
  pd->GetFieldData()->AddArray(version_info);

  // Output filename extension
  const std::string ext(itksys::SystemTools::GetFilenameExtension(filename));

  if (ext == ".vtp")
  {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    if (this->_writeBinary)
      writer->SetDataModeToBinary();
    else
      writer->SetDataModeToAscii();
  
    if (this->_writeCompressed)
      writer->SetCompressorTypeToZLib();
    else
      writer->SetCompressorTypeToNone();

#if (VTK_MAJOR_VERSION < 6)
    writer->SetInput(pd);
#else
    writer->SetInputData(pd);
#endif
    writer->SetFileName(filename);
    writer->Write();
  }
  else
  {
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    if (this->_writeBinary)
      writer->SetFileTypeToBinary();

#if (VTK_MAJOR_VERSION < 6)
    writer->SetInput(pd);
#else
    writer->SetInputData(pd);
#endif
    writer->SetFileName(filename);

    writer->Write();
  }
}

int VtkWriter::Write(const std::string &file_name,
                     const std::vector<UKFFiber> &fibers,
                     bool write_state,
                     bool store_glyphs)
{
  if (fibers.size() == 0)
  {
    std::cout << "No fiber exists." << std::endl;
    return EXIT_FAILURE;
  }

  //
  // Handle glyps
  if (store_glyphs)
  {
    std::stringstream ss;
    ss << file_name.substr(0, file_name.find_last_of(".")) << "_glyphs" << ".vtk";
    if (WriteGlyphs(ss.str(), fibers) == EXIT_FAILURE)
      return EXIT_FAILURE;
  }
  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  // handle fibers and tensors
  this->PopulateFibersAndTensors(polyData, fibers);

  // norm, fa etc hung as arrays on the point data for the polyData
  // object.
  vtkPointData *pointData = polyData->GetPointData();

  int num_fibers = fibers.size();
  int num_points = 0;
  for (int i = 0; i < num_fibers; ++i)
    num_points += fibers[i].position.size();

  // write norm
  {
    vtkSmartPointer<vtkFloatArray> norms = vtkSmartPointer<vtkFloatArray>::New();
    norms->SetNumberOfComponents(1);
    norms->Allocate(num_points);
    norms->SetName("EstimatedUncertainty");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        norms->InsertNextValue(fibers[i].norm[j]);
    }
    int idx = pointData->AddArray(norms);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  // write rtop1
  if (fibers[0].rtop1.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> fa = vtkSmartPointer<vtkFloatArray>::New();
    fa->SetNumberOfComponents(1);
    fa->Allocate(num_points);
    fa->SetName("RTOP1");

    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        fa->InsertNextValue(fibers[i].rtop1[j]);
    }
    int idx = pointData->AddArray(fa);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  // rtop2
  if (fibers[0].rtop2.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> fa2 = vtkSmartPointer<vtkFloatArray>::New();
    fa2->SetName("RTOP2");
    fa2->SetNumberOfComponents(1);
    fa2->Allocate(num_points);
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        fa2->InsertNextValue(fibers[i].rtop2[j]);
    }
    int idx = pointData->AddArray(fa2);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  // rtop3
  if (fibers[0].rtop3.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> fa3 = vtkSmartPointer<vtkFloatArray>::New();
    fa3->SetName("RTOP3");
    fa3->SetNumberOfComponents(1);
    fa3->Allocate(num_points);
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        fa3->InsertNextValue(fibers[i].rtop3[j]);
    }
    int idx = pointData->AddArray(fa3);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  // rtop model
  if (fibers[0].rtop_model.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> trace = vtkSmartPointer<vtkFloatArray>::New();
    trace->SetNumberOfComponents(1);
    trace->Allocate(num_points);
    trace->SetName("RTOP_model");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        trace->InsertNextValue(fibers[i].rtop_model[j]);
    }
    int idx = pointData->AddArray(trace);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  // rtop signal
  if (fibers[0].rtop_signal.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> trace2 = vtkSmartPointer<vtkFloatArray>::New();
    trace2->SetNumberOfComponents(1);
    trace2->Allocate(num_points);
    trace2->SetName("RTOP_signal");

    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        trace2->InsertNextValue(fibers[i].rtop_signal[j]);
    }
    int idx = pointData->AddArray(trace2);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].free_water.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> free_water = vtkSmartPointer<vtkFloatArray>::New();
    free_water->SetNumberOfComponents(1);
    free_water->Allocate(num_points);
    free_water->SetName("FreeWater");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        free_water->InsertNextValue(fibers[i].free_water[j]);
    }
    int idx = pointData->AddArray(free_water);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].w1.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> w1 = vtkSmartPointer<vtkFloatArray>::New();
    w1->SetNumberOfComponents(1);
    w1->Allocate(num_points);
    w1->SetName("w1");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        w1->InsertNextValue(fibers[i].w1[j]);
    }
    int idx = pointData->AddArray(w1);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].w2.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> w2 = vtkSmartPointer<vtkFloatArray>::New();
    w2->SetNumberOfComponents(1);
    w2->Allocate(num_points);
    w2->SetName("w2");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        w2->InsertNextValue(fibers[i].w2[j]);
    }
    int idx = pointData->AddArray(w2);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].w3.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> w3 = vtkSmartPointer<vtkFloatArray>::New();
    w3->SetNumberOfComponents(1);
    w3->Allocate(num_points);
    w3->SetName("w3");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        w3->InsertNextValue(fibers[i].w3[j]);
    }
    int idx = pointData->AddArray(w3);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].w1w2angle.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> w1w2angle = vtkSmartPointer<vtkFloatArray>::New();
    w1w2angle->SetNumberOfComponents(1);
    w1w2angle->Allocate(num_points);
    w1w2angle->SetName("d1_d2_angle");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        w1w2angle->InsertNextValue(fibers[i].w1w2angle[j]);
    }
    int idx = pointData->AddArray(w1w2angle);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].w1w3angle.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> w1w3angle = vtkSmartPointer<vtkFloatArray>::New();
    w1w3angle->SetNumberOfComponents(1);
    w1w3angle->Allocate(num_points);
    w1w3angle->SetName("d1_d3_angle");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        w1w3angle->InsertNextValue(fibers[i].w1w3angle[j]);
    }
    int idx = pointData->AddArray(w1w3angle);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].Fm1.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> Fm1 = vtkSmartPointer<vtkFloatArray>::New();
    Fm1->SetNumberOfComponents(1);
    Fm1->Allocate(num_points);
    Fm1->SetName("Frob dir 1");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        Fm1->InsertNextValue(fibers[i].Fm1[j]);
    }
    int idx = pointData->AddArray(Fm1);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].lmd1.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> lmd1 = vtkSmartPointer<vtkFloatArray>::New();
    lmd1->SetNumberOfComponents(1);
    lmd1->Allocate(num_points);
    lmd1->SetName("Frob lambda 1");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        lmd1->InsertNextValue(fibers[i].lmd1[j]);
    }
    int idx = pointData->AddArray(lmd1);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].Fm2.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> Fm2 = vtkSmartPointer<vtkFloatArray>::New();
    Fm2->SetNumberOfComponents(1);
    Fm2->Allocate(num_points);
    Fm2->SetName("Frob dir 2");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        Fm2->InsertNextValue(fibers[i].Fm2[j]);
    }
    int idx = pointData->AddArray(Fm2);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].lmd2.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> lmd2 = vtkSmartPointer<vtkFloatArray>::New();
    lmd2->SetNumberOfComponents(1);
    lmd2->Allocate(num_points);
    lmd2->SetName("Frob lambda 2");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        lmd2->InsertNextValue(fibers[i].lmd2[j]);
    }
    int idx = pointData->AddArray(lmd2);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].Fm3.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> Fm3 = vtkSmartPointer<vtkFloatArray>::New();
    Fm3->SetNumberOfComponents(1);
    Fm3->Allocate(num_points);
    Fm3->SetName("Frob dir 3");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        Fm3->InsertNextValue(fibers[i].Fm3[j]);
    }
    int idx = pointData->AddArray(Fm3);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].lmd3.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> lmd3 = vtkSmartPointer<vtkFloatArray>::New();
    lmd3->SetNumberOfComponents(1);
    lmd3->Allocate(num_points);
    lmd3->SetName("Frob lambda 3");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        lmd3->InsertNextValue(fibers[i].lmd3[j]);
    }
    int idx = pointData->AddArray(lmd3);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].varW1.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> varW1 = vtkSmartPointer<vtkFloatArray>::New();
    varW1->SetNumberOfComponents(1);
    varW1->Allocate(num_points);
    varW1->SetName("var w1");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        varW1->InsertNextValue(fibers[i].varW1[j]);
    }
    int idx = pointData->AddArray(varW1);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].varW2.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> varW2 = vtkSmartPointer<vtkFloatArray>::New();
    varW2->SetNumberOfComponents(1);
    varW2->Allocate(num_points);
    varW2->SetName("var w2");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        varW2->InsertNextValue(fibers[i].varW2[j]);
    }
    int idx = pointData->AddArray(varW2);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].varW3.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> varW3 = vtkSmartPointer<vtkFloatArray>::New();
    varW3->SetNumberOfComponents(1);
    varW3->Allocate(num_points);
    varW3->SetName("var w3");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        varW3->InsertNextValue(fibers[i].varW3[j]);
    }
    int idx = pointData->AddArray(varW3);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].varWiso.size() > 0)
  {
    vtkSmartPointer<vtkFloatArray> varWiso = vtkSmartPointer<vtkFloatArray>::New();
    varWiso->SetNumberOfComponents(1);
    varWiso->Allocate(num_points);
    varWiso->SetName("var wiso");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
        varWiso->InsertNextValue(fibers[i].varWiso[j]);
    }
    int idx = pointData->AddArray(varWiso);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (fibers[0].normMSE.size() > 0)
  {
    ukfPrecisionType nmse_sum(0);
    unsigned counter(0);
    vtkSmartPointer<vtkFloatArray> normMSE = vtkSmartPointer<vtkFloatArray>::New();
    normMSE->SetNumberOfComponents(1);
    normMSE->Allocate(num_points);
    normMSE->SetName("NormalizedSignalEstimationError");
    for (int i = 0; i < num_fibers; ++i)
    {
      int fiber_size = fibers[i].position.size();
      for (int j = 0; j < fiber_size; ++j)
      {
        normMSE->InsertNextValue(fibers[i].normMSE[j]);
        nmse_sum += fibers[i].normMSE[j];
        ++counter;
      }
    }
    int idx = pointData->AddArray(normMSE);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    std::cout << "nmse_avg=" << nmse_sum / counter << std::endl;
  }

  if (write_state)
  {
    int state_dim = fibers[0].state[0].size();
    vtkSmartPointer<vtkFloatArray> stateArray = vtkSmartPointer<vtkFloatArray>::New();
    stateArray->SetNumberOfComponents(state_dim);
    stateArray->Allocate(num_points);
    stateArray->SetName("state");

    float *tmpArray = new float[state_dim];

    for (int i = 0; i < num_fibers; i++)
    {
      const int fiber_size = static_cast<int>(fibers[i].position.size());
      for (int j = 0; j < fiber_size; ++j)
      {
        const ukfStateVector &state = fibers[i].state[j];
        for (int k = 0; k < state_dim; ++k)
          tmpArray[k] = state[k];
        stateArray->InsertNextTuple(tmpArray);
      }
    }
    int idx = pointData->AddArray(stateArray);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    delete[] tmpArray;
    tmpArray = NULL;
  }

  if (fibers[0].covariance.size() > 0)
  {
    int state_dim = fibers[0].state[0].size();

    int cov_dim = (state_dim * (state_dim + 1)) / 2;

    vtkSmartPointer<vtkFloatArray> covarianceArray = vtkSmartPointer<vtkFloatArray>::New();
    covarianceArray->SetNumberOfComponents(cov_dim);
    covarianceArray->Allocate(num_points);
    covarianceArray->SetName("covariance");

    float *tmpArray = new float[cov_dim];

    for (int i = 0; i < num_fibers; i++)
    {
      const int fiber_size = static_cast<int>(fibers[i].position.size());
      for (int j = 0; j < fiber_size; ++j)
      {
        int covIndex = 0;
        for (int a = 0; a < state_dim; ++a)
        {
          for (int b = a; b < state_dim; ++b)
          {
            tmpArray[covIndex] = fibers[i].covariance[j](a, b);
            ++covIndex;
          }
        }
        covarianceArray->InsertNextTuple(tmpArray);
      }
    }
    int idx = pointData->AddArray(covarianceArray);
    pointData->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    delete[] tmpArray;
    tmpArray = NULL;
  }

  WritePolyData(polyData, file_name.c_str());
  return EXIT_SUCCESS;
}

int VtkWriter::WriteWeight(const std::string &file_name, const std::vector<UKFFiber> &fibers)
{
  if (fibers.size() == 0)
  {
    std::cout << "No fiber exists." << std::endl;
    return EXIT_FAILURE;
  }

  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  // handle fibers and tensors
  this->PopulateFibersDirs(polyData, fibers);

  WritePolyData(polyData, file_name.c_str());
  return EXIT_SUCCESS;
}

int VtkWriter::WriteGlyphs(const std::string &file_name,
                           const std::vector<UKFFiber> &fibers)
{
  if (fibers.size() == 0)
  {
    std::cout << "No fibers existing." << std::endl;
    return EXIT_FAILURE;
  }
  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  int num_fibers = fibers.size();
  int num_points = 0;
  for (int i = 0; i < num_fibers; ++i)
    num_points += fibers[i].position.size();

  const ukfPrecisionType scale = ukfHalf * _scale_glyphs;

  for (int i = 0; i < num_fibers; ++i)
  {
    int fiber_size = fibers[i].position.size();
    for (int j = 0; j < fiber_size; ++j)
    {
      vec3_t point = fibers[i].position[j];
      const ukfStateVector &state = fibers[i].state[j];

      // Get the directions.
      vec3_t m1 = vec3_t::Zero();
      vec3_t m2 = vec3_t::Zero();
      vec3_t m3 = vec3_t::Zero();

      // For BiExp diffusion propagator (check if is correct)
      m1 << state[2], state[1], state[0];
      m2 << state[9], state[8], state[7];
      m3 << state[16], state[15], state[14];
      m1 = state[3] / 100.0 * m1;
      m2 = state[10] / 100.0 * m2;
      m3 = state[17] / 100.0 * m3;

      // Calculate the points. The glyphs are represented as two-point lines.
      vec3_t pos1, pos2;
      pos1 = point - scale * m1;
      pos2 = point + scale * m1;
      points->InsertNextPoint(pos1[0], pos1[1], pos1[2]);
      points->InsertNextPoint(pos2[0], pos2[1], pos2[2]);

      pos1 = point - scale * m2;
      pos2 = point + scale * m2;
      points->InsertNextPoint(pos1[0], pos1[1], pos1[2]);
      points->InsertNextPoint(pos2[0], pos2[1], pos2[2]);

      pos1 = point - scale * m3;
      pos2 = point + scale * m3;
      points->InsertNextPoint(pos1[0], pos1[1], pos1[2]);
      points->InsertNextPoint(pos2[0], pos2[1], pos2[2]);
    }
  }
  
  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  for (int i = 0; i < num_points * 3; ++i)
  {
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    vtkIdType ids[2];
    ids[0] = i * 2;
    ids[1] = ids[0] + 1;
    line->Initialize(2, ids, points);
    lines->InsertNextCell(line);
  }
  polyData->SetLines(lines);
  WritePolyData(polyData, file_name.c_str());
  return EXIT_SUCCESS;
}

vec3_t VtkWriter::PointConvert(const vec3_t &point)
{
  vec3_t rval;
  ukfVectorType p(4);
  p[0] = point[2]; // NOTICE the change of order here. Flips back to the original axis order
  p[1] = point[1];
  p[2] = point[0];

  if (_transform_position)
  {
    p[3] = ukfOne;
    ukfVectorType p_new(4);
    p_new = _signal_data->i2r() * p; // ijk->RAS transform
    rval[0] = p_new[0];
    rval[1] = p_new[1];
    rval[2] = p_new[2];
  }
  else
  {
    rval[0] = p[2];
    rval[1] = p[1];
    rval[2] = p[0];
  }
  return rval;
}

void VtkWriter::State2Tensor(const ukfStateVector &state, mat33_t &D, const int tensorNumber) const
{
  vec3_t eigenVec1;
  const int tensorIndex = (tensorNumber - 1);
  const size_t after_tensor_offset_index = _tensor_space * (tensorIndex);

  eigenVec1 << state[after_tensor_offset_index + _p_m1],
      state[after_tensor_offset_index + _p_m2],
      state[after_tensor_offset_index + _p_m3];

  // Perform ijk->RAS transform on eigen vectors
  eigenVec1 = _sizeFreeI2R * eigenVec1;

  // Renormalize eigenvectors
  eigenVec1.normalize();

  // Compute the diffusion matrix in RAS coordinate system
  // The transformed matrix is still positive-definite
  const float L1 = state[after_tensor_offset_index + _p_l1];
  const float L2 = state[after_tensor_offset_index + _p_l2];
#if 1
  D = diffusion_l2eql3(eigenVec1, L1, L2);
#else
  diagmat3_t lambdas;
  lambdas.diagonal()[0] = L1;
  lambdas.diagonal()[1] = L2;
  lambdas.diagonal()[2] = lambdas.diagonal()[1];
  D = diffusion(eigenVec1, lambdas);
#endif
}
