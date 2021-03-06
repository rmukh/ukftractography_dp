/**
 * \file vtk_writer.h
 * \brief UKFFiber writing functionality
 *
 * Contains Class definition of the VtkWriter, used for storing the fiber vector to a .vtk file
 * \todo Could be done more elegantly with VTK
*/

#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#include <string>
#include <vector>
#include "linalg.h"
#include "tractography.h"
#include "ukffiber.h"
#include "vtkByteSwap.h"
#include <fstream>
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

class ISignalData;
class vtkPointData;
/**
 * \class VtkWriter
 * \brief Class that allows to write a bunch of fibers to a .vtk file
*/
class VtkWriter
{
public:
  /**
   * \brief Constructor
  */
  VtkWriter(const ISignalData *signal_data, bool write_tensors);

  /** Destructor */
  virtual ~VtkWriter()
  {
  }

  /**
   * \brief Writes the fibers to the VTK file and attaches the selected values to the fiber
   * \param[in] file_name The path of the output fiber file (*.vtk)
   * \param[in] tractsWithSecondTensor File path for the fibers generated with the second tensor
   *                                   This one is optional.
   * \param[in] store_glyphs Write glyphs (i.e. main tensor directions) to a file named glyphs_{tracts}.
   * \return EXIT_FAILURE or EXIT_SUCCESS
  */
  int Write(const std::string &file_name,
            const std::vector<UKFFiber> &fibers,
            bool write_state, bool store_glyphs);

  int WriteWeight(const std::string &file_name, const std::vector<UKFFiber> &fibers);

  /** Write the glyphs (i.e. main tensor directions) to  a file named glyphs_{tracts}.
   * \return EXIT_FAILURE or EXIT_SUCCESS
   */
  int WriteGlyphs(const std::string &file_name, const std::vector<UKFFiber> &fibers);

  /** Sets the variable that toggles the transform from ijk to RAS before writing the fiber to VTK. */
  void set_transform_position(bool transform_position)
  {
    _transform_position = transform_position;
  }
  /** set the WriteBinary flag */
  void SetWriteBinary(bool wb) { this->_writeBinary = wb; }
  void SetWriteCompressed(bool wc) { this->_writeCompressed = wc; }

  /**
   * Writes the fibers and all values attached to them to a PolyData type
  */
  void PopulateFibersAndTensors(vtkPolyData *polyData,
                                const std::vector<UKFFiber> &fibers);

  void PopulateFibersDirs(vtkPolyData *polyData,
                          const std::vector<UKFFiber> &fibers);

protected:
  /**
   * Convert a point from the internal representation into what VTK expects
  */
  vec3_t PointConvert(const vec3_t &point);
  /**
   * Write a single scalar value out in binary.
   */
  template <typename TInput, typename TOutput>
  void WriteX(std::ofstream &output, const TInput toWrite)
  {
    TOutput tmp = toWrite;
    switch (sizeof(TOutput))
    {
    case 1:
      break;
    case 2:
      vtkByteSwap::Swap2BE(&tmp);
      break;
    case 4:
      vtkByteSwap::Swap4BE(&tmp);
      break;
    case 8:
      vtkByteSwap::Swap8BE(&tmp);
      break;
    }
    output.write(reinterpret_cast<const char *>(&tmp), sizeof(TOutput));
    if (output.fail())
    {
      throw;
    }
  }

  void WritePolyData(vtkSmartPointer<vtkPolyData> pd, const char *filename) const;

  /**
   * \brief Reconstructs the tensor from the state for each case
   * \param[out] D The calculated diffusion tensor
   * \todo I think there is something wrong with choosing a orthonormal basis for the tensor
  */
  void State2Tensor(const ukfStateVector &state, mat33_t &D, const int tensorNumber) const;

  /** The diffusion weighted signal data */
  const ISignalData *_signal_data;

  /** If set to true, WritePoint will convert the points from ijk to RAS space before writing */
  bool _transform_position;

  /** Scaling of the glyphs */
  const ukfPrecisionType _scale_glyphs;

  /** Whether to attach the tensors to the fiber */
  bool _write_tensors;

  // Positions of the eigenvalues, directions, or angles, respectively, in the state.
  int _p_l1, _p_l2, _p_l3;
  int _p_m1, _p_m2, _p_m3;

  /** How many indeces in the state are used up by 1 tensor */
  int _tensor_space;

  /** Additional scaling of the eigenvalues before writing */
  const ukfPrecisionType _eigenScaleFactor;

  /** Transformation matrix from ijk-RAS with voxel size normalized out */
  mat33_t _sizeFreeI2R;

  /** is the file to be written binary? */
  bool _writeBinary;
  /** is the file to be written compressed? */
  bool _writeCompressed;
};

#endif // VTK_WRITER_H_
