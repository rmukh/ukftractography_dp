#ifndef BRAINSThreadControl_h
#define BRAINSThreadControl_h

#include <itksys/SystemTools.hxx>
#include <sstream>

#include <itkMacro.h> // needed for ITK_VERSION_MAJOR
#if ITK_VERSION_MAJOR < 5
#include "itkMultiThreader.h"
#else
#include "itkMultiThreaderBase.h"
#endif

#include "ukf_exports.h"

namespace BRAINSUtils
{
/**
 * This class is designed so that
 * the ITK number of threads can be
 * adjusted to a different number of threads
 * during the running of this one call.
 * The desired number of threads
 * must be selected as part of the
 * construction process, and at
 * destruction, the original value is
 * restored.
 *
 * This type of functionality is needed
 * so that the shared libary version of
 * BRAINSFit does not globally change
 * the behavior of all other programs
 * in slicer.
 */
class UKFBASELIB_EXPORTS StackPushITKDefaultNumberOfThreads
{
public:
  explicit StackPushITKDefaultNumberOfThreads(const int desiredCount);
  ~StackPushITKDefaultNumberOfThreads();

protected:
  StackPushITKDefaultNumberOfThreads();                                                // Purposefully not implemented
  StackPushITKDefaultNumberOfThreads &operator=(StackPushITKDefaultNumberOfThreads &); // Purposefully not implemented

private:
  int m_originalThreadValue;
};
} // namespace BRAINSUtils

#endif // BRAINSThreadControl_h