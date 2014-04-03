#pragma once

#include <cstdlib> // exit()
#include "Definitions.h"
#include "Enums.h"

// include the shared infrastructure headers
#include "ApplyToAll.h"
#include "BlockSize.h"
#include "Tuple.h"
#include "TupleElements.h"
#include "TupleAlgorithms.h"
#include "Array.h"
#include "ArrayAlgorithms.h"
#include "JokerDataField.h"
#include "SwapDataField.h"
#include "ScalarStorage.h"
#include "DummyStorage.h"
#include "BlockStorage.h"
#include "ParameterTraits.h"
#include "UsageMeter.h"
#include "IJBoundary.h"
#include "IJKBoundary.h"
#include "KBoundary.h"

// define the field typed depending on the back end
#ifdef __CUDA_BACKEND__

#include "DataFieldCUDA.h"

typedef DataFieldAlignment<cDimI, 2 * cCacheLineSize / sizeof(Real)> CUDARealAlignment;
typedef DataFieldAlignment<cDimI, 2 * cCacheLineSize / sizeof(int)> CUDAIntAlignment;
typedef DataFieldIJBoundary<-cNumBoundaryLines, cNumBoundaryLines, -cNumBoundaryLines, cNumBoundaryLines> CUDAIJBoundary;

// default 3D data field type
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::KJI, CUDARealAlignment> > IJKRealField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::KJI, CUDARealAlignment> > IJKIntField;

// default 2D data field types
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::JI, CUDARealAlignment> > IJRealField;
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::KI, CUDARealAlignment> > IKRealField;
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::KJ, CUDARealAlignment> > JKRealField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::JI, CUDAIntAlignment> > IJIntField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::KI, CUDAIntAlignment> > IKIntField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::KJ, CUDAIntAlignment> > JKIntField;

// default 1D data field types
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::I, CUDARealAlignment> > IRealField;
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::J, CUDARealAlignment> > JRealField;
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::K, CUDARealAlignment> > KRealField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::I, CUDAIntAlignment> > IIntField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::J, CUDAIntAlignment> > JIntField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, StorageOrder::K, CUDAIntAlignment> > KIntField;

// scalar fields
typedef DataFieldCUDA<Real, DataFieldStorageFormat<CUDAIJBoundary, boost::mpl::vector0_c<Dimension>, CUDARealAlignment> > ScalarRealField;
typedef DataFieldCUDA<int, DataFieldStorageFormat<CUDAIJBoundary, boost::mpl::vector0_c<Dimension>, CUDAIntAlignment> > ScalarIntField;

// timer for fine grained timings which can be disabled
#ifdef ENABLE_PERFORMANCE_METERS
#include "TimerCUDA.h"
class PerformanceMeter : public TimerCUDA {};
#else
#include "TimerDummy.h"
class PerformanceMeter : public TimerDummy {};
#endif

// timer which is always enabled
#include "TimerCUDA.h"
class StopWatch : public TimerCUDA {};

#else

#include "DataFieldOpenMP.h"

typedef DataFieldAlignment<cDimK, 1> OpenMPAlignment;
typedef DataFieldIJBoundary<-cNumBoundaryLines, cNumBoundaryLines, -cNumBoundaryLines, cNumBoundaryLines> OpenMPIJBoundary;

// default 3D data field type
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::JIK, OpenMPAlignment> > IJKRealField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::JIK, OpenMPAlignment> > IJKIntField;

// default 2D data field types
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::JI, OpenMPAlignment> > IJRealField;
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::IK, OpenMPAlignment> > IKRealField;
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::JK, OpenMPAlignment> > JKRealField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::JI, OpenMPAlignment> > IJIntField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::IK, OpenMPAlignment> > IKIntField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::JK, OpenMPAlignment> > JKIntField;

// default 1D data field types
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::I, OpenMPAlignment> > IRealField;
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::J, OpenMPAlignment> > JRealField;
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::K, OpenMPAlignment> > KRealField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::I, OpenMPAlignment> > IIntField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::J, OpenMPAlignment> > JIntField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, StorageOrder::K, OpenMPAlignment> > KIntField;

// scalar fields
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<OpenMPIJBoundary, boost::mpl::vector0_c<Dimension>, OpenMPAlignment> > ScalarRealField;
typedef DataFieldOpenMP<int, DataFieldStorageFormat<OpenMPIJBoundary, boost::mpl::vector0_c<Dimension>, OpenMPAlignment> > ScalarIntField;

// timer for fine grained timings which can be disabled
#ifdef ENABLE_PERFORMANCE_METERS
#include "TimerOpenMP.h"
class PerformanceMeter : public TimerOpenMP {};
#else
#include "TimerDummy.h"
class PerformanceMeter : public TimerDummy {};
#endif

// portable timer which is always enabled
#include "TimerOpenMP.h"
class StopWatch : public TimerOpenMP {};

#endif
