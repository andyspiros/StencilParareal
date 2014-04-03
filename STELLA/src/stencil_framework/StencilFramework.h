#pragma once

// file including all headers necessary for stencil and update function definition
#include "SharedInfrastructure.h"

// includes needed for the definition of the update functions
#include "DataParameter.h"
#include "FunctionParameter.h"
#include "DataParameterEvaluator.h"
#include "FunctionParameterEvaluator.h"
#include "StencilStageMacros.h"
#include "StencilFunctionMacros.h"
#include "Call.h"

// includes needed for the stencil initialization
#include "StencilConfiguration.h"
#include "PackParameters.h"
#include "ParameterWrapper.h"
#include "DefineTemporaries.h"
#include "DefineLoops.h"
#include "DefineSwitch.h"
#include "DefineCase.h"
#include "DefineCaches.h"
#include "DefineSweep.h"
#include "DefineStages.h"
#include "StencilCompiler.h"



  
