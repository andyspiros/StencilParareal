#pragma once

/**
* @enum Dimension
* Enum listing all 3 dimensions I, J and K
*/
enum Dimension
{
    cDimI = 0, /**< Dimension I */
    cDimJ = 1, /**< Dimension J */
    cDimK = 2  /**< Dimension K */
};

/**
* @enum Direction
* Enum listing directions
*/
enum Direction
{
    cNoDir,   /**< No direction */
    cPlusDir, /**< Plus direction */
    cMinusDir /**< Minus direction */
};

/**
* @enum ParameterIntend
* Enum listing stencil parameter properties
*/
enum ParameterIntend
{
    cIn,        /**< Read only parameter */
    cInOut,     /**< Read write parameter */
    cScalar    /**< Scalar parameter */
};

/**
* @enum AccessPolicy
* Enum listing the iterator access policies
*/
enum AccessPolicy 
{
    cReadWrite, /**< Read and write accesses to the iterator data are possible */
    cReadOnly   /**< Read accesses to the iterator data are possible */
};

/**
* @enum CachePolicy
* Enum listing the context cache policies
*/
enum CachePolicy 
{
    cBypassCache, /**< Circumvent the cache when accessing context elements, used to fill software managed cache */
    cUseCache     /**< Standard context cache policy reading data through the cache */
};

/**
* @enum CacheIOPolicy
* Enum listing the cache IO policies
*/
enum CacheIOPolicy 
{
    cFillAndFlush,  /**< Read values from the cached field and write the result back */
    cFill,          /**< Read values form the cached field but do not write back */
    cFlush,         /**< Write values back the the cached field but do not read in */
    cLocal          /**< Local only cache, neither read nor write the the cached field */
};

/**
* @enum KLoopDirection
* Enum listing the 2 k loop directions
*/
enum KLoopDirection
{
    cKIncrement = 1,  /**< Loop forward over the column */
    cKDecrement = -1  /**< Loop backward over the column */
};

/**
* @enum KReferencePositions
* Enum listing the reference positions in k direction
* (Note that the numbers represent the order in k)
*/
enum KReferencePositions
{
    cKMinimumFlat = 0,      /**< k = 0  */
    cKMaximumFlat = 1,      /**< k = cFlatLimit-1 */
    cKMinimumTerrain = 2,   /**< k = cFlatLimit */
    cKMaximumTerrain = 3    /**< k = calculationDomain k size - 1 */ 
};

/**
* @enum CornerPolicy
* Enum listing the different policies of handling block boundaries. 
*/
enum CornerPolicy
{
    cComplete,   /**< update the corner regions */
    cIndented   /**< ignore the corners */
};

/**
* @enum BoundaryConditionType
* Enum listing the different boundary condition types
*/
enum BoundaryConditionType
{
    cVoidBoundaryCondition,         /**< Void boundary condition */
    cCopyBoundaryCondition,         /**< Copy boundary condition */ 
    cSwapBoundaryCondition,         /**< Swap boundary condition */ 
    cZeroGradientBoundaryCondition, /**< Zero gradient boundary condition */
    cValueBoundaryCondition,        /**< Value boundary condition */
    cInterpolationBoundaryCondition /**< Interpolation boundary condition */
};

/**
* @enum ImplementationPolicy
* Enum listing different method implementation types
*/
enum ImplementationPolicy
{
    cEmpty,         /**< no work needed implement an empty method */
    cCompileTime,   /**< implementation using compile-time parameters only */
    cRuntime        /**< implementation using runtime and compile-time parameters */
};