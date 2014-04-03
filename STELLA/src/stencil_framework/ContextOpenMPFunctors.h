#pragma once

#include "Enums.h"
#include "DataFieldOpenMPStorageIterator.h"

/**
* @struct IJOffset
* Structure storing an ij offset
*/
struct IJOffset
{
    int i;
    int j;
};

/**
* @struct IJKOffset
* Structure storing an ijk offset
*/
struct IJKOffset
{
    int i;
    int j;
    int k;
};

/**
* @struct AdvanceFunctor
* Functor advancing iterator positions by a constant offset
*/
template<int VI, int VJ, int VK>
struct AdvanceFunctor
{
    template<
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    static void Do(DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator)
    {
        iterator.Advance(VI, VJ, VK);
    }
};

/**
* @struct AdvanceFunctor
* Functor advancing iterator positions by a k offset
*/
template<int VI, int VJ>
struct AdvanceInKFunctor
{
    template<
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    static void Do(
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator,
        parameter_type<const int>::type kOffset)
    {
        iterator.Advance(VI, VJ, kOffset);
    }
};

/**
* @struct RestoreAndAdvanceInIJFunctor
* Functor restoring the iterator positions and advancing the iterators by an ij offset 
*/
template<int VK>
struct RestoreAndAdvanceInIJFunctor
{
    template<
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    static void Do(
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator, 
        parameter_type<const IJOffset>::type ijOffset)
    {
        iterator.RestoreMemorizedPosition();
        iterator.Advance(ijOffset.i, ijOffset.j, VK);
    }
};

/**
* @struct RestoreAndAdvanceInIJKFunctor
* Functor restoring the iterator positions and advancing the iterators by an ijk offset
*/
struct RestoreAndAdvanceInIJKFunctor
{
    template<
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    static void Do(
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator, 
        parameter_type<const IJKOffset>::type ijkOffset)
    {
        iterator.RestoreMemorizedPosition();
        iterator.Advance(ijkOffset.i, ijkOffset.j, ijkOffset.k);
    }
};

/**
* @struct AdvanceMemorizedPositionInIFunctor
* Functor advancing the memorized positions by an i offset
*/
struct AdvanceMemorizedPositionInIFunctor
{
    template<
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    static void Do(
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator, 
        parameter_type<const int>::type iOffset)
    {
        iterator.AdvanceMemorizedPosition(iOffset, 0, 0);
    }
};

/**
* @struct AdvanceMemorizedPositionInJFunctor
* Functor advancing the memorized positions by a j offset
*/
struct AdvanceMemorizedPositionInJFunctor
{
    template<
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    static void Do(
        DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>& iterator, 
        parameter_type<const int>::type jOffset)
    {
        iterator.AdvanceMemorizedPosition(0, jOffset, 0);
    }
};

  
