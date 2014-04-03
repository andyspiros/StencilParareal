#pragma once

/**
* @struct return_type
* Meta function calculating the return type of a stencil parameter
*/
template<typename T>
struct return_type;

/**
* @struct value_type
* Meta function calculating the value type of a stencil parameter
*/
template<typename T>
struct value_type;

/**
* @struct is_iterable
* Meta function returning true if ther parameter is iterable (no scalar parameter)
*/
template<typename T>
struct is_iterable;

