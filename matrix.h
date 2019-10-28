//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//  Copyright (c) 2016-2019 Leonardo Consoni <leonardojc@protonmail.com>            //
//                                                                                  //
//  This file is part of Simple Matrix.                                             //
//                                                                                  //
//  Simple Matrix is free software: you can redistribute it and/or modify           //
//  it under the terms of the GNU Lesser General Public License as published        //
//  by the Free Software Foundation, either version 3 of the License, or            //
//  (at your option) any later version.                                             //
//                                                                                  //
//  Simple Matrix is distributed in the hope that it will be useful,                //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of                  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                    //
//  GNU Lesser General Public License for more details.                             //
//                                                                                  //
//  You should have received a copy of the GNU Lesser General Public License        //
//  along with Simple Matrix. If not, see <http://www.gnu.org/licenses/>.           //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////


/// @file matrix.h
/// @brief Matrix data structure and operations abstractions

#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>

#define MATRIX_SIZE_MAX (50 * 50)   ///< Maximum allowed matrix number of elements (rows x columns)

#define MATRIX_IDENTITY 'I'         ///< Create square matrix as identity type (main diagonal filled with 1's)
#define MATRIX_ZERO '0'             ///< Create square matrix as zero type (completely zeroed)

#define MATRIX_TRANSPOSE 'T'        ///< Transpose matrix before multiplication
#define MATRIX_KEEP 'N'             ///< Keep matrix unadulterated before multiplication


typedef struct _MatrixData MatrixData;    ///< Matrix internal data structure
typedef MatrixData* Matrix;               ///< Opaque reference to Matrix data structure


/// @brief Creates matrix with specified values and dimensions                                               
/// @param[in] data array with values in row-major order to fill matrix data (NULL for filling with zeros)                                 
/// @param[in] rowsNumber number of rows                                         
/// @param[in] columnsNumber number of columns      
/// @return reference/pointer to allocated and filled matrix (NULL if dimensions are bigger than MATRIX_SIZE_MAX)
Matrix Mat_Create( double* data, size_t rowsNumber, size_t columnsNumber );     

/// @brief Creates square matrix of specified size and type                              
/// @param[in] size size/order of the square matrix (equal number of rows and cells)
/// @param[in] type defines if internal data is filled as zero (MATRIX_ZERO) or identity (MATRIX_IDENTITY) matrix       
/// @return reference/pointer to allocated and filled matrix (NULL if size is bigger than MATRIX_SIZE_MAX)
Matrix Mat_CreateSquare( size_t size, char type );

/// @brief Destroys/deallocates memory of matrix 
/// @param[in] matrix reference to matrix to be destroyed/deallocated
void Mat_Discard( Matrix matrix );
                                                                      
/// @brief Copies content from one matrix to another, previously allocated  
/// @param[in] source reference to matrix from which data will be copied
/// @param[in] destination matrix to which data will be copied
/// @return reference/pointer to destination matrix (NULL on errors)
Matrix Mat_Copy( Matrix source, Matrix destination );

/// @brief Sets all elements of given matrix to zero                             
/// @param[in] matrix reference to matrix to be cleared/zeroed
/// @return reference/pointer to cleared matrix
Matrix Mat_Clear( Matrix matrix );

/// @brief Gets columns number for given matrix                   
/// @param[in] matrix reference to matrix
/// @return number of columns for the matrix (0 on errors)
size_t Mat_GetWidth( Matrix matrix );

/// @brief Gets rows number for given matrix                        
/// @param[in] matrix reference to matrix
/// @return number of rows for the matrix (0 on errors)
size_t Mat_GetHeight( Matrix matrix );

/// @brief Gets value of given matrix element at specified position                              
/// @param[in] matrix reference to matrix
/// @param[in] row row position of accessed element                             
/// @param[in] column column position of accessed element
/// @return value of specified element (0 on error)
double Mat_GetElement( Matrix matrix, size_t row, size_t column );
                                                                            
/// @brief Sets value of given matrix element at specified position 
/// @param[in] matrix reference to matrix
/// @param[in] row row position of accessed element                            
/// @param[in] column column position of accessed element                                            
/// @param[in] value new value of updated element
void Mat_SetElement( Matrix matrix, size_t row, size_t column, double value );

/// @brief Gets value of given matrix element at specified position                         
/// @param[in] matrix reference to matrix
/// @param[out] buffer reference
/// @return pointer to filled buffer (NULL on errors)
double* Mat_GetData( Matrix matrix, double* buffer );

/// @brief Sets given matrix values from row-major order data array                          
/// @param[in] matrix reference to matrix
/// @param[in] data row-major order data array for filling the matrix
void Mat_SetData( Matrix matrix, double* data );

/// @brief Resizes/reallocates given matrix to specified dimensions. Fill new space with zeros when growing                    
/// @param[in] matrix reference to matrix to be resized
/// @param[in] rowsNumber new number of rows
/// @param[in] columnsNumber new number of columns
/// @return reference/pointer to resized/reallocated matrix (NULL on errors)
Matrix Mat_Resize( Matrix matrix, size_t rowsNumber, size_t columnsNumber );

/// @brief Multiply all given matrix elements by a specified value                              
/// @param[in] matrix reference to matrix to be scaled
/// @param[in] factor scalar value that multiplies the matrix
/// @param[in] result preallocated matrix to store the scaling result (can be the same as the input one)
/// @return reference/pointer to @a result scaled matrix (NULL on errors)
Matrix Mat_Scale( Matrix matrix, double factor, Matrix result );

/// @brief Defines desired internal count for given semaphore
/// @param[in] matrix_1 reference to first matrix
/// @param[in] weight_1 weight of first matrix on sum
/// @param[in] matrix_2 reference to second matrix
/// @param[in] weight_2 weight of second matrix on sum
/// @param[in] result preallocated matrix to store the sum result
/// @return reference/pointer to sum @a result matrix (NULL on errors)
Matrix Mat_Sum( Matrix matrix_1, double weight_1, Matrix matrix_2, double weight_2, Matrix result );

/// @brief Perform dot product/multiplication of 2 matrices
/// @param[in] matrix_1 reference to first matrix (nxk dimensions after transformation)
/// @param[in] trans_1 defines transformation applied to first matrix (MATRIX_TRANSPOSE or MATRIX_KEEP) 
/// @param[in] matrix_2 reference to second matrix (kxm dimensions after transformation) 
/// @param[in] trans_2 defines transformation applied to second matrix (MATRIX_TRANSPOSE or MATRIX_KEEP) 
/// @param[in] result preallocated matrix to store the dot product/multiplication result (nxm dimensions)
/// @return reference/pointer to multiplication @a result matrix (NULL on errors)
Matrix Mat_Dot( Matrix matrix_1, char trans_1, Matrix matrix_2, char trans_2, Matrix result );

/// @brief Calculates determinant of given matrix
/// @param[in] matrix reference to matrix
/// @return determinant value (0.0 on errors)
double Mat_Determinant( Matrix matrix );

/// @brief Transposes given matrix
/// @param[in] matrix reference to matrix to be trasnposed (nxm dimensions)
/// @param[in] result preallocated matrix (n*m data size) to store the transposition result (can be the same as the input one)
/// @return reference/pointer to transposed @a result matrix (NULL on errors)
Matrix Mat_Transpose( Matrix matrix, Matrix result );

/// @brief Calculate inverse of given square matrix                             
/// @param[in] matrix reference to matrix to be inverted (should be square)
/// @param[in] result preallocated matrix to store the inversion result (can be the same as the input one)
/// @return reference/pointer to inverted @a result matrix (NULL on errors)
Matrix Mat_Inverse( Matrix matrix, Matrix result );

/// @brief Print given matrix element values in a formatted way                             
/// @param[in] matrix reference to matrix to be displayed
void Mat_Print( Matrix matrix );

#endif // MATRICES_H
