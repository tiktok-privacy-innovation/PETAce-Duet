// Copyright 2023 TikTok Pte. Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <memory>

#include "gmp.h"
#include "gmpxx.h"

#include "duet/beaver_triple/boolean_triple.h"
#include "duet/beaver_triple/triple_factory.h"
#include "duet/paillier/paillier.h"
#include "duet/share_translation/st_generator.h"
#include "duet/util/defines.h"
#include "duet/util/io.h"
#include "duet/util/matrix.h"
#include "duet/util/prng.h"

namespace petace {
namespace duet {

/**
 * @brief Two party MPC framework.
 *
 * @see Refer to README.md for more details
 *
 * @par Example
 * Refer to example.cpp.
 */
class Duet {
public:
    Duet() = delete;

    Duet(const std::shared_ptr<network::Network>& net, std::size_t party_id);

    Duet(const Duet&) = delete;

    Duet& operator=(const Duet&) = delete;

    ~Duet() = default;

    /**
     * @brief Secret share a private matrix to both parties, where the private matrix contains double values.
     *
     * One party provides a private matrix, and the output is that each of both parties holds
     * additive share of the matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input private matrix. (non-sender party should provide an empty matrix)
     * @param[out] out The output additive share matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void share(const std::shared_ptr<network::Network>& net, const PrivateMatrix<double>& in, ArithMatrix& out) const;

    /**
     * @brief Secret share a private matrix to both parties, where the private matrix contains int64_t values.
     *
     * One party provides a private matrix, and the output is that each of both parties holds
     * additive share of the matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input private matrix. (non-sender party should provide an empty matrix)
     * @param[out] out The output additive share matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void share(const std::shared_ptr<network::Network>& net, const PrivateMatrixInt64& in, ArithMatrix& out) const;

    /**
     * @brief Secret share a private bool matrix to both parties, where the private matrix contains int64_t values.
     *
     * One party provides a private matrix, and the output is that each of both parties holds
     * xor share of the matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input private bool matrix. (non-sender party should provide an empty matrix)
     * @param[out] out The output xor share matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void share(const std::shared_ptr<network::Network>& net, const PrivateMatrixBool& in, BoolMatrix& out) const;

    /**
     * @brief Secret share a private bool matrix to both parties, where the private matrix contains int64_t values.
     *
     * One party provides a private matrix, and the output is that each of both parties holds
     * xor share of the matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input private bool matrix. (non-sender party should provide an empty matrix)
     * @param[out] out The output xor share matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void share(const std::shared_ptr<network::Network>& net, const PrivateMatrixBool& in, BoolMatrix& out);

    /**
     * @brief Reconstruct additive-secret-shared matrix to reveal the private matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party_id The party id of the private matrix receiver.
     * @param[in] in The input secret shared matrix.
     * @param[out] out The output private matrix. (non-receiver party should provide an empty matrix)
     * @throws std::invalid_argument if in.size() == 0.
     */
    void reveal(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PrivateMatrix<double>& out) const;

    /**
     * @brief Reconstruct additive-secret-shared matrix to reveal the public matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input secret shared matrix.
     * @param[out] out The output public matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void reveal(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PublicMatrix<double>& out) const;

    /**
     * @brief Reconstruct boolean-secret-shared matrix to reveal the private matrix.
     *
     * xor of two secret shared boolean values reveals the private boolean value.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party The party id of the private matrix receiver.
     * @param[in] in The input secret shared matrix.
     * @param[out] out The output private matrix. (non-receiver party should provide an empty matrix)
     * @throws std::invalid_argument if in.size() == 0.
     */
    void reveal(const std::shared_ptr<network::Network>& net, const BoolMatrix& in, PrivateMatrixBool& out) const;

    /**
     * @brief Reconstruct boolean-secret-shared matrix to reveal the public matrix.
     *
     * xor of two secret shared boolean values reveals the private boolean value.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party The party id of the private matrix receiver.
     * @param[in] in The input secret shared matrix.
     * @param[out] out The output private matrix. (non-receiver party should provide an empty matrix)
     * @throws std::invalid_argument if in.size() == 0.
     */
    void reveal(const std::shared_ptr<network::Network>& net, const BoolMatrix& in, PublicMatrixBool& out) const;

    /**
     * @brief Addition of two additively shared matrix
     *
     * z = x + y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input additively secret-shared matrix
     * @param[out] z The output matrix of the addition result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void add(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const;

    /**
     * @brief Addition of one additively shared matrix and public matrix (double)
     *
     * z = x + y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input double-precision public matrix
     * @param[out] z The output matrix of the addition result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void add(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const;

    /**
     * @brief Subtraction between an arith share matrix and a public double number.
     *
     * z = x + y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input public double number
     * @param[out] z The output matrix of the subtraction result
     */
    void add(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) const;

    /**
     * @brief Subtraction of two additively shared matrix
     *
     * z = x - y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input additively secret-shared matrix
     * @param[out] z The output matrix of the subtraction result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void sub(const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const;

    /**
     * @brief Subtraction of one additively shared matrix and public matrix (double)
     *
     * z = x - y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input double-precision public matrix
     * @param[out] z The output matrix of the subtraction result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void sub(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const;

    /**
     * @brief Subtraction between an arith share matrix and a public double number.
     *
     * z = x - y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input public double number
     * @param[out] z The output matrix of the subtraction result
     */
    void sub(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) const;

    /**
     * @brief Batched AND gate evaluation. Input matrices are boolean secret shared matrices.
     *
     * z = x AND y
     * Each element in matrix contains 64 boolean shares.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input boolean secret-shared matrix
     * @param[in] y The second input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing AND result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_bool_and(const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief Batched AND gate evaluation. Input matrices are boolean secret shared matrices.
     *
     * z = x AND y
     * Each element in matrix contains 64 boolean shares.
     *
     * @param[in] x The first input boolean public matrix
     * @param[in] y The second input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing AND result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_bool_and(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) const;

    /**
     * @brief OR gate evaluation. Input matrices are a boolean secret shared matrices and a public boolean matrix.
     *
     * z = x OR y
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input public boolean matrix
     * @param[in] y The second input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing OR result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_bool_or(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) const;

    /**
     * @brief OR gate evaluation. Input two matrices are a boolean secret shared matrices.
     *
     * z = x or y
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input boolean secret-shared matrix
     * @param[in] y The second input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing OR result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_bool_or(const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief XOR gate evaluation. Input two matrices are a boolean secret shared matrices.
     *
     * z = x XOR y
     *
     * @param[in] x The first input boolean secret-shared matrix
     * @param[in] y The second input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing XOR result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_bool_xor(const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z) const;

    /**
     * @brief XOR gate evaluation. Input two matrices are a boolean secret shared matrices.
     *
     * z = x xor y
     *
     * @param[in] x The first input boolean public matrix
     * @param[in] y The second input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing XOR result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_bool_xor(const PublicMatrixBool& x, const BoolMatrix& y, BoolMatrix& z) const;

    /**
     * @brief not gate evaluation. Input one matrices are a boolean secret shared matrices.
     *
     * z = not x
     *
     * @param[in] x The input boolean secret-shared matrix
     * @param[out] z The output boolean matrix storing NOT result
     * @throws std::invalid_argument if x.size() == 0.
     */
    void elementwise_bool_not(const BoolMatrix& x, BoolMatrix& z) const;

    /**
     * @brief Multiplication of two arithmetic share matrices.
     *
     * The operation is element-wise, e.g., z[i] = x[i] * y[i]
     * NOTE: this is not matrix multiplication.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[in] z The output multiplication result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_mul(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            ArithMatrix& z) const;

    /**
     * @brief Multiplication of an arithmetic share matrix and a public matrix.
     *
     * The operation is element-wise, e.g., z[i] = x[i] * y[i]
     * NOTE: this is not matrix multiplication.
     *
     * @param[in] x The first input public matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z The output of the multiplication result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_mul(const PublicMatrix<double>& x, const ArithMatrix& y, ArithMatrix& z) const;

    /**
     * @brief Multiplication of an arithmetic share matrix and a public number.
     *
     * The operation is element-wise, e.g., z[i] = x * y[i]
     * NOTE: this is not matrix multiplication.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface)
     * @param[in] x The first input public double number
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z The output of the multiplication result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_mul(PublicDouble x, const ArithMatrix& y, ArithMatrix& z) const;

    /**
     * @brief Division between an arith share matrix and a public double matrix.
     *
     * The operation is element-wise, e.g., z[i] = x[i] / y[i]
     *
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input public double matrix
     * @param[in] z The output division result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_div(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const;

    /**
     * @brief Subtraction between an arith share matrix and a public double number.
     *
     * z = x / y
     *
     * @param[in] x The first input additively secret-shared matrix
     * @param[in] y The second input public double number
     * @param[out] z The output matrix of the subtraction result
     */
    void elementwise_div(const ArithMatrix& x, PublicDouble y, ArithMatrix& z) const;

    /**
     * @brief Division of two arithmetic share matrices.
     *
     * The operation is element-wise, e.g., z[i] = x[i] / y[i]
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[in] z The output division result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void elementwise_div(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            ArithMatrix& z) const;

    /**
     * @brief Matrix multiplication of two arithmetic share matrices.
     *
     * The operation is element-wise, e.g., z = x * y
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[in] z The output multiplication result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void mat_mul(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            ArithMatrix& z) const;

    /**
     * @brief Matrix multiplication of public matrix arithmetic share matrices.
     *
     * The operation is element-wise, e.g., z = x * y
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input public matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[in] z The output multiplication result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void mat_mul(const PublicMatrix<double>& x, const ArithMatrix& y, ArithMatrix& z) const;

    /**
     * @brief Matrix multiplication of public matrix arithmetic share matrices.
     *
     * The operation is element-wise, e.g., z = x * y
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The second input arithmetic secret-shared matrix
     * @param[in] y The first input public matrix
     * @param[in] z The output multiplication result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void mat_mul(const ArithMatrix& x, const PublicMatrix<double>& y, ArithMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x > y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x > y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void greater(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x > y. results are stored in boolean secret shared matrix z.
     * y is public matrix.
     * ABY version of secure comparison is used.
     * z = (x > y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input public matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void greater(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x > y. results are stored in boolean secret shared matrix z.
     * y is public float.
     * ABY version of secure comparison is used.
     * z = (x > y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input public float
     * @param[out] z The output boolean share result
     */
    void greater(
            const std::shared_ptr<network::Network>& net, const ArithMatrix& x, PublicDouble y, BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x < y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x < y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z the output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x < y. results are stored in boolean secret shared matrix z.
     * y is public matrix.
     * ABY version of secure comparison is used.
     * z = (x < y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input public matrix
     * @param[out] z The output boolean share result
     */
    void less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, PublicDouble y, BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x < y. results are stored in boolean secret shared matrix z.
     * y is public matrix.
     * ABY version of secure comparison is used.
     * z = (x < y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input public float
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x >= y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x >= y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void greater_equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x >= y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x >= y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input public float
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void greater_equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x,
            const PublicMatrix<double>& y, BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x == y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x == y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x == y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x == y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input  public float
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PublicMatrix<double>& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure comparison.
     *
     * Compute if x != y. results are stored in boolean secret shared matrix z.
     * ABY version of secure comparison is used.
     * z = (x == y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input arithmetic secret-shared matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void not_equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y,
            BoolMatrix& z) const;

    /**
     * @brief Secure select.
     *
     * Accroding x, set z = 0 or y;
     * (e.g., z = (x == 0) ? 0 : y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The input the boolean shared matrix
     * @param[in] y The input of arithmetic secret-shared matrix
     * @param[out] z The output of arithmetic secret-shared matrix
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void multiplexer(const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const ArithMatrix& y,
            ArithMatrix& z) const;

    /**
     * @brief Secure select.
     *
     * Accroding x, set z = x or y;
     * (e.g., z = (alpha == 0) ? x : y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] alpha The input the boolean shared matrix.
     * @param[in] x The input of arithmetic secret-shared matrix.
     * @param[in] y The input of arithmetic secret-shared matrix.
     * @param[out] z The output of arithmetic secret-shared matrix
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void multiplexer(const std::shared_ptr<network::Network>& net, const BoolMatrix& alpha, const ArithMatrix& x,
            const ArithMatrix& y, ArithMatrix& z) const;

    /**
     * @brief Sum of all values column-wise in a arithmetic secret shared matrix.
     *
     *  Results are stored in a 1 * num_column secret shared matrix z.
     *
     * @param[in] x The input arithmetic secret-shared matrix
     * @param[out] z The output secret-shared matrix
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void sum(const ArithMatrix& x, ArithMatrix& z) const;

    /**
     * @brief Secert share based shuffling solution for shuffling a private matrix(for rows).
     *
     * One party has a private matrix, another contributes permutation to shuffling the private matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The private matrix.
     * @param[out] out Shuffled matrix.
     */
    void shuffle(const std::shared_ptr<network::Network>& net, const PrivateMatrix<double>& in, ArithMatrix& out) const;

    /**
     * @brief Secert share based shuffling solution for shuffling a secret shared matrix(for rows).
     *
     * Both parties have secret-shared matrix. Both parties contributes permutation to shuffling.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The matrix need to shuffle.
     * @param[out] out Shuffled matrix.
     */
    void shuffle(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out) const;

    /**
     * @brief Secert share based shuffling solution for shuffling a private matrix(for rows).
     *
     * One party has a private matrix, another contributes permutation to shuffling the private matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] perm The permutation which is belong to one party.
     * @param[in] in The private matrix.
     * @param[out] out Shuffled matrix.
     */
    void shuffle(const std::shared_ptr<network::Network>& net, const PrivatePermutation& perm, const ArithMatrix& in,
            ArithMatrix& out) const;

    /**
     * @brief Quick sort.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] matrix The input matrix to be sorted. The result will be in same matrix.
     */
    void quick_sort(const std::shared_ptr<network::Network>& net, ArithMatrix& matrix) const;

    /**
     * @brief Quick sort by column.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] col_index To be sorted by which column.
     * @param[in] in The input matrix to be sorted.
     * @param[out] out The output matrix which is sorted.
     */
    void quick_sort(const std::shared_ptr<network::Network>& net, std::size_t col_index, const ArithMatrix& in,
            ArithMatrix& out) const;

    /**
     * @brief Row major argmax and max interface.
     *
     * Return an arith matrix's row major argmax and max.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The arith matrix to compute.
     * @param[out] max_index The argmax arith share matrix of the input arith matrix.
     * @param[out] max_value The max arith share matrix of the input arith matrix.
     */
    void argmax_and_max(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& max_index,
            ArithMatrix& max_value) const;

    /**
     * @brief Row major argmin and min interface.
     *
     * Return an arith matrix's row major argmin and min.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The arith matrix to compute.
     * @param[out] min_index The argmin arith share matrix of the input arith matrix.
     * @param[out] min_value The min arith share matrix of the input arith matrix.
     */
    void argmin_and_min(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& max_index,
            ArithMatrix& max_value) const;

    /**
     * @brief Encrypt a private int64_t matrix with Paillier encryption using a party's public key.
     *
     * @param[in] input The input private matrix.
     * @param[out] output The paillier cipher matrix.
     * @param[in] using_self_pk If using_self_pk == 1, use self's public key, otherwise use the other's public key.
     * Default to 1.
     */
    void encrypt(const PrivateMatrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk = 1) const;

    /**
     * @brief Encrypt a private double matrix with Paillier encryption using a party's public key.
     *
     * @param[in] input The input private matrix.
     * @param[out] output The paillier cipher matrix.
     * @param[in] using_self_pk If using_self_pk == 1, use self's public key, otherwise use the other's public key.
     * Default to 1.
     */
    void encrypt(const PrivateMatrix<double>& input, PaillierMatrix& output, size_t using_self_pk = 1) const;

    /**
     * @brief Decrypt a Paillier cipher matrix and recover the plaintext, the plaintext is in the form of int64_t.
     * Only the party with the private key can decrypt the ciphertext.
     * NOTE: please make sure the output matrix is initialized with the correct size and party_id.
     *
     * @param[in] input The input Paillier matrix.
     * @param[out] output The decrypted plain matrix.
     */
    void decrypt(const PaillierMatrix& input, PrivateMatrix<std::int64_t>& output) const;

    /**
     * @brief Decrypt a Paillier cipher matrix and recover the plaintext, the plaintext is in the form of double.
     * Only the party with the private key can decrypt the ciphertext.
     * NOTE: please make sure the output matrix is initialized with the correct size and party_id.
     *
     * @param[in] input The input Paillier matrix.
     * @param[out] output The decrypted plain matrix.
     */
    void decrypt(const PaillierMatrix& input, PrivateMatrix<double>& output) const;

    /**
     * @brief Addition of two Paillier cipher matrices.
     * x and y should have the same size and party_id.
     * @param[in] x The first input Paillier matrix.
     * @param[in] y The second input Paillier matrix.
     * @param[out] z The output Paillier matrix of the addition result.
     */
    void add(const PaillierMatrix& x, const PaillierMatrix& y, PaillierMatrix& z, const bool self_pk = true) const;

    /**
     * @brief Addition of a Paillier cipher matrix and a Plaintext PrivateMatrix<std::int64_t>.
     * x and y should have the same size.
     * @param[in] x The first input PrivateMatrix<std::int64_t> matrix.
     * @param[in] y The second input Paillier matrix.
     * @param[out] z The output Paillier matrix of the addition result.
     */
    void add(const PrivateMatrix<std::int64_t>& x, const PaillierMatrix& y, PaillierMatrix& z,
            const bool self_pk = true) const;

    /**
     * @brief Addition of a Paillier cipher matrix and a Plaintext PrivateMatrix<std::double>.
     * x and y should have the same size.
     * @param[in] x The first input PrivateMatrix<std::double> matrix.
     * @param[in] y The second input Paillier matrix.
     * @param[out] z The output Paillier matrix of the addition result.
     */
    void add(const PrivateMatrix<double>& x, const PaillierMatrix& y, PaillierMatrix& z,
            const bool self_pk = true) const;

    /**
     * @brief Multiplication of a Paillier cipher matrix and a Plaintext PrivateMatrix<std::int64_t>.
     * x and y should have the same size.
     * We currently only support int64_t matrix as the plaintext.
     * @param[in] x The first input PrivateMatrix<std::int64_t> matrix.
     * @param[in] y The second input Paillier matrix.
     * @param[out] z The output Paillier matrix of the multiplication result.
     */
    void mul(const PrivateMatrix<std::int64_t>& x, const PaillierMatrix& y, PaillierMatrix& z,
            const bool self_pk = true) const;

    /**
     * @brief Convert a Paillier Matrix to a additive secret shared matrix ArithMatrix.
     * This function has very specific use case: party i has a Paillier matrix, which is encrypted by party (1 - i)'s
     * public key. Then we can use h2a to do the conversion. Otherwise, if the Paillier matrix is encrypted by party i's
     * public key, we should directly decrypt it then use share(). Besides, please make sure that the values in Paillier
     * cipher is less than 2^64, otherwise the values cannot be represented by int64_t share.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input Paillier matrix.
     * @param[in] out The output additive secret share.
     */
    void h2a(const std::shared_ptr<network::Network>& net, const PaillierMatrix& in, ArithMatrix& out) const;

    /**
     * @brief Convert an additive-shared  ArithMatrix to a Paillier Matrix.
     * This function has very specific use case: party 0 and party 1 have additvely shared input, then they invoke a2h
     * to convert it to a Paillier matrix.
     *
     * Before running this function, please initialize the out matrix with the party_id.
     * The cipher will be encrypted with party_id's pk, and the result will be sent to (1 - party_id) .
     * Otherwise, if you want the cipher to be sent to party_id, this is equivilant to reveal() as party_id can decrypt
     * it and know the plaintext values.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The input Paillier matrix.
     * @param[in] out The output additive secret share.
     */
    void a2h(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, PaillierMatrix& out) const;

    /**
     * @brief Group by with sum operation.
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in: The additively secret-shared input matrix.
     * @param[in] encoding: The one-hot coding of group by.
     * @param[out] out: The result of aggregation operation after group by.
     */
    void groupby_sum(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, const ArithMatrix& encoding,
            ArithMatrix& out) const;

    /**
     * @brief Group by with count operation.
     * @param[in] in: The additively secret-shared input matrix.
     * @param[in] encoding: The one-hot coding of group by.
     * @param[out] out: The result of aggregation operation after group by.
     */
    void groupby_count(const ArithMatrix& in, const ArithMatrix& encoding, ArithMatrix& out) const;

    /**
     * @brief Group by with max operation.
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in: The additively secret-shared input matrix.
     * @param[in] encoding: The one-hot coding of group by.
     * @param[out] out: The result of aggregation operation after group by.
     */
    void groupby_max(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, const ArithMatrix& encoding,
            ArithMatrix& out) const;

    /**
     * @brief Group by with min operation.
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in: The additively secret-shared input matrix.
     * @param[in] encoding: The one-hot coding of group by.
     * @param[out] out: The result of aggregation operation after group by.
     */
    void groupby_min(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, const ArithMatrix& encoding,
            ArithMatrix& out) const;

    /**
     * @brief Calculate sum based on count of group.
     *
     * @param[in] grouped_count The group count.
     * @param[in] in The additively secret-shared input matrix.
     * @param[out] out The additively secret-shared output matrix.
     */
    void group_then_sum_by_grouped_count(
            const PublicMatrix<double>& grouped_count, const ArithMatrix& in, ArithMatrix& out) const;

    /**
     * @brief Sigmoid f(x) = 1 / (1 + e^-x).
     *
     * We use a cubic polynomial to fit sigmoid.
     * f(x) = 0.4999 + 0.24955x -0.018715x^3;
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The additively secret-shared input matrix.
     * @param[out] out The additively secret-shared output matrix.
     */
    void sigmoid(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out) const;

    /**
     * @brief Secure comparison.
     *
     * Protocol from Cheetah paper. compute if input is less than zero.
     * (e.g., if the most significant bit of x is 1)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The input arithmetic secret-shared matrix
     * @param[out] y The output of the boolean share result
     * @param[in] bit_length The bit length of input matrix values, default is 64.
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void less_than_zero(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& y,
            std::size_t bit_length = 64) const;

    /**
     * @brief Reveal condition then split data by count.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] cond The condition.
     * @param[in] in The additively secret-shared input matrix.
     * @param[out] out0 The additively secret-shared output matrix corresponding to condition 0.
     * @param[out] out1 The additively secret-shared output matrix corresponding to condition 1.
     */
    void split_by_condition(const std::shared_ptr<network::Network>& net, const BoolMatrix& cond, const ArithMatrix& in,
            ArithMatrix& out0, ArithMatrix& out1) const;

    std::size_t party() const {
        return party_id_;
    }

    const std::shared_ptr<solo::ahepaillier::PublicKey>& get_pk() const {
        return paillier_engine_->get_pk();
    }

    const std::shared_ptr<solo::ahepaillier::PublicKey>& get_pk_other() const {
        return paillier_engine_->get_pk_other();
    }

private:
    void a2b(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& z) const;

    void millionaire(const std::shared_ptr<network::Network>& net, Matrix<std::int64_t>& x, BoolMatrix& y,
            std::size_t bit_length = 64) const;

    void kogge_stone_ppa(const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y,
            BoolMatrix& z) const;

    void pack_and_evaluate_and(
            const std::shared_ptr<network::Network>& net, BoolMatrix& x, BoolMatrix& y, BoolMatrix& z) const;

    void quick_sort(const std::shared_ptr<network::Network>& net, const std::vector<std::size_t>& index,
            ArithMatrix& matrix) const;

    void quick_sort(const std::shared_ptr<network::Network>& net, const std::vector<std::size_t>& index,
            std::size_t col_index, ArithMatrix& matrix) const;

    void encoding_parser(
            const std::shared_ptr<network::Network>& net, const std::vector<ArithMatrix>& in, ArithMatrix& out) const;

    void elementwise_mul(const std::shared_ptr<network::Network>& net, const PublicMatrix<std::int64_t>& precision,
            const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z) const;

    void sync_shape(const std::shared_ptr<network::Network>& net, const std::size_t party,
            const std::vector<std::size_t>& in, std::vector<std::size_t>& out) const;

    std::shared_ptr<PRNG> rand_generator_ = nullptr;
    std::shared_ptr<BooleanTriple> rand_bool_triple_generator_ = nullptr;
    std::shared_ptr<ArithmeticTriple> rand_arith_triple_generator_ = nullptr;
    std::shared_ptr<ObliviousTransfer> ot_generator_ = nullptr;
    std::shared_ptr<STGenerator> st_generator_ = nullptr;
    std::shared_ptr<Paillier> paillier_engine_ = nullptr;
    std::size_t party_id_ = 0;
};

}  // namespace duet
}  // namespace petace
