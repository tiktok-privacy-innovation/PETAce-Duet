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

#include "duet/beaver/arithmetic_triplet.h"
#include "duet/beaver/boolean_triplet.h"
#include "duet/st_generator/st_generator.h"
#include "duet/util/defines.h"
#include "duet/util/io.h"
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
     * @brief Secret share a plain matrix to both parties, where the plaintext contains double values.
     *
     * One party provides a plain matrix, and the output is that each of both parties holds
     * additive share of the matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party The party id of the plain matrix provider.
     * @param[in] in The input plain matrix. (non-sender party should provide an empty matrix)
     * @param[out] out The output additive share matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void share(const std::shared_ptr<network::Network>& net, std::size_t party, const PlainMatrix<double>& in,
            ArithMatrix& out);

    /**
     * @brief Secret share a plain matrix to both parties, where the plaintext contains int64_t values.
     *
     * One party provides a plain matrix, and the output is that each of both parties holds
     * additive share of the matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party The party id of the plain matrix provider.
     * @param[in] in The input plain matrix. (non-sender party should provide an empty matrix)
     * @param[out] out The output additive share matrix.
     * @throws std::invalid_argument if in.size() == 0.
     */
    void share(const std::shared_ptr<network::Network>& net, std::size_t party, const PlainMatrix<std::int64_t>& in,
            ArithMatrix& out);

    /**
     * @brief Reconstruct additive-secret-shared matrix to reveal the plain matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party The party id of the plain matrix receiver.
     * @param[in] in The input secret shared matrix.
     * @param[out] out The output plain matrix. (non-receiver party should provide an empty matrix)
     * @throws std::invalid_argument if in.size() == 0.
     */
    void reveal(const std::shared_ptr<network::Network>& net, std::size_t party, const ArithMatrix& in,
            PlainMatrix<double>& out);

    /**
     * @brief Reconstruct boolean-secret-shared matrix to reveal the plain matrix.
     *
     * xor of two secret shared boolean values reveals the plain boolean value.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party The party id of the plain matrix receiver.
     * @param[in] in The input secret shared matrix.
     * @param[out] out The output plain matrix. (non-receiver party should provide an empty matrix)
     * @throws std::invalid_argument if in.size() == 0.
     */
    void reveal_bool(const std::shared_ptr<network::Network>& net, std::size_t party, const BoolMatrix& in,
            PlainMatrix<std::int64_t>& out);

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
    void elementwise_bool_mul(
            const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z);

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
    void elementwise_mul(
            const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, ArithMatrix& z);

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
    void greater(
            const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z);

    /**
     * @brief Secure comparison.
     *
     * Compute if x > y. results are stored in boolean secret shared matrix z.
     * y is plaintext matrix.
     * ABY version of secure comparison is used.
     * z = (x > y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input plaintext secret-shared matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void greater(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PlainMatrix<double>& y,
            BoolMatrix& z);

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
    void less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z);

    /**
     * @brief Secure comparison.
     *
     * Compute if x < y. results are stored in boolean secret shared matrix z.
     * y is plaintext matrix.
     * ABY version of secure comparison is used.
     * z = (x < y)
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] x The first input arithmetic secret-shared matrix
     * @param[in] y The second input plaintext secret-shared matrix
     * @param[out] z The output boolean share result
     * @throws std::invalid_argument if x.size() != y.size().
     */
    void less(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const PlainMatrix<double>& y,
            BoolMatrix& z);

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
    void greater_equal(
            const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z);

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
    void equal(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, const ArithMatrix& y, BoolMatrix& z);

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
            std::size_t bit_length = 64);

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
    void multiplexer(
            const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const ArithMatrix& y, ArithMatrix& z);

    /**
     * @brief Secert share based shuffling solution for shuffling a plain matrix(for rows).
     *
     * One party has a plain matrix, another contributes permutation to shuffling the plain matrix.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] party Who has the plain matrix.
     * @param[in] in The plain matrix.
     * @param[out] out Shuffled matrix.
     */
    void shuffle(const std::shared_ptr<network::Network>& net, std::size_t party, const PlainMatrix<double>& in,
            ArithMatrix& out);

    /**
     * @brief Secert share based shuffling solution for shuffling a secret shared matrix(for rows).
     *
     * Both parties have secret-shared matrix. Both parties contributes permutation to shuffling.
     *
     * @param[in] net The network interface (e.g., PETAce-Network interface).
     * @param[in] in The matrix need to shuffle.
     * @param[out] out Shuffled matrix.
     */
    void shuffle(const std::shared_ptr<network::Network>& net, const ArithMatrix& in, ArithMatrix& out);

    std::size_t party() {
        return party_id_;
    }

private:
    void a2b(const std::shared_ptr<network::Network>& net, const ArithMatrix& x, BoolMatrix& z);

    void millionaire(const std::shared_ptr<network::Network>& net, PlainMatrix<std::int64_t>& x, BoolMatrix& y,
            std::size_t bit_length = 64);

    void kogge_stone_ppa(
            const std::shared_ptr<network::Network>& net, const BoolMatrix& x, const BoolMatrix& y, BoolMatrix& z);

    void pack_and_evaluate_and(
            const std::shared_ptr<network::Network>& net, BoolMatrix& x, BoolMatrix& y, BoolMatrix& z);

    std::shared_ptr<PRNG> rand_generator_ = nullptr;
    std::shared_ptr<BooleanTriplet> rand_bool_triplet_generator_ = nullptr;
    std::shared_ptr<ArithmeticTriplet> rand_arith_triplet_generator_ = nullptr;
    std::shared_ptr<OTGenerator> ot_generator_ = nullptr;
    std::shared_ptr<STGenerator> st_generator_ = nullptr;
    std::size_t party_id_ = 0;
};

}  // namespace duet
}  // namespace petace
