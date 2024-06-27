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
#include <vector>

#include "gmp.h"
#include "gmpxx.h"
#include "seal/seal.h"

#include "network/network.h"
#include "solo/prng.h"

#include "duet/beaver_triple/arithmetic_triple.h"
#include "duet/util/consts.h"

namespace petace {
namespace duet {

/**
 * @brief Implementation of Beaver's arithmetic triple generation protocol.
 *
 * This class implements the arithmetic triple generation protocol based on OT.
 */
class ArithmeticTripleFHE : public ArithmeticTriple {
public:
    /**
     * @brief Constructor of the class.
     *
     * @param[in] net The network instance (e.g., from PETAce-Network)
     * @param[in] party The party id
     * @param[in] prng The PRNG instance
     * @param[in] ot The OT instance
     */
    ArithmeticTripleFHE(const std::shared_ptr<network::Network>& net, std::size_t party_id, std::shared_ptr<PRNG> prng,
            const std::shared_ptr<ObliviousTransfer>& ot);

    ~ArithmeticTripleFHE() = default;

    std::vector<std::int64_t> get_rand_triple(const std::shared_ptr<network::Network>& net);

private:
    void refill_rand_triple_buffer(const std::shared_ptr<network::Network>& net);

    void gen_rand_triple(const std::shared_ptr<network::Network>& net);

    void crt_decoding(
            std::vector<mpz_class>& result, const std::vector<std::vector<std::uint64_t>>& v, const std::int64_t* m);

    void crt_encoding(const std::vector<std::int64_t>& values, const std::int64_t* primes,
            std::vector<std::vector<std::uint64_t>>& res);

    mpz_class get_random_mpz(std::size_t bits);

    void big_num_crt_encoding(std::vector<mpz_class>& values, std::size_t len_vec, const std::int64_t* primes,
            std::vector<std::vector<std::uint64_t>>& res);

    void crt_bfv_encrypt(const std::vector<std::shared_ptr<seal::BatchEncoder>>& batch_encoder,
            const std::vector<std::shared_ptr<seal::Encryptor>>& encryptor, const int64_t* primes,
            const std::vector<std::int64_t>& plain, std::vector<seal::Ciphertext>& ciphers);

    void send_fhe_cipher(const std::shared_ptr<network::Network>& net, seal::Ciphertext& cipher);

    void recv_fhe_cipher(
            const std::shared_ptr<network::Network>& net, seal::SEALContext& context, seal::Ciphertext& cipher);

    std::size_t rand_triple_idx_ = 0;
    std::vector<std::vector<std::int64_t>> rand_triple_buff_{};
    std::size_t party_id_ = 0;
    std::shared_ptr<ObliviousTransfer> ot_ = nullptr;
    std::shared_ptr<PRNG> prng_ = nullptr;

    std::vector<std::shared_ptr<seal::BatchEncoder>> batch_encoders_{};
    std::vector<std::shared_ptr<seal::Encryptor>> encryptors_{};
    std::vector<std::shared_ptr<seal::Decryptor>> decryptors_{};
    std::vector<std::shared_ptr<seal::Evaluator>> evaluators_{};
    std::vector<seal::SEALContext> vec_context_{};
    std::int64_t primes_[kCRTPrimeCount] = {};
    std::size_t num_iteration_ = 0;
};

inline std::shared_ptr<ArithmeticTriple> create_triple_engine_fhe(const std::shared_ptr<network::Network>& net,
        std::size_t party_id, std::shared_ptr<PRNG> prng, const std::shared_ptr<ObliviousTransfer>& ot) {
    return std::make_shared<ArithmeticTripleFHE>(net, party_id, prng, ot);
}

}  // namespace duet
}  // namespace petace
