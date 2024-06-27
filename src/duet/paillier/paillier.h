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

#include "gmp.h"
#include "gmpxx.h"

#include "solo/ahe_paillier.h"

#include "duet/util/consts.h"
#include "duet/util/defines.h"
#include "duet/util/matrix.h"

namespace petace {
namespace duet {

class Paillier {
public:
    explicit Paillier(std::size_t party);

    ~Paillier() = default;

    void encrypt(const PrivateMatrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk = 1) const;

    void encrypt(const Matrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk = 1) const;

    void decrypt(const PaillierMatrix& input, PrivateMatrix<std::int64_t>& output) const;

    void decrypt(const PaillierMatrix& input, Matrix<std::int64_t>& output) const;

    void encrypt_double(const PrivateMatrix<double>& input, PaillierMatrix& output, size_t using_self_pk = 1) const;

    void decrypt_double(const PaillierMatrix& input, PrivateMatrix<double>& output) const;

    void encode(const std::uint64_t in, solo::ahepaillier::Plaintext& out) const;

    void encode(const std::vector<std::uint64_t>& in, std::vector<solo::ahepaillier::Plaintext>& out) const;

    void encode(const Matrix<std::int64_t>& in, std::vector<solo::ahepaillier::Plaintext>& out) const;

    std::uint64_t decode(const solo::ahepaillier::Plaintext& in) const;

    void add(const std::vector<solo::ahepaillier::Ciphertext>& in_0,
            const std::vector<solo::ahepaillier::Ciphertext>& in_1, std::vector<solo::ahepaillier::Ciphertext>& out,
            const bool self_pk) const noexcept;

    void add(const std::vector<solo::ahepaillier::Ciphertext>& in_0,
            const std::vector<solo::ahepaillier::Plaintext>& in_1, std::vector<solo::ahepaillier::Ciphertext>& out,
            const bool self_pk) const noexcept;

    void mul(const std::vector<solo::ahepaillier::Ciphertext>& in_0,
            const std::vector<solo::ahepaillier::Plaintext>& in_1, std::vector<solo::ahepaillier::Ciphertext>& out,
            const bool self_pk) const noexcept;

    const std::shared_ptr<solo::ahepaillier::PublicKey>& get_pk() const;

    const std::shared_ptr<solo::ahepaillier::PublicKey>& get_pk_other() const;

    void set_pk_other(const std::shared_ptr<solo::ahepaillier::PublicKey>& in);

    void serialize_public_key_to_bytes(const std::shared_ptr<solo::ahepaillier::PublicKey>& pk, solo::Byte* out) const;

    void deserialize_public_key_from_bytes(
            const solo::Byte* in, std::shared_ptr<solo::ahepaillier::PublicKey>& pk) const;

    std::size_t get_public_key_byte_count() const;

private:
    Paillier(const Paillier&) = delete;
    std::size_t party_id_ = 0;
    std::shared_ptr<solo::ahepaillier::Encoder> encoder_ = nullptr;
    std::shared_ptr<solo::ahepaillier::Evaluator> evaluator_ = nullptr;
    std::shared_ptr<solo::ahepaillier::Evaluator> evaluator_the_other_ = nullptr;
    std::shared_ptr<solo::ahepaillier::KeyGenerator> key_generator_ = nullptr;
    std::shared_ptr<solo::ahepaillier::Encryptor> encryptor_ = nullptr;
    std::shared_ptr<solo::ahepaillier::Encryptor> encryptor_the_other_ = nullptr;
    std::shared_ptr<solo::ahepaillier::Decryptor> decryptor_ = nullptr;
    std::shared_ptr<solo::ahepaillier::SecretKey> sk_ = nullptr;
    std::shared_ptr<solo::ahepaillier::PublicKey> pk_ = nullptr;
    std::shared_ptr<solo::ahepaillier::PublicKey> pk_other_ = nullptr;
    std::shared_ptr<solo::ahepaillier::Serialization> serialization_ = nullptr;
};

}  // namespace duet
}  // namespace petace
