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
#include "duet/paillier/paillier.h"

#include "duet/util/common.h"

namespace petace {
namespace duet {

Paillier::Paillier(std::size_t party) : party_id_(party) {
    encoder_ = std::make_shared<solo::ahepaillier::Encoder>();
    key_generator_ = std::make_shared<solo::ahepaillier::KeyGenerator>(kPaillierKeySize);

    key_generator_->get_key_pair(sk_, pk_);
    evaluator_ = std::make_shared<solo::ahepaillier::Evaluator>(pk_, sk_);
    decryptor_ = std::make_shared<solo::ahepaillier::Decryptor>(sk_);
    encryptor_ = std::make_shared<solo::ahepaillier::Encryptor>(pk_, sk_);
}

void Paillier::encrypt(const PrivateMatrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk) const {
    std::vector<solo::ahepaillier::Plaintext> pt;
    output.resize(input.rows(), input.cols());
    pt.resize(input.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        encoder_->encode(input(i), pt[i]);
    }
    if (using_self_pk == 1) {
        encryptor_->encrypt_many(pt, output.ciphers(), kPaillierThreads);
    } else {
        encryptor_the_other_->encrypt_many(pt, output.ciphers(), kPaillierThreads);
    }
}

void Paillier::encrypt(const Matrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk) const {
    std::vector<solo::ahepaillier::Plaintext> pt;
    output.resize(input.rows(), input.cols());
    pt.resize(input.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        encoder_->encode(input(i), pt[i]);
    }
    if (using_self_pk == 1) {
        encryptor_->encrypt_many(pt, output.ciphers(), kPaillierThreads);
    } else {
        encryptor_the_other_->encrypt_many(pt, output.ciphers(), kPaillierThreads);
    }
}

void Paillier::decrypt(const PaillierMatrix& input, PrivateMatrix<std::int64_t>& output) const {
    if (input.party() != party_id_) {
        throw std::runtime_error("Cannot decrypt ciphertext from other party");
    }
    std::vector<solo::ahepaillier::Plaintext> pt;
    output.resize(input.rows(), input.cols());
    decryptor_->decrypt_many(input.ciphers(), pt, kPaillierThreads);
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        output(i) = encoder_->decode(pt[i]);
    }
}

void Paillier::decrypt(const PaillierMatrix& input, Matrix<std::int64_t>& output) const {
    if (input.party() != party_id_) {
        throw std::runtime_error("Cannot decrypt ciphertext from other party");
    }
    std::vector<solo::ahepaillier::Plaintext> pt;
    output.resize(input.rows(), input.cols());
    decryptor_->decrypt_many(input.ciphers(), pt, kPaillierThreads);
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        output(i) = encoder_->decode(pt[i]);
    }
}

void Paillier::encrypt_double(const PrivateMatrix<double>& input, PaillierMatrix& output, size_t using_self_pk) const {
    std::vector<solo::ahepaillier::Plaintext> pt;
    output.resize(input.rows(), input.cols());
    std::vector<std::uint64_t> vec;
    vec.resize(input.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        vec[i] = double_to_fixed(input(i));
    }
    encoder_->encode(vec, pt);
    if (using_self_pk == 1) {
        encryptor_->encrypt_many(pt, output.ciphers(), kPaillierThreads);
    } else {
        encryptor_the_other_->encrypt_many(pt, output.ciphers(), kPaillierThreads);
    }
}

void Paillier::decrypt_double(const PaillierMatrix& input, PrivateMatrix<double>& output) const {
    if (input.party() != party_id_) {
        throw std::runtime_error("Cannot decrypt ciphertext from other party");
    }
    std::vector<solo::ahepaillier::Plaintext> pt;
    output.resize(input.rows(), input.cols());
    decryptor_->decrypt_many(input.ciphers(), pt, kPaillierThreads);
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        output(i) = fixed_to_double(encoder_->decode(pt[i]));
    }
}

void Paillier::encode(const std::uint64_t in, solo::ahepaillier::Plaintext& out) const {
    encoder_->encode(in, out);
}

void Paillier::encode(const std::vector<std::uint64_t>& in, std::vector<solo::ahepaillier::Plaintext>& out) const {
    encoder_->encode(in, out);
}

std::uint64_t Paillier::decode(const solo::ahepaillier::Plaintext& in) const {
    return encoder_->decode(in);
}

const std::shared_ptr<solo::ahepaillier::PublicKey>& Paillier::get_pk() const {
    return pk_;
}

const std::shared_ptr<solo::ahepaillier::PublicKey>& Paillier::get_pk_other() const {
    return pk_other_;
}

void Paillier::set_pk_other(const std::shared_ptr<solo::ahepaillier::PublicKey>& in) {
    pk_other_ = in;
    encryptor_the_other_ = std::make_shared<solo::ahepaillier::Encryptor>(in);
    evaluator_the_other_ = std::make_shared<solo::ahepaillier::Evaluator>(in);
}

void Paillier::add(const std::vector<solo::ahepaillier::Ciphertext>& in_0,
        const std::vector<solo::ahepaillier::Ciphertext>& in_1, std::vector<solo::ahepaillier::Ciphertext>& out,
        const bool self_pk) const noexcept {
    if (self_pk) {
        evaluator_->add_many(in_0, in_1, out, kPaillierThreads);
    } else {
        evaluator_the_other_->add_many(in_0, in_1, out, kPaillierThreads);
    }
}

void Paillier::add(const std::vector<solo::ahepaillier::Ciphertext>& in_0,
        const std::vector<solo::ahepaillier::Plaintext>& in_1, std::vector<solo::ahepaillier::Ciphertext>& out,
        const bool self_pk) const noexcept {
    if (self_pk) {
        evaluator_->add_many(in_0, in_1, out, kPaillierThreads);
    } else {
        evaluator_the_other_->add_many(in_0, in_1, out, kPaillierThreads);
    }
}

void Paillier::mul(const std::vector<solo::ahepaillier::Ciphertext>& in_0,
        const std::vector<solo::ahepaillier::Plaintext>& in_1, std::vector<solo::ahepaillier::Ciphertext>& out,
        const bool self_pk) const noexcept {
    if (self_pk) {
        evaluator_->mul_many(in_0, in_1, out, kPaillierThreads);
    } else {
        evaluator_the_other_->mul_many(in_0, in_1, out, kPaillierThreads);
    }
}

void Paillier::serialize_public_key_to_bytes(
        const std::shared_ptr<solo::ahepaillier::PublicKey>& pk, solo::Byte* out) const {
    pk->serialize_to_bytes(out, pk->public_key_byte_count());
}

void Paillier::deserialize_public_key_from_bytes(
        const solo::Byte* in, std::shared_ptr<solo::ahepaillier::PublicKey>& pk) const {
    pk->deserialize_from_bytes(in, get_public_key_byte_count());
}

std::size_t Paillier::get_public_key_byte_count() const {
    return pk_->public_key_byte_count();
}

}  // namespace duet
}  // namespace petace
