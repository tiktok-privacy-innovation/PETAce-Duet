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
#include "duet/paillier_engine/paillier_engine.h"

#include "duet/util/common.h"

namespace petace {
namespace duet {

PaillierEngine::PaillierEngine(std::size_t party) : party_id_(party) {
    encoder_ = std::make_shared<solo::ahepaillier::Encoder>();
    key_generator_ = std::make_shared<solo::ahepaillier::KeyGenerator>(kPaillierKeySize);
    serialization_ = std::make_shared<solo::ahepaillier::Serialization>(kPaillierKeySize, true);
    key_generator_->get_key_pair(sk_, pk_);
    evaluator_ = std::make_shared<solo::ahepaillier::Evaluator>();
    decryptor_ = std::make_shared<solo::ahepaillier::Decryptor>(sk_);
    encryptor_ = std::make_shared<solo::ahepaillier::Encryptor>(pk_);
}

void PaillierEngine::encrypt(
        const PrivateMatrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk) const {
    solo::ahepaillier::Plaintext pt;
    output.resize(input.rows(), input.cols());
    encode(input.matrix(), pt);
    if (using_self_pk == 1) {
        encryptor_->encrypt(pt, output.ciphers());
    } else {
        encryptor_the_other_->encrypt(pt, output.ciphers());
    }
}

void PaillierEngine::encrypt(const Matrix<std::int64_t>& input, PaillierMatrix& output, size_t using_self_pk) const {
    solo::ahepaillier::Plaintext pt;
    output.resize(input.rows(), input.cols());
    encode(input, pt);
    if (using_self_pk == 1) {
        encryptor_->encrypt(pt, output.ciphers());
    } else {
        encryptor_the_other_->encrypt(pt, output.ciphers());
    }
}

void PaillierEngine::decrypt(const PaillierMatrix& input, PrivateMatrix<std::int64_t>& output) const {
    if (input.party() != party_id_) {
        throw std::runtime_error("Cannot decrypt ciphertext from other party");
    }
    solo::ahepaillier::Plaintext pt;
    std::vector<std::uint64_t> pt_decode;
    output.resize(input.rows(), input.cols());
    decryptor_->decrypt(input.ciphers(), pt);
    encoder_->decode(pt, pt_decode);
    for (std::size_t i = 0; i < static_cast<size_t>(input.size()); ++i) {
        output(i) = pt_decode[i];
    }
}

void PaillierEngine::decrypt(const PaillierMatrix& input, Matrix<std::int64_t>& output) const {
    if (input.party() != party_id_) {
        throw std::runtime_error("Cannot decrypt ciphertext from other party");
    }
    solo::ahepaillier::Plaintext pt;
    std::vector<std::uint64_t> pt_decode;
    output.resize(input.rows(), input.cols());
    decryptor_->decrypt(input.ciphers(), pt);
    encoder_->decode(pt, pt_decode);
    for (std::size_t i = 0; i < static_cast<size_t>(input.size()); ++i) {
        output(i) = pt_decode[i];
    }
}

void PaillierEngine::encrypt_double(
        const PrivateMatrix<double>& input, PaillierMatrix& output, size_t using_self_pk) const {
    solo::ahepaillier::Plaintext pt;
    output.resize(input.rows(), input.cols());
    std::vector<std::uint64_t> vec;
    vec.resize(input.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        vec[i] = double_to_fixed(input(i));
    }
    encoder_->encode(vec, pt);
    if (using_self_pk == 1) {
        encryptor_->encrypt(pt, output.ciphers());
    } else {
        encryptor_the_other_->encrypt(pt, output.ciphers());
    }
}

void PaillierEngine::decrypt_double(const PaillierMatrix& input, PrivateMatrix<double>& output) const {
    if (input.party() != party_id_) {
        throw std::runtime_error("Cannot decrypt ciphertext from other party");
    }
    solo::ahepaillier::Plaintext pt;
    std::vector<std::uint64_t> pt_decode;
    output.resize(input.rows(), input.cols());
    decryptor_->decrypt(input.ciphers(), pt);
    encoder_->decode(pt, pt_decode);
    for (std::size_t i = 0; i < static_cast<std::size_t>(input.size()); ++i) {
        output(i) = fixed_to_double(pt_decode[i]);
    }
}

void PaillierEngine::encode(const std::uint64_t in, solo::ahepaillier::Plaintext& out) const {
    encoder_->encode(in, out);
}

void PaillierEngine::encode(const std::vector<std::uint64_t>& in, solo::ahepaillier::Plaintext& out) const {
    encoder_->encode(in, out);
}

void PaillierEngine::encode(const Matrix<std::int64_t>& input, solo::ahepaillier::Plaintext& out) const {
    std::vector<uint64_t> vec(input.data(), input.data() + input.rows() * input.cols());
    encoder_->encode(vec, out);
}

std::uint64_t PaillierEngine::decode(const solo::ahepaillier::Plaintext& in) const {
    return encoder_->decode(in);
}

// Helper functions for serialization
ByteVector PaillierEngine::encode(const solo::ahepaillier::BigNum& bn, bool is_n_square) const {
    std::size_t bytes_len = (kPaillierKeySize + 7) / 8;
    bytes_len = bytes_len * (1 + is_n_square);
    ByteVector out(bytes_len);
    serialization_->bn_to_bytes(bn, out.data(), out.size());
    return out;
}

solo::ahepaillier::BigNum PaillierEngine::decode(const ByteVector& in) const {
    solo::ahepaillier::BigNum out;
    serialization_->bn_from_bytes(in.data(), in.size(), out);
    return out;
}

const std::shared_ptr<solo::ahepaillier::PublicKey>& PaillierEngine::get_pk() const {
    return pk_;
}

const std::shared_ptr<solo::ahepaillier::PublicKey>& PaillierEngine::get_pk_other() const {
    return pk_other_;
}

void PaillierEngine::set_pk_other(const std::shared_ptr<solo::ahepaillier::PublicKey>& in) {
    pk_other_ = in;
    encryptor_the_other_ = std::make_shared<solo::ahepaillier::Encryptor>(in);
}

void PaillierEngine::bn_to_pt(
        const std::vector<solo::ahepaillier::BigNum>& in, solo::ahepaillier::Plaintext& out) const {
    out = solo::ahepaillier::Plaintext(in);
}

void PaillierEngine::bn_to_ct(const std::vector<solo::ahepaillier::BigNum>& in,
        const std::shared_ptr<solo::ahepaillier::PublicKey>& pk, solo::ahepaillier::Ciphertext& out) const noexcept {
    out = solo::ahepaillier::Ciphertext(*pk, in);
}

void PaillierEngine::add(const solo::ahepaillier::Ciphertext& in_0, const solo::ahepaillier::Ciphertext& in_1,
        solo::ahepaillier::Ciphertext& out) const noexcept {
    evaluator_->add(in_0, in_1, out);
}

void PaillierEngine::add(const solo::ahepaillier::Ciphertext& in_0, const solo::ahepaillier::Plaintext& in_1,
        solo::ahepaillier::Ciphertext& out) const noexcept {
    evaluator_->add(in_0, in_1, out);
}

void PaillierEngine::mul(const solo::ahepaillier::Ciphertext& in_0, const solo::ahepaillier::Plaintext& in_1,
        solo::ahepaillier::Ciphertext& out) const noexcept {
    evaluator_->mul(in_0, in_1, out);
}

void PaillierEngine::serialize_public_key_to_bytes(
        const std::shared_ptr<solo::ahepaillier::PublicKey>& pk, solo::Byte* out) const {
    serialization_->serialize_public_key_to_bytes(pk, out, serialization_->public_key_byte_count());
}

void PaillierEngine::deserialize_public_key_from_bytes(
        const solo::Byte* in, std::shared_ptr<solo::ahepaillier::PublicKey>& pk) const {
    serialization_->deserialize_public_key_from_bytes(in, serialization_->public_key_byte_count(), pk);
}

std::size_t PaillierEngine::get_public_key_byte_count() const {
    return serialization_->public_key_byte_count();
}

}  // namespace duet
}  // namespace petace
