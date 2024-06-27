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

#include "duet/beaver_triple/arithmetic_triple_fhe.h"

namespace petace {
namespace duet {

ArithmeticTripleFHE::ArithmeticTripleFHE(const std::shared_ptr<network::Network>& net, std::size_t party_id,
        std::shared_ptr<PRNG> prng, const std::shared_ptr<ObliviousTransfer>& ot)
        : party_id_(party_id), ot_(ot), prng_(prng) {
    num_iteration_ = kDefaulArithmeticTripleBufferSize / kFHEBatchSize +
                     (kDefaulArithmeticTripleBufferSize % kFHEBatchSize != 0);
    rand_triple_buff_.resize(num_iteration_ * kFHEBatchSize);
    // Only first four primes are used, the reason to generate 8 values here is to avoid prime collision caused by SEAL.
    std::vector<seal::Modulus> vec_pi =
            seal::PlainModulus::Batching(kPolyModulusDegree, {44, 44, 44, 44, 44, 44, 44, 44});
    std::vector<seal::PublicKey> vec_pk;

    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        primes_[i] = vec_pi[i].value();
        seal::EncryptionParameters parms(seal::scheme_type::bfv);
        parms.set_poly_modulus_degree(kPolyModulusDegree);
        parms.set_coeff_modulus(seal::CoeffModulus::BFVDefault(kPolyModulusDegree));
        parms.set_plain_modulus(vec_pi[i]);
        vec_context_.emplace_back(seal::SEALContext(parms));
    }

    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        // key generation
        if (party_id_ == 0) {
            seal::KeyGenerator keygen(vec_context_[i]);
            seal::SecretKey secret_key = keygen.secret_key();
            seal::Serializable<seal::PublicKey> public_key_serialized = keygen.create_public_key();
            std::vector<seal::seal_byte> byte_buffer(static_cast<size_t>(public_key_serialized.save_size()));
            std::size_t s = byte_buffer.size();
            public_key_serialized.save(reinterpret_cast<seal::seal_byte*>(byte_buffer.data()), s);
            seal::PublicKey public_key;
            public_key.load(vec_context_[i], reinterpret_cast<const seal::seal_byte*>(byte_buffer.data()), s);
            encryptors_.emplace_back(std::make_shared<seal::Encryptor>(vec_context_[i], public_key));
            batch_encoders_.emplace_back(std::make_shared<seal::BatchEncoder>(vec_context_[i]));
            decryptors_.emplace_back(std::make_shared<seal::Decryptor>(vec_context_[i], secret_key));
        } else {
            batch_encoders_.emplace_back(std::make_shared<seal::BatchEncoder>(vec_context_[i]));
            evaluators_.emplace_back(std::make_shared<seal::Evaluator>(vec_context_[i]));
        }
    }
    refill_rand_triple_buffer(net);
}

std::vector<std::int64_t> ArithmeticTripleFHE::get_rand_triple(const std::shared_ptr<network::Network>& net) {
    if (rand_triple_idx_ >= rand_triple_buff_.size()) {
        refill_rand_triple_buffer(net);
    }
    std::vector<std::int64_t> ret;
    ret.emplace_back(rand_triple_buff_[rand_triple_idx_][0]);
    ret.emplace_back(rand_triple_buff_[rand_triple_idx_][1]);
    ret.emplace_back(rand_triple_buff_[rand_triple_idx_][2]);
    rand_triple_idx_++;
    return ret;
}

void ArithmeticTripleFHE::refill_rand_triple_buffer(const std::shared_ptr<network::Network>& net) {
    gen_rand_triple(net);
    rand_triple_idx_ = 0;
    return;
}

void ArithmeticTripleFHE::crt_decoding(
        std::vector<mpz_class>& result, const std::vector<std::vector<std::uint64_t>>& v, const std::int64_t* m) {
    std::size_t size_input = v[0].size();
    mpz_class N = 1;
    std::vector<mpz_class> p(kCRTPrimeCount);
    std::vector<mpz_class> b(kCRTPrimeCount);
    std::vector<mpz_class> b_inverse(kCRTPrimeCount);
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        N *= static_cast<std::uint64_t>(m[i]);
    }
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        p[i] = static_cast<std::uint64_t>(m[i]);
        mpz_divexact(b[i].get_mpz_t(), N.get_mpz_t(), p[i].get_mpz_t());
        mpz_invert(b_inverse[i].get_mpz_t(), b[i].get_mpz_t(), p[i].get_mpz_t());
    }
    for (std::size_t i = 0; i < size_input; i++) {
        mpz_class res = 0;
        for (std::size_t j = 0; j < kCRTPrimeCount; j++) {
            res += mpz_class(v[j][i]) * b[j] * b_inverse[j];
        }
        mpz_mod(res.get_mpz_t(), res.get_mpz_t(), N.get_mpz_t());
        result[i] = res;
    }
}

void ArithmeticTripleFHE::big_num_crt_encoding(std::vector<mpz_class>& values, std::size_t len_vec,
        const std::int64_t* primes, std::vector<std::vector<std::uint64_t>>& res) {
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        res.emplace_back(std::vector<std::uint64_t>());
    }
    mpz_class r, current_prime;
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        current_prime = primes[i];
        for (std::size_t j = 0; j < len_vec; j++) {
            mpz_mod(r.get_mpz_t(), values[j].get_mpz_t(), current_prime.get_mpz_t());
            res[i].emplace_back(mpz_get_ui(r.get_mpz_t()));
        }
    }
}

mpz_class ArithmeticTripleFHE::get_random_mpz(std::size_t bits) {
    mpz_class out = 0;
    do {
        std::size_t byte_count = (bits + 7 / 8);
        std::vector<solo::Byte> in(byte_count);
        prng_->get_unique_rand_gen()->generate(byte_count, in.data());
        mpz_import(out.get_mpz_t(), in.size(), 1, 1, 0, 0, in.data());
        out >>= (byte_count * 8 - bits);
    } while (out == 0 || mpz_sizeinbase(out.get_mpz_t(), 2) != bits);
    return out;
}

void ArithmeticTripleFHE::crt_encoding(const std::vector<std::int64_t>& values, const std::int64_t* primes,
        std::vector<std::vector<std::uint64_t>>& res) {
    std::size_t len_vec = values.size();
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        res.emplace_back(std::vector<std::uint64_t>());
    }
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        std::uint64_t current_prime = static_cast<uint64_t>(primes[i]);
        for (std::size_t j = 0; j < len_vec; j++) {
            res[i].emplace_back(static_cast<uint64_t>(values[j]) % current_prime);
        }
    }
}

void ArithmeticTripleFHE::crt_bfv_encrypt(const std::vector<std::shared_ptr<seal::BatchEncoder>>& batch_encoder,
        const std::vector<std::shared_ptr<seal::Encryptor>>& encryptor, const int64_t* primes,
        const std::vector<std::int64_t>& plain, std::vector<seal::Ciphertext>& ciphers) {
    std::vector<std::vector<std::uint64_t>> crt_plaintext;
    crt_encoding(plain, primes, crt_plaintext);
    seal::Plaintext plain_matrix;
    for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
        batch_encoder[i]->encode(crt_plaintext[i], plain_matrix);
        encryptor[i]->encrypt(plain_matrix, ciphers[i]);
    }
}

void ArithmeticTripleFHE::send_fhe_cipher(const std::shared_ptr<network::Network>& net, seal::Ciphertext& cipher) {
    std::vector<seal::seal_byte> byte_buffer(static_cast<size_t>(cipher.save_size()));
    std::size_t buffer_size = byte_buffer.size();
    cipher.save(reinterpret_cast<seal::seal_byte*>(byte_buffer.data()), buffer_size, seal::compr_mode_type::zstd);
    net->send_data(&buffer_size, sizeof(std::size_t));
    net->send_data(byte_buffer.data(), byte_buffer.size() * sizeof(seal::seal_byte));
}

void ArithmeticTripleFHE::recv_fhe_cipher(
        const std::shared_ptr<network::Network>& net, seal::SEALContext& context, seal::Ciphertext& cipher) {
    std::size_t buffer_size;
    net->recv_data(&buffer_size, sizeof(std::size_t));
    std::vector<seal::seal_byte> byte_buffer(buffer_size);
    net->recv_data(byte_buffer.data(), buffer_size * sizeof(seal::seal_byte));
    cipher.load(context, reinterpret_cast<const seal::seal_byte*>(byte_buffer.data()), buffer_size);
}

void ArithmeticTripleFHE::gen_rand_triple(const std::shared_ptr<network::Network>& net) {
    std::vector<seal::Ciphertext> a_0_encrypted(kCRTPrimeCount);
    std::vector<seal::Ciphertext> b_0_encrypted(kCRTPrimeCount);
    for (std::size_t iter = 0; iter < num_iteration_; iter++) {
        if (party_id_ == 0) {
            std::vector<std::int64_t> a_0(kFHEBatchSize);
            std::vector<std::int64_t> b_0(kFHEBatchSize);
            for (std::size_t i = 0; i < kFHEBatchSize; i++) {
                a_0[i] = prng_->get_unique_rand<std::int64_t>();
                b_0[i] = prng_->get_unique_rand<std::int64_t>();
            }
            crt_bfv_encrypt(batch_encoders_, encryptors_, primes_, a_0, a_0_encrypted);
            crt_bfv_encrypt(batch_encoders_, encryptors_, primes_, b_0, b_0_encrypted);

            for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
                send_fhe_cipher(net, a_0_encrypted[i]);
                send_fhe_cipher(net, b_0_encrypted[i]);
            }
            std::vector<seal::Ciphertext> c_0_encrypted(kCRTPrimeCount);
            std::vector<seal::Plaintext> c_0_plaintext(kCRTPrimeCount);
            std::vector<std::vector<std::uint64_t>> c_0_crt(kCRTPrimeCount);
            std::vector<std::int64_t> c_0(kFHEBatchSize);
            std::vector<mpz_class> c_0_mpz(kFHEBatchSize);
            for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
                recv_fhe_cipher(net, vec_context_[i], c_0_encrypted[i]);
                decryptors_[i]->decrypt(c_0_encrypted[i], c_0_plaintext[i]);
                batch_encoders_[i]->decode(c_0_plaintext[i], c_0_crt[i]);
            }
            crt_decoding(c_0_mpz, c_0_crt, primes_);
            for (std::size_t i = 0; i < kFHEBatchSize; i++) {
                c_0_mpz[i] += mpz_class(a_0[i] * b_0[i]);
                c_0[i] = mpz_get_ui(c_0_mpz[i].get_mpz_t());
            }
            for (std::size_t i = 0; i < kFHEBatchSize; i++) {
                rand_triple_buff_[iter * kFHEBatchSize + i].emplace_back(a_0[i]);
                rand_triple_buff_[iter * kFHEBatchSize + i].emplace_back(b_0[i]);
                rand_triple_buff_[iter * kFHEBatchSize + i].emplace_back(c_0[i]);
            }

        } else {
            // Party 1
            for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
                recv_fhe_cipher(net, vec_context_[i], a_0_encrypted[i]);
                recv_fhe_cipher(net, vec_context_[i], b_0_encrypted[i]);
            }
            std::vector<std::int64_t> a_1(kFHEBatchSize);
            std::vector<std::int64_t> b_1(kFHEBatchSize);
            std::vector<std::int64_t> c_1(kFHEBatchSize);
            std::vector<mpz_class> r(kFHEBatchSize);

            for (std::size_t i = 0; i < kFHEBatchSize; i++) {
                a_1[i] = prng_->get_unique_rand<std::int64_t>();
                b_1[i] = prng_->get_unique_rand<std::int64_t>();
                r[i] = get_random_mpz(kFHERandomBitLength);
            }

            std::vector<std::vector<std::uint64_t>> crt_a_1;
            std::vector<std::vector<std::uint64_t>> crt_b_1;
            std::vector<std::vector<std::uint64_t>> crt_r;
            crt_encoding(a_1, primes_, crt_a_1);
            crt_encoding(b_1, primes_, crt_b_1);
            big_num_crt_encoding(r, kFHEBatchSize, primes_, crt_r);

            std::vector<seal::Plaintext> crt_a_1_plaintext;
            std::vector<seal::Plaintext> crt_b_1_plaintext;
            std::vector<seal::Plaintext> crt_r_plaintext;

            for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
                crt_a_1_plaintext.emplace_back(seal::Plaintext());
                batch_encoders_[i]->encode(crt_a_1[i], crt_a_1_plaintext[i]);
                crt_b_1_plaintext.emplace_back(seal::Plaintext());
                batch_encoders_[i]->encode(crt_b_1[i], crt_b_1_plaintext[i]);
                crt_r_plaintext.emplace_back(seal::Plaintext());
                batch_encoders_[i]->encode(crt_r[i], crt_r_plaintext[i]);
                evaluators_[i]->multiply_plain_inplace(a_0_encrypted[i], crt_b_1_plaintext[i]);
                evaluators_[i]->multiply_plain_inplace(b_0_encrypted[i], crt_a_1_plaintext[i]);
                evaluators_[i]->add_inplace(a_0_encrypted[i], b_0_encrypted[i]);
                evaluators_[i]->add_plain_inplace(a_0_encrypted[i], crt_r_plaintext[i]);
            }
            for (std::size_t i = 0; i < kCRTPrimeCount; i++) {
                send_fhe_cipher(net, a_0_encrypted[i]);
            }
            for (std::size_t i = 0; i < kFHEBatchSize; i++) {
                r[i] = mpz_class(a_1[i] * b_1[i]) - r[i];
                mpz_mod(r[i].get_mpz_t(), r[i].get_mpz_t(), mpz_class(kTwoPowerSixtyFour).get_mpz_t());
                c_1[i] = mpz_get_ui(r[i].get_mpz_t());
            }
            for (std::size_t i = 0; i < kFHEBatchSize; i++) {
                rand_triple_buff_[iter * kFHEBatchSize + i].emplace_back(a_1[i]);
                rand_triple_buff_[iter * kFHEBatchSize + i].emplace_back(b_1[i]);
                rand_triple_buff_[iter * kFHEBatchSize + i].emplace_back(c_1[i]);
            }
        }
    }
    return;
}

}  // namespace duet
}  // namespace petace
