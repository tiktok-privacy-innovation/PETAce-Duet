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

#include "Eigen/Dense"

#include "network/network.h"

#include "duet/util/common.h"
#include "duet/util/consts.h"
#include "duet/util/defines.h"
#include "duet/util/matrix.h"

namespace petace {
namespace duet {

void send_block(const std::shared_ptr<network::Network>& net, const block* data, std::size_t nblock);

void recv_block(const std::shared_ptr<network::Network>& net, block* data, std::size_t nblock);

void send_bool(const std::shared_ptr<network::Network>& net, bool* data, std::size_t length);

void recv_bool(const std::shared_ptr<network::Network>& net, bool* data, std::size_t length);

void send_matrix(const std::shared_ptr<network::Network>& net, const Matrix<std::int64_t>* data, std::size_t nmatrix);

void recv_matrix(const std::shared_ptr<network::Network>& net, Matrix<std::int64_t>* data, std::size_t nmatrix);

void send_cipher(const std::shared_ptr<network::Network>& net, const PaillierMatrix& in, std::size_t paillier_key_size);

void recv_cipher(const std::shared_ptr<network::Network>& net, const std::shared_ptr<solo::ahepaillier::PublicKey>& pk,
        PaillierMatrix& out, std::size_t paillier_key_size);

}  // namespace duet
}  // namespace petace
