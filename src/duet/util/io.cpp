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

#include "duet/util/io.h"

#include <vector>

namespace petace {
namespace duet {

void send_block(const std::shared_ptr<network::Network>& net, const block* data, std::size_t nblock) {
    net->send_data(data, nblock * sizeof(block));
}

void recv_block(const std::shared_ptr<network::Network>& net, block* data, std::size_t nblock) {
    net->recv_data(data, nblock * sizeof(block));
}

void send_bool(const std::shared_ptr<network::Network>& net, bool* data, std::size_t length) {
    std::size_t byte_len = (length + 7) / 8;
    std::vector<uint8_t> byte_buf(byte_len, 0);
    for (std::size_t i = 0; i < length; i++) {
        byte_buf[i / 8] |= static_cast<uint8_t>(static_cast<uint8_t>(data[i]) << (7 - (i % 8)));
    }
    net->send_data(byte_buf.data(), byte_len);
}

void recv_bool(const std::shared_ptr<network::Network>& net, bool* data, std::size_t length) {
    std::size_t byte_len = (length + 7) / 8;
    std::vector<uint8_t> byte_buf(byte_len);
    net->recv_data(byte_buf.data(), byte_len);
    for (std::size_t i = 0; i < length; i++) {
        data[i] = static_cast<bool>((byte_buf[i / 8] >> (7 - (i % 8))) & 0x1);
    }
}

void send_matrix(const std::shared_ptr<network::Network>& net, const Matrix<std::int64_t>* data, std::size_t nmatrix) {
    net->send_data(data->data(), nmatrix * data->size() * sizeof(std::int64_t));
}

void recv_matrix(const std::shared_ptr<network::Network>& net, Matrix<std::int64_t>* data, std::size_t nmatrix) {
    std::vector<std::int64_t> buf(data->size() * nmatrix);
    net->recv_data(&buf[0], nmatrix * data->size() * sizeof(std::int64_t));
    *data = Eigen::Map<Matrix<std::int64_t>>(buf.data(), data->rows(), data->cols());
}

void send_cipher(
        const std::shared_ptr<network::Network>& net, const PaillierMatrix& in, std::size_t paillier_key_size) {
    std::size_t data_size;
    std::size_t party = in.party();
    std::size_t row = in.rows();
    std::size_t col = in.cols();
    net->send_data(&row, sizeof(std::size_t));
    net->send_data(&col, sizeof(std::size_t));
    net->send_data(&party, sizeof(std::size_t));

    ByteVector send_buffer(row * col * kPaillierCipherSize);
    for (std::size_t i = 0; i < in.size(); i++) {
        in(i).serialize_to_bytes(send_buffer.data() + i * kPaillierCipherSize, kPaillierCipherSize);
    }
    data_size = send_buffer.size();
    net->send_data(&data_size, sizeof(std::size_t));
    net->send_data(send_buffer.data(), data_size * sizeof(solo::Byte));

    return;
}

void recv_cipher(const std::shared_ptr<network::Network>& net, const std::shared_ptr<solo::ahepaillier::PublicKey>& pk,
        PaillierMatrix& out, std::size_t paillier_key_size) {
    std::size_t party;
    std::size_t row;
    std::size_t col;
    std::size_t data_size;
    petace::solo::ahepaillier::Encoder encoder;
    net->recv_data(&row, sizeof(std::size_t));
    net->recv_data(&col, sizeof(std::size_t));
    net->recv_data(&party, sizeof(std::size_t));

    out.resize(row, col);
    out.set_party(party);

    ByteVector receive_buffer;
    net->recv_data(&data_size, sizeof(std::size_t));
    receive_buffer.resize(data_size);
    net->recv_data(receive_buffer.data(), data_size * sizeof(solo::Byte));

    for (std::size_t i = 0; i < row * col; i++) {
        out(i).deserialize_from_bytes(pk, receive_buffer.data() + i * kPaillierCipherSize, kPaillierCipherSize);
    }
}

}  // namespace duet
}  // namespace petace
