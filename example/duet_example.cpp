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

#include "duet_example.h"

#include <iostream>

void add_test(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 2;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PrivateMatrix<double> private_0(rows, cols, 0);
    petace::duet::PrivateMatrix<double> private_1(rows, cols, 1);
    petace::duet::PrivateMatrix<double> ret(rows, cols, 0);
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < size; i++) {
            private_0(i) = static_cast<double>(i);
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            private_1(i) = static_cast<double>(i * 99 + 1);
        }
    }

    std::cout << "party " << duet->party() << " input is: " << std::endl;
    if (duet->party() == 0) {
        std::cout << private_0.matrix() << std::endl;
    } else {
        std::cout << private_1.matrix() << std::endl;
    }

    petace::duet::ArithMatrix cipher0(rows, cols);
    petace::duet::ArithMatrix cipher1(rows, rows);
    petace::duet::ArithMatrix cipher2(rows, rows);

    duet->share(net, private_0, cipher0);
    duet->share(net, private_1, cipher1);
    duet->add(cipher0, cipher1, cipher2);

    duet->reveal(net, cipher2, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party " << ret.party_id() << ", ret is: " << std::endl;
        std::cout << ret.matrix() << std::endl;
    }
}

void mul_test(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 2;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PrivateMatrix<double> private_0(rows, cols, 0);
    petace::duet::PrivateMatrix<double> private_1(rows, cols, 1);
    petace::duet::PrivateMatrix<double> ret(rows, cols, 0);
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < size; i++) {
            private_0(i) = static_cast<double>(i);
        }
    } else {
        for (std::size_t i = 0; i < size; i++) {
            private_1(i) = static_cast<double>(i * 99 + 1);
        }
    }

    std::cout << "party " << duet->party() << " input is: " << std::endl;
    if (duet->party() == 0) {
        std::cout << private_0.matrix() << std::endl;
    } else {
        std::cout << private_1.matrix() << std::endl;
    }

    petace::duet::ArithMatrix cipher0(rows, cols);
    petace::duet::ArithMatrix cipher1(rows, rows);
    petace::duet::ArithMatrix cipher2(rows, rows);

    duet->share(net, private_0, cipher0);
    duet->share(net, private_1, cipher1);
    duet->elementwise_mul(net, cipher0, cipher1, cipher2);

    duet->reveal(net, cipher2, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party " << ret.party_id() << ", ret is: " << std::endl;
        std::cout << ret.matrix() << std::endl;
    }
}

void private_shuffle_test(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 4;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PrivateMatrix<double> private_0(rows, cols, 0);
    petace::duet::PrivateMatrix<double> ret(rows, cols, 0);
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < size; i++) {
            private_0(i) = static_cast<double>(i);
        }
        std::cout << "party " << duet->party() << " input is: " << std::endl;
        std::cout << private_0.matrix() << std::endl;
    }

    petace::duet::ArithMatrix cipher(rows, cols);
    duet->shuffle(net, private_0, cipher);
    duet->reveal(net, cipher, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party " << ret.party_id() << ", ret is: " << std::endl;
        std::cout << ret.matrix() << std::endl;
    }
}

void less_than_zero_test(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 4;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PrivateMatrix<double> private_0(rows, cols, 0);
    petace::duet::PrivateMatrixBool ret(rows, cols, 0);
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < size; i++) {
            private_0(i) = static_cast<double>(i);
        }
        std::cout << "party " << duet->party() << " input is: " << std::endl;
        std::cout << private_0.matrix() << std::endl;
    }

    petace::duet::ArithMatrix share(rows, cols);
    petace::duet::BoolMatrix boolen_ret(rows, cols);

    duet->share(net, private_0, share);
    duet->less_than_zero(net, share, boolen_ret);
    duet->reveal_bool(net, boolen_ret, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party " << ret.party_id() << ", ret is: " << std::endl;
        std::cout << ret.matrix() << std::endl;
    }
}

void millionaires(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 1;
    std::size_t cols = 1;
    petace::duet::PrivateMatrix<double> private_w_alice(rows, cols, 0);
    petace::duet::PrivateMatrix<double> private_w_bob(rows, cols, 1);
    petace::duet::ArithMatrix share_w_alice(rows, cols);
    petace::duet::ArithMatrix share_w_bob(rows, cols);

    double input = 0;
    if (duet->party() == 0) {
        std::cout << "I'm Alice, input your value:" << std::endl;
        std::cin >> input;
        private_w_alice(0) = input;
    } else {
        std::cout << "I'm Bob, input your value:" << std::endl;
        std::cin >> input;
        private_w_bob(0) = input;
    }

    // secret share W_Alice and W_Bob
    duet->share(net, private_w_alice, share_w_alice);
    duet->share(net, private_w_bob, share_w_bob);

    // compute a secure comparison and get protocol result
    petace::duet::BoolMatrix boolen_ret(rows, cols);
    petace::duet::PrivateMatrixBool ret;
    petace::duet::ArithMatrix a_minus_b(rows, cols);

    duet->sub(share_w_alice, share_w_bob, a_minus_b);
    duet->less_than_zero(net, a_minus_b, boolen_ret);
    // reveal result to Alice
    duet->reveal_bool(net, boolen_ret, ret);
    // print the result on Alice side
    if (duet->party() == 0) {
        if (ret(0) == 1) {
            std::cout << "Bob win!" << std::endl;
        } else {
            std::cout << "Alice win!" << std::endl;
        }
    }
}

void ppml(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 0;
    std::size_t cols = 1;
    std::cout << "Please set model parmas!" << std::endl;
    std::cout << "Firstly, input parmas number:" << std::endl;
    std::cin >> rows;

    petace::duet::PrivateMatrix<double> private_model(rows, cols, 0);
    petace::duet::PrivateMatrix<double> private_input(rows, cols, 1);
    petace::duet::PrivateMatrix<double> private_ret;
    petace::duet::ArithMatrix share_model(rows, cols);
    petace::duet::ArithMatrix share_input(rows, cols);

    if (duet->party() == 0) {
        // Server
        // Server puts private model here
        std::cout << "Next, I'm server, input model params:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th params:" << std::endl;
            std::cin >> private_model(i);
        }
    } else {
        // Client
        // Client puts private input here
        std::cout << "Next, I'm client, input sample params:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th params:" << std::endl;
            std::cin >> private_input(i);
        }
    }
    // secret share the model and client input
    duet->share(net, private_model, share_model);
    duet->share(net, private_input, share_input);

    // inner product computation
    petace::duet::ArithMatrix element_mul_res(rows, cols);
    petace::duet::ArithMatrix res(1, 1);
    duet->elementwise_mul(net, share_model, share_input, element_mul_res);
    // column-wise sum, then the result secret share will be in res(0)
    duet->sum(element_mul_res, res);

    duet->reveal(net, res, private_ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << private_ret.matrix() << std::endl;
    }
}

void distance(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 0;
    std::size_t cols = 1;
    std::cout << "Firstly, input vector length:" << std::endl;
    std::cin >> rows;
    petace::duet::PrivateMatrix<double> private_x(rows, cols, 0);
    petace::duet::PrivateMatrix<double> private_y(rows, cols, 1);
    petace::duet::PrivateMatrix<double> private_ret;
    petace::duet::ArithMatrix share_x(rows, cols);
    petace::duet::ArithMatrix share_y(rows, cols);

    if (duet->party() == 0) {
        // party_0 puts his private coordinate here
        std::cout << "Next, I'm Alice, input vector:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th value:" << std::endl;
            std::cin >> private_x(i);
        }
    } else {
        // party_1 puts his private coordinate here
        std::cout << "Next, I'm Bob, input vector:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th value:" << std::endl;
            std::cin >> private_y(i);
        }
    }
    // secret share the model and client input
    duet->share(net, private_x, share_x);
    duet->share(net, private_y, share_y);

    petace::duet::ArithMatrix x_minus_y(rows, cols);
    petace::duet::ArithMatrix x_minus_y_square(rows, cols);
    petace::duet::ArithMatrix res(1, 1);
    duet->sub(share_x, share_y, x_minus_y);
    duet->elementwise_mul(net, x_minus_y, x_minus_y, x_minus_y_square);
    // column-wise sum, then the result secret share will be in res(0)
    duet->sum(x_minus_y_square, res);
    duet->reveal(net, res, private_ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << private_ret.matrix() << std::endl;
    }
}

void dot_product_paillier(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 0;
    std::size_t cols = 1;
    std::cout << "Firstly, input vector length:" << std::endl;
    std::cin >> rows;
    petace::duet::PrivateMatrix<int64_t> private_x(rows, cols, 0);
    petace::duet::PrivateMatrix<int64_t> private_y(rows, cols, 1);
    petace::duet::PrivateMatrix<int64_t> private_ret(rows, cols, 0);
    petace::duet::PaillierMatrix cipher_x(rows, cols);
    petace::duet::PaillierMatrix cipher_x_received(rows, cols);
    petace::duet::PaillierMatrix cipher_result(rows, cols);
    if (duet->party() == 0) {
        // party_0 puts his private coordinate here
        std::cout << "Next, I'm Alice, input vector: (please use integer)" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th value:" << std::endl;
            std::cin >> private_x(i);
        }
        duet->encrypt(private_x, cipher_x);
        petace::duet::send_cipher(net, cipher_x, petace::duet::kPaillierKeySize);
        petace::duet::recv_cipher(net, duet->get_pk(), cipher_result, petace::duet::kPaillierKeySize);
        duet->decrypt(cipher_result, private_ret);
        std::int64_t result = 0;
        for (std::size_t i = 0; i < rows; i++) {
            result += private_ret(i);
        }
        std::cout << "result is available to party 0, result is: " << result << std::endl;
    } else {
        // party_1 puts his private coordinate here
        std::cout << "Next, I'm Bob, input vector: (please use integer)" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th value:" << std::endl;
            std::cin >> private_y(i);
        }
        petace::duet::recv_cipher(net, duet->get_pk_other(), cipher_x_received, petace::duet::kPaillierKeySize);
        duet->mul(private_y, cipher_x_received, cipher_result);
        petace::duet::send_cipher(net, cipher_result, petace::duet::kPaillierKeySize);
    }
}
