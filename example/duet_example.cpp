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

    petace::duet::PlainMatrix<double> plain(rows, cols);
    petace::duet::PlainMatrix<double> ret(rows, cols);
    for (std::size_t i = 0; i < size; i++) {
        plain(i) = static_cast<double>(duet->party() == 0 ? i : (i * 99 + 1));
    }
    std::cout << "party " << duet->party() << " input is: " << std::endl;
    std::cout << plain << std::endl;
    petace::duet::ArithMatrix cipher0(rows, cols);
    petace::duet::ArithMatrix cipher1(rows, rows);
    petace::duet::ArithMatrix cipher2(rows, rows);

    duet->share(net, 0, plain, cipher0);
    duet->share(net, 1, plain, cipher1);
    duet->add(cipher0, cipher1, cipher2);

    duet->reveal(net, 0, cipher2, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << ret << std::endl;
    }
}

void mul_test(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 2;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PlainMatrix<double> plain(rows, cols);
    petace::duet::PlainMatrix<double> ret(rows, cols);
    for (std::size_t i = 0; i < size; i++) {
        plain(i) = static_cast<double>(duet->party() == 0 ? i : (i * 99 + 1));
    }
    std::cout << "party " << duet->party() << " input is: " << std::endl;
    std::cout << plain << std::endl;
    petace::duet::ArithMatrix cipher0(rows, cols);
    petace::duet::ArithMatrix cipher1(rows, rows);
    petace::duet::ArithMatrix cipher2(rows, rows);

    duet->share(net, 0, plain, cipher0);
    duet->share(net, 1, plain, cipher1);
    duet->elementwise_mul(net, cipher0, cipher1, cipher2);

    duet->reveal(net, 0, cipher2, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << ret << std::endl;
    }
}

void plain_shuffle_test(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 4;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PlainMatrix<double> plain(rows, cols);
    petace::duet::PlainMatrix<double> ret(rows, cols);
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < size; i++) {
            plain(i) = static_cast<double>(i);
        }
        std::cout << "party " << duet->party() << " input is: " << std::endl;
        std::cout << plain << std::endl;
    }

    petace::duet::ArithMatrix cipher(rows, cols);
    duet->shuffle(net, 0, plain, cipher);
    duet->reveal(net, 0, cipher, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << ret << std::endl;
    }
}

void less_than_zero_test(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 4;
    std::size_t cols = 2;
    std::size_t size = rows * cols;

    petace::duet::PlainMatrix<double> plain(rows, cols);
    petace::duet::PlainMatrix<std::int64_t> ret(rows, cols);
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < size; i++) {
            plain(i) = static_cast<double>(i);
        }
        std::cout << "party " << duet->party() << " input is: " << std::endl;
        std::cout << plain << std::endl;
    }

    petace::duet::ArithMatrix share(rows, cols);
    petace::duet::BoolMatrix boolen_ret(rows, cols);

    duet->share(net, 0, plain, share);
    duet->less_than_zero(net, share, boolen_ret);
    duet->reveal_bool(net, 0, boolen_ret, ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << ret << std::endl;
    }
}

void millionaires(
        const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 1;
    std::size_t cols = 1;
    petace::duet::PlainMatrix<double> plain_W_Alice(rows, cols);
    petace::duet::PlainMatrix<double> plain_W_Bob(rows, cols);
    petace::duet::ArithMatrix share_W_Alice(rows, cols);
    petace::duet::ArithMatrix share_W_Bob(rows, cols);

    double input = 0;
    if (duet->party() == 0) {
        std::cout << "I'm Alice, input your value:" << std::endl;
        std::cin >> input;
        plain_W_Alice(0) = input;
    } else {
        std::cout << "I'm Bob, input your value:" << std::endl;
        std::cin >> input;
        plain_W_Bob(0) = input;
    }

    // secret share W_Alice and W_Bob
    duet->share(net, 0, plain_W_Alice, share_W_Alice);
    duet->share(net, 1, plain_W_Bob, share_W_Bob);

    // compute a secure comparison and get protocol result
    petace::duet::BoolMatrix boolen_ret(rows, cols);
    petace::duet::PlainMatrix<std::int64_t> ret(rows, cols);
    petace::duet::ArithMatrix A_minus_B(rows, cols);

    duet->sub(share_W_Alice, share_W_Bob, A_minus_B);
    duet->less_than_zero(net, A_minus_B, boolen_ret);
    // reveal result to Alice
    duet->reveal_bool(net, 0, boolen_ret, ret);
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

    petace::duet::PlainMatrix<double> plain_model(rows, cols);
    petace::duet::PlainMatrix<double> plain_input(rows, cols);
    petace::duet::PlainMatrix<double> plain_ret(1, 1);
    petace::duet::ArithMatrix share_model(rows, cols);
    petace::duet::ArithMatrix share_input(rows, cols);

    if (duet->party() == 0) {
        // Server
        // Server puts private model here
        std::cout << "Next, I'm server, input model params:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th params:" << std::endl;
            std::cin >> plain_model(i);
        }
    } else {
        // Client
        // Client puts private input here
        std::cout << "Next, I'm client, input sample params:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th params:" << std::endl;
            std::cin >> plain_input(i);
        }
    }
    // secret share the model and client input
    duet->share(net, 0, plain_model, share_model);
    duet->share(net, 1, plain_input, share_input);

    // inner product computation
    petace::duet::ArithMatrix element_mul_res(rows, cols);
    petace::duet::ArithMatrix res(1, 1);
    duet->elementwise_mul(net, share_model, share_input, element_mul_res);
    // column-wise sum, then the result secret share will be in res(0)
    duet->sum(element_mul_res, res);

    duet->reveal(net, 0, res, plain_ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << plain_ret << std::endl;
    }
}

void distance(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet) {
    std::size_t rows = 0;
    std::size_t cols = 1;
    std::cout << "Firstly, input vector length:" << std::endl;
    std::cin >> rows;
    petace::duet::PlainMatrix<double> plain_X(rows, cols);
    petace::duet::PlainMatrix<double> plain_Y(rows, cols);
    petace::duet::PlainMatrix<double> plain_ret(1, 1);
    petace::duet::ArithMatrix share_X(rows, cols);
    petace::duet::ArithMatrix share_Y(rows, cols);

    if (duet->party() == 0) {
        // party_0 puts his private coordinate here
        std::cout << "Next, I'm Alice, input vector:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th value:" << std::endl;
            std::cin >> plain_X(i);
        }
    } else {
        // party_1 puts his private coordinate here
        std::cout << "Next, I'm Bob, input vector:" << std::endl;
        for (std::size_t i = 0; i < rows; i++) {
            std::cout << "input " << i << "th value:" << std::endl;
            std::cin >> plain_Y(i);
        }
    }
    // secret share the model and client input
    duet->share(net, 0, plain_X, share_X);
    duet->share(net, 1, plain_Y, share_Y);

    petace::duet::ArithMatrix X_minus_Y(rows, cols);
    petace::duet::ArithMatrix X_minus_Y_square(rows, cols);
    petace::duet::ArithMatrix res(1, 1);
    duet->sub(share_X, share_Y, X_minus_Y);
    duet->elementwise_mul(net, X_minus_Y, X_minus_Y, X_minus_Y_square);
    // column-wise sum, then the result secret share will be in res(0)
    duet->sum(X_minus_Y_square, res);
    duet->reveal(net, 0, res, plain_ret);
    if (duet->party() == 0) {
        std::cout << "reveal to party 0, ret is: " << std::endl;
        std::cout << plain_ret << std::endl;
    }
}
