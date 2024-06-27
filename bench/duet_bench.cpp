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

#include "duet_bench.h"

#include <chrono>
#include <iomanip>
#include <stdexcept>

#include "glog/logging.h"

double get_unix_timestamp() {
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();

    std::chrono::duration<double> duration_since_epoch =
            std::chrono::duration_cast<std::chrono::duration<double>>(now.time_since_epoch());

    return duration_since_epoch.count();
}

void mul_bench(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet,
        std::size_t test_number, std::size_t rows, std::size_t cols) {
    try {
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

        petace::duet::ArithMatrix cipher0(rows, cols);
        petace::duet::ArithMatrix cipher1(rows, cols);
        petace::duet::ArithMatrix cipher2(rows, cols);

        duet->share(net, private_0, cipher0);
        duet->share(net, private_1, cipher1);

        double begin = get_unix_timestamp();
        LOG(INFO) << std::fixed << "case mul_" << rows << "_" << cols << "_bench"
                  << " begin " << begin << " " << test_number << " " << rows << " " << cols;
        for (size_t i = 0; i < test_number; i++) {
            duet->elementwise_mul(net, cipher0, cipher1, cipher2);
        }

        double end = get_unix_timestamp();

        LOG(INFO) << std::fixed << "case mul_" << rows << "_" << cols << "_bench"
                  << " end " << end << " " << end - begin << "s " << net->get_bytes_sent() << " "
                  << net->get_bytes_received();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

void inner_product_bench(const std::shared_ptr<petace::network::Network>& net,
        const std::shared_ptr<petace::duet::Duet>& duet, std::size_t test_number, std::size_t vector_size) {
    petace::duet::PrivateMatrix<double> private_0(vector_size, 1, 0);
    petace::duet::PrivateMatrix<double> private_1(vector_size, 1, 1);
    petace::duet::PrivateMatrix<double> ret(vector_size, 1, 0);
    double begin = 0;
    double end = 0;
    double case_begin = 0;
    std::uint64_t bytes_sent_now = 0;
    std::uint64_t bytes_received_now = 0;
    if (duet->party() == 0) {
        for (std::size_t i = 0; i < vector_size; i++) {
            private_0(i) = static_cast<double>(i);
        }
    } else {
        for (std::size_t i = 0; i < vector_size; i++) {
            private_1(i) = static_cast<double>(i * 99 + 1);
        }
    }

    petace::duet::ArithMatrix cipher0(vector_size, 1);
    petace::duet::ArithMatrix cipher1(vector_size, 1);
    petace::duet::ArithMatrix cipher2(vector_size, 1);
    petace::duet::ArithMatrix res(1, 1);

    duet->share(net, private_0, cipher0);
    duet->share(net, private_1, cipher1);

    begin = get_unix_timestamp();
    case_begin = begin;
    LOG(INFO) << std::fixed << "case inner_product_" << vector_size << "_bench"
              << " begin " << begin << " " << test_number << " " << vector_size;
    LOG(INFO) << std::fixed << "stage mul begin " << begin;
    // mul
    for (size_t i = 0; i < test_number; i++) {
        duet->elementwise_mul(net, cipher0, cipher1, cipher2);
    }
    end = get_unix_timestamp();
    bytes_sent_now = net->get_bytes_sent();
    bytes_received_now = net->get_bytes_received();
    LOG(INFO) << std::fixed << "stage mul end " << end << " " << end - begin << "s " << bytes_sent_now << " "
              << bytes_received_now;

    begin = get_unix_timestamp();
    LOG(INFO) << std::fixed << "stage sum begin " << begin;
    // sum
    for (size_t i = 0; i < test_number; i++) {
        duet->sum(cipher2, res);
    }
    end = get_unix_timestamp();
    bytes_sent_now = net->get_bytes_sent();
    bytes_received_now = net->get_bytes_received();
    LOG(INFO) << std::fixed << "stage sum end " << end << " " << end - begin << "s " << bytes_sent_now << " "
              << bytes_received_now;

    LOG(INFO) << std::fixed << "case inner_product_" << vector_size << "_bench"
              << " end " << end << " " << end - case_begin << "s " << bytes_sent_now << " " << bytes_received_now;
}

void equal_bench(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet,
        std::size_t test_number, std::size_t rows, std::size_t cols) {
    try {
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
                private_1(i) = static_cast<double>(i + i % 2);
            }
        }

        petace::duet::ArithMatrix cipher0(rows, cols);
        petace::duet::ArithMatrix cipher1(rows, cols);
        petace::duet::BoolMatrix cipher2(rows, cols);

        duet->share(net, private_0, cipher0);
        duet->share(net, private_1, cipher1);

        double begin = get_unix_timestamp();
        LOG(INFO) << std::fixed << "case equal_" << rows << "_" << cols << "_bench"
                  << " begin " << begin << " " << test_number << " " << rows << " " << cols;
        for (size_t i = 0; i < test_number; i++) {
            duet->equal(net, cipher0, cipher1, cipher2);
        }

        double end = get_unix_timestamp();

        LOG(INFO) << std::fixed << "case equal_" << rows << "_" << cols << "_bench"
                  << " end " << end << " " << end - begin << "s " << net->get_bytes_sent() << " "
                  << net->get_bytes_received();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

void greater_bench(const std::shared_ptr<petace::network::Network>& net,
        const std::shared_ptr<petace::duet::Duet>& duet, std::size_t test_number, std::size_t rows, std::size_t cols) {
    try {
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
                private_1(i) = static_cast<double>(i + i % 2);
            }
        }

        petace::duet::ArithMatrix cipher0(rows, cols);
        petace::duet::ArithMatrix cipher1(rows, cols);
        petace::duet::BoolMatrix cipher2(rows, cols);

        duet->share(net, private_0, cipher0);
        duet->share(net, private_1, cipher1);

        double begin = get_unix_timestamp();
        LOG(INFO) << std::fixed << "case greater_" << rows << "_" << cols << "_bench"
                  << " begin " << begin << " " << test_number << " " << rows << " " << cols;
        for (size_t i = 0; i < test_number; i++) {
            duet->greater(net, cipher0, cipher1, cipher2);
        }

        double end = get_unix_timestamp();

        LOG(INFO) << std::fixed << "case greater_" << rows << "_" << cols << "_bench"
                  << " end " << end << " " << end - begin << "s " << net->get_bytes_sent() << " "
                  << net->get_bytes_received();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

void less_bench(const std::shared_ptr<petace::network::Network>& net, const std::shared_ptr<petace::duet::Duet>& duet,
        std::size_t test_number, std::size_t rows, std::size_t cols) {
    try {
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
                private_1(i) = static_cast<double>(i + i % 2);
            }
        }

        petace::duet::ArithMatrix cipher0(rows, cols);
        petace::duet::ArithMatrix cipher1(rows, cols);
        petace::duet::BoolMatrix cipher2(rows, cols);

        duet->share(net, private_0, cipher0);
        duet->share(net, private_1, cipher1);

        double begin = get_unix_timestamp();
        LOG(INFO) << std::fixed << "case less_" << rows << "_" << cols << "_bench"
                  << " begin " << begin << " " << test_number << " " << rows << " " << cols;
        for (size_t i = 0; i < test_number; i++) {
            duet->less(net, cipher0, cipher1, cipher2);
        }

        double end = get_unix_timestamp();

        LOG(INFO) << std::fixed << "case less_" << rows << "_" << cols << "_bench"
                  << " end " << end << " " << end - begin << "s " << net->get_bytes_sent() << " "
                  << net->get_bytes_received();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

void less_than_zero_bench(const std::shared_ptr<petace::network::Network>& net,
        const std::shared_ptr<petace::duet::Duet>& duet, std::size_t test_number, std::size_t rows, std::size_t cols) {
    try {
        std::size_t size = rows * cols;

        petace::duet::PrivateMatrix<double> private_0(rows, cols, 0);
        petace::duet::PrivateMatrix<double> ret(rows, cols, 0);
        if (duet->party() == 0) {
            for (std::size_t i = 0; i < size; i++) {
                private_0(i) = static_cast<double>(i);
            }
        }

        petace::duet::ArithMatrix cipher0(rows, cols);
        petace::duet::BoolMatrix cipher2(rows, cols);

        duet->share(net, private_0, cipher0);

        double begin = get_unix_timestamp();
        LOG(INFO) << std::fixed << "case less_than_zero_" << rows << "_" << cols << "_bench"
                  << " begin " << begin << " " << test_number << " " << rows << " " << cols;
        for (size_t i = 0; i < test_number; i++) {
            duet->less_than_zero(net, cipher0, cipher2);
        }

        double end = get_unix_timestamp();

        LOG(INFO) << std::fixed << "case less_than_zero_" << rows << "_" << cols << "_bench"
                  << " end " << end << " " << end - begin << "s " << net->get_bytes_sent() << " "
                  << net->get_bytes_received();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}
