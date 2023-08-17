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

#include "duet/duet.h"

#include <stdlib.h>
#include <unistd.h>

#include <memory>
#include <thread>
#include <utility>

#include "gtest/gtest.h"

#include "network/net_factory.h"
#include "network/net_socket.h"
#include "network/network.h"
#include "solo/prng.h"
#include "verse/verse_factory.h"

#include "duet/util/common.h"

namespace petace {
namespace duet {

class AbyTest : public ::testing::Test {
public:
    void SetUp() {
        net_params0.remote_addr = "127.0.0.1";
        net_params0.remote_port = 8890;
        net_params0.local_port = 8891;

        net_params1.remote_addr = "127.0.0.1";
        net_params1.remote_port = 8891;
        net_params1.local_port = 8890;
    }

    bool is_equal_plain_matrix(const PlainMatrix<double>& m1, const PlainMatrix<double>& m2, double e) {
        if ((m1.rows() != m2.rows()) || (m1.cols() != m2.cols())) {
            return false;
        }

        for (std::size_t i = 0; i < static_cast<std::size_t>(m1.size()); i++) {
            double d = m1(i) - m2(i);
            if ((d < -e) || (d > e)) {
                return false;
            }
        }

        return true;
    }

    bool is_equal_plain_matrix(const PlainMatrix<std::int64_t>& m1, const PlainMatrix<std::int64_t>& m2) {
        if ((m1.rows() != m2.rows()) || (m1.cols() != m2.cols())) {
            return false;
        }

        for (std::size_t i = 0; i < static_cast<std::size_t>(m1.size()); i++) {
            int64_t d = m1(i) - m2(i);
            if (d != 0) {
                return false;
            }
        }

        return true;
    }

    double get_rand_double() {
        double res;
        double threshold = 1.0 * (1 << 16);

        do {
            res = 1.0 * static_cast<double>(random()) / static_cast<double>(random());
            if (random() & 1) {
                res = -res;
            }
        } while ((res > threshold) || (res < -threshold));

        return res;
    }

    void get_rand_plain_matrix(PlainMatrix<double>& plain) {
        for (std::size_t i = 0; i < static_cast<std::size_t>(plain.size()); i++) {
            plain(i) = get_rand_double();
        }
    }

    std::int64_t get_rand_int64() {
        std::int64_t res;
        std::int64_t threshold = (1 << 16);
        do {
            res = static_cast<std::int64_t>(random());
            if (random() & 1) {
                res = -res;
            }
        } while ((res > threshold) || (res < -threshold));
        return res;
    }

    void get_rand_plain_matrix(PlainMatrix<std::int64_t>& plain) {
        for (std::size_t i = 0; i < static_cast<std::size_t>(plain.size()); i++) {
            plain(i) = get_rand_int64();
        }
    }

    void add_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<double> plain_c(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);

        aby_test->add(cipher_a, cipher_b, cipher_c);

        aby_test->reveal(net, 0, cipher_c, reveal_plain_);
        plain_ = plain_a + plain_b;

        return;
    }

    void sub_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<double> plain_c(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);

        aby_test->sub(cipher_a, cipher_b, cipher_c);

        aby_test->reveal(net, 0, cipher_c, reveal_plain_);
        plain_ = plain_a - plain_b;

        return;
    }

    void mul_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<double> plain_c(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);

        aby_test->elementwise_mul(net, cipher_a, cipher_b, cipher_c);

        aby_test->reveal(net, 0, cipher_c, reveal_plain_);
        plain_ = plain_a.cwiseProduct(plain_b);

        return;
    }

    void greater_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<double> plain_c(4, 3);
        PlainMatrix<double> plain_d(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);
        ArithMatrix cipher_d(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);
        get_rand_plain_matrix(plain_d);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);
        aby_test->share(net, 0, plain_d, cipher_d);

        aby_test->greater(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);

        reveal_plain_.resize(4, 3);
        aby_test->reveal(net, 0, cipher_b, reveal_plain_);

        plain_.resize(4, 3);
        for (std::size_t i = 0; i < static_cast<std::size_t>(plain_d.size()); i++) {
            plain_(i) = (plain_a(i) > plain_b(i)) ? plain_d(i) : 0;
        }

        return;
    }

    void greater_plain_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<std::int64_t> plain_c(4, 3);
        PlainMatrix<double> plain_d(4, 3);
        ArithMatrix cipher_a(4, 3);
        BoolMatrix cipher_c(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->greater(net, cipher_a, plain_b, cipher_c);
        aby_test->reveal_bool(net, 0, cipher_c, plain_c);

        plain_bool_.resize(4, 3);
        reveal_plain_bool_.resize(4, 3);
        for (std::size_t i = 0; i < static_cast<std::size_t>(plain_d.size()); i++) {
            plain_bool_(i) = (plain_a(i) > plain_b(i)) ? 1 : 0;
        }
        if (party_id == 0) {
            for (std::size_t i = 0; i < static_cast<std::size_t>(plain_c.size()); i++) {
                reveal_plain_bool_(i) = plain_c(i);
            }
        }

        return;
    }

    void less_plain_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<std::int64_t> plain_c(4, 3);
        PlainMatrix<double> plain_d(4, 3);
        ArithMatrix cipher_a(4, 3);
        BoolMatrix cipher_c(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->less(net, cipher_a, plain_b, cipher_c);
        aby_test->reveal_bool(net, 0, cipher_c, plain_c);

        plain_bool_.resize(4, 3);
        reveal_plain_bool_.resize(4, 3);
        for (std::size_t i = 0; i < static_cast<std::size_t>(plain_d.size()); i++) {
            plain_bool_(i) = (plain_a(i) < plain_b(i)) ? 1 : 0;
        }
        if (party_id == 0) {
            for (std::size_t i = 0; i < static_cast<std::size_t>(plain_c.size()); i++) {
                reveal_plain_bool_(i) = plain_c(i);
            }
        }

        return;
    }

    void less_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<double> plain_c(4, 3);
        PlainMatrix<double> plain_d(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);
        ArithMatrix cipher_d(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);
        get_rand_plain_matrix(plain_d);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);
        aby_test->share(net, 0, plain_d, cipher_d);

        aby_test->less(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);

        reveal_plain_.resize(4, 3);
        aby_test->reveal(net, 0, cipher_b, reveal_plain_);

        plain_.resize(4, 3);
        for (std::size_t i = 0; i < static_cast<std::size_t>(plain_d.size()); i++) {
            plain_(i) = (plain_a(i) < plain_b(i)) ? plain_d(i) : 0;
        }

        return;
    }

    void equal_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        PlainMatrix<double> plain_c(4, 3);
        PlainMatrix<double> plain_d(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);
        ArithMatrix cipher_d(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);
        get_rand_plain_matrix(plain_d);

        plain_b(0) = plain_a(0);
        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);
        aby_test->share(net, 0, plain_d, cipher_d);

        aby_test->equal(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);

        reveal_plain_.resize(4, 3);
        aby_test->reveal(net, 0, cipher_b, reveal_plain_);

        plain_.resize(4, 3);

        for (std::size_t i = 0; i < static_cast<std::size_t>(plain_d.size()); i++) {
            plain_(i) = (plain_a(i) == plain_b(i)) ? plain_d(i) : 0;
        }

        return;
    }

    void less_than_zero_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        const std::size_t num_element = 40;

        PlainMatrix<std::int64_t> plain_x(1, num_element);
        PlainMatrix<std::int64_t> reveal_plain_x(1, num_element);
        ArithMatrix share_x(1, num_element);
        BoolMatrix res(1, num_element);

        get_rand_plain_matrix(plain_x);
        plain_bool_.resize(1, num_element);
        reveal_plain_bool_.resize(1, num_element);

        std::size_t plain_x_size = plain_x.size();
        for (std::size_t i = 0; i < plain_x_size; i++) {
            plain_bool_(i) = (plain_x(i) < 0) ? 1 : 0;
        }

        aby_test->share(net, 0, plain_x, share_x);
        aby_test->less_than_zero(net, share_x, res);
        aby_test->reveal_bool(net, 0, res, reveal_plain_x);
        if (party_id == 0) {
            for (std::size_t i = 0; i < plain_x_size; i++) {
                reveal_plain_bool_(i) = reveal_plain_x(i);
            }
        }

        return;
    }

    void sum_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PlainMatrix<double> plain_a(4, 3);
        PlainMatrix<double> plain_b(4, 3);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);

        get_rand_plain_matrix(plain_a);
        get_rand_plain_matrix(plain_b);

        aby_test->share(net, 0, plain_a, cipher_a);
        aby_test->share(net, 0, plain_b, cipher_b);

        aby_test->sum(cipher_a, cipher_b);

        reveal_plain_.resize(1, 3);
        aby_test->reveal(net, 0, cipher_b, reveal_plain_);

        plain_.resize(1, 3);
        plain_ = plain_a.colwise().sum();

        return;
    }

    void plain_shuffle_test(std::size_t party_id, PlainMatrix<std::int64_t>& plain, PlainMatrix<std::int64_t>& share0,
            PlainMatrix<std::int64_t>& share1) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        std::size_t rows = 10;
        std::size_t cols = 3;
        std::size_t matrix_size = rows * cols;

        PlainMatrix<double> plain_a(rows, cols);
        ArithMatrix cipher_a(rows, cols);
        if (party_id == 0) {
            get_rand_plain_matrix(plain_a);
            plain.resize(rows, cols);
            for (std::size_t i = 0; i < matrix_size; ++i) {
                plain(i) = float_to_fixed(plain_a(i));
            }
        }
        aby_test->shuffle(net, 0, plain_a, cipher_a);
        if (party_id == 0) {
            share0 = cipher_a.shares;
        } else {
            share1 = cipher_a.shares;
        }
        return;
    }

    void share_shuffle_test(std::size_t party_id, PlainMatrix<std::int64_t>& in0, PlainMatrix<std::int64_t>& in1,
            PlainMatrix<std::int64_t>& out0, PlainMatrix<std::int64_t>& out1) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        std::size_t rows = 10;
        std::size_t cols = 3;
        std::size_t matrix_size = rows * cols;
        ArithMatrix cipher_a(rows, cols);
        ArithMatrix cipher_b;
        for (std::size_t i = 0; i < matrix_size; ++i) {
            cipher_a.shares(i) = random();
        }

        aby_test->shuffle(net, cipher_a, cipher_b);
        if (party_id == 0) {
            in0 = cipher_a.shares;
            out0 = cipher_b.shares;
        } else {
            in1 = cipher_a.shares;
            out1 = cipher_b.shares;
        }
    }

public:
    PlainMatrix<double> reveal_plain_;
    PlainMatrix<double> plain_;
    PlainMatrix<std::int64_t> reveal_plain_bool_;
    PlainMatrix<std::int64_t> plain_bool_;
    petace::network::NetParams net_params0;
    petace::network::NetParams net_params1;
};

TEST_F(AbyTest, add_test) {
    pid_t pid;
    int status;

    pid = fork();
    if (pid < 0) {
        status = -1;
        exit(EXIT_FAILURE);
    } else if (pid == 0) {
        add_test(1);
        exit(EXIT_SUCCESS);
    } else {
        add_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
        return;
    }
}

TEST_F(AbyTest, sub_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        sub_test(1);
        exit(EXIT_SUCCESS);
    } else {
        sub_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
    }
}

TEST_F(AbyTest, mul_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        mul_test(1);
        exit(EXIT_SUCCESS);
    } else {
        mul_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
    }
}

TEST_F(AbyTest, greater_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        greater_test(1);
        exit(EXIT_SUCCESS);
    } else {
        greater_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
    }
}

TEST_F(AbyTest, greater_plain_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        greater_plain_test(1);
        exit(EXIT_SUCCESS);
    } else {
        greater_plain_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_bool_, plain_bool_));
    }
}

TEST_F(AbyTest, less_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        less_test(1);
        exit(EXIT_SUCCESS);
    } else {
        less_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
    }
}

TEST_F(AbyTest, less_plain_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        less_plain_test(1);
        exit(EXIT_SUCCESS);
    } else {
        less_plain_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_bool_, plain_bool_));
    }
}

TEST_F(AbyTest, equal_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        equal_test(1);
        exit(EXIT_SUCCESS);
    } else {
        equal_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
    }
}

TEST_F(AbyTest, less_than_zero_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        less_than_zero_test(1);
        exit(EXIT_SUCCESS);
    } else {
        less_than_zero_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_bool_, plain_bool_));
    }
}

TEST_F(AbyTest, sum_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        sum_test(1);
        exit(EXIT_SUCCESS);
    } else {
        sum_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_plain_matrix(reveal_plain_, plain_, 0.001));
    }
}

TEST_F(AbyTest, plain_shuffle_test) {
    std::vector<std::thread> threads;
    PlainMatrix<std::int64_t> plain;
    PlainMatrix<std::int64_t> share0;
    PlainMatrix<std::int64_t> share1;
    PlainMatrix<std::int64_t> shuffled_plain;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { plain_shuffle_test(i, plain, share0, share1); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    shuffled_plain = share0 + share1;

    // row = 10, col = 3;
    // mirror probability equal
    EXPECT_NE(shuffled_plain, plain);
    std::sort(plain.data(), plain.data() + 30);
    std::sort(shuffled_plain.data(), shuffled_plain.data() + 30);
    EXPECT_EQ(shuffled_plain, plain);
}

TEST_F(AbyTest, share_shuffle_test) {
    std::vector<std::thread> threads;
    PlainMatrix<std::int64_t> in0;
    PlainMatrix<std::int64_t> out0;
    PlainMatrix<std::int64_t> in1;
    PlainMatrix<std::int64_t> out1;
    PlainMatrix<std::int64_t> plain;
    PlainMatrix<std::int64_t> shuffled_plain;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { share_shuffle_test(i, in0, in1, out0, out1); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    shuffled_plain = out0 + out1;
    plain = in0 + in1;

    // row = 10, col = 3;
    // mirror probability equal
    EXPECT_NE(shuffled_plain, plain);
    std::sort(plain.data(), plain.data() + 30);
    std::sort(shuffled_plain.data(), shuffled_plain.data() + 30);
    EXPECT_EQ(shuffled_plain, plain);
}

}  // namespace duet
}  // namespace petace
