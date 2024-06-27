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

#include <cmath>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <thread>
#include <utility>

#include "gtest/gtest.h"
#include "print_utils.h"
#include "test_utils.h"

#include "network/net_factory.h"
#include "network/net_socket.h"
#include "network/network.h"
#include "solo/prng.h"
#include "verse/verse_factory.h"

#include "duet/util/common.h"
#include "duet/util/consts.h"
#include "duet/util/io.h"
#include "duet/util/matrix.h"

namespace petace {
namespace duet {

enum class AirthOP { Add, Sub, Mul, Div };

enum class BooleanhOP { And, Xor, Or, Not };

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

    void add_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        PublicMatrix<double> reveal_public(4, 3);
        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);

        aby_test->add(cipher_a, cipher_b, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_public);
        if (party_id == 0) {
            private_.matrix() = private_a.matrix() + private_b.matrix();
            reveal_private_.matrix() = reveal_public.matrix();
        }

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

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->sub(cipher_a, cipher_b, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);
        if (party_id == 0) {
            private_.matrix() = private_a.matrix() - private_b.matrix();
        }

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

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);

        aby_test->elementwise_mul(net, cipher_a, cipher_b, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);
        if (party_id == 0) {
            private_.matrix() = private_a.matrix().cwiseProduct(private_b.matrix());
        }

        return;
    }

    void public_mul_share_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicMatrix<double> public_a(4, 3);
        PrivateMatrix<double> private_b(4, 3, 0);

        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_plain_matrix(public_a.matrix());
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);
        aby_test->elementwise_mul(public_a, cipher_b, cipher_c);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);

        private_.resize(4, 3);
        if (party_id == 0) {
            private_.matrix() = public_a.matrix().cwiseProduct(private_b.matrix());
        }

        return;
    }

    void bool_public_share_test(std::size_t party_id, BooleanhOP op) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);
        PublicMatrixBool public_a(4, 3);
        PrivateMatrixBool private_b(4, 3, 0);

        BoolMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);

        get_rand_private_bool_matrix(private_b, 0);
        get_rand_public_bool_matrix(public_a);

        duet_test->share(net, private_b, cipher_b);
        duet_test->reveal(net, cipher_b, private_b);

        if (op == BooleanhOP::And) {
            duet_test->elementwise_bool_and(public_a, cipher_b, cipher_c);
        } else if (op == BooleanhOP::Xor) {
            duet_test->elementwise_bool_xor(public_a, cipher_b, cipher_c);
        } else if (op == BooleanhOP::Or) {
            duet_test->elementwise_bool_or(public_a, cipher_b, cipher_c);
        }

        duet_test->reveal(net, cipher_c, reveal_private_bool_);

        private_bool_.resize(4, 3);
        if (party_id == 0) {
            if (op == BooleanhOP::And) {
                for (size_t i = 0; i < public_a.size(); ++i) {
                    private_bool_(i) = public_a(i) & private_b(i);
                }
            } else if (op == BooleanhOP::Xor) {
                for (size_t i = 0; i < public_a.size(); ++i) {
                    private_bool_(i) = public_a(i) ^ private_b(i);
                }
            } else if (op == BooleanhOP::Or) {
                for (size_t i = 0; i < public_a.size(); ++i) {
                    private_bool_(i) = public_a(i) | private_b(i);
                }
            }
        }
    }

    void share_div_public_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicMatrix<double> public_a(4, 3);
        PrivateMatrix<double> private_b(4, 3, 0);

        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_plain_matrix(public_a.matrix());
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);
        aby_test->elementwise_div(cipher_b, public_a, cipher_c);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);

        private_.resize(4, 3);
        if (party_id == 0) {
            private_.matrix() = private_b.matrix().cwiseQuotient(public_a.matrix());
        }

        return;
    }

    void public_matrix_test(std::size_t party_id, AirthOP op) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicMatrix<double> public_a(4, 3);
        PrivateMatrix<double> private_b(4, 3, 0);
        PrivateMatrix<double> aim_ret(4, 3, 0);
        PrivateMatrix<double> real_ret(4, 3, 0);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);

        public_a.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699, 0.390775,
                -0.0008, 0.61128, -0.464325;

        get_rand_private_matrix(private_b, party_id);
        aby_test->share(net, private_b, cipher_b);

        if (op == AirthOP::Add) {
            aby_test->add(cipher_b, public_a, cipher_c);
        } else if (op == AirthOP::Sub) {
            aby_test->sub(cipher_b, public_a, cipher_c);
        } else if (op == AirthOP::Mul) {
            aby_test->elementwise_mul(public_a, cipher_b, cipher_c);
        } else if (op == AirthOP::Div) {
            aby_test->elementwise_div(cipher_b, public_a, cipher_c);
        }

        aby_test->reveal(net, cipher_c, real_ret);

        if (party_id == 0) {
            if (op == AirthOP::Add) {
                aim_ret.matrix() = private_b.matrix() + public_a.matrix();
            } else if (op == AirthOP::Sub) {
                aim_ret.matrix() = private_b.matrix() - public_a.matrix();
            } else if (op == AirthOP::Mul) {
                aim_ret.matrix() = private_b.matrix().cwiseProduct(public_a.matrix());
            } else if (op == AirthOP::Div) {
                for (std::size_t i = 0; i < private_b.size(); i++) {
                    aim_ret(i) = private_b(i) / public_a(i);
                }
            }
            ASSERT_TRUE(is_equal_private_matrix(real_ret, aim_ret, 0.1));
        }
        return;
    }

    void public_scalar_test(std::size_t party_id, AirthOP op) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicDouble public_a = 0.25;
        PrivateMatrix<double> private_b(4, 3, 0);
        PrivateMatrix<double> aim_ret(4, 3, 0);
        PrivateMatrix<double> real_ret(4, 3, 0);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);
        if (op == AirthOP::Add) {
            aby_test->add(cipher_b, public_a, cipher_c);
        } else if (op == AirthOP::Sub) {
            aby_test->sub(cipher_b, public_a, cipher_c);
        } else if (op == AirthOP::Mul) {
            aby_test->elementwise_mul(public_a, cipher_b, cipher_c);
        } else if (op == AirthOP::Div) {
            aby_test->elementwise_div(cipher_b, public_a, cipher_c);
        }

        aby_test->reveal(net, cipher_c, real_ret);

        if (party_id == 0) {
            if (op == AirthOP::Add) {
                aim_ret.matrix() = private_b.matrix().array() + public_a;
            } else if (op == AirthOP::Sub) {
                aim_ret.matrix() = private_b.matrix().array() - public_a;
            } else if (op == AirthOP::Mul) {
                aim_ret.matrix() = private_b.matrix().array() * public_a;
            } else if (op == AirthOP::Div) {
                aim_ret.matrix() = private_b.matrix().array() / public_a;
            }
            ASSERT_TRUE(is_equal_private_matrix(real_ret, aim_ret, 0.001));
        }

        return;
    }

    void mat_mul_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(3, 5, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(3, 5);
        ArithMatrix cipher_c(4, 5);
        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);

        aby_test->mat_mul(net, cipher_a, cipher_b, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);
        if (party_id == 0) {
            private_.matrix() = private_a.matrix() * private_b.matrix();
        }
        return;
    }

    void public_mat_mul_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicMatrix<double> public_a(4, 3);
        PrivateMatrix<double> private_b(3, 5, 0);
        public_a.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699, 0.390775,
                -0.0008, 0.61128, -0.464325;

        ArithMatrix cipher_b(3, 5);
        ArithMatrix cipher_c(4, 5);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);

        aby_test->mat_mul(public_a, cipher_b, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);
        private_.resize(4, 5);
        if (party_id == 0) {
            private_.matrix() = public_a.matrix() * private_b.matrix();
        }
        return;
    }

    void mat_mul_public_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicMatrix<double> public_a(4, 3);
        PrivateMatrix<double> private_b(5, 4, 0);
        public_a.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699, 0.390775,
                -0.0008, 0.61128, -0.464325;

        ArithMatrix cipher_b(5, 4);
        ArithMatrix cipher_c(5, 3);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);

        aby_test->mat_mul(cipher_b, public_a, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);
        private_.resize(4, 5);
        if (party_id == 0) {
            private_.matrix() = private_b.matrix() * public_a.matrix();
        }
        return;
    }

    void bool_ss_test(std::size_t party_id, BooleanhOP op) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrixBool private_a(4, 3, 0);
        PrivateMatrixBool private_b(4, 3, 0);
        PrivateMatrixBool aim_ret(4, 3, 0);
        PrivateMatrixBool real_ret(4, 3, 0);

        BoolMatrix cipher_a(4, 3);
        BoolMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);

        get_rand_plain_bool_matrix(cipher_a.matrix());
        get_rand_plain_bool_matrix(cipher_b.matrix());

        duet_test->reveal(net, cipher_a, private_a);
        duet_test->reveal(net, cipher_b, private_b);

        if (op == BooleanhOP::And) {
            duet_test->elementwise_bool_and(net, cipher_a, cipher_b, cipher_c);
        } else if (op == BooleanhOP::Xor) {
            duet_test->elementwise_bool_xor(cipher_a, cipher_b, cipher_c);
        } else if (op == BooleanhOP::Or) {
            duet_test->elementwise_bool_or(net, cipher_a, cipher_b, cipher_c);
        } else if (op == BooleanhOP::Not) {
            duet_test->elementwise_bool_not(cipher_a, cipher_c);
        }

        duet_test->reveal(net, cipher_c, real_ret);

        if (party_id == 0) {
            if (op == BooleanhOP::And) {
                for (size_t i = 0; i < private_a.size(); ++i) {
                    aim_ret(i) = private_a(i) & private_b(i);
                }
            } else if (op == BooleanhOP::Xor) {
                for (size_t i = 0; i < private_a.size(); ++i) {
                    aim_ret(i) = private_a(i) ^ private_b(i);
                }
            } else if (op == BooleanhOP::Or) {
                for (size_t i = 0; i < private_a.size(); ++i) {
                    aim_ret(i) = private_a(i) | private_b(i);
                }
            } else if (op == BooleanhOP::Not) {
                for (size_t i = 0; i < private_a.size(); ++i) {
                    aim_ret(i) = ~private_a(i);
                }
            }
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, real_ret));
        }

        return;
    }

    void div_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        private_a.matrix() << -2.13039, -0.875821, -0.436358, -1.1604, -0.710526, 1.44121, -0.233295, -1.7699, 0.390775,
                -0.108926, 0.61128, -0.464325;
        private_b.matrix() << -0.507397, -1.46439, 0.448924, 0.43636, -0.0734787, 0.447709, 2.55159, -0.0437469,
                0.245489, -3.19096, -0.493532, -1.25605;
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);

        aby_test->elementwise_div(net, cipher_a, cipher_b, cipher_c);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);

        if (party_id == 0) {
            private_.resize(4, 3);
            private_.set_party_id(0);
            for (std::size_t i = 0; i < private_a.size(); i++) {
                private_(i) = private_a(i) / private_b(i);
            }
        }

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

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        PrivateMatrix<double> private_c(4, 3, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);
        ArithMatrix cipher_d(4, 3);

        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);
        get_rand_private_matrix(private_c, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);

        aby_test->greater(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(4, 3);

        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) > private_b(i)) ? private_c(i) : 0;
            }
        }

        return;
    }

    void multiplexer_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        std::size_t row = 4;
        std::size_t col = 3;

        PrivateMatrix<double> private_a(row, col, 0);
        PrivateMatrix<double> private_b(row, col, 0);
        PrivateMatrix<double> private_c(row, col, 0);
        PrivateMatrix<double> private_d(row, col, 0);

        ArithMatrix cipher_a(row, col);
        ArithMatrix cipher_b(row, col);
        BoolMatrix cipher_c(row, col);
        ArithMatrix cipher_d(row, col);
        ArithMatrix cipher_e(row, col);

        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);
        get_rand_private_matrix(private_c, party_id);
        get_rand_private_matrix(private_d, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);
        aby_test->share(net, private_d, cipher_e);

        aby_test->greater(net, cipher_a, cipher_b, cipher_c);

        std::cout << net->get_bytes_sent() << " " << net->get_bytes_received() << std::endl;
        aby_test->multiplexer(net, cipher_c, cipher_e, cipher_d, cipher_b);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(row, col);

        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) > private_b(i)) ? private_c(i) : private_d(i);
            }
        }
        std::cout << net->get_bytes_sent() << " " << net->get_bytes_received() << std::endl;

        return;
    }

    void greater_public_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(4, 3, 0);
        PublicMatrix<double> public_b(4, 3);
        ArithMatrix cipher_a(4, 3);
        BoolMatrix cipher_c(4, 3);

        get_rand_private_matrix(private_a, party_id);
        get_rand_plain_matrix(public_b.matrix());

        aby_test->share(net, private_a, cipher_a);
        aby_test->greater(net, cipher_a, public_b, cipher_c);
        reveal_private_bool_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_bool_);

        private_bool_.resize(4, 3);
        if (party_id == 0) {
            for (std::size_t i = 0; i < static_cast<std::size_t>(public_b.size()); i++) {
                private_bool_(i) = (private_a(i) > public_b(i)) ? 1 : 0;
            }
        }

        return;
    }

    void less_public_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(4, 3, 0);
        PublicMatrix<double> public_b(4, 3);
        ArithMatrix cipher_a(4, 3);
        BoolMatrix cipher_c(4, 3);

        get_rand_private_matrix(private_a, party_id);
        get_rand_plain_matrix(public_b.matrix());

        aby_test->share(net, private_a, cipher_a);
        aby_test->less(net, cipher_a, public_b, cipher_c);
        reveal_private_bool_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_bool_);

        private_bool_.resize(4, 3);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_a.size(); i++) {
                private_bool_(i) = (private_a(i) < public_b(i)) ? 1 : 0;
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

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        PrivateMatrix<double> private_c(4, 3, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);
        ArithMatrix cipher_d(4, 3);

        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);
        get_rand_private_matrix(private_c, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);

        aby_test->less(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(4, 3);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) < private_b(i)) ? private_c(i) : 0;
            }
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

        std::size_t row = 4;
        std::size_t col = 3;

        PrivateMatrix<double> private_a(row, col, 0);
        PrivateMatrix<double> private_b(row, col, 0);
        PrivateMatrix<double> private_c(row, col, 0);
        ArithMatrix cipher_a(row, col);
        ArithMatrix cipher_b(row, col);
        BoolMatrix cipher_c(row, col);
        ArithMatrix cipher_d(row, col);

        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);
        get_rand_private_matrix(private_c, party_id);

        private_b(0) = private_a(0);
        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);

        aby_test->equal(net, cipher_a, cipher_b, cipher_c);
        std::cout << net->get_bytes_sent() << " " << net->get_bytes_received() << std::endl;
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(row, col);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) == private_b(i)) ? private_c(i) : 0;
            }
        }

        return;
    }

    void not_equal_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        std::size_t row = 4;
        std::size_t col = 3;

        PrivateMatrix<double> private_a(row, col, 0);
        PrivateMatrix<double> private_b(row, col, 0);
        PrivateMatrix<double> private_c(row, col, 0);
        ArithMatrix cipher_a(row, col);
        ArithMatrix cipher_b(row, col);
        BoolMatrix cipher_c(row, col);
        ArithMatrix cipher_d(row, col);

        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);
        get_rand_private_matrix(private_c, party_id);

        private_b(0) = private_a(0);
        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);

        aby_test->not_equal(net, cipher_a, cipher_b, cipher_c);
        std::cout << net->get_bytes_sent() << " " << net->get_bytes_received() << std::endl;
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(row, col);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) != private_b(i)) ? private_c(i) : 0;
            }
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

        PrivateMatrixBool private_x(1, num_element, 0);
        ArithMatrix share_x(1, num_element);
        BoolMatrix res(1, num_element);

        get_rand_private_matrix(private_x, party_id);

        private_bool_.resize(1, num_element);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_x.size(); i++) {
                private_bool_(i) = (private_x(i) < 0) ? 1 : 0;
            }
        }

        aby_test->share(net, private_x, share_x);
        aby_test->less_than_zero(net, share_x, res);
        reveal_private_bool_.set_party_id(0);
        aby_test->reveal(net, res, reveal_private_bool_);

        return;
    }

    void paillier_h2a_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        std::size_t test_size = 100;
        PaillierMatrix cipher_he_1(1, test_size, 0);
        ArithMatrix shared_cipher(1, test_size);
        private_.set_party_id(1);
        private_.resize(1, test_size);
        if (party_id == 0) {
            PrivateMatrix<double> plain_he_1(1, test_size, 0);
            for (std::size_t i = 0; i < test_size; ++i) {
                plain_he_1(i) = static_cast<double>(i);
                private_(i) = static_cast<double>(i);
            }
            aby_test->encrypt(plain_he_1, cipher_he_1);
            send_cipher(net, cipher_he_1, kPaillierKeySize);
        } else {
            recv_cipher(net, aby_test->get_pk_other(), cipher_he_1, kPaillierKeySize);
        }
        aby_test->h2a(net, cipher_he_1, shared_cipher);
        reveal_private_.set_party_id(0);
        reveal_private_.resize(1, test_size);
        aby_test->reveal(net, shared_cipher, reveal_private_);

        return;
    }

    void thread_h2a_test(std::size_t party_id, PrivateMatrix<double>& plain, PaillierMatrix& cipher_he_0,
            PaillierMatrix& cipher_he_1, ArithMatrix& shared_cipher_0, ArithMatrix& shared_cipher_1,
            PrivateMatrix<double>& private_revealed) {
        std::shared_ptr<petace::network::Network> net;
        std::size_t test_size = 100;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        if (party_id == 0) {
            PrivateMatrix<double> plain_he_1(1, test_size, 0);
            for (std::size_t i = 0; i < test_size; ++i) {
                plain_he_1(i) = plain(i);
            }
            aby_test->encrypt(plain_he_1, cipher_he_0);
            send_cipher(net, cipher_he_0, kPaillierKeySize);
            aby_test->h2a(net, cipher_he_0, shared_cipher_0);
            aby_test->reveal(net, shared_cipher_0, private_revealed);
        } else {
            recv_cipher(net, aby_test->get_pk_other(), cipher_he_1, kPaillierKeySize);
            aby_test->h2a(net, cipher_he_1, shared_cipher_1);
            aby_test->reveal(net, shared_cipher_1, private_revealed);
        }
        return;
    }

    void paillier_a2h_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        std::size_t test_size = 100;
        auto aby_test = std::make_shared<Duet>(net, party_id);
        PrivateMatrix<double> plain_1_a2h(1, test_size, 0);
        ArithMatrix shared_1_a2h(1, test_size);
        PaillierMatrix cipher_1_a2h(1, test_size, 0);
        private_.set_party_id(0);
        private_.resize(1, test_size);
        if (party_id == 0) {
            for (std::size_t i = 0; i < test_size; ++i) {
                plain_1_a2h(i) = static_cast<double>(i);
                private_(i) = static_cast<double>(i);
            }
        }

        aby_test->share(net, plain_1_a2h, shared_1_a2h);
        aby_test->a2h(net, shared_1_a2h, cipher_1_a2h);
        if (party_id == 1) {
            send_cipher(net, cipher_1_a2h, kPaillierKeySize);
        } else {
            recv_cipher(net, aby_test->get_pk(), cipher_1_a2h, kPaillierKeySize);
        }
        reveal_private_.set_party_id(0);
        reveal_private_.resize(1, test_size);
        if (party_id == 0) {
            aby_test->decrypt(cipher_1_a2h, reveal_private_);
        }

        return;
    }

    void paillier_add_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        PrivateMatrix<double> plain_1(1, 5, 0);
        PrivateMatrix<double> plain_2(1, 5, 0);
        PrivateMatrix<double> plain_3(1, 5, 0);
        PaillierMatrix cipher_1(1, 5);
        PaillierMatrix cipher_2(1, 5);
        PaillierMatrix cipher_3(1, 5);
        PaillierMatrix cipher_sum(1, 5);
        PaillierMatrix cipher_sum2(1, 5);
        private_.set_party_id(0);
        private_.resize(1, 5);
        reveal_private_.set_party_id(0);
        reveal_private_.resize(1, 5);
        if (aby_test->party() == 0) {
            for (std::size_t i = 0; i < 5; ++i) {
                plain_1(i) = static_cast<double>(i);
                plain_2(i) = static_cast<double>(i * 3);
                plain_3(i) = static_cast<double>(i * 5);
                private_(i) = static_cast<double>(i * 9);
            }
            aby_test->encrypt(plain_1, cipher_1);
            aby_test->encrypt(plain_2, cipher_2);
            aby_test->add(cipher_1, cipher_2, cipher_sum);
            aby_test->add(plain_3, cipher_sum, cipher_sum2);
            aby_test->decrypt(cipher_sum2, reveal_private_);
        }
        return;
    }

    void paillier_matrix_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        petace::solo::ahepaillier::Encoder encoder;
        PrivateMatrix<double> plain_1(1, 10, 0);
        PrivateMatrix<double> plain_2(1, 1, 0);
        PaillierMatrix cipher_1(1, 10);
        PaillierMatrix cipher_2(1, 1, 0);
        private_.set_party_id(0);
        private_.resize(1, 10);
        reveal_private_.set_party_id(0);
        reveal_private_.resize(1, 10);
        if (aby_test->party() == 0) {
            for (std::size_t i = 0; i < 10; ++i) {
                plain_1(i) = static_cast<double>(i);
                private_(i) = static_cast<double>(i);
            }
            aby_test->encrypt(plain_1, cipher_1);
            for (std::size_t i = 0; i < 5; ++i) {
                cipher_2(0) = petace::solo::ahepaillier::Ciphertext(cipher_1(i));
                aby_test->decrypt(cipher_2, plain_2);
                reveal_private_(i) = plain_2(0);
            }
            for (std::size_t i = 5; i < 10; ++i) {
                cipher_2(0) = petace::solo::ahepaillier::Ciphertext(cipher_1(0, i));
                aby_test->decrypt(cipher_2, plain_2);
                reveal_private_(i) = plain_2(0);
            }
        }
        return;
    }

    void paillier_add_int_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        PrivateMatrix<std::int64_t> plain_1(1, 5, 0);
        PrivateMatrix<std::int64_t> plain_2(1, 5, 0);
        PaillierMatrix cipher_1(1, 5);
        PaillierMatrix cipher_2(1, 5);
        PaillierMatrix cipher_sum(1, 5);
        private_bool_.set_party_id(0);
        private_bool_.resize(1, 5);
        reveal_private_bool_.set_party_id(0);
        reveal_private_bool_.resize(1, 5);
        if (aby_test->party() == 0) {
            for (std::size_t i = 0; i < 5; ++i) {
                plain_1(i) = static_cast<std::int64_t>(i);
                plain_2(i) = static_cast<std::int64_t>(i * 4);
                private_bool_(i) = static_cast<std::int64_t>(i * 5);
            }
            aby_test->encrypt(plain_1, cipher_1);
            aby_test->add(plain_2, cipher_1, cipher_sum);
            aby_test->decrypt(cipher_sum, reveal_private_bool_);
        }
        return;
    }

    void paillier_mul_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        PrivateMatrix<std::int64_t> plain_1(1, 5, 0);
        PrivateMatrix<std::int64_t> plain_2(1, 5, 0);
        PaillierMatrix cipher_1(1, 5);
        PaillierMatrix cipher_mul(1, 5);
        private_bool_.set_party_id(0);
        private_bool_.resize(1, 5);
        reveal_private_bool_.set_party_id(0);
        reveal_private_bool_.resize(1, 5);
        if (party_id == 0) {
            for (std::size_t i = 0; i < 5; ++i) {
                plain_1(i) = static_cast<std::int64_t>(i);
                plain_2(i) = static_cast<std::int64_t>(i * 3);
                private_bool_(i) = static_cast<std::int64_t>(i * i * 3);
            }
            aby_test->encrypt(plain_1, cipher_1);
            aby_test->mul(plain_2, cipher_1, cipher_mul);
            aby_test->decrypt(cipher_mul, reveal_private_bool_);
        }
        return;
    }

    void paillier_encryption_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);
        PrivateMatrix<std::int64_t> plain_1(1, 5, 0);
        PaillierMatrix cipher_1(1, 5);
        if (party_id == 0) {
            private_bool_.set_party_id(0);
            private_bool_.resize(1, 5);
            reveal_private_bool_.set_party_id(0);
            reveal_private_bool_.resize(1, 5);
            for (std::size_t i = 0; i < 5; ++i) {
                plain_1(i) = static_cast<std::int64_t>(i);
                private_bool_(i) = static_cast<std::int64_t>(i);
            }
            aby_test->encrypt(plain_1, cipher_1);
            aby_test->decrypt(cipher_1, reveal_private_bool_);
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

        PrivateMatrix<double> private_a(4, 3, 0);
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);

        get_rand_private_matrix(private_a, party_id);

        aby_test->share(net, private_a, cipher_a);

        aby_test->sum(cipher_a, cipher_b);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);
        if (party_id == 0) {
            private_.matrix() = private_a.matrix().colwise().sum();
        }
        return;
    }

    void argmax_and_max_test(std::size_t party_id, Matrix<std::int64_t>& aim_index, Matrix<std::int64_t>& real_index,
            Matrix<double>& aim_value, Matrix<double>& real_value) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(5, 3, 0);
        PrivateMatrix<double> private_b(0);
        PrivateMatrix<double> private_c(0);
        ArithMatrix cipher_a(5, 3);
        ArithMatrix cipher_b(5, 3);
        ArithMatrix cipher_c(5, 3);

        get_rand_private_matrix(private_a, party_id);

        aby_test->share(net, private_a, cipher_a);

        aby_test->argmax_and_max(net, cipher_a, cipher_b, cipher_c);

        aby_test->reveal(net, cipher_b, private_b);
        aby_test->reveal(net, cipher_c, private_c);
        if (party_id == 0) {
            std::vector<Matrix<double>::Index> max_indices(private_a.cols());
            for (std::size_t i = 0; i < private_a.cols(); ++i) {
                private_a.matrix().col(i).maxCoeff(&max_indices[i]);
            }
            aim_value = private_a.matrix().colwise().maxCoeff();
            real_value = private_c.matrix();
            aim_index.resize(3, 1);
            real_index.resize(3, 1);
            for (std::size_t i = 0; i < 3; ++i) {
                aim_index(i) = static_cast<std::int64_t>(max_indices[i]);
                real_index(i) = static_cast<std::int64_t>(private_b(i));
            }
        }
        return;
    }

    void argmin_and_min_test(std::size_t party_id, Matrix<std::int64_t>& aim_index, Matrix<std::int64_t>& real_index,
            Matrix<double>& aim_value, Matrix<double>& real_value) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PrivateMatrix<double> private_a(5, 3, 0);
        PrivateMatrix<double> private_b(0);
        PrivateMatrix<double> private_c(0);
        ArithMatrix cipher_a(5, 3);
        ArithMatrix cipher_b(5, 3);
        ArithMatrix cipher_c(5, 3);

        get_rand_private_matrix(private_a, party_id);

        aby_test->share(net, private_a, cipher_a);

        aby_test->argmin_and_min(net, cipher_a, cipher_b, cipher_c);

        aby_test->reveal(net, cipher_b, private_b);
        aby_test->reveal(net, cipher_c, private_c);
        if (party_id == 0) {
            std::vector<Matrix<double>::Index> min_indices(private_a.cols());
            for (std::size_t i = 0; i < private_a.cols(); ++i) {
                private_a.matrix().col(i).minCoeff(&min_indices[i]);
            }
            aim_value = private_a.matrix().colwise().minCoeff();
            real_value = private_c.matrix();
            aim_index.resize(3, 1);
            real_index.resize(3, 1);
            for (std::size_t i = 0; i < 3; ++i) {
                aim_index(i) = static_cast<std::int64_t>(min_indices[i]);
                real_index(i) = static_cast<std::int64_t>(private_b(i));
            }
        }
        return;
    }

    void private_shuffle_test(std::size_t party_id, Matrix<std::int64_t>& plain, Matrix<std::int64_t>& share0,
            Matrix<std::int64_t>& share1) {
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

        PrivateMatrix<double> private_a(rows, cols, 0);
        ArithMatrix cipher_a(rows, cols);
        get_rand_private_matrix(private_a, party_id);
        if (party_id == 0) {
            plain.resize(rows, cols);
            for (std::size_t i = 0; i < matrix_size; ++i) {
                plain(i) = double_to_fixed(private_a(i));
            }
        }
        aby_test->shuffle(net, private_a, cipher_a);
        if (party_id == 0) {
            share0 = cipher_a.shares();
        } else {
            share1 = cipher_a.shares();
        }
        return;
    }

    void share_shuffle_test(std::size_t party_id, Matrix<std::int64_t>& in0, Matrix<std::int64_t>& in1,
            Matrix<std::int64_t>& out0, Matrix<std::int64_t>& out1) {
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
            cipher_a(i) = random();
        }

        aby_test->shuffle(net, cipher_a, cipher_b);
        if (party_id == 0) {
            in0 = cipher_a.shares();
            out0 = cipher_b.shares();
        } else {
            in1 = cipher_a.shares();
            out1 = cipher_b.shares();
        }
    }

    void he_share_shuffle_test(std::size_t party_id, Matrix<std::int64_t>& in0, Matrix<std::int64_t>& in1,
            Matrix<std::int64_t>& out0, Matrix<std::int64_t>& out1, std::vector<std::size_t>& permutation) {
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
            cipher_a(i) = random();
        }

        PrivatePermutation perm(10);
        perm.set_party_id(0);
        if (party_id == 0) {
            permutation = perm.data();
        }

        aby_test->shuffle(net, perm, cipher_a, cipher_b);
        if (party_id == 0) {
            in0 = cipher_a.shares();
            out0 = cipher_b.shares();
        } else {
            in1 = cipher_a.shares();
            out1 = cipher_b.shares();
        }
    }

    void quick_sort_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        std::size_t matrix_size = 20;

        PrivateMatrix<double> private_a(matrix_size, 1, 0);
        PrivateMatrix<double> aim_ret(matrix_size, 1, 0);
        PrivateMatrix<double> real_ret(0);
        ArithMatrix cipher_a(matrix_size, 1);

        auto duet_test = std::make_shared<Duet>(net, party_id);

        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator(seed);

        std::vector<double> arr(matrix_size, 0);
        for (std::size_t i = 0; i < arr.size(); ++i) {
            arr[i] = static_cast<double>(i);
        }
        std::shuffle(arr.begin(), arr.end(), generator);

        if (party_id == 0) {
            for (std::size_t i = 0; i < arr.size(); ++i) {
                private_a(i) = arr[i];
                aim_ret(i) = static_cast<double>(i);
            }
        }
        duet_test->share(net, private_a, cipher_a);
        duet_test->quick_sort(net, cipher_a);
        duet_test->reveal(net, cipher_a, real_ret);
        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, real_ret, 0.001));
        }
    }

    void quick_sort_by_column_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        std::size_t rows = 4;
        std::size_t cols = 2;

        PrivateMatrix<double> private_a(rows, cols, 0);
        PrivateMatrix<double> aim_ret0(rows, cols, 0);
        PrivateMatrix<double> aim_ret1(rows, cols, 0);
        PrivateMatrix<double> real_ret0(0);
        PrivateMatrix<double> real_ret1(0);
        ArithMatrix cipher_a;
        ArithMatrix cipher_b;
        ArithMatrix cipher_c;

        auto duet_test = std::make_shared<Duet>(net, party_id);

        std::vector<double> arr = {1, 3, 4, 1, 2, 2, 3, 4};
        std::vector<double> aim_arr0 = {1, 3, 2, 2, 3, 4, 4, 1};  // 0
        std::vector<double> aim_arr1 = {4, 1, 2, 2, 1, 3, 3, 4};  // 1

        if (party_id == 0) {
            for (std::size_t i = 0; i < arr.size(); ++i) {
                private_a(i) = arr[i];
                aim_ret0(i) = aim_arr0[i];
                aim_ret1(i) = aim_arr1[i];
            }
        }
        duet_test->share(net, private_a, cipher_a);
        duet_test->quick_sort(net, 0, cipher_a, cipher_b);
        duet_test->quick_sort(net, 1, cipher_a, cipher_c);
        duet_test->reveal(net, cipher_b, real_ret0);
        duet_test->reveal(net, cipher_c, real_ret1);
        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret0, real_ret0, 0.001));
            ASSERT_TRUE(is_equal_private_matrix(aim_ret1, real_ret1, 0.001));
        }
    }

    void split_by_condition_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);
        BoolMatrix cond(4, 1);
        PrivateMatrix<double> aim_ret_0(2, 3, 0);
        PrivateMatrix<double> aim_ret_1(2, 3, 0);
        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> real_ret_0(0);
        PrivateMatrix<double> real_ret_1(0);
        get_rand_private_matrix(private_a, party_id);
        std::size_t index_0 = 0;
        std::size_t index_1 = 0;

        for (std::size_t i = 0; i < 4; i++) {
            if (party_id == 0) {
                if (i % 2 == 0) {
                    cond(i) = 0;
                    aim_ret_0.matrix().row(index_0) = private_a.matrix().row(i);
                    index_0++;
                } else {
                    cond(i) = 1;
                    aim_ret_1.matrix().row(index_1) = private_a.matrix().row(i);
                    index_1++;
                }
            } else {
                cond(i) = 0;
            }
        }
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b;
        ArithMatrix cipher_c;

        duet_test->share(net, private_a, cipher_a);
        duet_test->split_by_condition(net, cond, cipher_a, cipher_b, cipher_c);
        duet_test->reveal(net, cipher_b, real_ret_0);
        duet_test->reveal(net, cipher_c, real_ret_1);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret_0, real_ret_0, 0.001));
            ASSERT_TRUE(is_equal_private_matrix(aim_ret_1, real_ret_1, 0.001));
        }
    }

    void sigmoid_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);
        PrivateMatrix<double> private_a(4, 5, 0);
        PrivateMatrix<double> aim_ret(4, 5, 0);
        PrivateMatrix<double> real_ret(0);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_a.size(); ++i) {
                private_a(i) = -1.0 + static_cast<double>(i) * 0.1;
                aim_ret(i) = 1.0 / (1.0 + std::exp(-private_a(i)));
            }
        }

        ArithMatrix cipher_a(4, 5);
        ArithMatrix cipher_b(4, 5);
        duet_test->share(net, private_a, cipher_a);
        duet_test->sigmoid(net, cipher_a, cipher_b);
        duet_test->reveal(net, cipher_b, real_ret);
        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(aim_ret, real_ret, 0.001));
        }
    }

    void groupby_sum_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);

        const std::size_t row = 6;
        const std::size_t col = 2;

        PrivateMatrix<double> plain0(0), plain1(0);
        ArithMatrix cipher0, cipher1;
        PrivateMatrix<double> one_hot_matrix(0);
        ArithMatrix one_hot_cipher_matrix(6, 2);

        plain0.resize(row, col);
        plain0.matrix() << -15.812, -14.7387, 7.8120, -9.7387, -2.8120, 6.7387, 1.8120, 1.7387, 5.8120, 0.7387, -7.8120,
                -9.7387;

        one_hot_matrix.resize(row, col);
        one_hot_matrix.matrix() << 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1;

        plain1.resize(2, 2);
        plain1.matrix() << -14, -8.188, -32.4774, -31.7387;

        duet_test->share(net, plain0, cipher0);
        duet_test->share(net, one_hot_matrix, one_hot_cipher_matrix);

        duet_test->groupby_sum(net, cipher0, one_hot_cipher_matrix, cipher1);
        duet_test->reveal(net, cipher1, plain0);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(plain0, plain1, 0.001));
        }
    }

    void groupby_count_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);

        const std::size_t row = 6;
        const std::size_t col = 2;

        PrivateMatrix<double> plain0(0), plain1(0);
        ArithMatrix cipher0, cipher1;
        PrivateMatrix<double> one_hot_matrix(0);
        ArithMatrix one_hot_cipher_matrix;

        plain0.resize(row, col);
        plain0.matrix() << -15.812, -14.7387, 7.8120, -9.7387, -2.8120, 6.7387, 1.8120, 1.7387, 5.8120, 0.7387, -7.8120,
                -9.7387;

        one_hot_matrix.resize(row, col);
        one_hot_matrix.matrix() << 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1;

        plain1.resize(2, 2);
        plain1.matrix() << 4, 5, 4, 5;

        duet_test->share(net, plain0, cipher0);
        duet_test->share(net, one_hot_matrix, one_hot_cipher_matrix);

        duet_test->groupby_count(cipher0, one_hot_cipher_matrix, cipher1);
        duet_test->reveal(net, cipher1, plain0);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(plain0, plain1, 0.001));
        }
    }

    void groupby_max_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);

        const std::size_t row = 6;
        const std::size_t col = 2;

        PrivateMatrix<double> plain0(0), plain1(0);
        ArithMatrix cipher0, cipher1;
        PrivateMatrix<double> one_hot_matrix(0);
        ArithMatrix one_hot_cipher_matrix(6, 2);

        plain0.resize(row, col);
        plain0.matrix() << -15.812, -14.7387, 7.8120, -9.7387, -2.8120, 6.7387, 1.8120, 1.7387, 5.8120, 0.7387, -7.8120,
                -9.7387;

        one_hot_matrix.resize(row, col);
        one_hot_matrix.matrix() << 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1;

        plain1.resize(2, 2);
        plain1.matrix() << 7.8120, 7.8120, 1.73869, 1.7387;

        duet_test->share(net, plain0, cipher0);
        duet_test->share(net, one_hot_matrix, one_hot_cipher_matrix);

        duet_test->groupby_max(net, cipher0, one_hot_cipher_matrix, cipher1);
        duet_test->reveal(net, cipher1, plain0);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(plain0, plain1, 0.001));
        }
    }

    void groupby_min_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);

        const std::size_t row = 6;
        const std::size_t col = 2;

        PrivateMatrix<double> plain0(0), plain1(0);
        ArithMatrix cipher0, cipher1;
        PrivateMatrix<double> one_hot_matrix(0);
        ArithMatrix one_hot_cipher_matrix(6, 2);

        plain0.resize(row, col);
        plain0.matrix() << -15.812, -14.7387, 7.8120, -9.7387, -2.8120, 6.7387, 1.8120, 1.7387, 5.8120, 0.7387, -7.8120,
                -9.7387;

        one_hot_matrix.resize(row, col);
        one_hot_matrix.matrix() << 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1;

        plain1.resize(2, 2);
        plain1.matrix() << -15.812, -15.812, -14.7387, -14.7387;

        duet_test->share(net, plain0, cipher0);
        duet_test->share(net, one_hot_matrix, one_hot_cipher_matrix);

        duet_test->groupby_min(net, cipher0, one_hot_cipher_matrix, cipher1);
        duet_test->reveal(net, cipher1, plain0);

        if (party_id == 0) {
            ASSERT_TRUE(is_equal_private_matrix(plain0, plain1, 0.001));
        }
    }

public:
    PrivateMatrix<double> reveal_private_;
    PrivateMatrix<double> private_;
    PrivateMatrixBool reveal_private_bool_;
    PrivateMatrixBool private_bool_;
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, mat_mul_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        mat_mul_test(1);
        exit(EXIT_SUCCESS);
    } else {
        mat_mul_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.01));
    }
}

TEST_F(AbyTest, public_mat_mul_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        public_mat_mul_test(1);
        exit(EXIT_SUCCESS);
    } else {
        public_mat_mul_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.01));
    }
}

TEST_F(AbyTest, mat_mul_public_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        mat_mul_public_test(1);
        exit(EXIT_SUCCESS);
    } else {
        mat_mul_public_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.01));
    }
}

TEST_F(AbyTest, public_mul_share_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        public_mul_share_test(1);
        exit(EXIT_SUCCESS);
    } else {
        public_mul_share_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, share_div_public_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        share_div_public_test(1);
        exit(EXIT_SUCCESS);
    } else {
        share_div_public_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, bool_public_share_test) {
    pid_t pid;
    int status;

    for (int op = static_cast<int>(BooleanhOP::And); op <= static_cast<int>(BooleanhOP::Or); ++op) {
        if ((pid = fork()) < 0) {
            status = -1;
        } else if (pid == 0) {
            bool_public_share_test(1, static_cast<BooleanhOP>(op));
            exit(EXIT_SUCCESS);
        } else {
            bool_public_share_test(0, static_cast<BooleanhOP>(op));
            while (waitpid(pid, &status, 0) < 0) {
                if (errno != EINTR) {
                    status = -1;
                    break;
                }
            }
            ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
        }
    }
}

TEST_F(AbyTest, public_matrix_test) {
    for (int op = static_cast<int>(AirthOP::Add); op <= static_cast<int>(AirthOP::Div); ++op) {
        std::vector<std::thread> threads;

        for (std::size_t i = 0; i < 2; ++i) {
            threads.emplace_back([&, i]() { public_matrix_test(i, static_cast<AirthOP>(op)); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }
}

TEST_F(AbyTest, public_scalar_test) {
    for (int op = static_cast<int>(AirthOP::Add); op <= static_cast<int>(AirthOP::Div); ++op) {
        std::vector<std::thread> threads;

        for (std::size_t i = 0; i < 2; ++i) {
            threads.emplace_back([&, i]() { public_scalar_test(i, static_cast<AirthOP>(op)); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }
}

TEST_F(AbyTest, multiplexer_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        multiplexer_test(1);
        exit(EXIT_SUCCESS);
    } else {
        multiplexer_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, div_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        div_test(1);
        exit(EXIT_SUCCESS);
    } else {
        div_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.1));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, greater_public_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        greater_public_test(1);
        exit(EXIT_SUCCESS);
    } else {
        greater_public_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, less_public_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        less_public_test(1);
        exit(EXIT_SUCCESS);
    } else {
        less_public_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, not_equal_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        not_equal_test(1);
        exit(EXIT_SUCCESS);
    } else {
        not_equal_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, paillier_a2h_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_encryption_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_encryption_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
    }
}

TEST_F(AbyTest, paillier_encryption_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_encryption_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_encryption_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
    }
}

TEST_F(AbyTest, paillier_add_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_add_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_add_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, paillier_add_int_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_add_int_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_add_int_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
    }
}

TEST_F(AbyTest, paillier_mul_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_mul_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_mul_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
    }
}

TEST_F(AbyTest, paillier_h2a_test) {
    std::vector<std::thread> threads;
    std::size_t test_size = 100;
    PrivateMatrix<double> plain(1, test_size, 0);
    PaillierMatrix cipher_he_0, cipher_he_1;
    ArithMatrix shared_cipher_0(1, test_size);
    ArithMatrix shared_cipher_1(1, test_size);
    PrivateMatrix<double> private_revealed(1, test_size, 0);

    for (std::size_t i = 0; i < test_size; ++i) {
        plain(i) = static_cast<double>(i);
    }
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() {
            thread_h2a_test(i, plain, cipher_he_0, cipher_he_1, shared_cipher_0, shared_cipher_1, private_revealed);
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    for (std::size_t i = 0; i < test_size; ++i) {
        EXPECT_EQ(plain(i), private_revealed(i));
    }
}

TEST_F(AbyTest, paillier_matrix_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_matrix_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_matrix_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
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
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, argmax_and_max_test) {
    std::vector<std::thread> threads;
    Matrix<std::int64_t> aim_index;
    Matrix<std::int64_t> real_index;
    Matrix<double> aim_value;
    Matrix<double> real_value;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { argmax_and_max_test(i, aim_index, real_index, aim_value, real_value); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    EXPECT_EQ(is_equal_plain_matrix(aim_value, real_value, 0.001), true);
    EXPECT_EQ(is_equal_plain_matrix(aim_index, real_index), true);
}

TEST_F(AbyTest, argmin_and_min_test) {
    std::vector<std::thread> threads;
    Matrix<std::int64_t> aim_index;
    Matrix<std::int64_t> real_index;
    Matrix<double> aim_value;
    Matrix<double> real_value;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { argmin_and_min_test(i, aim_index, real_index, aim_value, real_value); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    EXPECT_EQ(is_equal_plain_matrix(aim_value, real_value, 0.001), true);
    EXPECT_EQ(is_equal_plain_matrix(aim_index, real_index), true);
}

TEST_F(AbyTest, private_shuffle_test) {
    std::vector<std::thread> threads;
    Matrix<std::int64_t> plain;
    Matrix<std::int64_t> share0;
    Matrix<std::int64_t> share1;
    Matrix<std::int64_t> shuffled_plain;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { private_shuffle_test(i, plain, share0, share1); });
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
    Matrix<std::int64_t> in0;
    Matrix<std::int64_t> out0;
    Matrix<std::int64_t> in1;
    Matrix<std::int64_t> out1;
    Matrix<std::int64_t> plain;
    Matrix<std::int64_t> shuffled_plain;

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

TEST_F(AbyTest, he_share_shuffle_test) {
    std::vector<std::thread> threads;
    Matrix<std::int64_t> in0;
    Matrix<std::int64_t> out0;
    Matrix<std::int64_t> in1;
    Matrix<std::int64_t> out1;
    Matrix<std::int64_t> plain;
    Matrix<std::int64_t> shuffled_plain;
    Matrix<std::int64_t> shuffled_plain_compare(10, 3);
    std::vector<std::size_t> perm;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { he_share_shuffle_test(i, in0, in1, out0, out1, perm); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    shuffled_plain = out0 + out1;
    plain = in0 + in1;

    for (std::size_t i = 0; i < static_cast<std::size_t>(shuffled_plain_compare.rows()); i++) {
        shuffled_plain_compare.row(i) = plain.row(perm[i]);
    }

    EXPECT_EQ(shuffled_plain, shuffled_plain_compare);
}

TEST_F(AbyTest, quick_sort_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { quick_sort_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, quick_sort_by_column_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { quick_sort_by_column_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, bool_ss_test) {
    for (int op = static_cast<int>(BooleanhOP::And); op <= static_cast<int>(BooleanhOP::Not); ++op) {
        std::vector<std::thread> threads;

        for (std::size_t i = 0; i < 2; ++i) {
            threads.emplace_back([&, i]() { bool_ss_test(i, static_cast<BooleanhOP>(op)); });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }
}

TEST_F(AbyTest, split_by_condition_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { split_by_condition_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, sigmoid_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { sigmoid_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, groupby_sum_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { groupby_sum_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, groupby_count_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { groupby_count_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, groupby_max_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { groupby_max_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

TEST_F(AbyTest, groupby_min_test) {
    std::vector<std::thread> threads;

    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back([&, i]() { groupby_min_test(i); });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}

}  // namespace duet
}  // namespace petace
