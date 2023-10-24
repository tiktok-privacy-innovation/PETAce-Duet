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

        if (party_id == 0) {
            private_.matrix() = public_a.matrix().cwiseProduct(private_b.matrix());
        }

        return;
    }

    void public_or_share_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto duet_test = std::make_shared<Duet>(net, party_id);
        PublicMatrixBool public_a(4, 3);
        PrivateMatrixBool private_a(4, 3, 0);

        BoolMatrix cipher_a(4, 3);
        BoolMatrix cipher_b(4, 3);

        get_rand_public_bool_matrix(public_a);

        get_rand_plain_bool_matrix(cipher_a.matrix());

        duet_test->reveal_bool(net, cipher_a, private_a);
        duet_test->elementwise_bool_or(public_a, cipher_a, cipher_b);
        duet_test->reveal_bool(net, cipher_b, reveal_private_bool_);
        private_bool_.resize(4, 3);
        if (party_id == 0) {
            for (size_t i = 0; i < private_a.size(); ++i) {
                private_bool_(i) = public_a(i) | private_a(i);
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
        aby_test->elementwise_share_div_public(cipher_b, public_a, cipher_c);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);

        if (party_id == 0) {
            private_.matrix() = private_b.matrix().cwiseQuotient(public_a.matrix());
        }

        return;
    }

    void public_scalar_mul_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicDouble public_a = 0.25;
        PrivateMatrix<double> private_b(4, 3, 0);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);
        aby_test->scalar_mul(public_a, cipher_b, cipher_c);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);

        if (party_id == 0) {
            private_.matrix() = public_a * private_b.matrix();
        }
        return;
    }

    void public_scalar_sub_test(std::size_t party_id) {
        std::shared_ptr<petace::network::Network> net;
        if (party_id == 0) {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params0);
        } else {
            net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params1);
        }
        auto aby_test = std::make_shared<Duet>(net, party_id);

        PublicDouble public_a = 0.25;
        PrivateMatrix<double> private_b(4, 3, 0);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_private_matrix(private_b, party_id);

        aby_test->share(net, private_b, cipher_b);
        aby_test->share_sub_public_double(cipher_b, public_a, cipher_c);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_c, reveal_private_);

        if (party_id == 0) {
            private_.matrix() = private_b.matrix().array() - public_a;
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
        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        ArithMatrix cipher_c(4, 3);
        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);

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

        PrivateMatrix<double> private_a(4, 3, 0);
        PrivateMatrix<double> private_b(4, 3, 0);
        PrivateMatrix<double> private_c(4, 3, 0);
        PrivateMatrix<double> private_d(4, 3, 0);

        ArithMatrix cipher_a(4, 3);
        ArithMatrix cipher_b(4, 3);
        BoolMatrix cipher_c(4, 3);
        ArithMatrix cipher_d(4, 3);
        ArithMatrix cipher_e(4, 3);

        get_rand_private_matrix(private_a, party_id);
        get_rand_private_matrix(private_b, party_id);
        get_rand_private_matrix(private_c, party_id);
        get_rand_private_matrix(private_d, party_id);

        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);
        aby_test->share(net, private_d, cipher_e);

        aby_test->greater(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_e, cipher_d, cipher_b);

        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(4, 3);

        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) > private_b(i)) ? private_c(i) : private_d(i);
            }
        }

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
        aby_test->reveal_bool(net, cipher_c, reveal_private_bool_);

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
        aby_test->reveal_bool(net, cipher_c, reveal_private_bool_);

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

        private_b(0) = private_a(0);
        aby_test->share(net, private_a, cipher_a);
        aby_test->share(net, private_b, cipher_b);
        aby_test->share(net, private_c, cipher_d);

        aby_test->equal(net, cipher_a, cipher_b, cipher_c);
        aby_test->multiplexer(net, cipher_c, cipher_d, cipher_b);
        reveal_private_.set_party_id(0);
        aby_test->reveal(net, cipher_b, reveal_private_);

        private_.resize(4, 3);
        if (party_id == 0) {
            for (std::size_t i = 0; i < private_c.size(); i++) {
                private_(i) = (private_a(i) == private_b(i)) ? private_c(i) : 0;
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
        aby_test->reveal_bool(net, res, reveal_private_bool_);

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
        std::size_t test_size = 10;
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
        std::vector<petace::solo::ahepaillier::BigNum> test;
        test.resize(1);
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
                test[0] = cipher_1(i);
                cipher_2.ciphers() = petace::solo::ahepaillier::Ciphertext(*aby_test->get_pk(), test);
                aby_test->decrypt(cipher_2, plain_2);
                reveal_private_(i) = plain_2(0);
            }
            for (std::size_t i = 5; i < 10; ++i) {
                test[0] = cipher_1(0, i);
                cipher_2.ciphers() = petace::solo::ahepaillier::Ciphertext(*aby_test->get_pk(), test);
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
        private_bool_.set_party_id(0);
        private_bool_.resize(1, 5);
        reveal_private_bool_.set_party_id(0);
        reveal_private_bool_.resize(1, 5);
        if (party_id == 0) {
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

    void row_major_argmax_and_max_test(std::size_t party_id, Matrix<std::int64_t>& aim_index,
            Matrix<std::int64_t>& real_index, Matrix<double>& aim_value, Matrix<double>& real_value) {
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

        aby_test->row_major_argmax_and_max(net, cipher_a, cipher_b, cipher_c);

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

TEST_F(AbyTest, public_or_share_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        public_or_share_test(1);
        exit(EXIT_SUCCESS);
    } else {
        public_or_share_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_bool_, private_bool_));
    }
}

TEST_F(AbyTest, public_scalar_mul_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        public_scalar_mul_test(1);
        exit(EXIT_SUCCESS);
    } else {
        public_scalar_mul_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
    }
}

TEST_F(AbyTest, public_scalar_sub_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        public_scalar_sub_test(1);
        exit(EXIT_SUCCESS);
    } else {
        public_scalar_sub_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
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

TEST_F(AbyTest, paillier_a2h_test) {
    pid_t pid;
    int status;

    if ((pid = fork()) < 0) {
        status = -1;
    } else if (pid == 0) {
        paillier_a2h_test(1);
        exit(EXIT_SUCCESS);
    } else {
        paillier_a2h_test(0);
        while (waitpid(pid, &status, 0) < 0) {
            if (errno != EINTR) {
                status = -1;
                break;
            }
        }
        ASSERT_TRUE(is_equal_private_matrix(reveal_private_, private_, 0.001));
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

TEST_F(AbyTest, row_major_argmax_and_max_test) {
    std::vector<std::thread> threads;
    Matrix<std::int64_t> aim_index;
    Matrix<std::int64_t> real_index;
    Matrix<double> aim_value;
    Matrix<double> real_value;
    for (std::size_t i = 0; i < 2; ++i) {
        threads.emplace_back(
                [&, i]() { row_major_argmax_and_max_test(i, aim_index, real_index, aim_value, real_value); });
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

}  // namespace duet
}  // namespace petace
