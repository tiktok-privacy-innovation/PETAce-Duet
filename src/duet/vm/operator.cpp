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

#include "duet/vm/operator.h"

namespace petace {
namespace duet {

void ShareDoublePA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PrivateMatrix<double>>& in = reg->get_data<PrivateMatrix<double>>(addr[0]);
    const std::shared_ptr<ArithMatrix> out = reg->get_data<ArithMatrix>(addr[1]);
    duet->share(net, *in, *out);
}

void RevealDoubleAP::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PrivateMatrix<double>>& out = reg->get_data<PrivateMatrix<double>>(addr[1]);
    duet->reveal(net, *in, *out);
}

void RevealDoubleAC::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& out = reg->get_data<PublicMatrix<double>>(addr[1]);
    duet->reveal(net, *in, *out);
}

void RevealBoolBP::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& in = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<PrivateMatrixBool>& out = reg->get_data<PrivateMatrixBool>(addr[1]);
    duet->reveal_bool(net, *in, *out);
}

void AddAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->add(*x, *y, *z);
}

void SubAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->sub(*x, *y, *z);
}

void MulBBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_mul(net, *x, *y, *z);
}

void MulAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_mul(net, *x, *y, *z);
}

void SubACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& y = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->share_sub_public_double(*x, *y, *z);
}

void MulACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& x = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->scalar_mul(*x, *y, *z);
}

void DivACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_share_div_public(*x, *y, *z);
}

void DivAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_div(net, *x, *y, *z);
}

void OrCBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PublicMatrixBool>& x = reg->get_data<PublicMatrixBool>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_or(*x, *y, *z);
}

void GreaterAAB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->greater(net, *x, *y, *z);
}

void LessAAB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->less(net, *x, *y, *z);
}

void GreatEqualAAB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->greater_equal(net, *x, *y, *z);
}

void EqualAAB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->equal(net, *x, *y, *z);
}

void SumAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    duet->sum(*x, *y);
}

void LessThanZeroAB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    duet->less_than_zero(net, *x, *y);
}

void MultiplexerBAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->multiplexer(net, *x, *y, *z);
}

void MultiplexerBAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    const std::shared_ptr<ArithMatrix>& a = reg->get_data<ArithMatrix>(addr[3]);
    duet->multiplexer(net, *x, *y, *z, *a);
}

void ShufflePA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PrivateMatrix<double>>& in = reg->get_data<PrivateMatrix<double>>(addr[0]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[1]);
    duet->shuffle(net, *in, *out);
}

void ShuffleAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    duet->shuffle(net, *x, *y);
}

void RowMajorArgmaxAndMaxAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& max_index = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& max_value = reg->get_data<ArithMatrix>(addr[2]);
    duet->row_major_argmax_and_max(net, *x, *max_index, *max_value);
}

}  // namespace duet
}  // namespace petace
