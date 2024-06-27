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

void ShareBoolPA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PrivateMatrixBool>& in = reg->get_data<PrivateMatrixBool>(addr[0]);
    const std::shared_ptr<BoolMatrix> out = reg->get_data<BoolMatrix>(addr[1]);
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
    duet->reveal(net, *in, *out);
}

void AddAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->add(*x, *y, *z);
}

void AddACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& y = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->add(*x, *y, *z);
}

void AddACdmA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
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

void SubACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& y = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->sub(*x, *y, *z);
}

void SubACdmA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->sub(*x, *y, *z);
}

void MulBBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_and(net, *x, *y, *z);
}

void MulAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_mul(net, *x, *y, *z);
}

void MulACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& x = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_mul(*x, *y, *z);
}

void MulACdmA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_mul(*y, *x, *z);
}

void MatMulAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->mat_mul(net, *x, *y, *z);
}

void MatMulCdmAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PublicMatrix<double>>& x = reg->get_data<PublicMatrix<double>>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->mat_mul(*x, *y, *z);
}

void MatMulACdmA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->mat_mul(*x, *y, *z);
}

void DivACdmA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_div(*x, *y, *z);
}

void DivACdA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& y = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[2]);
    duet->elementwise_div(*x, *y, *z);
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

void OrBBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_or(net, *x, *y, *z);
}

void AndCBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PublicMatrixBool>& x = reg->get_data<PublicMatrixBool>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_and(*x, *y, *z);
}

void AndBBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_and(net, *x, *y, *z);
}

void XorCBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PublicMatrixBool>& x = reg->get_data<PublicMatrixBool>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_xor(*x, *y, *z);
}

void XorBBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& y = reg->get_data<BoolMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->elementwise_bool_xor(*x, *y, *z);
}

void NotBB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& x = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[1]);
    duet->elementwise_bool_not(*x, *z);
    for (std::size_t i = 0; i < z->size(); i++) {
        z->shares()(i) = z->shares()(i) & 0x1;
    }
}

void GreaterAAB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->greater(net, *x, *y, *z);
}

void GreaterACdB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& y = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->greater(net, *x, *y, *z);
}

void GreaterACdmB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
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

void LessACdB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicDouble>& y = reg->get_data<PublicDouble>(addr[1]);
    const std::shared_ptr<BoolMatrix>& z = reg->get_data<BoolMatrix>(addr[2]);
    duet->less(net, *x, *y, *z);
}

void LessACdmB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
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

void GreatEqualACdmB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
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

void EqualACdmB::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicMatrix<double>>& y = reg->get_data<PublicMatrix<double>>(addr[1]);
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
    duet->argmax_and_max(net, *x, *max_index, *max_value);
}

void GroupThenSumByGroupedCountCdAA::run(const std::shared_ptr<Duet>& duet,
        const std::shared_ptr<network::Network>& /*net*/, const std::shared_ptr<Register>& reg,
        const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<PublicMatrix<double>>& group_count = reg->get_data<PublicMatrix<double>>(addr[0]);
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[2]);
    duet->group_then_sum_by_grouped_count(*group_count, *in, *out);
}

void QuikSortA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    duet->quick_sort(net, *in);
}

void SplitByConditionBAAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<BoolMatrix>& cond = reg->get_data<BoolMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out_0 = reg->get_data<ArithMatrix>(addr[2]);
    const std::shared_ptr<ArithMatrix>& out_1 = reg->get_data<ArithMatrix>(addr[3]);
    duet->split_by_condition(net, *cond, *in, *out_0, *out_1);
}

void SigmoidAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[1]);
    duet->sigmoid(net, *in, *out);
}

void QuikSortByColumnAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<PublicIndex>& col_index = reg->get_data<PublicIndex>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[2]);
    duet->quick_sort(net, *col_index, *in, *out);
}

void GroupBySumAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& encoding = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[2]);

    duet->groupby_sum(net, *in, *encoding, *out);
}

void GroupByCountAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& encoding = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[2]);

    duet->groupby_count(*in, *encoding, *out);
}

void GroupByMaxAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& encoding = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[2]);

    duet->groupby_max(net, *in, *encoding, *out);
}

void GroupByMinAA::run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& in = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& encoding = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<ArithMatrix>& out = reg->get_data<ArithMatrix>(addr[2]);

    duet->groupby_min(net, *in, *encoding, *out);
}

void TransposeAA::run(const std::shared_ptr<Duet>& /*duet*/, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& y = reg->get_data<ArithMatrix>(addr[1]);
    y->resize(x->cols(), x->rows());
    y->shares() = x->shares().transpose();
}

void SetItemA::run(const std::shared_ptr<Duet>& /*duet*/, const std::shared_ptr<network::Network>& /*net*/,
        const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
    const std::shared_ptr<ArithMatrix>& x = reg->get_data<ArithMatrix>(addr[0]);
    const std::shared_ptr<ArithMatrix>& z = reg->get_data<ArithMatrix>(addr[1]);
    const std::shared_ptr<PublicIndex>& start_row = reg->get_data<PublicIndex>(addr[2]);
    const std::shared_ptr<PublicIndex>& start_col = reg->get_data<PublicIndex>(addr[3]);
    const std::shared_ptr<PublicIndex>& block_row = reg->get_data<PublicIndex>(addr[4]);
    const std::shared_ptr<PublicIndex>& block_col = reg->get_data<PublicIndex>(addr[5]);
    z->shares().block(static_cast<Eigen::Index>(*start_row), static_cast<Eigen::Index>(*start_col),
            static_cast<Eigen::Index>(*block_row), static_cast<Eigen::Index>(*block_col)) = x->shares();
}

}  // namespace duet
}  // namespace petace
