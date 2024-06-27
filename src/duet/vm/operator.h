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
#include <vector>

#include "network/network.h"

#include "duet/duet.h"
#include "duet/util/defines.h"
#include "duet/vm/register.h"

namespace petace {
namespace duet {

class OperatorBase {
public:
    virtual ~OperatorBase() = default;
    virtual void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) = 0;
};

class ShareDoublePA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) override;
};

class ShareBoolPA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) override;
};

class RevealDoubleAP : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class RevealDoubleAC : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class RevealBoolBP : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class AddAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class AddACdA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class AddACdmA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class SubAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class SubACdA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class SubACdmA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MulBBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MulAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MulACdA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MulACdmA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MatMulAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MatMulCdmAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MatMulACdmA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class DivACdmA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class DivACdA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class DivAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class OrCBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class OrBBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class AndCBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class AndBBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class XorCBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class XorBBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class NotBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreaterAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreaterACdB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreaterACdmB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class LessAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class LessACdB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class LessACdmB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreatEqualAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreatEqualACdmB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class EqualAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class EqualACdmB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class SumAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class LessThanZeroAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MultiplexerBAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class MultiplexerBAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class ShufflePA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class ShuffleAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class RowMajorArgmaxAndMaxAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GroupThenSumByGroupedCountCdAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class QuikSortA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class SplitByConditionBAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class SigmoidAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class QuikSortByColumnAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GroupBySumAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GroupByCountAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GroupByMaxAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GroupByMinAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class TransposeAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

template <typename T>
class Reshape : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& /*duet*/, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
        const std::shared_ptr<T>& x = reg->get_data<T>(addr[0]);
        const std::shared_ptr<PublicIndex>& r = reg->get_data<PublicIndex>(addr[1]);
        const std::shared_ptr<PublicIndex>& c = reg->get_data<PublicIndex>(addr[2]);
        const std::shared_ptr<T>& z = reg->get_data<T>(addr[3]);
        z->resize(*r, *c);
        for (std::size_t i = 0; i < x->size(); i++) {
            z->shares()(i) = x->shares()(i);
        }
    }
};

class SetItemA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& /*duet*/, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

template <typename T>
class Resize : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& /*duet*/, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr) {
        const std::shared_ptr<T>& x = reg->get_data<T>(addr[0]);
        const std::shared_ptr<PublicIndex>& r = reg->get_data<PublicIndex>(addr[1]);
        const std::shared_ptr<PublicIndex>& c = reg->get_data<PublicIndex>(addr[2]);
        const std::shared_ptr<T>& z = reg->get_data<T>(addr[3]);
        z->resize(*r, *c);
        for (std::size_t i = 0; i < z->size(); i++) {
            z->shares()(i) = x->shares()(i % (x->size()));
        }
    }
};

}  // namespace duet
}  // namespace petace
