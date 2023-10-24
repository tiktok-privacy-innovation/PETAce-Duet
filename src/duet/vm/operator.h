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

class DivACdA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class OrCBB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& /*net*/,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class DivAAA : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreaterAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class LessAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class GreatEqualAAB : public OperatorBase {
public:
    void run(const std::shared_ptr<Duet>& duet, const std::shared_ptr<network::Network>& net,
            const std::shared_ptr<Register>& reg, const std::vector<RegisterAddress>& addr);
};

class EqualAAB : public OperatorBase {
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

}  // namespace duet
}  // namespace petace
