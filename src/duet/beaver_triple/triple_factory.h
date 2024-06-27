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

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "duet/beaver_triple/arithmetic_triple.h"
#include "duet/beaver_triple/arithmetic_triple_fhe.h"
#include "duet/beaver_triple/arithmetic_triple_ot.h"
namespace petace {
namespace duet {

enum class TripleScheme : std::uint32_t { OT = 0, FHE = 1 };

using TripleEngineCreator =
        std::function<std::shared_ptr<ArithmeticTriple>(const std::shared_ptr<network::Network>& net,
                std::size_t party_id, std::shared_ptr<PRNG> prng, const std::shared_ptr<ObliviousTransfer>& ot)>;

/**
 * @brief Beaver Triple Generator factory.
 *
 * @see Refer to README.md for more details.
 *
 * @par Example.
 * Refer to example.cpp.
 */
class TripleFactory {
public:
    /**
     * @brief Get a factory singleton.
     *
     * @return Return factory instance.
     */
    static TripleFactory& get_instance() {
        static TripleFactory factory;
        return factory;
    }

    /**
     * @brief Get triple generator object based on net scheme.
     *
     * @param[in] scheme: Name of the triple generation scheme
     */
    std::shared_ptr<ArithmeticTriple> build(const TripleScheme& scheme, const std::shared_ptr<network::Network>& net,
            std::size_t party_id, std::shared_ptr<PRNG> prng, const std::shared_ptr<ObliviousTransfer>& ot) {
        auto where = creator_map_.find(scheme);
        if (where == creator_map_.end()) {
            throw std::invalid_argument("net creator is not registered.");
        }
        return where->second(net, party_id, prng, ot);
    }

    /**
     * @brief Register triple schemes
     *
     * @param[in] scheme: Name of the triple generation scheme
     * @param[in] creator: creator of the triple generation scheme
     */
    void register_triple(const TripleScheme& scheme, TripleEngineCreator creator) {
        creator_map_.insert(std::make_pair(scheme, creator));
    }

protected:
    TripleFactory() {
        register_triple(TripleScheme::OT, create_triple_engine_ot);
        register_triple(TripleScheme::FHE, create_triple_engine_fhe);
    }
    ~TripleFactory() {
    }
    TripleFactory(const TripleFactory&) = delete;
    TripleFactory& operator=(const TripleFactory&) = delete;
    TripleFactory(TripleFactory&&) = delete;
    TripleFactory& operator=(TripleFactory&&) = delete;

private:
    std::map<TripleScheme, TripleEngineCreator> creator_map_;
};

}  // namespace duet
}  // namespace petace
