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

#include <string>
#include <vector>

namespace petace {
namespace duet {

class Instruction {
public:
    explicit Instruction(const std::vector<std::string>& in);

    bool operator==(const Instruction& other) const;

    const std::vector<std::string>& data() const;

private:
    std::vector<std::string> code_{};
};

}  // namespace duet
}  // namespace petace

namespace std {

// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template <>
struct hash<petace::duet::Instruction> {
    std::size_t operator()(const petace::duet::Instruction& k) const {
        std::size_t seed = k.data().size();
        std::hash<std::string> hasher;
        for (auto& i : k.data()) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

}  // namespace std
