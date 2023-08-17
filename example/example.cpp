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

#include <iostream>

#include "duet_example.h"
#include "tclap/CmdLine.h"

#include "network/net_factory.h"
#include "network/network.h"

int main(int argc, char** argv) {
    try {
        TCLAP::CmdLine cmd("duet demo", ' ', "0.1");
        // add cmd params
        TCLAP::ValueArg<std::size_t> party_arg("p", "party", "which party", false, 0, "std::size_t");
        TCLAP::ValueArg<std::string> host_arg("", "host", "host of party 0", false, "127.0.0.1", "string");
        TCLAP::ValueArg<std::uint16_t> port0_arg(
                "", "port0", "port of party 0, defalut 8089", false, 8089, "std::uint16_t");
        TCLAP::ValueArg<std::uint16_t> port1_arg(
                "", "port1", "port of party 1, defalut 8090", false, 8090, "std::uint16_t");
        TCLAP::ValueArg<std::string> which_case_arg("e", "example", "which example to run", false, "all", "string");

        cmd.add(party_arg);
        cmd.add(host_arg);
        cmd.add(port0_arg);
        cmd.add(port1_arg);
        cmd.add(which_case_arg);

        // parse
        cmd.parse(argc, argv);
        std::size_t party = party_arg.getValue();
        std::string host = host_arg.getValue();
        std::uint16_t port0 = port0_arg.getValue();
        std::uint16_t port1 = port1_arg.getValue();
        std::string which_case = which_case_arg.getValue();

        // print info
        std::cout << "party: " << party << std::endl;
        std::cout << "Running: " << which_case << std::endl;
        std::cout << "initialize network" << std::endl;

        // init net and duet
        petace::network::NetParams net_params;
        if (party == 0) {
            net_params.remote_addr = "127.0.0.1";
            net_params.remote_port = port1;
            net_params.local_port = port0;
        } else {
            net_params.remote_addr = "127.0.0.1";
            net_params.remote_port = port0;
            net_params.local_port = port1;
        }
        auto net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params);
        auto duet = std::make_shared<petace::duet::Duet>(net, party);

        if (which_case == "add") {
            add_test(net, duet);
        } else if (which_case == "mul") {
            mul_test(net, duet);
        } else if (which_case == "shuffle") {
            plain_shuffle_test(net, duet);
        } else if (which_case == "less than zero") {
            less_than_zero_test(net, duet);
        } else if (which_case == "millionaires") {
            millionaires(net, duet);
        } else if (which_case == "ppml") {
            ppml(net, duet);
        } else if (which_case == "distance") {
            distance(net, duet);
        } else if (which_case == "all") {
            std::cout << "Running: "
                      << "add" << std::endl;
            add_test(net, duet);
            std::cout << "Running: "
                      << "mul" << std::endl;
            mul_test(net, duet);
            std::cout << "Running: "
                      << "shuffle" << std::endl;
            plain_shuffle_test(net, duet);
            std::cout << "Running: "
                      << "less than zero" << std::endl;
            less_than_zero_test(net, duet);
            std::cout << "Running: "
                      << "millionaires" << std::endl;
            millionaires(net, duet);
            std::cout << "Running: "
                      << "ppml" << std::endl;
            ppml(net, duet);
            std::cout << "Running: "
                      << "distance" << std::endl;
            distance(net, duet);
        } else {
            std::cout << "unknown case" << std::endl;
        }
    } catch (TCLAP::ArgException& e) {
        std::cerr << "err: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }
    return 0;
}
