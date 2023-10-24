# Copyright 2023 TikTok Pte. Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
import petace.duet.securenumpy as snp
from petace.duet.pyduet import NetParams, DuetVM, NetScheme


def numpy_linear_regression(model, sample):
    return np.sum(model * sample)

def securenumpy_linear_regression(model, sample):
    return snp.sum(model * sample)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='PETAce-Duet demo.')
    parser.add_argument("-p", "--party", type=int, help="which party")
    parser.add_argument("--port0", type=int,
                        help="port of party 0, defalut 8089", default=8089)
    parser.add_argument("--port1", type=int,
                        help="port of party 1, defalut 8090", default=8090)
    parser.add_argument("--host", type=str,
                        help="host of this party", default="127.0.0.1")

    args = parser.parse_args()
    party = args.party
    port0 = args.port0
    port1 = args.port1
    host = args.host

    net_params = NetParams()
    if party == 0:
        net_params.remote_addr = host
        net_params.remote_port = port1
        net_params.local_port = port0
    else:
        net_params.remote_addr = host
        net_params.remote_port = port0
        net_params.local_port = port1

    duet = DuetVM(net_params, NetScheme.SOCKET, party)

    data0 = np.empty((0, 0))
    data1 = np.empty((0, 0))
    if duet.party_id() == 0:
        data0 = np.array([[2], [2], [3], [5]], dtype=np.float64)
        print("Party 0's data:")
        print(data0)
    else:
        data1 = np.array([[3], [5], [7], [9]], dtype=np.float64)
        print("Party 1's data:")
        print(data1)

    sample = snp.Private(data0, 0, duet)
    model = snp.Private(data1, 1, duet)

    ret = securenumpy_linear_regression(model, sample)

    p0 = ret.reveal_to(0)
    if duet.party_id() == 0:
        print("Model inference's ret is:")
        print(p0.to_numpy())
