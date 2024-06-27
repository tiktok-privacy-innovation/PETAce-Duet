# PETAce-Duet Benchmark

We have prepared a series of performance tests for PETAce-Duet.

## How to run the benchmark

The benchmark codes of all test cases are available in `duet_bench.cpp`.
You should set the option `DUET_BUILD_BENCH` to be on when building Duet, and the binary file `duet_bench` will be found in `/build/bin`.
To run the benchmark cases, open two terminal sessions and run the following command in the first terminal:

```bash
./build/bin/duet_bench -p 0 --log_path ./duet0.log
```
This terminal will act as party 0 and execute the benchmark program.
Similarly, run the following command in the second terminal.

```bash
./build/bin/duet_bench -p 1 --log_path ./duet1.log
```

Then, you can find logs files in the `/build` directory, named `duet0.log` and `duet1.log`, which represent the performance metrics of Party 0 and Party 1, respectively.

## Logging Format
Logging information follows a specific format as follows:

```bash
case <case name> begin <timestamp> <test_iteration> <parameter>
case <case name> end <timestamp> <cost time> <bytes_send> <bytes_received>
```

```bash
case <case name> begin <timestamp> <test_iteration> <parameter>
stage <stage name> begin <timestamp>
stage <stage name> end <timestamp> <cost time> <bytes_send> <bytes_received>
case <case name> end <timestamp> <cost time> <bytes_send> <bytes_received>
```

`case * begin` and `case * end` indicates the start and end of the test case, respectively. A test case may include multiple stages. For example, when calculating the dot product, one must first compute the element-wise multiplication of two vectors, followed by the summation of these products. We denote the start and end of a stage with the labels `stage begin` and `stage end`.

The meanings of other fields are as follows:

- `case name`: the name of the test case
- `timestamp`: the timestamp of the start or end of the test case or stage
- `test_iteration`: the number of the test case executions
- `parameter`: the parameters of the test case
- `stage name`: the name of the stage
- `cost time`: the time cost of the test case or stage
- `bytes_send`: the number of bytes sent
- `bytes_received`: the number of bytes received

## Benchmark with Various Network Conditions

For MPC (Multi-Party Computation), network overhead is a critical metric. We offer a straightforward method to simulate various network conditions.

The following command can be used to simulate a network with a given latency and bandwidth.

```bash
tc qdisc add dev lo root netem delay {latency}ms rate {bandwidth}mbit
```

For example, to simulate a network with 10ms latency and 100Mbps bandwidth, run the following command.

```bash
tc qdisc add dev lo root netem delay 10ms rate 100mbit
```

To remove this limitation, you can use the following command:

```bash
tc qdisc del dev lo root
```
