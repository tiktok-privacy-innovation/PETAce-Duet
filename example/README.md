# PETAce-Duet Examples

We illustrate three examples here to show you how to use PETAce-Duet.

## How to run the following examples
The example codes of all three use cases are available in `duet_example.cpp`.
You should set the option `DUET_BUILD_EXAMPLE` to be on when building Duet, and the binary file `duet_example` will be found in `/build/bin`.
To run the examples, open two terminal sessions and run this command in the first terminal:

```bash
./build/bin/duet_example -p 0 -e millionaires
```

Then this terminal will act as party 0 and execute the millionaire example.
Similarly, run the following command in the second terminal.

```bash
./build/bin/duet_example -p 1 -e millionaires
```

Then you just need to follow the prompts to provide the inputs, and the result will be revealed to party 0.

If you want to run the other two examples, simply replace `millionaires` with `ppml` for Example 2, or `distance` for Example 3.

## Use case 1: The Millionaires' problem

### Background

We have two millionaires negotiating collaborations.
They are curious about who is richer.
However, they do not want to reveal their actual wealth to each other, as this information is treated as a commercial secret.

### How to solve this with MPC

MPC parties in this case are two millionaires (let's name them Alice and Bob), and the private inputs are their actual wealth $W_{Alice}$ and $W_{Bob}$.
The output of the protocol is a bit: $b=1$, if $W_{Alice} < W_{Bob}$; and $b=0$, otherwise.
Before building protocols, we should initialize the Duet instance first.
You can follow the instructions below to initialize both parties on your local machine:

```c++
petace::network::NetParams net_params;
if (party == 0) {
    net_params.remote_addr = "127.0.0.1";
    net_params.remote_port = 8089;
    net_params.local_port = 8090;
} else {
    net_params.remote_addr = "127.0.0.1";
    net_params.remote_port = 8090;
    net_params.local_port = 8089;
}
auto net = petace::network::NetFactory::get_instance().build(petace::network::NetScheme::SOCKET, net_params);
auto duet = std::make_shared<petace::duet::Duet>(net, party);
```

Now we have Duet instance ready, we can start building the protocol.

```c++
std::size_t num_rows = 1;
std::size_t num_cols = 1;
petace::duet::PlainMatrix<double> plain_W_Alice(num_rows, num_cols);
petace::duet::PlainMatrix<double> plain_W_Bob(num_rows, num_cols);
petace::duet::ArithMatrix share_W_Alice(num_rows, num_cols);
petace::duet::ArithMatrix share_W_Bob(num_rows, num_cols);

if (duet->party() == 0) {
    // Alice
    // assume Alice has 1000 million dollar
    plain_W_Alice(0) = 1000;
} else {
    // assume Bob has 2000 million dollar
    plain_W_Bob(0) = 2000;
}

// secret share W_Alice and W_Bob
duet->share(net, 0, plain_W_Alice, share_W_Alice);
duet->share(net, 1, plain_W_Bob, share_W_Bob);

// compute a secure comparison and get protocol result
petace::duet::BoolMatrix boolen_ret(num_rows, num_cols);
petace::duet::PlainMatrix<std::int64_t> ret(num_rows, num_cols);
petace::duet::ArithMatrix A_minus_B(num_rows, num_cols);

duet->sub(share_W_Alice, share_W_Bob, A_minus_B);
duet->less_than_zero(net, A_minus_B, boolen_ret);
// reveal result to Alice
duet->reveal(net, 0, boolen_ret, ret);
// print the result on Alice side
if (duet->party() == 0) {
    std::cout << "reveal to party 0, ret is: " << std::endl;
    std::cout << ret(0) << std::endl;
}
```

## Use case 2: privacy-preserving linear machine learning model inference

### Background

Let's assume that the server has a machine learning model and it can give a recommendation score whether you should buy some company's stock.
The input of the model includes factors such as the company's profitability, the company's costs, the company's competition, and the company's assets.
Investors may want to use this model for investment suggestions, but they do not want to reveal their interests to these companies and stocks publicly.
Meanwhile, the server does not want to make the model public as well so that it can make a profit by providing inference service.
For simplicity, let's assume a linear model, which can be represented by a vector $M = \{m_1, ..., m_n\}$.
Client input can also be represented as a vector $X = \{x_1, x_2, ..., x_n\}$.
The inference result can be computed through the inner product of $M$ and $X$: $x_0m_0+x_1m_1+\cdots+x_nm_n$.

### How to solve this with MPC

First, we initialize the Duet instance with the same instructions as use case 1, then we can build the inference protocol as follows:

```c++
// let's assume the model has 10 parameters.
std::size_t num_rows = 10;
std::size_t num_cols = 1;
petace::duet::PlainMatrix<double> plain_model(num_rows, num_cols);
petace::duet::PlainMatrix<double> plain_input(num_rows, num_cols);
petace::duet::ArithMatrix share_model(num_rows, num_cols);
petace::duet::ArithMatrix share_input(num_rows, num_cols);

if (duet->party() == 0) {
    // Server
    // Server puts private model here
    plain_model = ...;
} else {
    // Client
    // Client puts private input here
    plain_input = ...;
}
// secret share the model and client input
duet->share(net, 0, plain_model, share_model);
duet->share(net, 1, plain_input, share_input);

// inner product computation
petace::duet::ArithMatrix element_mul_res(num_rows, num_cols);
petace::duet::ArithMatrix res(1,1);
duet->elementwise_mul(net, share_model, share_input, element_mul_res);
// column-wise sum, then the result secret share will be in res(0)
duet->sum(element_mul_res, res);
```

If you want to reveal the result to designated parties, simply use `reveal()` to reconstruct a result from secret shares.

# Use case 3:  distance of two vectors
## Background

Computing the distance of two vectors has a wide range of applications such as clustering and recommendation system.
One of the classic definitions of distance is Euclidian distance.
Given two vectors $X$ and $Y$ of the same length $n$, the Euclidian distance is defined as follows.

$$d(X, Y) = \sqrt{\sum_{i = 1}^{n}(x_i - y_i)^2}$$

Let's consider a scenario of COVID-19 contact tracing: Alice wants to figure out if she has been closely contacted with someone with COVID-19.
There is a database storing information about the patients and the locations that they visited.
However, this database is not available to Alice as the information is highly sensitive.
We can use MPC to evaluate this problem, the corresponding computation is to compute the distance between the locations that Alice provides and the locations in the databases.
If one of the distances is close enough, we confirm that Alice has closely contacted with someone with COVID-19.

The location is in the form of coordinates, which can be represented by a pair $(x, y)$.
### How to solve this with MPC

Computing square root is quite expensive in MPC.
However, we can compute the square of the Euclidian distance and reveal it, then we can compute square root in plaintext.
Therefore, the goal of MPC is just to compute $\sum(x_i - y_i)^2$, which only involves addition (subtraction) and multiplication (square).

```c++
// let's assume the the length of vectors is 10.
std::size_t num_rows = 2;
std::size_t num_cols = 1;
petace::duet::PlainMatrix<double> plain_X(num_rows, num_cols);
petace::duet::PlainMatrix<double> plain_Y(num_rows, num_cols);
petace::duet::ArithMatrix share_X(num_rows, num_cols);
petace::duet::ArithMatrix share_Y(num_rows, num_cols);

if (duet->party() == 0) {
    // party_0 puts his private coordinate here
    plain_X = ...;
} else {
    // party_1 puts his private coordinate here
    plain_Y = ...;
}
// secret share the model and client input
duet->share(net, 0, plain_X, share_X);
duet->share(net, 1, plain_Y, share_Y);

petace::duet::ArithMatrix X_minus_Y(num_rows, num_cols);
petace::duet::ArithMatrix X_minus_Y_square(num_rows, num_cols);
petace::duet::ArithMatrix res(1,1);
duet->sub(share_X, share_Y, X_minus_Y);
duet->mul(net, X_minus_Y, X_minus_Y, X_minus_Y_square);
// column-wise sum, then the result secret share will be in res(0)
duet->sum(X_minus_Y_square, res);
```
