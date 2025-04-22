# Helix: Scalable Multi-Party Machine Learning Inference against Malicious Adversaries
This is the repo of the paper `Helix: Scalable Multi-Party Machine Learning Inference against Malicious Adversaries`. It provides a scalable framework for maliciously secure privacy-preserving machine learning in the honest majority setting. This work builds on [hmmpc-public](https://github.com/f7ed/hmmpc-public) and the coin-tossing protocol uses the backend of [MP-SPDZ](https://github.com/data61/MP-SPDZ), which relies on the commitment scheme.

### Table of Contents
1. [Requirement](#requirement)
2. [Local Inference on a Single Machine](#local-inference-on-a-single-machine)
3. [Remote Inference over Multiple Servers](#remote-inference-over-multiple-servers)
    1. [Traffic Control to Simulate LAN and WAN](#traffic-control-to-simulate-lan-and-wan)
    2. [Separate Execution with Scripts](#separate-execution-with-scripts)

## Requirement 

This code should work on macOS and most Linux distributions.

- The experiments in the paper are conducted with Ubuntu 22.04 LTS.

Requirement packages:

- g++ or clang
- cmake
- make
- OpenSSL
- libsodium library
- Eigen library

On Linux:

```shell
sudo apt update

# install cmake openssl libsodium tc 
# tc (traffic control) is used to simulate LAN & WAN
sudo apt install -y build-essential cmake libssl-dev libsodium-dev iproute2 
```
We have added the Eigen's header files in our artifact, located in `Eigen/` subdirectory, which are the only required files to compile with Eigen.



## Local Inference on a Single Machine

It will make an inference on a 3-layer DNN in the setting of 3PC with the threshold $t=1$.

```shell
# Usage: ./Scripts/inference.sh <npc>
./Scripts/inference.sh 3
```

This script, located in `Scripts/inference.sh`, is only recommended for testing the functionality of the inference since it only simulates $n$ parties locally.

Each party occupies different port numbers in the local IP address `127.0.0.1`. They execute the MPC protocol described in `ML/inference.cpp`. 

The output of the terminal is defaulted to P0's output.


## Remote Inference over Multiple Servers

We conducted our real experiments on several servers to simulate a large number of parties, giving us reliable results.

In the setting of 3PC, 7PC, and 11PC, a single party controls a single server.

In the settings of 21PC, 31PC, and 63PC, we run multi-parties on a single server simultaneously with different ports. We guarantee the number of parties on different servers differs by up to one for load balancing. The assignment is reflected in the IP file. All IP files are located in `Inference/IP_HOSTS`.

### Traffic Control to Simulate LAN and WAN

We use the Linux `tc` command to set the latency and the bandwidth.

For the **LAN setting**, we need to set the latency on the NIC **with regard to the IP address for parties on the different servers** and the **local loopback NIC for parties on the same server.** Besides, we need to set the bandwidth for the loopback NIC. And in our experimental setup, the bandwidth between different servers is inherently 10 Gbps.

```shell
# Set latency and bandwidth on the local NIC
sudo tc qdisc add dev lo root netem delay 0.13ms rate 10000mbit
# Set latency on the NIC with regard to the IP address
sudo tc qdisc add dev eth0 root netem delay 0.13ms 
```

For the WAN setting, we need to set the latency on the NIC **with regard to the IP address for parties on the different servers** and the **local loopback NIC for parties on the same server.** You can run the mini script to check the ping time for each pair of parties. Besides, we need to set the bandwidth for the corresponding NIC.

```shell
# Set latency and bandwidth on the local NIC
sudo tc qdisc add dev lo root netem delay 20ms rate 100mbit 
# Set latency and bandwidth on the NIC with regard to the IP address
sudo tc qdisc add dev eth0 root netem delay 20ms rate 100mbit 
```

### Separate Execution with Scripts

We also support various scripts to facilitate the separate execution with different settings.

All scripts are located in `Scripts/infer-<npc>`  for $n=3,7,11,21,31,63$ .

Let me show the details of how to conduct the experiment for `7pc` with 7 separate servers.

1. Prepare 7 servers with requirements installed, and separately compile `inference.x`. Then run the following command to perform the inference locally to ensure the environment is ready.

   ```shell
   make -j8 inference.x
   
   ./Scripts/inference.sh 3
   ```

2. Use `tc` to set latency and bandwidth to simulate the LAN or WAN.

3. Modify the IP file `Inference/IP_HOSTS/IP_7` to assign the IP of each server to the corresponding party. You can also assign the port numbers in the form of `IP:port number`.

4. For each separate server, run the script `./Scripts/infer-7pc/inference_7pc.sh <id>`. With such a script, it allows us to run the same program with different options by just modifying the corresponding script while only compiling once. Note that `<id>` stats from `0`  to `6`, which must correspond to the assignment specified in the IP file.

   ```shell
   ./Scripts/infer-7pc/inference_7pc.sh 0 # Server 0
   ./Scripts/infer-7pc/inference_7pc.sh 1 # Server 1
   ./Scripts/infer-7pc/inference_7pc.sh 2 # Server 2
   ./Scripts/infer-7pc/inference_7pc.sh 3 # Server 3
   ./Scripts/infer-7pc/inference_7pc.sh 4 # Server 4
   ./Scripts/infer-7pc/inference_7pc.sh 5 # Server 5
   ./Scripts/infer-7pc/inference_7pc.sh 6 # Server 6
   ```

5. Modify the script to execute the protocol with a different setting. The options include


   -  `PRIME=PR61`: Specified arithmetic field over a Mersenne prime. In the malicious setting, to achieve the statistical security, we can only use `PR61`.

      - `PR61`: Mersenne prime $2^{61}-1$.

   - `NETWORK=SecureML`: Only support SecureML, Sarda and MiniONN now.

   - `DATASET=MNIST`: We only support MNIST dataset right now.

   - `NET=WAN`: Conduct the experiment in the WAN setting.

   - `CORES=8`: Control the number of threads for Eigen's algorithms with the use of multi-threading.

   - `MULT_COMPRESSION=32`: The compression parameter utilized in the multiplication verification protocol. It is set to 4 for LAN and 32 for WAN in our experimental evaluations in the paper.

   - `VERIFY_PROTOCOL=GS20`: The protocol for multiplication verification includes [GS20](https://eprint.iacr.org/2020/134.pdf), [LN17](https://eprint.iacr.org/2017/816.pdf) and [Falcon](https://arxiv.org/pdf/2004.02229) now. The experiments in the paper are conducted using the GS20 protocol.

   - `IP_FILE=Inference/IP_HOSTS/IP_$NPC` : The default location in this example is `Inference/IP_HOSTS/IP_7`. 

   - `OFFLINE_ARG=`: The location of the file, which specifies the number of various random sharings required to be generated in the preprocessing phase.  We precompute it for some settings, stored in the following location.

     ```bash
     OFFLINE_ARG=Inference/$NETWORK/offline/${PRIME}_offline_b${TEST_DATA_SIZE}_c${MULT_COMPRESSION}_${VERIFY_PROTOCOL}.txt
     # In this example, the default location is Inference/SecureML/offline/PR61_offline_b1_c32_GS20.txt
     ```

   - `Output_file=` : The location of the output file. In order to ensure reliable results, a single execution contains 10 independent runs. For simplicity, we only collect the output of P0's. The default location is set to 

     ```bash
     Output_file=Inference/$NETWORK/test${PRIME}/${NET}/${NETWORK}.b${TEST_DATA_SIZE}_c${MULT_COMPRESSION}_${NPC}PC
     # In this example, the default location is Inference/SecureML/testPR61/WAN/SecureML.b1_c32_7PC
     ```