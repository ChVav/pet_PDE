
/*Example how to bechmark, probably need to restructure code to make use this for testing a dense vs scarce implementation.*/

#include <benchmark/benchmark.h>

#include <vector>

using namespace std;

// Function 1 to benchmark: forward Euler, calculate one timestep
vector<double> one_timestep(double dt, const vector<double>& B, const vector<double>& u_n) {
	//returns u_n+1=u_n+dt*B*u_n
	vector<double> u_n1(u_n.size());
	for (int i = 0; i < u_n.size(); i++) {
		u_n1[i] = u_n[i];
		for (int j = 0; j < u_n.size(); j++) {
			u_n1[i] += dt * B[j + i * u_n.size()] * u_n[j];
		}
	}
	return u_n1;
}

// wrapper used by google benchmarking
//testing how fast the forward Euler implemetnation is for a range of dt
static void one_timestepBench(benchmark::State& s) {
	vector<double> B(25, 3.6);  //invent a B matrix
	vector<double> u_n(4, 0.2); // invent a u_n
	vector<double> steps(10);
	steps[0] = 0.01;
	for (int i = 1; i < 10; i++) {
		steps[i] = steps[i - 1] + 0.01;
	}
	for (auto _ : s) {
		//code that gets timed
		double dt = steps[s.range(0)];
		one_timestep(dt, B, u_n);
	}
}

// Register the benchmark, test a range of dt, usual input arguments are integers, here indices to the dt vector
BENCHMARK(one_timestepBench)->DenseRange(0, 9);
//BENCHMARK(one_timestepBench)->DenseRange(0, 9)->Unit(benchmark::kMillisecond);


// replaces normal main function
BENCHMARK_MAIN();
