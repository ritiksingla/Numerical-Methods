#include <bits/stdc++.h>
using namespace std;

template<typename T, class F = function<T(const T&)>>
class TranscendalEquations {
  public:
	TranscendalEquations(const F& f_): f(f_) {
		epsilon = 0.001;
	}
	TranscendalEquations(const F& f_, T epsilon_): f(f_), epsilon(epsilon_) {}
	T bisection(T l, T r, T tolerance, int threshold = 15, bool verbose = true) {
		assert(threshold >= 1);
		T fl = f(l);
		T fr = f(r);
		assert(fl * fr < 0);
		if (abs(fr - fl) < tolerance) {
			printf("endpoints are within tolerance error\n");
			return l;
		}
		T eval = fl;
		for (int i = 1; i <= threshold; ++i) {
			T m = (l + r) / 2;
			eval = f(m);
			if (fl < 0 && eval < 0) {
				l = m;
				fl = eval;
			} else {
				r = m;
				fr = eval;
			}
			if (verbose) {
				debug(i, m, eval);
			}
			if (is_approx_zero(eval)) {
				return m;
			}
			if (abs(fr - fl) < tolerance) {
				printf("endpoints are within tolerance error\n");
				return m;
			}
		}
		return l;
	}

	T regula_falsi(T x0, T x1, T tolerance, int threshold = 10, bool verbose = true) {
		T f0 = f(x0);
		T f1 = f(x1);
		assert(f0 * f1 < 0);
		if (abs(f1 - f0) < tolerance) {
			printf("endpoints are within tolerance error\n");
			return x0;
		}
		if (f1 < 0) {
			std::swap(x0, x1);
			std::swap(f0, f1);
		}
		// fix x1
		for (int i = 1; i <= threshold; ++i) {
			x0 = x1 - ((x1 - x0) / (f1 - f0)) * f1;
			f0 = f(x0);
			if (verbose) {
				debug(i, x0, f0);
			}
			if (is_approx_zero(f0)) {
				return x0;
			}
			if (abs(f1 - f0) < tolerance) {
				printf("endpoints are within tolerance error\n");
				return x0;
			}
		}
		return x0;
	}

	T secant(T x0, T x1, T tolerance, int threshold = 10, bool verbose = true) {
		T f0 = f(x0);
		T f1 = f(x1);
		if (abs(f1 - f0) < tolerance) {
			printf("endpoints are within tolerance error\n");
			return x0;
		}
		if (f1 > f0) {
			std::swap(x0, x1);
			std::swap(f0, f1);
		}
		// fk < fk-1 < fk-2
		for (int i = 1; i <= threshold; ++i) {
			T x_i = x1 - ((x1 - x0) / (f1 - f0)) * f1;
			T f_i = f(x_i);

			// update
			x0 = x1;
			f0 = f1;

			x1 = x_i;
			f1 = f_i;

			if (verbose) {
				debug(i, x1, f1);
			}
			if (is_approx_zero(f_i)) {
				return x1;
			}
			if (abs(f1 - f0) < tolerance) {
				printf("endpoints are within tolerance error\n");
				return x1;
			}
		}
		return x1;
	}

	T newton_raphson(T x0, F df, int threshold = 10, bool verbose = true) {
		T f0 = f(x0);
		T df0 = df(x0);
		for (int i = 1; i <= threshold; ++i) {
			x0 = x0 - f0 / df0;
			f0 = f(x0);
			df0 = df(x0);

			if (verbose) {
				debug(i, x0, f0);
			}
			if (is_approx_zero(f0)) {
				return x0;
			}
		}
		return x0;
	}
  private:
	F f; // equation to be solved
	T epsilon; // negative power of 10 to be considered as small enough

	bool is_approx_zero(T x) {
		if (abs(x) <= epsilon) {
			printf("evaluated function reached epsilon\n");
			return true;
		}
		return false;
	}
	void debug(int iteration, T x, T eval) {
		printf(
		    "approx. root after %d iterations is %.6f and evaluates to %.6f\n",
		    iteration,
		    x,
		    eval
		);
	}
};
int main() {
	auto f = [](double x) {return (x - 1) * sin(x) - x - 1;};
	auto df = [](double x) {return (x - 1) * cos(x) + sin(x) - 1;};
	TranscendalEquations<double>T(f, 0.001);
	double approx_root = T.newton_raphson(0, df);
	return 0;
}