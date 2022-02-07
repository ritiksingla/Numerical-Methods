#include <bits/stdc++.h>
#include "segtree.h"
#include "VectorOperations.h"
#include "io.h"
using namespace std;

using ld = long double;

namespace InterpolationUnequal {
	vvt<ld> get_divided_diff(const vt<ld>& x, const vt<ld>& y, int n) {
		vvt<ld>divided_diff(n);
		divided_diff[0] = y;
		for (int i = 1; i < n; ++i) {
			divided_diff[i].resize(n - i);
			for (int j = 0; j < n - i; ++j) {
				divided_diff[i][j] = (divided_diff[i - 1][j + 1] - divided_diff[i - 1][j]);
				divided_diff[i][j] /= (x[j + i] - x[j]);
			}
		}
		return divided_diff;
	}
	ld lagrange(const vt<ld>& x, const vt<ld>& y, ld qx) {
		assert(x.size() == y.size());
		int n = (int) x.size();
		assert(n >= 1);
		Segtree D(n), N(n);
		for (int i = 0; i < n; ++i) {
			D.modify(i, i, -x[i]);
			N.modify(i, i, qx - x[i]);
		}
		ld res = 0;
		const ld eps = 1e-4;
		for (int i = 0; i < n; ++i) {
			D.modify(0, n - 1, x[i]);
			ld den = 1;
			ld num = 1;
			if (i > 0) {
				den *= D.get(0, i - 1);
				num *= N.get(0, i - 1);
			}
			if (i + 1 < n) {
				den *= D.get(i + 1, n - 1);
				num *= N.get(i + 1, n - 1);
			}
			assert(abs(den) > eps);
			res += y[i] * (num / den);
			D.modify(0, n - 1, -x[i]);
		}
		return res;
	}
	// return interpolated polynomial
	vt<ld> newton_divided_difference(const vt<ld>& x, const vt<ld>& y) {
		assert(x.size() == y.size());
		int n = (int) x.size();
		assert(n >= 1);
		auto divided_diff = get_divided_diff(x, y, n);
		vt<ld>res(n);
		res[0] = divided_diff[0][0];
		vt<ld>poly(n);
		poly[n - 1] = 1.0;
		for (int i = 1; i < n; ++i) {
			for (int j = n - i - 1; j + 1 < n; ++j) {
				poly[j] -= poly[j + 1] * x[i - 1];
			}
			VectorOperations::add(res, VectorOperations::mul(poly, divided_diff[i][0], n - i - 1), i + 1, 0, n - i - 1);
		}
		return res;
	}
	// find f(x)
	ld eval(const vt<ld>& f, ld x) {
		ld res = 0;
		ld power = 1;
		int n = (int) f.size();
		for (int i = 0; i < n; ++i) {
			res += f[i] * power;
			if (i + 1 == n) break;
			power *= x;
		}
		return res;
	}
} // namespace InterpolationUnequal
using namespace InterpolationUnequal;
int main() {
	std::ios::sync_with_stdio(false); std::cin.tie(nullptr);
	int n; cin >> n;
	vt<ld>x(n), y(n);
	cin >> x >> y;
	cout << newton_divided_difference(x, y);
	return 0;
}