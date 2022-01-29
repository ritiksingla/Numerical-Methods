#include <bits/stdc++.h>
using namespace std;
namespace Interpolation {
	vector<vector<double>> get_difference_table(const vector<double>& y, int n) {
		vector<vector<double>>diff_table(n);
		diff_table[0] = y;
		for (int i = 1; i < n; ++i) {
			diff_table[i].resize(n - i);
			for (int j = 0; j < n - i; ++j) {
				diff_table[i][j] = (diff_table[i - 1][j + 1] - diff_table[i - 1][j]);
			}
		}
		return diff_table;
	}
	// type = 1 for forward, type = 2 for backward
	double newton_interpolation(const vector<double>& x, const vector<double>& y, double qx, int type) {
		assert(x.size() == y.size() && x.size() > 1);
		assert(type == 1 || type == 2);
		int n = (int) x.size();
		double diff = x[1] - x[0];
		double q = qx;
		if (type == 1)q = (q - x[0]) / diff;
		else q = (q - x[n - 1]) / diff;
		auto diff_table = get_difference_table(y, n);
		double res = 0;
		double fact = 1;
		double num = 1;
		for (int i = 0; i < n; ++i) {
			if (i != 0) {
				fact *= i;
				if (type == 1) num *= (q - i + 1);
				else num *= (q + i - 1);
			}
			res += (num * (type == 1 ? diff_table[i][0] : diff_table[i].back())) / fact;
		}
		return res;
	}

	// central differences
	// type = 1 for forward, type = 2 for backward, x0 = x[n/2], nx = x0 + offset * diff
	double gauss_interpolation(const vector<double>& x, const vector<double>& y, double qx, int type) {
		assert(x.size() == y.size() && x.size() > 1);
		assert(type == 1 || type == 2);
		int n = (int) x.size();
		double diff = x[1] - x[0];
		double x0 = x[n / 2];
		int j = n / 2;
		if (n % 2 == 0 && type == 1) { // for even case type-1 x0 can be x[n/2-1] also
			// change according to question need
			x0 = x[n / 2 - 1];
			j = n / 2 - 1;
		}
		double offset = (qx - x0) / diff;
		auto diff_table = get_difference_table(y, n);
		double res = 0;
		double fact = 1;
		double num = 1;
		for (int i = 0; i < n; ++i) {
			if (i != 0) {
				int delta = i / 2;
				fact *= i;
				if (i & 1) {
					if (type == 1)
						num *= (offset + delta);
					else {
						num *= (offset - delta);
						// decrease the ptr to difference table after every two steps in type-2
						j--;
					}
				} else {
					if (type == 1) {
						num *= (offset - delta);
						// decrease the ptr to difference table after every two steps in type-1
						j--;
					} else
						num *= (offset + delta);
				}
			}
			// if (type == 2 && (i % 2 == 1)) j--;
			res += (num * diff_table[i][j]) / fact;
		}
		return res;
	}

	double stirling(const vector<double>& x, const vector<double>& y, double qx) {
		assert(y.size() > 1);
		int n = (int) y.size();
		assert(n % 2 == 1);

		double diff = x[1] - x[0];
		double x0 = x[n / 2];
		double offset = (qx - x0) / diff;
		auto diff_table = get_difference_table(y, n);
		double res = 0;
		double fact = 1;
		int j = n / 2;
		double num1 = 1, num2 = offset;
		for (int i = 0; i < n; ++i) {
			double median = diff_table[i][j];
			if (i != 0) {
				fact *= i;
				int delta = i / 2;
				if (i & 1) {
					if (delta > 0)
						num2 *= (offset - delta) * (offset + delta);

					median = (median + diff_table[i][j - 1]) / 2;
					j--;
				} else {
					delta--;
					num1 *= (offset - delta) * (offset + delta);
				}
			}
			if (i & 1) {
				res += (num2 * median) / fact;
			} else {
				res += (num1 * median) / fact;
			}
		}
		return res;
	}

	double bessel(const vector<double>& x, const vector<double>& y, double qx) {
		assert(y.size() > 1);
		int n = (int) y.size();
		assert(n % 2 == 0);

		double diff = x[1] - x[0];
		double x0 = x[n / 2 - 1];
		double offset = (qx - x0) / diff;
		auto diff_table = get_difference_table(y, n);
		double res = (diff_table[0][n / 2] + diff_table[0][n / 2 - 1]) / 2;
		double fact = 1;
		double num1 = offset - 0.5, num2 = 1;
		for (int i = 1; i < n; ++i) {
			int j = (n - i) / 2;
			double median = diff_table[i][j];
			if (i & 1 ^ 1) {
				median = (median + diff_table[i][j - 1]) / 2;
			}
			fact *= i;
			int delta = i / 2;
			if (i & 1) {
				if (delta > 0)
					num1 *= (offset - delta) * (offset + delta - 1);
			} else {
				num2 *= (offset - delta) * (offset + delta - 1);
			}
			if (i & 1) {
				res += (num1 * median) / fact;
			} else {
				res += (num2 * median) / fact;
			}
		}
		return res;
	}
} // namespace Interpolation
using namespace Interpolation;
int main() {
	int n; cin >> n;
	vector<double>x(n), y(n);
	for (auto& it : x)cin >> it;
	for (auto& it : y)cin >> it;
	double nx; cin >> nx;
	cout << fixed << setprecision(5) << stirling(x, y, nx);
	return 0;
}