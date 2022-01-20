#include <bits/stdc++.h>
using namespace std;

namespace SLAE {
	vector<int>idxs;
	vector<bool>used;

	template<typename T>
	bool is_zero(const T& x, const T& epsilon) {
		return abs(x) <= epsilon;
	}

	template<typename T>
	int partial_pivot(const vector<vector<T>>& a, int n, int rank, int col) {
		int sel_row = rank;
		for (int i = rank; i < n; ++i) {
			if (abs(a[i][col]) > abs(a[sel_row][col])) {
				sel_row = i;
			}
		}
		return sel_row;
	}

	template<typename T>
	int non_zero_pivot(const vector<vector<T>>& a, int n, int rank, int col, T epsilon) {
		int sel_row = rank;
		for (int i = rank; i < n; ++i) {
			if (!is_zero(a[i][col], epsilon)) {
				sel_row = i;
				break;
			}
		}
		return sel_row;
	}

	template <typename T>
	pair<int, T> gauss_elimination(vector<vector<T>>& a, T epsilon, int pivot_end = -1, bool diagonalize = false) {
		const int n = (int)a.size();
		assert(n >= 1);
		const int m = (int)a[0].size();
		int rank = 0;
		T determinant = static_cast<T>(1);
		if (pivot_end == -1) {
			pivot_end = m;
		}
		used.resize(pivot_end);
		std::fill(used.begin(), used.end(), false);
		idxs.clear();

		for (int col = 0; col < pivot_end; ++col) {
			// Partial Pivoting
			// int sel_row = partial_pivot(a, n, rank, col);
			int sel_row = non_zero_pivot(a, n, rank, col, epsilon);
			if (is_zero(a[sel_row][col], epsilon)) {
				determinant = 0;
				continue;
			}
			if (sel_row != rank) {
				determinant *= -1;
				a[sel_row].swap(a[rank]);
			}
			determinant *= a[rank][col];

			// R_rank = R_rank / R_rank[rank][col]
			if (diagonalize && a[rank][col] != 1) {
				T coef = 1 / a[rank][col];
				for (int k = col; k < m; ++k) {
					a[rank][k] *= coef;
				}
			}

			int is = diagonalize ? 0 : rank + 1;
			// R_i = R_i - R_rank * (R_i[col] / R_rank[col]) (for all i != rank)
			for (int i = is; i < n; ++i) {
				if (is_zero(a[i][col], epsilon) || i == rank) {
					continue;
				}
				T freq = a[i][col] / a[rank][col];
				for (int k = col; k < m; k++) {
					a[i][k] -= freq * a[rank][k];
				}
			}
			rank++;
			used[col] = true;
			idxs.push_back(col);
			if (rank == n) {
				break;
			}
		}
		return make_pair(rank, determinant);
	}

	template <typename T>
	vector<vector<T>> inverse_matrix(const vector<vector<T>>& a, T epsilon) {
		assert(a.size() == a[0].size());
		const int n = (int) a.size();

		// augment input matrix with identity matrix
		vector<vector<T>>b(n, vector<T>(2 * n));
		for (int i = 0; i < n; ++i) {
			std::copy(a[i].cbegin(), a[i].cend(), b[i].begin());
			b[i][i + n] = 1;
		}
		auto [rank, determinant] = gauss_elimination(b, epsilon, n, true);
		if (rank != n) {
			return vector<vector<T>>();
		}
		vector<vector<T>>ans(n, n);
		for (int i = 0; i < n; ++i) {
			std::copy(b[i].cbegin() + n, b[i].cend(), ans[i].begin());
		}
		return ans;
	}

	// solve ax = b
	template <typename T>
	vector<T> solve_gauss_elimination(const vector<vector<T>>& a, const vector<T>& b, T epsilon) {
		assert(a.size() == b.size());
		const int n = (int)a.size();
		const int m = (int)a[0].size();
		assert(n == m);

		// augment matrix [a | b]
		vector<vector<T>> aug(n, vector<T>(m + 1));
		for (int i = 0; i < n; ++i) {
			std::copy(a[i].cbegin(), a[i].cend(), aug[i].begin());
			aug[i][m] = b[i];
		}
		// solve for m variables in case (n != m)
		auto [rank, determinant] = gauss_elimination(aug, epsilon, m);
		vector<T>solution;
		for (int i = rank; i < n; ++i) {
			if (!is_zero(aug[i][m], epsilon)) {
				return solution;
			}
		}
		/*
		'aug' has now upper triangular matrix in it's rank size minor.
		Find solution using back substitution.
		*/
		solution.resize(m);
		for (int i = rank - 1; i >= 0; --i) {
			int col = idxs[i];
			solution[col] = aug[i][m];
			for (int j = col + 1; j < m; ++j) {
				solution[col] -= aug[i][j] * solution[j];
			}
			solution[col] /= aug[i][col];
		}
		return solution;
	}

	// solve ax = b
	template <typename T>
	vector<T> solve_gauss_jordan(const vector<vector<T>>& a, const vector<T>& b, T epsilon) {
		assert(a.size() == b.size());
		const int n = (int)a.size();
		const int m = (int)a[0].size();
		assert(n == m);

		// augment matrix [a | b]
		vector<vector<T>> aug(n, vector<T>(m + 1));
		for (int i = 0; i < n; ++i) {
			std::copy(a[i].cbegin(), a[i].cend(), aug[i].begin());
			aug[i][m] = b[i];
		}
		// solve for m variables in case (n != m)
		auto [rank, determinant] = gauss_elimination(aug, epsilon, m);

		vector<T>solution;
		for (int i = rank; i < n; ++i) {
			if (!is_zero(aug[i][m], epsilon)) {
				return solution;
			}
		}
		for (int i = 0; i < n; ++i) {
			std::reverse(aug[i].begin(), aug[i].begin() + m);
		}
		for (int i = 0; i < n - i - 1; ++i) {
			aug[i].swap(aug[n - i - 1]);
		}
		// Note: must use non_zero_pivot here
		gauss_elimination(aug, epsilon, m);
		/*
		'aug' has now diagonal matrix in it's rank size minor.
		Find solution using direct assignment.
		*/
		solution.resize(m);
		for (int i = 0; i < rank; ++i) {
			int col = idxs[i];
			solution[m - col - 1] = aug[i][m] / aug[i][col];
		}
		return solution;
	}

	// solve ax = b
	template <typename T>
	vector<T> solve_crout(const vector<vector<T>>& a, const vector<T>& b, T epsilon) {
		assert(a.size() == b.size());
		const int n = (int)a.size();
		const int m = (int)a[0].size();
		assert(n == m);

		// augment matrix [a | b]
		vector<vector<T>> aug(n, vector<T>(m + 1));
		for (int i = 0; i < n; ++i) {
			std::copy(a[i].cbegin(), a[i].cend(), aug[i].begin());
			aug[i][m] = b[i];
		}
		// LU = A
		// UX = Y
		// LY = B
		// auxiliary matrix [LU | Y]
		// [row >= col] if for L
		// [row < col] is for U and Y
		vector<vector<T>> aux(n, vector<T>(m + 1));
		int row, col;
		for (int i = 0; i < n; ++i) {
			// find ith column of aux
			// [row >= col]
			col = i;
			for (row = i; row < n; ++row) {
				T inner_prod = 0;
				for (int j = 0; j < col; ++j) inner_prod += (aux[row][j] * aux[j][col]);
				aux[row][col] = aug[row][col] - inner_prod;
			}
			assert(!is_zero(aux[i][i], epsilon));
			// find ith row of aux
			// [row < col]
			row = i;
			for (col = i + 1; col <= m; ++col) {
				T inner_prod = 0;
				for (int j = 0; j < row; ++j) inner_prod += (aux[row][j] * aux[j][col]);
				aux[row][col] = (aug[row][col] - inner_prod) / aux[row][row];
			}
		}
		vector<T>solution(n);
		for (row = n - 1; row >= 0; --row) {
			solution[row] = aux[row][m];
			for (col = row + 1; col < m; ++col) {
				solution[row] -= aux[row][col] * solution[col];
			}
		}
		return solution;
	}
}
int main() {
	std::ios::sync_with_stdio(false); std::cin.tie(nullptr);
	vector<vector<double>>a = {
		{1, 2, 1},
		{2 , -3, -1},
		{3 , 1 , 2}
	};
	vector<double>b = {4, -3, 3};
	auto ans = SLAE::solve_crout(a, b, 0.0);
	return 0;
}