#include <bits/stdc++.h>
using namespace std;
namespace VectorOperations {
	template<typename T>
	bool is_valid(const vector<T>& a, int l, int r) {
		return (l >= 0 && r >= l && r < (int) a.size());
	}
	// a[ptr1:ptr2 + len] += b[ptr2:ptr2 + len]
	template<typename T>
	void add(vector<T>& a, vector<T> b, int len, int ptr1 = 0, int ptr2 = 0) {
		assert(is_valid(a, ptr1, ptr1 + len - 1));
		assert(is_valid(b, ptr2, ptr2 + len - 1));
		int idx = 0;
		while (idx < len) {
			a[ptr1 + idx] += b[ptr2 + idx];
			idx++;
		}
	}

	template<typename T>
	vector<T> mul(const vector<T>& a, const T& coef, int l = 0, int r = -1) {
		if (r == -1) r = (int) a.size() - 1;
		assert(is_valid(a, l, r));
		vector<T>b = a;
		while (l <= r) {
			b[l++] *= coef;
		}
		return b;
	}
} // namespace VectorOperations