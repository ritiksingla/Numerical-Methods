#include <bits/stdc++.h>
using namespace std;
using ld = long double;
class Segtree {
  public:
	struct Node {
		ld sum;
		ld mul;
		Node() {
			sum = 0;
			mul = 1;
		}
		void apply(int l, int r, ld v) {
			sum += (r - l + 1) * v;
			mul = sum;
		}
	}; // Node

	int n;
	vector<Node>tree;

	static Node unite(const Node& a, const Node& b) {
		Node res;
		res.sum = a.sum + b.sum;
		res.mul = a.mul * b.mul;
		return res;
	}
	void init(int _n) {
		n = _n;
		tree.resize(2 * n - 1);
		build(0, 0, n - 1);
	}

	Segtree(int n) {
		init(n);
	};

	template<typename T>
	void modify(int l, int r, const T& v) {
		modify(0, 0, n - 1, l, r, v);
	}
	ld get(int l, int r) {
		return get(0, 0, n - 1, l, r);
	}
  private:

	void pull(int x, int y) {
		tree[x] = unite(tree[x + 1], tree[y]);
	}

	template<typename T>
	void build(int x, int l, int r, const vector<T>& v) {
		if (l == r) {
			tree[x].apply(l, r, v[l]);
			return;
		}
		int m = (l + r) >> 1;
		int y = x + ((m - l + 1) << 1);
		build(x + 1, l, m, v);
		build(y, m + 1, r, v);
		pull(x, y);
	}
	void build(int x, int l, int r) {
		if (l == r) {
			tree[x].apply(l, r, 0.0);
			return;
		}
		int m = (l + r) >> 1;
		int y = x + ((m - l + 1) << 1);
		build(x + 1, l, m);
		build(y, m + 1, r);
		pull(x, y);
	}
	ld get(int x, int l, int r, int ql, int qr) {
		if (ql <= l && r <= qr) {
			return tree[x].mul;
		}
		int m = (l + r) >> 1;
		int y = x + ((m - l + 1) << 1);
		if (qr <= m) {
			return get(x + 1, l, m, ql, qr);
		} else if (ql > m) {
			return get(y, m + 1, r, ql, qr);
		} else {
			return (get(x + 1, l, m, ql, qr) * get(y, m + 1, r, ql, qr));
		}
	}
	template<typename T>
	void modify(int x, int l, int r, int ql, int qr, const T& v) {
		if (l == r) {
			tree[x].apply(l, r, v);
			return;
		}
		int m = (l + r) >> 1;
		int y = x + ((m - l + 1) << 1);
		if (ql <= m) {
			modify(x + 1, l, m, ql, qr, v);
		}
		if (qr > m) {
			modify(y, m + 1, r, ql, qr, v);
		}
		pull(x, y);
	}
}; // Segtree