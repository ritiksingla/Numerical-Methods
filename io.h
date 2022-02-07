#include <bits/stdc++.h>
using namespace std;
#define trav(it, x) for(auto& it : x)

template<class T> using vt = vector<T>;
template<class T> using vvt = vt<vt<T>>;
template<class T> using pr = pair<T, T>;
namespace IO {
	template<class T> istream& operator >>(istream& is, vt<T>& v) {
		trav(x, v) is >> x; return is;
	}
	template<class T> ostream& operator <<(ostream& os, const vt<T>& v) {
		bool first = true;
		for (const T& x : v) {if (!first) os << ' '; first = false; os << x; } return os;
	}
	template<class T> istream& operator >> (istream& is, vvt<T>& vv) {
		trav(v, vv) is >> v; return is;
	}
	template<class T> ostream& operator << (ostream& os, const vvt<T>& vv) {
		bool first = true;
		for (const vt<T>& v : vv) {if (!first) os << '\n'; first = false; os << v; } return os;
	}
	template<class T, class U> istream& operator >>(istream& is, pair<T, U>& p) {
		is >> p.first >> p.second; return is;
	}
	template<class T, class U> ostream& operator <<(ostream& os, const pair<T, U>& p) {
		os << p.first << ' ' << p.second; return os;
	}
} // namespace IO
using namespace IO;