#include <iostream>
#include <vector>
#include <string>
using namespace std;

struct Fenwick {
    vector<long long> tree;
    int n;
    Fenwick(int sz) : n(sz), tree(sz, 0) {}
    void add(int i, long long x) {
        for (; i < n; i |= (i + 1)) tree[i] += x;
    }
    long long sum(int r) {
        long long res = 0;
        for (--r; r >= 0; r = (r & (r + 1)) - 1) res += tree[r];
        return res;
    }
    long long range(int l, int r) { 
        return sum(r) - sum(l);
    }
};

int main() {
    int n;
    cin >> n;
    vector<long long> a(n);
    for (int i = 0; i < n; ++i) cin >> a[i];
    Fenwick b(n);
    for (int i = 0; i < n; ++i) b.add(i, a[i]);
    int q;
    cin >> q;
    while (q--) {
        string cmd;
        cin >> cmd;
        if (cmd == "Add") {
            int i;
            long long x;
            cin >> i >> x;
            b.add(i, x);
        }
        else if (cmd == "FindSum") {
            int l, r;
            cin >> l >> r;
            cout << b.range(l, r) << "\n";
        }
    }
    return 0;
}

