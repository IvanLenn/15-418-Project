#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int, int> pi;
typedef pair<ll, ll> pl;
typedef vector<int> vi;
typedef vector<ll> vl;
#define eb emplace_back
mt19937_64 rnd(chrono::steady_clock::now().time_since_epoch().count());

int rand(int N, int M) {
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::uniform_int_distribution<> distr(N, M); // Define the range
    return distr(gen); // Generate a random number within the range [N, M]
}

struct FlowEdge {
    int v, u;
    long long cap, flow = 0;
    FlowEdge(int v, int u, long long cap) : v(v), u(u), cap(cap) {}
};

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    // freopen("path", "w", stdout);
    const int N = 300;
    const int M = 4000;
    int s = rand(0, N - 1);
    int t = rand(0, N - 1);
    while (t == s) t = rand(0, N - 1);
    cout << N << ' ' << M << ' ' << s << ' ' << t << '\n';
    vi a(N);
    for (int i = 0; i < M; i++) {
        int u = rand(0, N - 1);
        int v = rand(0, N - 1);
        a[u]++; a[v]++;
        int cap = rand(0, 1e9 + 7);
        cout << u << ' ' << v << ' ' << cap << '\n';
    }

    int mn = M, sum = 0;
    for (auto i : a) {
        sum += i;
        mn = min(mn, i);
    }
    cout << sum << ' ' << mn;
}