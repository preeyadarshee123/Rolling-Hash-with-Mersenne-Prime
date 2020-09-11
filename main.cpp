//  Rolling Hash

constexpr uint64_t mod = (1ull<<61) - 1;
int gen_base(const int before, const int after) {
    seed_seq seq{
        (uint64_t) chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count(),
        (uint64_t) __builtin_ia32_rdtsc(),
        (uint64_t) (uintptr_t) make_unique<char>().get()
    };
    mt19937 rng(seq);
    int  base = uniform_int_distribution<int>(before+1, after)(rng);
    return base % 2 == 0 ? base-1 : base;
}
uint64_t mul(uint64_t a, uint64_t b){
    uint64_t l1 = (uint32_t)a, h1 = a>>32, l2 = (uint32_t)b, h2 = b>>32;
    uint64_t l = l1*l2, m = l1*h2 + l2*h1, h = h1*h2;
    uint64_t ret = (l&mod) + (l>>61) + (h << 3) + (m >> 29) + (m << 35 >> 3) + 1;
    ret = (ret & mod) + (ret>>61);
    ret = (ret & mod) + (ret>>61);
    return ret-1;
}
uint64_t sub(uint64_t a, uint64_t b) {
    return (a -= b) >= mod ? a + mod : a;
}
uint64_t add(uint64_t a, uint64_t b) {
    return (a += b) < mod ? a : a - mod;
}
uint64_t modInverse(uint64_t x,uint64_t n){
    uint64_t res = 1;
    while(n){
        if(n%2){
            res = mul(res , x);
        }
        x = mul(x , x);
        n = n/2;
    }
    return res;
}
class PolyHash{
    public :
    vector<uint64_t> pre,pow,invPow,suf;
    uint64_t base;
    uint64_t invBase;
    void calc(string s){
        int n = (int)s.length();
        pre.resize(n+1,0);
        pow.resize(n+1,1);
        invPow.resize(n+1,1);
        suf.resize(n+2,0);
        for(int i = 0;i < n;i++){
            pow[i+1] = mul(pow[i],base);
            invPow[i+1] = mul(invPow[i],invBase);
            pre[i+1] = add(mul(s[i],pow[i]), pre[i]);
            suf[n-i] = add(mul(s[n-i-1],pow[i]),suf[n-i+1]);
        }
    }
    // Compares two string starting from index id1 and id2 of length len
    int cmp(int id1,int id2,int len){
        int lo = 0,hi = len+1;
        while(hi - lo > 1){
            int mid = (hi + lo)/2;
            ull p1 = mul(sub(pre[id1+mid],pre[id1]) , invPow[id1]);
            ull p2 = mul(sub(pre[id2+mid],pre[id2]) , invPow[id2]);
            if(p1 == p2){
                lo = mid;
            }
            else{
                hi = mid;
            }
        }
        return lo;
    }
};

