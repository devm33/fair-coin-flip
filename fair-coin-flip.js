
/* testing */
module.exports = square_root_modp;
console.log('sqr(9) mod 17 = '+square_root_modp(9,17));

function square_root_modpq(x, p, q) {
    /* Finds the square roots of x in Z_n, where n = p * q
     * Returns an array of the two "positive" (< n/2) square roots
     * assumes p <> q */
    var rp, rq;
    /* Compute square roots in the fields of the two factors */
    rp = square_root_modp(x, p);
    rq = square_root_modp(x, q);
    /* Use Chinese Remainder Theorem (and Bezout's Identity) map roots
     * to Z_n */
    //TODO
    return [-1,1];
}

function square_root_modp(x, p) {
    /* Returns a square root of x in mod p, p an odd prime, could be any
     * either root of the forms: a, p-a
     * (assumes x is indeed a square)
     * Uses Cipolla's algorithm:
     * http://people.math.gatech.edu/~mbaker/pdf/cipolla2011.pdf */
    var t, w, i, p2, a, b, a0, b0;
    do {
        t = random_integer_between(1, p);
        w = t * t - x;
    } while(legendre(w, p) == 1);
    p12 = (p+1)/2;
    a = t; /* use a,b to track squaring: (t + sq(w)) ^ k = a + b sq(w) */ 
    b = 1; /* as such: t + sq(w)) ^ k + 1 = at + bw + (a + bt)sq(w) */
    for(i = 1; i < p12; i++) {
        a0 = a * t + b * w;
        b0 = a + b * t;
        a = a0 % p;
        b = b0 % p;
    }
    while(a < 0) {
        a = (a + p) % p;
    }
    return a;
}

function legendre(a, p) {
    /* Returns the Legendre symbol of a/p where p is an odd prime
     * Will be of the set {1, 0, p-1} */
    var l = modular_power(a, (p-1)/2, p);
    while(l != 1 && l != 0 &&  l != p-1) {
        l = (l + p) % p;
    }
    return l;
};

/* TODO get bezouts coefs:
    private static Pair extended_euclidean(int a, int b) {
        Pair s = new Pair(0, 1);
        Pair t = new Pair(1, 0);
        Pair r = new Pair(max(a,b), min(a,b));
        int q;
        while (r.b != 0) {
            q = r.a / r.b;
            r = new Pair(r.b, r.a - q * r.b);
            s = new Pair(s.b, s.a - q * s.b);
            t = new Pair(t.b, t.a - q * t.b);
        }
        return new Pair(s.a, t.a);
    }
*/

function generate_prime_of_length(n) {
    /* Returns a random prime' number that has n digits
     * where prime' means probably prime, uses Miller-Rabin
     */
    var p;
    while(true) {
         p = random_integer_of_length(n);
         if(!has_small_prime_factor(p) && miller_rabin_test(p, 100)) {
             return p;
         }
    }
}

function random_integer_of_length(n) {
    /* Returns a random integer with n digits */
    return random_integer_between(Math.pow(10,n-1), Math.pow(10,n));
    /* note: this doesnt support n > 16 and may not be the fastest way
     * to do this */
}

function random_integer_between(min, max) {
    /* Returns a random integer in [min, max) */
    return Math.floor(Math.random() * (max - min) + min);
}

function miller_rabin_test(n, k) {
    /* Returns false if n is composite and true if n is probably prime
     * with probably 1/4^k
     */
    var dk, d, k, a, x, i, j;
    dk = separate_two_factor(n); d = dk.d; k = dk.k;
    witness: for(i = 0; i < k; i++) {
        a = random_integer_between(2, n-2);
        x = modular_power(a, d, n);
        if(x == 1 || x == n - 1) {
            continue witness;
        }
        for(j = 0; j < k - 1; j++) {
            x = (x * x) % n;
            if(x == 1) {
                return false;
            }
            if(x == n -1) {
                continue witness;
            }
        }
        return false;
    }
    return true;
}

function separate_two_factor(n) {
    /* Returns an object {'d':d,'k':k} where n = d * 2^k, for maximal k
     * helper to Miller-Rabin test
     */
    var d = n, k = 0;
    while(d % 2 == 0) {
        d = d / 2;
        k++;
    }
    return {'d': d, 'k': k};
}

function modular_power(b, k, n) {
    /* Returns b raised to the k in mod n */
    var ans = 1;
    for(var i = 0; i < k; i++) {
        ans = ans * b;
        if(ans > n) {
            ans = ans % n;
        }
    }
    return ans;
}

var has_small_prime_factor = (function(){
    /* anonymous function to trap list in scope so is only made once */
    var small_prime = find_all_primes_less_than(1e4);
    var len = small_prime.length;
    return function(n) {
        /* Returns true if n has a factor in list of small primes */
        for(var i = 0; i < len; i++) {
            if(n % small_prime[i] == 0) {
                if(n == small_prime[i]) { /* quick check in case n is a small prime -- could handle small n better, but not worth building out right now */
                    return false;
                }
                return true;
            }
        }
        return false;
    }
})();

function find_all_primes_less_than(n) {
    /* Uses Sieve of Eratosthenes to return an array of primes < n */
    var i, j, primes = [], candidate = {};
    /* Initialize map to true for all 1 < i < max */
    for(i = 2; i < n; i++) {
        candidate[i] = true;
    }
    /* Sieve out numbers that are multiples of primes */
    for(i = 2; i < n; i++) {
        if(candidate[i]) {
            for(j = i+i; j < n; j += i) {
                candidate[j] = false;
            }
        }
    }
    /* Filter through candidates to return array of primes */
    for(i = 2; i < n; i++) {
        if(candidate[i]) {
            primes.push(i);
        }
    }
    return primes;
}
