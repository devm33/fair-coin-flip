
function square_root_modpq(x, p, q) {
    /* Finds the square roots of x in Z_n, where n = p * q
     * Returns an array of the two "positive" (< n/2) square roots
     * assumes p <> q */
    /* Compute square roots in the fields of the two factors */
    rp = square_root_modp(x, p);
    rq = square_root_modp(x, q);
    //TODO CRT to find answer, find bezouts: http://en.wikipedia.org/wiki/Chinese_remainder_theorem#Case_of_two_equations
    return [-1,1];
}

function square_root_modp(x, p) {
    /* Returns the positive (< p/2) square of x in mod p, p an odd prime
     * (assumes x is indeed a square)
     * Uses Cipolla's algorithm: http://en.wikipedia.org/wiki/Cipolla%27s_algorithm */
    var a, a2n;
    do {
        a = random_integer_between(1, p);
        a2n = a * a - n;
    } while(legendre(a2n, p) != p-1);
    //TODO something about the square root of a2n - even though that makes no sense, im tired
    return 1;
}

function legendre(a, p) {
    /* Returns the Legendre symbol of a/p where p is an odd prime */
    return modular_power(a, (p-1)/2, p);
};

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
    /* note: there could be some performance gains from doing this differently */
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
