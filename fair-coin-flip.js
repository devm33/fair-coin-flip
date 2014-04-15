
function square_root_modpq(x, p, q) {
    /* Finds the square roots of x in Z_n, where n = p * q
     * Returns an array of the two "positive" (< n/2) square roots
     * assumes p <> q */
    var roots, inverses, partials, n = p * q;
    /* Compute square roots in the fields of the two factors */
    roots = [square_root_modp(x, p), square_root_modp(x, q)];
    /* Use Chinese Remainder Theorem (and Bezout's Identity) map roots
     * to Z_n */
    inverses = bezouts_coefs(p, q);
    partials = [inverses[0] * p * roots[1], inverses[1] * q * roots[0]];
    return [
        (partials[0] + partials[1]) % n,
        -(partials[0] + partials[1]) % n,
        (partials[0] - partials[1]) % n,
        (partials[1] - partials[0]) % n
    ].map(function(v) {
        return v > 0 ? v : v + n; // return all roots, but greater than zero
    });
}

function square_root_modp(a, p) { //console.log('root '+a+' mod '+p);
    /* Returns a square root of a in mod p, p an odd prime, could be any
     * either root of the forms: x, p-x
     * (assumes a is indeed a quadratic residue)
     * First checks easy case of p=3(4) then uses Cipolla's algorithm:
     * http://people.math.gatech.edu/~mbaker/pdf/cipolla2011.pdf */
    var t, w, i, p2, x, y, x0, y0;
    if(p % 4 == 3) {
        return modular_power(a, Math.floor(p / 4) + 1, p);
    } 
    do {
        t = random_integer_between(1, p);
        w = t * t - a;
    } while(jacobi(w, p) == 1);
    p12 = (p+1)/2;
    x = t; /* use x,y to track squaring: (t + root(w)) ^ k = x + y root(w) */ 
    y = 1; /* as such: t + root(w)) ^ k + 1 = xt + yw + (x + yt) root(w) */
    for(i = 1; i < p12; i++) {
        x0 = x * t + y * w;
        y0 = x + y * t;
        x = x0 % p;
        y = y0 % p;
    }
    if(x < 0) {
        x = x + p;
    }
    return x;
}

function legendre(a, p) {
    /* Returns the Legendre symbol of a/p where p is an odd prime
     * Will be of the set {1, 0, p-1} */
    var l = modular_power(a, (p-1)/2, p);
    if(l < 0) {
        l = l + p;
    }
    return l;
}

function jacobi(A, B) { //console.log('jacobi '+A+' '+B);
    /* Returns the Jacobi symbol of a by b for b odd */
    var a = A, b = B, c, s, sign = 1, t;
    while(b > 1) { //console.log('> jacobi '+a+' '+b);
        if(a >= b) {
            a = a % b;
        }
        else if(a < 0) {
            a = (b - a) % b; // -a, negated and mod simplified
            a = b - a; // --a = a, double negated
        }
        if(a === 0) {
            return 0;
        }
        if((a & 1) == 1) { //a odd
            if((a & 3) == 3 && (b & 3) == 3) { // both a, b are 3 mod 4 - cor of quadratic reciprocity
                sign = sign * -1;
            }
            t = a; a = b; b = t; // swap a and b
        }
        else { // a even
            c = a, s = 0; // factor out 2 from a = c * 2 ^ s
            while((c & 1) === 0) {
                c = c >> 1;
                s++;
            }
            if((s & 1) == 1 && ((b & 7) == 5 || (b & 7) == 3)) {
                // s (power of factor of 2) must be odd and b must be 5 or 3 mod 8 - cor of Gauss lemma
                sign = sign * -1;
            }
            a = c;
        }
    }
    return sign;
    /* In an effort to make this faster, much of this was stolen
     * gratuitously from: http://yacas.sourceforge.net/ */ 
}

function bezouts_coefs(a, b) {
    var s,t,r,q,flag;
    s = [0, 1];
    t = [1, 0];
    r = [Math.max(a,b), Math.min(a,b)];
    flag = (a > b);
    while(r[1] !== 0) {
            q = Math.floor(r[0] / r[1]);
            r = [r[1], r[0] - q * r[1]];
            s = [s[1], s[0] - q * s[1]];
            t = [t[1], t[0] - q * t[1]];
    }
    if(flag) {
        return [t[0], s[0]];
    }
    return [s[0], t[0]];
}

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
    var ds, d, s, a, x, i, j;
    ds = separate_two_factor(n-1); d = ds[0]; s = ds[1];
    witness: for(i = 0; i < k; i++) {
        a = random_integer_between(2, n-2);
        x = modular_power(a, d, n);
        if(x == 1 || x == n - 1) {
            continue witness;
        }
        for(j = 0; j < s - 1; j++) {
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
    /* Returns an array [d, k] where n = d * 2^k, for maximal k
     * helper to Miller-Rabin test
     */
    var d = n, k = 0;
    while((d & 1) === 0) {
        d = d >> 1;
        k++;
    }
    return [d, k];
}

function modular_power(b, k, n) {
    /* Returns b raised to the k in mod n */
    var p = b, j = k, r = 1;
    while(j > 0) {
        if((j & 1) == 1) { // j odd
            r = (r * p) % n;
        }
        p = (p * p) % n;
        j = j >> 1;
    }
    return r;
    /* In an effort to make this faster, this was stolen gratuitously
     * from: http://yacas.sourceforge.net/ */ 
}

function modular_power_slower(b, k, n) {
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
    var small_prime = find_all_primes_less_than(1e6);
    var len = small_prime.length;
    return function(n) {
        /* Returns true if n has a factor in list of small primes */
        for(var i = 0; i < len; i++) {
            if(n % small_prime[i] === 0) {
                if(n == small_prime[i]) { /* quick check in case n is a small prime */
                    return false; /* could handle small n better, but small dont need to be handled better */
                }
                return true;
            }
        }
        return false;
    };
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
