#include <string.h>
#include "gghlite-internals.h"
#include "gghlite.h"
#include "oz/oz.h"
#include "oz/flint-addons.h"

void
gghlite_sk_set_D_g(gghlite_sk_t self)
{
    assert(self);
    assert(self->params);
    assert(mpfr_cmp_d(self->params->sigma_p, 0) > 0);

    const oz_flag_t flags = (self->params->flags & GGHLITE_FLAGS_QUIET) ? 0 : OZ_VERBOSE;
    self->t_D_g = ggh_walltime(0);
    self->D_g = _gghlite_dgsl_from_poly(self->g, self->params->sigma_p, NULL, DGSL_INLATTICE, flags);
    self->t_D_g = ggh_walltime(self->t_D_g);
}

static void
_gghlite_sk_sample_g(gghlite_sk_t self, aes_randstate_t randstate)
{
    assert(self->params);
    assert(self->params->n);

    fmpz_poly_init(self->g);
    fmpq_poly_init(self->g_inv);

    mpfr_t g_inv_norm;
    mpfr_init2(g_inv_norm, fmpz_sizeinbase(self->params->q,2));

    mpfr_t norm;
    mpfr_init2(norm, mpfr_get_prec(self->params->sigma));

    mpfr_t sqrtn_sigma;
    mpfr_init2(sqrtn_sigma, mpfr_get_prec(self->params->sigma));
    mpfr_set_si(sqrtn_sigma, self->params->n, MPFR_RNDN);
    mpfr_sqrt(sqrtn_sigma, sqrtn_sigma, MPFR_RNDN);
    mpfr_mul(sqrtn_sigma, sqrtn_sigma, self->params->sigma, MPFR_RNDN);

    const oz_flag_t flags = (self->params->flags & GGHLITE_FLAGS_QUIET) ? 0 : OZ_VERBOSE;
    dgsl_rot_mp_t *D = _gghlite_dgsl_from_n(self->params->n, self->params->sigma, flags);

    long fail[4] = {0,0,0,0};

    fmpq_poly_t g_q;
    fmpq_poly_init(g_q);


    const int check_prime = self->params->flags & GGHLITE_FLAGS_PRIME_G;

    int prime_pass = 0;
    mp_limb_t *primes_s = NULL, *primes_p;

    const int nsp = _gghlite_nsmall_primes(self->params);
    primes_p = _fmpz_poly_oz_ideal_probable_prime_factors(self->params->n, nsp);

    if (!check_prime) {
        primes_s = _fmpz_poly_oz_ideal_small_prime_factors(self->params->n, 2*(self->params->kappa+1));
    }

    fmpz_t N;
    fmpz_init(N);

    while(1) {
        ggh_fprintf(stderr, self->params, "\r      Computing g:: !n: %4ld, !p: %4ld, !i: %4ld, !N: %4ld",
                    fail[0], fail[1], fail[2], fail[3]);

        uint64_t t = ggh_walltime(0);
        fmpz_poly_sample_D(self->g, D, randstate);
        self->t_sample += ggh_walltime(t);

        fmpz_poly_2norm_mpfr(norm, self->g, MPFR_RNDN);
        if(mpfr_cmp(norm, sqrtn_sigma)>0) {
            fail[0]++;
            continue;
        }

        /* 1. check if prime */
        t = ggh_walltime(0);
        if (check_prime)
            prime_pass = fmpz_poly_oz_ideal_is_probaprime(self->g, self->params->n, 0, primes_p);
        else {
            /** we first check for probable prime factors */
            prime_pass = fmpz_poly_oz_ideal_not_prime_factors(self->g, self->params->n, primes_p);
            if (prime_pass) {
                /* if that passes we exclude small prime factors, regardless of how 
                 * probable they are */
                prime_pass = fmpz_poly_oz_ideal_not_prime_factors(self->g, self->params->n, primes_s);
            }
        }
        self->t_is_prime += ggh_walltime(t);
        if (!prime_pass) {
            fail[1]++;
            continue;
        }

        /* 2. check norm of inverse */
        fmpq_poly_set_fmpz_poly(g_q, self->g);
        _fmpq_poly_oz_invert_approx(self->g_inv, g_q, self->params->n, 2*self->params->lambda);
        if (!_gghlite_g_inv_check(self->params, self->g_inv)) {
            fail[2]++;
            continue;
        }

        fmpz_poly_oz_ideal_norm(N, self->g, self->params->n, 2);
        if (fmpz_sizeinbase(N, 2) < (size_t)self->params->n) {
            fail[3]++;
            continue;
        }

        break;
    }
    //4096 seems like a good choice
    const long prec = (self->params->n/4 < 8192) ? 8192 : self->params->n/4;
    if (self->params->flags & GGHLITE_FLAGS_GOOD_G_INV) {
        /** we compute the inverse in high precision for gghlite_enc_set_gghlite_clr **/
        _fmpq_poly_oz_invert_approx(self->g_inv, g_q, self->params->n, prec);
    }

    fmpz_clear(N);
    free(primes_p);
    if (!check_prime)
        free(primes_s);

    ggh_fprintf(stderr, self->params, "\n");
  
    mpfr_clear(norm);
    mpfr_clear(sqrtn_sigma);
    mpfr_clear(g_inv_norm);
    fmpq_poly_clear(g_q);
    dgsl_rot_mp_clear(D);
}


void
_gghlite_sk_set_pzt(gghlite_sk_t self)
{
    assert(self->params);
    assert(self->params->n);
    assert(fmpz_cmp_ui(self->params->q,0)>0);
    assert(!fmpz_poly_is_zero(self->h));

    fmpz_mod_poly_t z_kappa;  fmpz_mod_poly_init(z_kappa, self->params->q);

    if (gghlite_sk_is_symmetric(self)) {
        assert(!fmpz_mod_poly_is_zero(self->z[0]));
        fmpz_mod_poly_set(z_kappa, self->z[0]);
        fmpz_mod_poly_oz_ntt_pow_ui(z_kappa, z_kappa, self->params->kappa, self->params->n);
    } else {
        fmpz_mod_poly_oz_ntt_set_ui(z_kappa, 1, self->params->n);
        uint64_t t = ggh_walltime(0);
        for(size_t i=0; i<self->params->gamma; i++) {
            assert(!fmpz_mod_poly_is_zero(self->z[i]));
            fmpz_mod_poly_oz_ntt_mul(z_kappa, z_kappa, self->z[i], self->params->n);
            timer_printf("\r    Progress: [%lu / %lu] %8.2fs", i+1,
                         self->params->gamma, ggh_seconds(ggh_walltime(t)));
            fflush(stdout);
        }
        timer_printf("\n");
    }

    fmpz_mod_poly_t g_inv;  fmpz_mod_poly_init(g_inv, self->params->q);
    fmpz_mod_poly_oz_ntt_enc_fmpz_poly(g_inv, self->g, self->params->ntt);
    fmpz_mod_poly_oz_ntt_inv(g_inv, g_inv, self->params->n);

    fmpz_mod_poly_t pzt;  fmpz_mod_poly_init(pzt, self->params->q);
    fmpz_mod_poly_oz_ntt_mul(pzt, z_kappa, g_inv, self->params->n);

    fmpz_mod_poly_t h;  fmpz_mod_poly_init(h, self->params->q);
    fmpz_mod_poly_oz_ntt_enc_fmpz_poly(h, self->h, self->params->ntt);

    fmpz_mod_poly_oz_ntt_mul(pzt, pzt, h, self->params->n);

    fmpz_mod_poly_init(self->params->pzt, self->params->q);
    fmpz_mod_poly_set(self->params->pzt, pzt);

    fmpz_mod_poly_clear(h);
    fmpz_mod_poly_clear(pzt);
    fmpz_mod_poly_clear(z_kappa);
    fmpz_mod_poly_clear(g_inv);
}

void
_gghlite_sk_sample_h(gghlite_sk_t self, aes_randstate_t randstate)
{
    assert(self->params);
    assert(self->params->n);
    assert(fmpz_cmp_ui(self->params->q,0)>0);

    mpfr_t sqrt_q;
    mpfr_init2(sqrt_q, fmpz_sizeinbase(self->params->q, 2));
    _gghlite_params_get_q_mpfr(sqrt_q, self->params, MPFR_RNDN);
    mpfr_sqrt(sqrt_q, sqrt_q, MPFR_RNDN);

    /* dgsl samples proportionally to `\exp(-(x-c)²/(2σ²))` but GGHLite is
       specifiied with respect to `\exp(-π(x-c)²/σ²)`. So we divide by \sqrt{2π} */
    mpfr_mul_d(sqrt_q, sqrt_q, S_TO_SIGMA, MPFR_RNDN);

    fmpz_poly_init(self->h);

    /* we already ruled out probable prime factors when sampling <g> */
    mp_limb_t *primes = _fmpz_poly_oz_ideal_probable_prime_factors(self->params->n, 2);

    int coprime = 0;
    while(!coprime) {
        uint64_t t = ggh_walltime(0);
        fmpz_poly_sample_sigma(self->h, self->params->n, sqrt_q, randstate);
        self->t_sample += ggh_walltime(t);
        t = ggh_walltime(0);

        coprime = fmpz_poly_oz_coprime(self->g, self->h, self->params->n, 0, primes);
        self->t_coprime +=  ggh_walltime(t);
    }


    free(primes);
    mpfr_clear(sqrt_q);
}

void
_gghlite_sk_sample_z(gghlite_sk_t self, aes_randstate_t randstate)
{
    assert(self->params);
    assert(self->params->n);
    assert(fmpz_cmp_ui(self->params->q, 0)>0);

    const size_t bound = (gghlite_sk_is_symmetric(self)) ? 1 : self->params->gamma;
  
    /* do not parallelize calls to randomness generation! */  
    uint64_t t_init = ggh_walltime(0);
    for(size_t i=0; i<bound; i++) {
        fmpz_mod_poly_init(self->z[i], self->params->q);
        fmpz_mod_poly_randtest_aes(self->z[i], randstate, self->params->n);
        timer_printf("\r    Init Progress: [%lu / %lu] %8.2fs", i+1,
                     bound, ggh_seconds(ggh_walltime(t_init)));
    }
    timer_printf("\n");

    int progress_count_approx = 0;
    uint64_t t = ggh_walltime(0);
#pragma omp parallel for
    for(size_t i = 0; i < bound; i++) {
        fmpz_mod_poly_oz_ntt_enc(self->z[i], self->z[i], self->params->ntt);
        fmpz_mod_poly_init(self->z_inv[i], self->params->q);
        fmpz_mod_poly_oz_ntt_inv(self->z_inv[i], self->z[i], self->params->n);
#pragma omp critical
        {
            progress_count_approx++;
            timer_printf("\r    Computation Progress (Parallel): [%lu / %lu] %8.2fs",
                         progress_count_approx, bound, ggh_seconds(ggh_walltime(t)));
        }
    }
    timer_printf("\n");
}

void
gghlite_sk_init(gghlite_sk_t self, aes_randstate_t randstate)
{
    assert(self->params->lambda);
    assert(self->params->kappa);
    assert(self->params->gamma);

    {
        size_t nbytes;
        unsigned char *buf = random_aes(randstate, 128, &nbytes);
        aes_randinit_seedn(self->rng, (char *) buf, nbytes, NULL, 0);
        free(buf);
    }

    self->t_coprime = 0;
    self->t_is_prime = 0;
    self->t_sample = 0;

    self->z     = calloc(self->params->gamma, sizeof(gghlite_enc_t));
    self->z_inv = calloc(self->params->gamma, sizeof(gghlite_enc_t));
  
    start_timer();
    timer_printf("Starting precomp init...\n");
    fmpz_mod_poly_oz_ntt_precomp_init(self->params->ntt, self->params->n, self->params->q);
    timer_printf("Finished precomp init");
    print_timer();
    timer_printf("\n");

    start_timer();
    timer_printf("Starting sampling g...\n");
    _gghlite_sk_sample_g(self, randstate);
    timer_printf("Finished sampling g");
    print_timer();
    timer_printf("\n");
  
    start_timer();
    timer_printf("Starting sampling z...\n");
    _gghlite_sk_sample_z(self, randstate);
    timer_printf("Finished sampling z");
    print_timer();
    timer_printf("\n");
  
    start_timer();
    timer_printf("Starting sampling h...\n");
    _gghlite_sk_sample_h(self, randstate);
    timer_printf("Finished sampling h");
    print_timer();
    timer_printf("\n");

    start_timer();
    timer_printf("Starting setting D_g...\n");
    gghlite_sk_set_D_g(self);
    timer_printf("Finished setting D_g");
    print_timer();
    timer_printf("\n");

    start_timer();
    timer_printf("Starting setting pzt...\n");
    _gghlite_sk_set_pzt(self);
    timer_printf("Finished setting pzt");
    print_timer();
    timer_printf("\n");
}

void
gghlite_init(gghlite_sk_t self, const size_t lambda, const size_t kappa,
             const size_t gamma, const uint64_t rerand_mask,
             const gghlite_flag_t flags, aes_randstate_t randstate)
{
    memset(self, 0, sizeof(struct _gghlite_sk_struct));
    gghlite_params_init_gamma(self->params, lambda, kappa, gamma, rerand_mask, flags);
    gghlite_sk_init(self, randstate);
}

void
gghlite_sk_clear(gghlite_sk_t self, int clear_params)
{
    const size_t bound = (gghlite_sk_is_symmetric(self)) ? 1 : self->params->gamma;

    for(size_t i=0; i<bound; i++) {
        fmpz_mod_poly_clear(self->z[i]);
        fmpz_mod_poly_clear(self->z_inv[i]);
    }

    fmpz_poly_clear(self->h);
    fmpz_poly_clear(self->g);
    fmpq_poly_clear(self->g_inv);
    dgsl_rot_mp_clear(self->D_g);

    free(self->z);
    free(self->z_inv);

    if (clear_params)
        gghlite_params_clear(self->params);
}


void
gghlite_sk_print_norms(const gghlite_sk_t self)
{
    assert(self->params->n);
    assert(!fmpz_poly_is_zero(self->g));

    mpfr_t norm;
    mpfr_init2(norm, _gghlite_prec(self->params));

    mpfr_t bound;
    mpfr_init2(bound, _gghlite_prec(self->params));

    mpfr_t sqrt_n;
    mpfr_init2(sqrt_n, _gghlite_prec(self->params));
    mpfr_set_ui(sqrt_n, self->params->n, MPFR_RNDN);
    mpfr_sqrt(sqrt_n, sqrt_n, MPFR_RNDN);
    mpfr_log2(sqrt_n, sqrt_n, MPFR_RNDN);

    fmpz_poly_2norm_mpfr(norm, self->g, MPFR_RNDN);
    mpfr_log2(norm, norm, MPFR_RNDN);

    mpfr_set(bound, self->params->sigma, MPFR_RNDN);
    mpfr_log2(bound, bound, MPFR_RNDN);
    mpfr_add(bound, bound, sqrt_n, MPFR_RNDN);

    printf("  log(|g|): %6.1f ?< %6.1f\n", mpfr_get_d(norm, MPFR_RNDN), mpfr_get_d(bound, MPFR_RNDN));

    mpfr_set(bound, self->params->sigma_p, MPFR_RNDN);
    mpfr_log2(bound, bound, MPFR_RNDN);
    mpfr_add(bound, bound, sqrt_n, MPFR_RNDN);

    _gghlite_params_get_q_mpfr(bound, self->params, MPFR_RNDN);
    mpfr_sqrt(bound, bound, MPFR_RNDN);
    mpfr_log2(bound, bound, MPFR_RNDN);
    mpfr_add(bound, bound, sqrt_n, MPFR_RNDN);

    fmpz_poly_2norm_mpfr(norm, self->h, MPFR_RNDN);
    mpfr_log2(norm, norm, MPFR_RNDN);
    printf("  log(|h|): %6.1f ?< %6.1f\n", mpfr_get_d(norm, MPFR_RNDN), mpfr_get_d(bound, MPFR_RNDN));
}

void
gghlite_sk_print_times(const gghlite_sk_t self)
{
    printf("           sampling: %7.1fs\n", ggh_seconds(self->t_sample));
    printf("     primality test: %7.1fs\n", ggh_seconds(self->t_is_prime));
    printf("                D_g: %7.1fs\n", ggh_seconds(self->t_D_g));
    printf("gcd(N(g),N(h)) == 1: %7.1fs\n", ggh_seconds(self->t_coprime));
    /* printf("   <b_0,b_1> == <g>: %7.1fs\n", ggh_seconds(self->t_is_subideal)); */
}

dgsl_rot_mp_t *
_gghlite_dgsl_from_poly(fmpz_poly_t g, mpfr_t sigma, fmpq_poly_t c,
                        const dgsl_alg_t algorithm, const oz_flag_t flags)
{
    mpfr_t sigma_;
    mpfr_init2(sigma_, mpfr_get_prec(sigma));
    mpfr_set(sigma_, sigma, MPFR_RNDN);

    /* dgsl samples proportionally to $\exp(-(x-c)²/(2σ²))$ but GGHLite is
       specified with respect to $\exp(-π(x-c)²/σ²)$. So we divide by $\sqrt{2π}$
    */

    mpfr_mul_d(sigma_, sigma_, S_TO_SIGMA, MPFR_RNDN);

    dgsl_rot_mp_t *D = dgsl_rot_mp_init(fmpz_poly_length(g), g, sigma_, c, algorithm, flags);
    mpfr_clear(sigma_);
    return D;
}

dgsl_rot_mp_t *
_gghlite_dgsl_from_n(const long n, mpfr_t sigma, const oz_flag_t flags)
{
    fmpz_poly_t I;
    fmpz_poly_init(I);
    fmpz_poly_one(I);

    mpfr_t sigma_;
    mpfr_init2(sigma_, mpfr_get_prec(sigma));
    mpfr_set(sigma_, sigma, MPFR_RNDN);

    /* dgsl samples proportionally to $\exp(-(x-c)²/(2σ²))$ but GGHLite is
       specifiied with respect to $\exp(-π(x-c)²/σ²)$. So we divide by \sqrt{2π}
    */

    mpfr_mul_d(sigma_, sigma_, S_TO_SIGMA, MPFR_RNDN);

    dgsl_rot_mp_t *D = dgsl_rot_mp_init(n, I, sigma_, NULL, DGSL_IDENTITY, flags);
    fmpz_poly_clear(I);
    mpfr_clear(sigma_);
    return D;
}

/* Compute || [pzt \cdot op]_q ||_\infty */
void
_gghlite_enc_extract_raw(gghlite_clr_t rop, const gghlite_params_t self,
                         const gghlite_enc_t op)
{
    gghlite_enc_t t;
    gghlite_enc_init(t, self);
    fmpz_mod_poly_oz_ntt_mul(t, self->pzt, op, self->n);
    fmpz_mod_poly_oz_ntt_dec(t, t, self->ntt);
    fmpz_poly_set_fmpz_mod_poly(rop, t);
    gghlite_enc_clear(t);
}
