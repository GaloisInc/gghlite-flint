#include <stdarg.h>
#include <stdio.h>
#include "gghlite.h"
#include "gghlite-internals.h"

void
gghlite_enc_init(gghlite_enc_t op, const gghlite_params_t self)
{
    assert(!fmpz_is_zero(self->q));
    fmpz_mod_poly_init2(op, self->q, self->n);
}

void
gghlite_enc_set_gghlite_clr(gghlite_enc_t rop, const gghlite_sk_t self,
                            const gghlite_clr_t f, const size_t k, int *group,
                            const int rerand)
{
    fmpz_poly_t t_o; fmpz_poly_init(t_o);
    const mp_bitcnt_t prec = (self->params->n/4 < 8192) ? 8192 : self->params->n/4;
    const oz_flag_t flags = (self->params->flags & GGHLITE_FLAGS_VERBOSE) ? OZ_VERBOSE : 0;
    _fmpz_poly_oz_rem_small_iter(t_o, f, self->g, self->params->n, self->g_inv, prec, flags);

    if (rerand)
        dgsl_rot_mp_call_plus_fmpz_poly(t_o, self->D_g, t_o, *self->rng);

    // encode at level zero
    fmpz_mod_poly_oz_ntt_enc_fmpz_poly(rop, t_o, self->params->ntt);

    fmpz_poly_clear(t_o);

    // encode at level k
    if(k > 0) {
        if(!gghlite_sk_is_symmetric(self) && (k>1))
            ggh_die("Raising to higher levels than 1 not supported. Instead, multiply by the right combination of y_i.");

        for (unsigned int r = 0; r < self->params->gamma; r++) {
            if (group[r]) {
                if(gghlite_sk_is_symmetric(self)) {
                    for(size_t j=0; j<k; j++) // divide by z_i^k
                        fmpz_mod_poly_oz_ntt_mul(rop, rop, self->z_inv[r], self->params->n);
                    break;
                } else {
                    fmpz_mod_poly_oz_ntt_mul(rop, rop, self->z_inv[r], self->params->n);
                }
            }
        }
    }
}

int
gghlite_enc_is_zero(const gghlite_params_t self, const fmpz_mod_poly_t op)
{
    gghlite_clr_t t;
    mpfr_t norm, bound;
    int r;

    gghlite_clr_init(t);
    _gghlite_enc_extract_raw(t, self, op);

    mpfr_init2(norm, _gghlite_prec(self));
    mpfr_init2(bound, _gghlite_prec(self));
    fmpz_poly_2norm_mpfr(norm, t, MPFR_RNDN);

    {
        /* Set bound = q */
        mpz_t q_z;
        mpz_init(q_z);
        fmpz_get_mpz(q_z, self->q);
        mpfr_set_z(bound, q_z, MPFR_RNDN);
        mpz_clear(q_z);
    }

    {
        /* compute q^{1-xi} */
        mpfr_t ex;
        mpfr_init2(ex, _gghlite_prec(self));
        mpfr_set_ui(ex, 1, MPFR_RNDN);
        mpfr_sub(ex, ex, self->xi, MPFR_RNDN);
        mpfr_pow(bound, bound, ex, MPFR_RNDN);
        mpfr_clear(ex);
    }

    r = mpfr_cmp(norm, bound);

    mpfr_clear(bound);
    mpfr_clear(norm);
    gghlite_clr_clear(t);

    if (r <= 0)
        return 1;
    else
        return 0;
}
