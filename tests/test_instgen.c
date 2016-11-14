#include <gghlite/gghlite.h>

int
test_instgen_asymm(const size_t lambda, const size_t kappa, const uint64_t rerand, aes_randstate_t randstate)
{
    printf("symm: 0, λ: %4zu, κ: %2zu, rerand: 0x%016zx", lambda, kappa, rerand);

    gghlite_sk_t self;
    const gghlite_flag_t flags = GGHLITE_FLAGS_GDDH_HARD | GGHLITE_FLAGS_QUIET | GGHLITE_FLAGS_ASYMMETRIC;
    gghlite_init(self, lambda, kappa, kappa /* gamma */, rerand, flags, randstate);

    int status = 0;

    for(size_t i=0; i< kappa; i++) {
        /* we want all z[i], z_inv[i] to be != 0 */
        if (fmpz_mod_poly_is_zero(self->z[i]))         status++;
        if (fmpz_mod_poly_is_zero(self->z_inv[i]))     status++;
    }

    if (status == 0)
        printf(" (%d) PASS\n", status);
    else
        printf(" (%d) FAIL\n", status);

    gghlite_sk_clear(self, 1);
    return status;
}


int
main(int argc, char *argv[])
{
    aes_randstate_t randstate;
    aes_randinit(randstate);

    int status = 0;

    status += test_instgen_asymm(20, 2, 0x0, randstate);
    status += test_instgen_asymm(20, 4, 0x0, randstate);

    aes_randclear(randstate);
    flint_cleanup();
    mpfr_free_cache();
    return status;
}
