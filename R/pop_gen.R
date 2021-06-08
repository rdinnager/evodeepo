evd_fst <- function(x, pops = NULL, method = c("WG17")) {

  if(!inherits(snps, "torch_tensor")) {
    snps <- as_tensor(x, NA.method = "asis")
  }

  if(is.null(pops) & inherits(x, "genlight")) {
    pops <- adegenet::pop(x)
  }

  n_snps <- ncol(snps)

  max_count <- torch::torch_max(snps[!torch::torch_isnan(snps)])

  al_1 <- snps
  al_2 <- max_count - snps

  ## nans
  ind_snps <- (!torch::torch_isnan(snps)) %>%
    torch::torch_sum(2L)

  al_1 <- torch::torch_where(torch::torch_isnan(al_1),
                             torch::torch_zeros_like(al_1),
                             al_1)

  al_2 <- torch::torch_where(torch::torch_isnan(al_2),
                             torch::torch_zeros_like(al_2),
                             al_2)

  snp_tens <- torch::torch_stack(list(al_1, al_2), 3L)

  dos_sums <- torch::torch_tensordot(snp_tens, snp_tens$permute(c(3L, 2L, 1L)),
                                     list(c(-1L, -2L), c(1L, 2L))) / 4

  dos_sums / torch::torch_unsqueeze(ind_snps, 1L) /
    torch::torch_unsqueeze(ind_snps, 2L) / (n_snps^2)

}
