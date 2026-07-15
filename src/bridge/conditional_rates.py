"""
Algorithm 1, line 11: ConditionalReferenceRates R^{0|T1}_t.

The Schrodinger-bridge training target for the mutation head. We match the
learned per-site jump distribution R_theta to the reference process P^0 *Doob
h-transformed to hit the observed terminal AA* x1 — i.e. bridge matching, not
terminal cross-entropy.

"""

import torch
import torch.nn.functional as F

PAD_IDX = 20


def conditional_bridge_log_target(
    log_R0_mut: torch.Tensor,  
    x1_idx: torch.Tensor,       
    t: float,
    c: float = 1.0,
) -> torch.Tensor:
    """
    Return log pi_cond, the log conditional-bridge target distribution over the
    20 amino acids, shape [..., 20], normalized along the last dim.
    """
    log_q = F.log_softmax(log_R0_mut, dim=-1)            
    log_kappa = -c * (1.0 - t)
    kappa = torch.exp(torch.tensor(log_kappa, device=log_R0_mut.device))
    log_1m_kappa = torch.log1p(-kappa).item()          

    x1_safe = x1_idx.clamp(0, 19).unsqueeze(-1)        
    log_q_x1 = log_q.gather(-1, x1_safe).squeeze(-1)   

    log_h_x1 = torch.logaddexp(
        torch.full_like(log_q_x1, log_kappa),
        log_1m_kappa + log_q_x1,
    )                                             
    log_h_other = log_1m_kappa + log_q_x1             
    boost = log_h_x1 - log_h_other                

    log_target = log_q.clone()
    log_target.scatter_add_(-1, x1_safe, boost.unsqueeze(-1))
    log_target = F.log_softmax(log_target, dim=-1)      
    return log_target


def conditional_bridge_kl(
    log_R_theta_mut: torch.Tensor, 
    log_R0_mut: torch.Tensor,     
    x1_idx: torch.Tensor,        
    t: float,
    c: float = 1.0,
) -> torch.Tensor:
    """
    Per-position KL( pi_cond || p_theta ), shape [...].  D_KL(R^{0|T1} || R_theta)
    localized to the site's destination distribution.
    """
    log_target = conditional_bridge_log_target(log_R0_mut, x1_idx, t, c)
    log_p_theta = F.log_softmax(log_R_theta_mut, dim=-1)               
    target = log_target.exp()
    return (target * (log_target - log_p_theta)).sum(-1)         


if __name__ == "__main__":
    torch.manual_seed(0)
    n, L = 2, 5
    log_R0 = torch.randn(n, L, 20)
    x1 = torch.randint(0, 20, (n, L))

    # 1. target is a valid distribution (sums to 1, nonneg)
    lt = conditional_bridge_log_target(log_R0, x1, t=0.5, c=1.0)
    p = lt.exp()
    assert torch.allclose(p.sum(-1), torch.ones(n, L), atol=1e-5), "target not normalized"
    assert (p >= 0).all()

    # 2. t -> 1 : target concentrates on x1, and KL -> cross-entropy -log p_theta(x1)
    log_theta = torch.randn(n, L, 20)
    gaps = []
    for t in (0.9, 0.99, 0.999, 0.9999):
        lt = conditional_bridge_log_target(log_R0, x1, t=t, c=1.0)
        mass_on_x1 = lt.exp().gather(-1, x1.unsqueeze(-1)).squeeze(-1)
        kl = conditional_bridge_kl(log_theta, log_R0, x1, t=t, c=1.0)
        ce = F.cross_entropy(log_theta.reshape(-1, 20), x1.reshape(-1), reduction="none").reshape(n, L)
        gap = (kl - ce).abs().max().item()
        gaps.append(gap)
        print(f"t={t}: min mass_on_x1={mass_on_x1.min():.5f}  max|KL-CE|={gap:.5f}")
    assert mass_on_x1.min() > 0.99, "target should concentrate on x1 as t->1"
    assert gaps == sorted(gaps, reverse=True), "KL-CE gap should shrink monotonically as t->1"
    assert gaps[-1] < 2e-3, "KL should approach CE as t->1"

    # 3. conserved site (x1 == current a): target favors staying at a
    a = torch.tensor([[3]])
    log_R0_c = torch.randn(1, 1, 20)
    lt_c = conditional_bridge_log_target(log_R0_c, a, t=0.5, c=1.0)
    q = F.softmax(log_R0_c, -1)
    # mass on x1 under target should exceed reference q(x1) (up-weighted)
    assert lt_c.exp()[0, 0, 3] > q[0, 0, 3], "conserved target should up-weight staying"

    # 4. reference sensitivity: boosting q at an off-target AA raises its target mass
    base = torch.zeros(1, 1, 20)
    x1_t = torch.tensor([[5]])
    lt0 = conditional_bridge_log_target(base, x1_t, t=0.5, c=1.0)
    ref = base.clone(); ref[0, 0, 9] += 2.0  # reference favors AA 9 (an off-target site)
    lt1 = conditional_bridge_log_target(ref, x1_t, t=0.5, c=1.0)
    assert lt1.exp()[0, 0, 9] > lt0.exp()[0, 0, 9], "reference should raise off-target target mass"

    print("conditional_rates self-tests passed")
