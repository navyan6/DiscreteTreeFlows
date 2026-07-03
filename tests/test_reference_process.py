"""Tests for ReferenceProcess and branching intensity."""
import torch
import numpy as np

from src.reference_process import BranchingIntensityMLP


class TestBranchingIntensityMLP:
    """Test 1-layer branching intensity MLP."""

    def test_forward_shape(self):
        """MLP outputs correct shape."""
        mlp = BranchingIntensityMLP(esm_c_dim=1024)

        # Batch of embeddings
        embeddings = torch.randn(4, 1024)
        lambdas = mlp(embeddings)

        assert lambdas.shape == (4,)
        assert torch.all(lambdas > 0)  # Softplus ensures positive

    def test_forward_single(self):
        """MLP works on single embedding."""
        mlp = BranchingIntensityMLP(esm_c_dim=1024)

        embedding = torch.randn(1024)
        lambda_val = mlp(embedding.unsqueeze(0))[0].item()

        assert lambda_val > 0
        assert isinstance(lambda_val, float)

    def test_parameters(self):
        """MLP has correct number of parameters."""
        mlp = BranchingIntensityMLP(esm_c_dim=1024)

        # Should have: 1024*1 weights + 1 bias = 1025 parameters
        total_params = sum(p.numel() for p in mlp.parameters())
        assert total_params == 1025

    def test_softplus_ensures_positive(self):
        """Output is always positive."""
        mlp = BranchingIntensityMLP(esm_c_dim=10)

        # Test with extreme values
        embeddings = torch.randn(100, 10) * 10
        lambdas = mlp(embeddings)

        assert torch.all(lambdas > 0)
        assert torch.all(torch.isfinite(lambdas))


class TestReferenceProcess:
    """Test reference process (requires ESM models, so basic only)."""

    def test_branching_mlp_trainable(self):
        """MLP can be optimized with Poisson loss."""
        mlp = BranchingIntensityMLP(esm_c_dim=1024)
        optimizer = torch.optim.Adam(mlp.parameters(), lr=0.01)

        # Simulate training
        embeddings = torch.randn(10, 1024)
        observed_offspring = torch.tensor([0.0, 1.0, 2.0, 1.0, 3.0, 0.0, 1.0, 2.0, 1.0, 0.0])

        initial_loss = None
        for _ in range(50):
            lambdas = mlp(embeddings)

            # Poisson loss: -log p(b | λ) = -b*log(λ) + λ + log(b!)
            loss = -(
                observed_offspring * torch.log(lambdas + 1e-8)
                - lambdas
                - torch.lgamma(observed_offspring + 1)
            )
            loss = loss.mean()

            if initial_loss is None:
                initial_loss = loss.item()

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        final_loss = loss.item()
        assert final_loss < initial_loss  # Loss should decrease


if __name__ == "__main__":
    import pytest

    pytest.main([__file__, "-v"])
