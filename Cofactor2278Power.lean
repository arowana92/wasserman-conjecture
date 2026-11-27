import Lean

namespace Wasserman

set_option exponentiation.threshold 3000

private theorem powerResidue (x e m : Nat) :
    x^e % m = Nat.powMod (x % m) e m := by
  rw [Nat.powMod_def, ← Nat.pow_mod]

/-- A unit raised to the index 52 * 53 has residue one modulo 53 squared. -/
theorem unitPower2756 (x : Nat) (unit : x % 53 ≠ 0) : x^2756 % 2809 = 1 := by
  have fermat : ∀ i : Fin 53, i.val ≠ 0 → Nat.powMod i.val 52 53 = 1 := by decide
  have low : x^52 % 53 = 1 := by
    rw [powerResidue]
    exact fermat ⟨x % 53, Nat.mod_lt _ (by decide)⟩ unit
  have lift : ∀ i : Fin 53, Nat.powMod (1+53*i.val) 53 2809 = 1 := by decide
  have high : (x^52)^53 % 2809 = 1 := by
    let i : Fin 53 := ⟨(x^52 % 2809)/53, by omega⟩
    have h : x^52 % 2809 = 1+53*i.val := by dsimp [i]; omega
    rw [powerResidue, h]
    exact lift i
  simpa only [← Nat.pow_mul] using high

/-- The residue at the second 53-adic layer forces a carry in V, 2V, 4V. -/
theorem cofactor2278_auxiliaryCarry (x V : Nat) (unit : x % 53 ≠ 0)
    (weight : 7*V = 2278*x^2756+1) :
    V % 2809 = 2332 ∧ 2809 ≤ V % 2809 + (2*V) % 2809 + (4*V) % 2809 := by
  have h := unitPower2756 x unit
  have coeff : (2278*x^2756+1) % 2809 = 2279 := by
    rw [Nat.add_mod, Nat.mul_mod, h]
  have eqn : (7*V) % 2809 = 2279 := by rw [weight, coeff]
  have inv : V % 2809 = 2332 := by
    have scaled : (1204*(7*V)) % 2809 = 2332 := by rw [Nat.mul_mod, eqn]
    have identity : 1204*(7*V) = V+2809*(3*V) := by omega
    simpa [identity] using scaled
  exact ⟨inv, by omega⟩

end Wasserman
