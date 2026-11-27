import Lean

namespace Wasserman.UniformCofactorPowers

/-- A reference weight residue persists at every exponent of a base congruent to one. -/
theorem referenceResidue (d s q e m : Nat) (unit : q % m = 1) :
    (1+d*s*q^e) % m = (1+d*s) % m := by
  have power : q^e % m = 1 % m := by
    rw [Nat.pow_mod, unit, Nat.one_pow]
  calc
    (1+d*s*q^e) % m = (1 % m + ((d*s) % m * (q^e % m)) % m) % m := by
      rw [Nat.add_mod 1, Nat.mul_mod (d*s)]
    _ = (1+d*s*1) % m := by
      rw [power, ← Nat.mul_mod, ← Nat.add_mod]
    _ = _ := by simp

/-- A common factor at two doubled exponents already divides the reference weight. -/
theorem doubledRestriction (m d x : Int) (first : m ∣ 1+d*x) (second : m ∣ 1+d*x^2) :
    m ∣ 1+d := by
  obtain ⟨a, ha⟩ := first
  obtain ⟨b, hb⟩ := second
  have diff : x-1 = m*(x*a-b) := by grind
  exact ⟨a-d*(x*a-b), by grind⟩

private theorem inverseThree (V : Nat) (eqn : (3*V) % 65536 = 64768) :
    V % 65536 = 65280 := by
  have scaled : (43691*(3*V)) % 65536 = 65280 := by rw [Nat.mul_mod, eqn]
  have id : 43691*(3*V) = V+65536*(2*V) := by omega
  simpa [id] using scaled

def s : Nat := 1306232213153877917593933367823536683066159885737114363997573716868895663450681
def q : Nat := 322880541893073397389196821002188355954113056549133075089295162936641374385930241

/-- For the recorded pair, the weight prime 251 certifies every exponent and the
    three identical companion parts have their first binary one at position 8. -/
theorem exampleLayers (e V : Nat) (weight : 3*V = 1+247*s*q^e) :
    (1+247*s*q^e)/251 % 251 = 1 ∧ (1+247*s*q^e) % 251 = 0 ∧ V % 65536 = 65280 := by
  have high := referenceResidue 247 s q e 63001 (by decide)
  have low := referenceResidue 247 s q e 65536 (by decide)
  have a : (1+247*s) % 63001 = 251 := by decide
  have b : (1+247*s) % 65536 = 64768 := by decide
  rw [a] at high
  rw [b] at low
  constructor
  · rw [← Nat.mod_mul_right_div_self, show 251*251 = 63001 by decide, high]
  constructor
  · rw [← Nat.mod_mul_right_mod _ 251 251, show 251*251 = 63001 by decide, high]
  · exact inverseThree V (by rw [weight, low])

end Wasserman.UniformCofactorPowers
