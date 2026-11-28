import Lean

namespace Wasserman.PairDivisor

private def fac : Nat → Nat
  | 0 => 1
  | n + 1 => (n+1) * fac n

private theorem fac_pos : ∀ n, 0 < fac n
  | 0 => by decide
  | n + 1 => Nat.mul_pos (Nat.zero_lt_succ n) (fac_pos n)

/-- Number of shuffles of two lists of specified lengths. -/
def shuffles : Nat → Nat → Nat
  | 0, _ => 1
  | _ + 1, 0 => 1
  | a + 1, b + 1 => shuffles a (b+1) + shuffles (a+1) b

private theorem shuffle_factorial : ∀ a b, shuffles a b * (fac a * fac b) = fac (a+b)
  | 0, b => by simp [shuffles, fac]
  | a+1, 0 => by simp [shuffles, fac]
  | a+1, b+1 => by
    have h₁ := shuffle_factorial a (b+1)
    have h₂ := shuffle_factorial (a+1) b
    have same : a+1+b = a+(b+1) := by omega
    rw [same] at h₂
    calc
      _ = (a+1) * (shuffles a (b+1) * (fac a * fac (b+1))) +
          (b+1) * (shuffles (a+1) b * (fac (a+1) * fac b)) := by
            rw [shuffles, Nat.add_mul]
            congr 1 <;> simp only [fac] <;> ac_rfl
      _ = (a+1+(b+1)) * fac (a+(b+1)) := by rw [h₁, h₂, ← Nat.add_mul]
      _ = _ := by
        change _ = (a+1+b+1) * fac (a+1+b)
        rw [same]
        congr 1 <;> omega

/-- The positive-row theorem also uses this coefficient after deleting units. -/
def trinomial (x y z : Nat) : Nat := shuffles x y * shuffles (x+y) z

private theorem coefficient_factorial (x y z : Nat) :
    trinomial x y z * (fac x * fac y * fac z) = fac (x+y+z) := by
  unfold trinomial
  calc
    _ = shuffles (x+y) z * ((shuffles x y * (fac x * fac y)) * fac z) := by ac_rfl
    _ = _ := by rw [shuffle_factorial, shuffle_factorial]

private theorem delete_left (x y z : Nat) :
    trinomial (x+1) y z * (x+1) = (x+1+y+z) * trinomial x y z := by
  apply Nat.eq_of_mul_eq_mul_right (Nat.mul_pos (Nat.mul_pos (fac_pos x) (fac_pos y)) (fac_pos z))
  calc
    _ = trinomial (x+1) y z * (fac (x+1) * fac y * fac z) := by simp only [fac]; ac_rfl
    _ = fac (x+1+y+z) := coefficient_factorial ..
    _ = (x+1+y+z) * fac (x+y+z) := by
      have e : x+1+y+z = (x+y+z)+1 := by omega
      rw [e, fac]
    _ = _ := by rw [← coefficient_factorial x y z]; ac_rfl

private theorem swap (x y z : Nat) : trinomial x y z = trinomial y x z := by
  apply Nat.eq_of_mul_eq_mul_right (Nat.mul_pos (Nat.mul_pos (fac_pos x) (fac_pos y)) (fac_pos z))
  rw [coefficient_factorial, Nat.mul_comm (fac x) (fac y), coefficient_factorial, Nat.add_comm x y]

private theorem delete_right (x y z : Nat) :
    trinomial x (y+1) z * (y+1) = (x+(y+1)+z) * trinomial x y z := by
  rw [swap x (y+1) z, delete_left, swap y x z, Nat.add_comm (y+1) x]

private theorem quotient_dvd {e k t : Nat} (pos : 0 < e) (h : e ∣ k*t) :
    e / e.gcd k ∣ t := by
  apply Nat.dvd_of_mul_dvd_mul_left (Nat.gcd_pos_of_pos_left k pos)
  rw [Nat.mul_div_cancel' (Nat.gcd_dvd_left e k), ← Nat.gcd_mul_right]
  exact Nat.dvd_gcd (Nat.dvd_mul_right e t) h

private theorem pair_bound (x y : Nat) :
    ((x*y).gcd (x*(x-1))).gcd (y*(y-1)) ∣
      x.gcd y * x.gcd (y-1) * (x-1).gcd y := by
  let w := y.gcd (x-1)
  obtain ⟨r, hr⟩ := Nat.gcd_dvd_left y (x-1)
  change y = w*r at hr
  have h₁ : r ∣ y := ⟨w, by rw [hr, Nat.mul_comm]⟩
  have h₂ := Nat.gcd_dvd_gcd_of_dvd_right x h₁
  have h₃ : x.gcd (r*(y-1)) ∣ x.gcd r * x.gcd (y-1) := by
    first | exact Nat.gcd_mul_right_dvd_mul_gcd .. | exact Nat.gcd_mul_dvd_mul_gcd ..
  rw [Nat.gcd_mul_left, ← Nat.gcd_comm y (x-1)]
  change (x*w).gcd (y*(y-1)) ∣ x.gcd y * x.gcd (y-1) * w
  have h₄ : y*(y-1) = (r*(y-1))*w := calc
    _ = (w*r)*(y-1) := congrArg (fun t => t*(y-1)) hr
    _ = _ := by ac_rfl
  rw [h₄, Nat.gcd_mul_right]
  exact Nat.mul_dvd_mul (Nat.dvd_trans h₃ (Nat.mul_dvd_mul h₂ (Nat.dvd_refl _))) (Nat.dvd_refl _)

private theorem delete_two_left (x y z : Nat) :
    (x+y+z)*((x+y+z)-1) ∣ trinomial x y z * (x*(x-1)) := by
  cases x with
  | zero => simp [Nat.dvd_zero]
  | succ x =>
    cases x with
    | zero => simp [Nat.dvd_zero]
    | succ x =>
      refine ⟨trinomial x y z, ?_⟩
      have prev : x+2+y+z-1 = x+1+y+z := by omega
      change trinomial (x+2) y z * ((x+2)*(x+1)) =
        (x+2+y+z)*(x+2+y+z-1)*trinomial x y z
      rw [← Nat.mul_assoc, delete_left, Nat.mul_assoc, delete_left, prev, Nat.mul_assoc]

private theorem delete_separate (x y z : Nat) :
    (x+y+z)*((x+y+z)-1) ∣ trinomial x y z * (x*y) := by
  cases x with
  | zero => simp [Nat.dvd_zero]
  | succ x =>
    cases y with
    | zero => simp [Nat.dvd_zero]
    | succ y =>
      refine ⟨trinomial x y z, ?_⟩
      have prev : (x+1+(y+1)+z)-1 = x+(y+1)+z := by omega
      rw [← Nat.mul_assoc, delete_left, Nat.mul_assoc, delete_right, prev, Nat.mul_assoc]

/-- Two selected parts bound the portion of the endpoint coefficient left after
taking its gcd with the trinomial coefficient. -/
theorem endpoint_quotient_dvd_pair (x y z : Nat) (hx : 0 < x) (hy : 0 < y) :
    let e := (x+y+z)*((x+y+z)-1)
    e / e.gcd (trinomial x y z) ∣ x.gcd y * x.gcd (y-1) * (x-1).gcd y := by
  let e := (x+y+z)*((x+y+z)-1)
  have he : 0 < e := Nat.mul_pos (by omega) (by omega)
  have yy := delete_two_left y x z
  rw [Nat.add_comm y x, swap y x z] at yy
  exact Nat.dvd_trans (Nat.dvd_gcd (Nat.dvd_gcd
    (quotient_dvd he (delete_separate x y z))
    (quotient_dvd he (delete_two_left x y z))) (quotient_dvd he yy)) (pair_bound x y)

/-- Two rows sharing a part h have exact endpoint gcd Q if the complement
h*(h-1) is coprime to their common coefficient divisor. -/
theorem exact_shared_part_gcd (h a b c d n Q : Nat) (hh : 2 ≤ h)
    (w₁ : h+a+b = n) (w₂ : h+c+d = n)
    (factor : n*(n-1) = (h*(h-1))*Q)
    (avoid : (h*(h-1)).Coprime ((trinomial h a b).gcd (trinomial h c d))) :
    (n*(n-1)).gcd ((trinomial h a b).gcd (trinomial h c d)) = Q := by
  have pos : 0 < h*(h-1) := Nat.mul_pos (by omega) (by omega)
  have h₁ := delete_two_left h a b
  have h₂ := delete_two_left h c d
  rw [w₁, factor, Nat.mul_comm (trinomial h a b)] at h₁
  rw [w₂, factor, Nat.mul_comm (trinomial h c d)] at h₂
  have common := Nat.dvd_gcd ((Nat.mul_dvd_mul_iff_left pos).mp h₁)
                              ((Nat.mul_dvd_mul_iff_left pos).mp h₂)
  obtain ⟨t, ht⟩ := common
  have htdiv : t ∣ (trinomial h a b).gcd (trinomial h c d) :=
    ⟨Q, by rw [ht, Nat.mul_comm]⟩
  have cop := Nat.Coprime.coprime_dvd_left htdiv avoid.symm
  rw [factor, ht, Nat.mul_comm (h*(h-1)) Q, Nat.gcd_mul_left, cop.symm.gcd_eq_one,
      Nat.mul_one]

end Wasserman.PairDivisor
