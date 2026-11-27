import Lean

namespace Wasserman.CanonicalCofactorCut

/-- The scale times the residue equals the seed times the normalized part. -/
theorem companion_projection (B U V f n r : Int)
    (weight : B ∣ U*V-n) (part : B ∣ V*f-r)
    (lo : -B < U*r-n*f) (hi : U*r-n*f < B) : U*r = n*f := by
  obtain ⟨x, hx⟩ := weight
  obtain ⟨y, hy⟩ := part
  have h₁ := congrArg (fun z : Int => f*z) hx
  have h₂ := congrArg (fun z : Int => U*z) hy
  have eqn : U*r-n*f = B*(f*x-U*y) := by grind
  have small : (U*r-n*f).natAbs < B.natAbs := by omega
  have zero := Int.eq_zero_of_dvd_of_natAbs_lt_natAbs
    (show B ∣ U*r-n*f from ⟨_, eqn⟩) small
  omega

/-- Bézout coefficients for a primitive triple recover its denominator. -/
theorem denominator (u n f g h r s t x y z : Int)
    (primitive : x*f+y*g+z*h = 1)
    (e₁ : u*r = n*f) (e₂ : u*s = n*g) (e₃ : u*t = n*h) : u ∣ n := by
  refine ⟨x*r+y*s+z*t, ?_⟩
  have h₀ := congrArg (fun v : Int => n*v) primitive
  have h₁ := congrArg (fun v : Int => x*v) e₁
  have h₂ := congrArg (fun v : Int => y*v) e₂
  have h₃ := congrArg (fun v : Int => z*v) e₃
  grind

end Wasserman.CanonicalCofactorCut
