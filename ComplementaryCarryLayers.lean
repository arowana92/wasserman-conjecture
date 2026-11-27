import Lean

namespace Wasserman.ComplementaryCarryLayers

/-- A perturbation smaller than M cannot move this nonintegral rational across a floor. -/
private theorem centered_floor (M d x q r u : Int) (hm : 0 < M) (hd : 0 < d)
    (hr : 0 < r ∧ r < d) (small : -M < u ∧ u < M)
    (center : d*x = M*(d*q+r)+u) : x/M = q := by
  apply (Int.ediv_eq_iff_of_pos hm).2
  have low := Int.mul_le_mul_of_nonneg_left (show 1 ≤ r by omega) (Int.le_of_lt hm)
  have high := Int.mul_le_mul_of_nonneg_left (show r ≤ d-1 by omega) (Int.le_of_lt hm)
  constructor
  · apply Int.le_of_mul_le_mul_left (a := d) _ hd
    grind
  · apply Int.lt_of_mul_lt_mul_left (a := d) _ (Int.le_of_lt hd)
    grind

/-- The two quotients for a single nonintegral center have a constant defect of one. -/
theorem floor_pair (M B d h x u z : Int) (hm : 0 < M) (hb : 1 ≤ B) (hd : 0 < d)
    (half : B+1=d*h) (small : -M < u ∧ u < M)
    (center : d*x=M*B*z+u) (nonintegral : ¬ d ∣ u) :
    x/M + x/(M*B) = h*z-1 := by
  let r := z%d
  let q := z/d
  have decomp : z = d*q+r := by
    have := Int.emod_add_mul_ediv z d
    dsimp [r,q]
    omega
  have hrem : 0 < r ∧ r < d := by
    have lo := Int.emod_nonneg z (by omega : d ≠ 0)
    have hi := Int.emod_lt_of_pos z hd
    have ne : r ≠ 0 := by
      intro zero
      apply nonintegral
      refine ⟨x-M*B*q, ?_⟩
      grind
    dsimp [r] at *
    omega
  have big := Int.mul_le_mul_of_nonneg_left hb (Int.le_of_lt hm)
  simp only [Int.mul_one] at big
  have second : x/(M*B) = q := centered_floor (M*B) d x q r u
    (by omega) hd hrem (by omega) (by grind)
  have first : x/M = h*z-q-1 := centered_floor M d x (h*z-q-1) (d-r) u
    hm hd (by omega) small (by grind)
  omega

/-- The formula with integral centers included. The implication on divisibility
    holds whenever M and B are units modulo d. -/
theorem signed_floor_pair (M B d h x u z : Int) (hm : 0 < M) (hb : 1 ≤ B)
    (hd : 0 < d) (half : B+1=d*h) (small : -M < u ∧ u < M)
    (center : d*x=M*B*z+u) (unit : d ∣ u → d ∣ z) :
    x/M + x/(M*B) = h*z - (if d ∣ u then (if u < 0 then 2 else 0) else 1) := by
  by_cases hu : d ∣ u
  · obtain ⟨v, hv⟩ := hu
    obtain ⟨r, hr⟩ := unit ⟨v,hv⟩
    have exactX : x = M*B*r+v := by
      apply Int.eq_of_mul_eq_mul_left (a := d) (by omega)
      grind
    have big := Int.mul_le_mul_of_nonneg_left hb (Int.le_of_lt hm)
    simp only [Int.mul_one] at big
    have vs : -M < v ∧ v < M := by
      by_cases sign : 0 ≤ v
      · have := Int.mul_le_mul_of_nonneg_right (show 1 ≤ d by omega) sign
        grind
      · have := Int.mul_le_mul_of_nonneg_right (show 1 ≤ d by omega) (show 0 ≤ -v by omega)
        grind
    have integral : d ∣ u := ⟨v,hv⟩
    simp only [ite_eq_left integral]
    by_cases neg : u < 0
    · have vn : v < 0 := by
        by_cases hn : v < 0
        · exact hn
        have := Int.mul_nonneg (Int.le_of_lt hd) (show 0 ≤ v by omega)
        omega
      have a : x/M = B*r-1 := (Int.ediv_eq_iff_of_pos hm).2 (by grind)
      have b : x/(M*B) = r-1 := (Int.ediv_eq_iff_of_pos (by omega)).2 (by grind)
      simp only [ite_eq_left neg]
      grind
    · have vpos : 0 ≤ v := by
        by_cases hn : 0 ≤ v
        · exact hn
        have := Int.mul_neg_of_pos_of_neg hd (show v < 0 by omega)
        omega
      have a : x/M = B*r := (Int.ediv_eq_iff_of_pos hm).2 (by grind)
      have b : x/(M*B) = r := (Int.ediv_eq_iff_of_pos (by omega)).2 (by grind)
      simp only [ite_eq_right neg]
      grind
  · simpa only [ite_eq_right hu] using floor_pair M B d h x u z hm hb hd half small center hu

/-- For a binomial, the two Legendre layer contributions sum to exactly one. -/
theorem binomial_pair (M B d h x y u v z t : Int) (hm : 0 < M) (hb : 1 ≤ B)
    (hd : 0 < d) (half : B+1=d*h)
    (su : -M < u ∧ u < M) (sv : -M < v ∧ v < M) (ss : -M < u+v ∧ u+v < M)
    (cx : d*x=M*B*z+u) (cy : d*y=M*B*t+v)
    (nu : ¬ d ∣ u) (nv : ¬ d ∣ v) (ns : ¬ d ∣ u+v) :
    ((x+y)/M-x/M-y/M) + ((x+y)/(M*B)-x/(M*B)-y/(M*B)) = 1 := by
  have a := floor_pair M B d h x u z hm hb hd half su cx nu
  have b := floor_pair M B d h y v t hm hb hd half sv cy nv
  have c := floor_pair M B d h (x+y) (u+v) (z+t) hm hb hd half ss (by grind) ns
  grind

end Wasserman.ComplementaryCarryLayers
