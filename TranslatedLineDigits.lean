import Lean

namespace Wasserman.TranslatedLineDigits

/-- A direct certificate of carry-free addition, reading the low digits last. -/
inductive Free (p : Nat) : Nat → Nat → Prop where
  | zero : Free p 0 0
  | digit {a b i j : Nat} : Free p a b → i+j<p → Free p (p*a+i) (p*b+j)

def value (p : Nat) : List Nat → Nat
  | [] => 0
  | d :: ds => p * value p ds + d

private theorem zero_left (p : Nat) (ds : List Nat) (valid : ∀ d ∈ ds, d<p) :
    Free p 0 (value p ds) := by
  induction ds with
  | nil => exact .zero
  | cons d ds ih =>
    have h := ih (fun e he => valid e (by simp [he]))
    simpa [value] using Free.digit (i:=0) (j:=d) h (by simpa using valid d (by simp))

/-- Allocate the u lowest available digit tokens to the shifted first coordinate. -/
theorem allocate (p u : Nat) (ds : List Nat) (valid : ∀ d ∈ ds, d<p)
    (small : u<p) (enough : u≤ds.sum) :
    ∃ T Y, p*T+u+Y = T+value p ds ∧ Free p (p*T+u) Y ∧ Free p T (value p ds) := by
  induction ds generalizing u with
  | nil =>
    simp at enough
    subst u
    exact ⟨0,0,by simp [value],by simpa using Free.zero (p:=p),by exact Free.zero⟩
  | cons d ds ih =>
    have hd : d<p := valid d (by simp)
    have tail : ∀ e ∈ ds, e<p := fun e he => valid e (by simp [he])
    by_cases hu : u≤d
    · refine ⟨0,p*value p ds+(d-u),?_,?_,zero_left p (d::ds) valid⟩
      · simp [value]; omega
      · have h := Free.digit (i:=u) (j:=d-u) (zero_left p ds tail) (by omega)
        simpa using h
    · have rest : u-d ≤ ds.sum := by simp at enough; omega
      obtain ⟨t,y,eq,hxy,ht⟩ := ih (u-d) tail (by omega) rest
      let x := p*t+(u-d)
      refine ⟨x,p*y,?_,?_,?_⟩
      · have h : x+y=t+value p ds := eq
        have scaled := congrArg (p*·) h
        simp [Nat.mul_add] at scaled
        simp only [value]
        change p*x+u+p*y = x+(p*value p ds+d)
        have hx : x+d=p*t+u := by dsimp [x]; omega
        omega
      · simpa [x] using Free.digit (i:=u) (j:=0) hxy (by omega)
      · simpa [x,value] using Free.digit (i:=u-d) (j:=d) ht (by omega)

end Wasserman.TranslatedLineDigits
