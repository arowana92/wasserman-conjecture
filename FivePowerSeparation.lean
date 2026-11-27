import Lean

namespace Wasserman.FivePowerSeparation

/-- The low binary digits needed by the cofactor-five obstruction. -/
def passes (q e m : Nat) : Bool :=
  let v := Nat.powMod q e m
  ((3*v % m) &&& (2*v % m)) == 0

def atMostOne (q m : Nat) : Bool :=
  let a := passes q 1 m
  let b := passes q 121 m
  let c := passes q 251 m
  !(a && b || a && c || b && c)

set_option maxRecDepth 20000
set_option maxHeartbeats 0

private theorem table :
    ∀ lo hi : Fin 64, atMostOne (3+8*(lo.val+64*hi.val)) 32768 = true := by decide

/-- The binary tests for the first, 121st and 251st powers cannot hold in pairs. -/
theorem threePowers (q : Nat) (residue : q % 8 = 3) :
    atMostOne (q % 32768) 32768 = true := by
  let v := q % 32768
  have hv : v < 32768 := Nat.mod_lt _ (by decide)
  have low : v % 8 = 3 := by dsimp [v]; omega
  let lo : Fin 64 := ⟨v/8 % 64, Nat.mod_lt _ (by decide)⟩
  let hi : Fin 64 := ⟨v/512, by omega⟩
  have eqn : 3+8*(lo.val+64*hi.val) = v := by dsimp [lo,hi]; omega
  simpa only [eqn] using table lo hi

/-- The nineteenth power is already incompatible at nine bits. -/
theorem nineteenTable :
    ∀ i : Fin 64, (passes (3+8*i.val) 1 512 && passes (3+8*i.val) 19 512) = false := by decide

#print axioms threePowers
#print axioms nineteenTable

end Wasserman.FivePowerSeparation
