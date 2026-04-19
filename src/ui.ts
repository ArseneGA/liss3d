import type { Params } from "./types";

type NumericId = "p" | "q" | "r" | "A" | "B" | "C" | "phix" | "phiy" | "phiz" | "L" | "R" | "k";

const NUMERIC_DECIMALS: Record<NumericId, number> = {
  p: 2, q: 2, r: 2,
  A: 2, B: 2, C: 2,
  phix: 2, phiy: 2, phiz: 2,
  L: 0,
  R: 3, k: 3,
};

function el<T extends HTMLElement>(id: string): T {
  const node = document.getElementById(id);
  if (!node) throw new Error(`#${id} missing`);
  return node as T;
}

export function readParams(): Params {
  const num = (id: NumericId): number => parseFloat(el<HTMLInputElement>(id).value);
  const modeInput = document.querySelector<HTMLInputElement>('input[name="mode"]:checked');
  return {
    p: num("p"),
    q: num("q"),
    r: num("r"),
    A: num("A"),
    B: num("B"),
    C: num("C"),
    phix: num("phix"),
    phiy: num("phiy"),
    phiz: num("phiz"),
    L: num("L"),
    R: num("R"),
    k: num("k"),
    mode: (modeInput?.value === "B" ? "B" : "A"),
    resolution: parseInt(el<HTMLSelectElement>("res").value, 10),
    scaleMm: parseFloat(el<HTMLInputElement>("scaleMm").value),
  };
}

function clamp(v: number, lo: number, hi: number): number {
  return v < lo ? lo : v > hi ? hi : v;
}

function bindNumeric(id: NumericId, onChange: () => void): void {
  const range = el<HTMLInputElement>(id);
  const num = el<HTMLInputElement>(`${id}-num`);
  const dec = NUMERIC_DECIMALS[id];

  // Init : mirror range → number
  num.value = parseFloat(range.value).toFixed(dec);

  // Slider → number input + callback
  range.addEventListener("input", () => {
    num.value = parseFloat(range.value).toFixed(dec);
    onChange();
  });

  // Number input → slider (parse, clamp, snap au step)
  // On ne réécrit pas num.value en plein `input` pour préserver le curseur
  // et permettre de taper des valeurs intermédiaires (ex. "3.").
  num.addEventListener("input", () => {
    const v = parseFloat(num.value);
    if (!Number.isFinite(v)) return;
    const lo = parseFloat(range.min);
    const hi = parseFloat(range.max);
    const step = parseFloat(range.step);
    const clamped = clamp(v, lo, hi);
    const snapped = Number.isFinite(step) && step > 0
      ? lo + Math.round((clamped - lo) / step) * step
      : clamped;
    range.value = String(snapped);
    onChange();
  });

  // Au blur, reformater le number input avec la valeur effective du slider
  num.addEventListener("blur", () => {
    num.value = parseFloat(range.value).toFixed(dec);
  });

  // Enter = blur + commit
  num.addEventListener("keydown", (ev) => {
    if (ev.key === "Enter") {
      num.blur();
    }
  });
}

export function bindUI(onChange: (p: Params) => void): void {
  const numericIds: NumericId[] = ["p", "q", "r", "A", "B", "C", "phix", "phiy", "phiz", "L", "R", "k"];

  const fire = (): void => onChange(readParams());

  for (const id of numericIds) {
    bindNumeric(id, fire);
  }

  for (const id of ["res", "scaleMm"]) {
    el<HTMLInputElement>(id).addEventListener("input", fire);
  }

  for (const radio of document.querySelectorAll<HTMLInputElement>('input[name="mode"]')) {
    radio.addEventListener("change", fire);
  }
}
