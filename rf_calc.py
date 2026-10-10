"""TLC 1 枚の Rf と溶媒比から、同じ化合物で目標の Rf になる溶媒比を見積もる（シリカプレート）。

モード 1: 同じ溶媒系で比を変える
  R_M = log10((1 - Rf) / Rf)
  R_M = a - n * log10(X_B)        (X_B: 強溶媒 B のモル分率, n は 1〜2 くらい)
  → X_B(新) = X_B(今) * 10 ** ((R_M(今) - R_M(目標)) / n)
  1 枚だけだと n が決まらないので n = 1 と n = 2 の範囲で出す。2 枚目を入力すると n を実測から決める。

モード 2: 別の溶媒系に移す (ε_mix が同じなら似た保持になる、という Snyder の考え方)
  2 枚目（今の溶媒系の別の比）を入力すると、s = (R_M1 - R_M2) / (ε2 - ε1) を実測から決める。
  ε_AB = ε_A + log10( N_B * 10^(α n_b (ε_B - ε_A)) + (1 - N_B) ) / (α n_b)   (Snyder-Soczewinski)
  目標の Rf に必要な ε の変化: Δε = (R_M(今) - R_M(目標)) / s
  s = Δlog10 k' / Δε。経験則「ε が 0.05 変わると k' は約 2〜3 倍」から s = 6〜9.5 とし、範囲で出す。
  n_b (溶媒 B の分子が占める吸着面積) は分子体積からの近似 n_b = 6.0 * (V/89.4)^(2/3)、α = 0.57。
  どちらも近似で、溶媒の選択性（水素結合など）の違いは入っていない。新しい系で 1 枚測って確かめること。

出典:
  ε° (シリカ) : Snyder, in Horvath (ed.), High-Performance Liquid Chromatography Vol.3 (1983) /
                Burdick & Jackson Solvent Guide 3rd ed. (1990) の表（教科書に引用されたもの）
  モル体積    : Zeng et al., Polymer Handbook Ch.16 (cm3/mol)
"""
import math

# モル体積 (cm3/mol)。番号は、TLC でよく使う順（若い番号ほどよく使う）。
MOLAR_VOLUME = {"Hexane": 131.6, "AcOEt": 98.5, "CHCl3": 80.7, "MeOH": 40.7,
                "CH2Cl2": 63.9, "Toluene": 106.8, "Et2O": 104.8, "Acetone": 74.0,
                "THF": 81.7, "MeCN": 52.6, "EtOH": 58.5, "2-Propanol": 76.8,
                "Benzene": 89.4, "Dioxane": 85.7}
# Snyder の溶離力 ε° (シリカ)。表にない溶媒（EtOH, Benzene など）は、モード 2 では使えない。
EPS_SILICA = {"Hexane": 0.01, "Toluene": 0.22, "CHCl3": 0.26, "CH2Cl2": 0.30, "Et2O": 0.43,
              "AcOEt": 0.48, "Dioxane": 0.51, "MeCN": 0.52, "Acetone": 0.53, "THF": 0.53,
              "2-Propanol": 0.60, "MeOH": 0.70}
ALPHA = 0.57
S_RANGE = (math.log10(2) / 0.05, math.log10(3) / 0.05)   # k' が 0.05 の ε で 2 倍〜3 倍


def r_m(rf):
    return math.log10((1 - rf) / rf)


def parse_ratio(text):
    weak, strong = (float(v) for v in text.split(":"))
    return weak, strong


def mole_fraction_strong(weak, strong, ratio):
    """体積比 弱:強 → 強溶媒のモル分率。"""
    vw, vs = parse_ratio(ratio)
    nw = vw / MOLAR_VOLUME[weak]
    ns = vs / MOLAR_VOLUME[strong]
    return ns / (nw + ns)


def volume_ratio_from_x(weak, strong, x):
    """強溶媒のモル分率 → (弱:強 を強=1 にした体積比, 強の体積分率)。"""
    vol_strong = x * MOLAR_VOLUME[strong]
    vol_weak = (1 - x) * MOLAR_VOLUME[weak]
    return vol_weak / vol_strong, vol_strong / (vol_strong + vol_weak)


def eps_mix(weak, strong, x_strong):
    e_a, e_b = EPS_SILICA[weak], EPS_SILICA[strong]
    k = ALPHA * 6.0 * (MOLAR_VOLUME[strong] / 89.4) ** (2 / 3)
    return e_a + math.log10(x_strong * 10 ** (k * (e_b - e_a)) + (1 - x_strong)) / k


def x_for_eps(weak, strong, target):
    """ε_mix = target になる強溶媒のモル分率（二分法）。範囲外なら None。"""
    lo, hi = EPS_SILICA[weak], EPS_SILICA[strong]
    if not lo < target < hi:
        return None
    a, b = 0.0, 1.0
    for _ in range(80):
        mid = (a + b) / 2
        if eps_mix(weak, strong, mid) < target:
            a = mid
        else:
            b = mid
    return (a + b) / 2


# ---- 入力 ------------------------------------------------------------------

def ask(label, default=None):
    text = input(f"{label}" + (f" [{default}]" if default is not None else "") + ": ").strip()
    return text if text or default is None else default


def print_solvents(names):
    for i, name in enumerate(names, 1):
        print(f"  {i:2d}: {name}")


def ask_solvent(label, default, names):
    """番号でも名前（大文字小文字は区別しない）でも選べる。"""
    while True:
        text = ask(label, default)
        if text.isdigit() and 1 <= int(text) <= len(names):
            return names[int(text) - 1]
        for known in names:
            if known.lower() == text.lower():
                return known
        print("  番号か名前で選んでください。")
        print_solvents(names)


def ask_ratio(label, optional=False):
    while True:
        text = input(f"{label}: ").strip()
        if optional and not text:
            return None
        try:
            if all(v > 0 for v in parse_ratio(text)):
                return text
        except ValueError:
            pass
        print("  例 3:1 のように、正の数で入力してください。")


def ask_rf(label, optional=False):
    while True:
        text = input(f"{label}: ").strip()
        if optional and not text:
            return None
        try:
            value = float(text)
            if 0 < value < 1:
                return value
        except ValueError:
            pass
        print("  0 と 1 の間の数で入力してください。")


def ask_system(title, names, default_weak, default_strong):
    print(title)
    weak = ask_solvent("  弱い溶媒", default_weak, names)
    strong = ask_solvent("  強い溶媒", default_strong, names)
    return weak, strong


def order_by_strength(weak, strong):
    """ε° の小さいほうを弱、大きいほうを強にそろえる（モード 2 用）。"""
    return (weak, strong) if EPS_SILICA[weak] <= EPS_SILICA[strong] else (strong, weak)


# ---- モード 1 ---------------------------------------------------------------

def same_system():
    names = list(MOLAR_VOLUME)
    print("使える溶媒（番号でも名前でも入力できます）:")
    print_solvents(names)
    weak, strong = ask_system("溶媒系", names, "例:Hexane, 1", "例:AcOEt, 2")
    ratio = ask_ratio("今の比 (弱:強, 体積比) 例 3:1")
    rf = ask_rf("今の Rf")
    target = ask_rf("目標の Rf")
    ratio2 = ask_ratio("2 枚目の比 (なければ空 Enter)", optional=True)
    rf2 = ask_rf("2 枚目の Rf") if ratio2 else None

    x1 = mole_fraction_strong(weak, strong, ratio)
    rm1, rm_t = r_m(rf), r_m(target)
    print()
    print(f"今: {weak}:{strong} = {ratio}  X_B={x1:.3f}  Rf={rf}  R_M={rm1:.3f}")
    print(f"目標: Rf={target}  R_M={rm_t:.3f}  (ΔR_M={rm1 - rm_t:+.3f})")

    if ratio2:
        x2 = mole_fraction_strong(weak, strong, ratio2)
        if x2 == x1 or rf2 == rf:
            print("2 枚目が 1 枚目と同じ比・Rf なので、n は決まりません。")
            return
        n = (r_m(rf) - r_m(rf2)) / (math.log10(x2) - math.log10(x1))
        print(f"2 枚目: {ratio2}  X_B={x2:.3f}  Rf={rf2}  → 実測の n = {n:.2f}")
        if n <= 0:
            print("n が 0 以下です（強溶媒を増やしても Rf が上がっていない）。測定を確認してください。")
            return
        ns = [n]
    else:
        ns = [1.0, 2.0]
        print("2 枚目なし: n=1 と n=2 の範囲で見積もる")

    for n in ns:
        x_new = x1 * 10 ** ((rm1 - rm_t) / n)
        if x_new >= 1:
            print(f"n={n:.2f}: 強溶媒 100% でも届かない（X_B={x_new:.2f}）")
            continue
        weak_per_strong, phi = volume_ratio_from_x(weak, strong, x_new)
        print(f"n={n:.2f}: X_B={x_new:.3f}  → {weak}:{strong} ≒ {weak_per_strong:.2f}:1"
              f"  (強溶媒の体積分率 {phi * 100:.0f}%)")


# ---- モード 2 ---------------------------------------------------------------

def other_system():
    names = [s for s in MOLAR_VOLUME if s in EPS_SILICA]
    print("使える溶媒（シリカの ε° とモル体積が両方ある。番号でも名前でも入力できます）:")
    print_solvents(names)
    weak, strong = ask_system("今の溶媒系", names, "例:Hexane, 1", "例:AcOEt, 2")
    weak, strong = order_by_strength(weak, strong)
    ratio = ask_ratio(f"今の比 ({weak}:{strong}, 体積比) 例 3:1")
    rf = ask_rf("今の Rf")
    target = ask_rf("目標の Rf")
    new_weak, new_strong = ask_system("移したい溶媒系", names, "例:CHCl3, 3", "例:MeOH, 4")
    new_weak, new_strong = order_by_strength(new_weak, new_strong)
    if (new_weak, new_strong) == (weak, strong):
        print("同じ溶媒系です。モード 1 を使ってください。")
        return
    ratio2 = ask_ratio(f"2 枚目の比 (今の溶媒系 {weak}:{strong} の別の比。なければ空 Enter)", optional=True)
    rf2 = ask_rf("2 枚目の Rf") if ratio2 else None

    x_old =mole_fraction_strong(weak, strong, ratio)
    e_old = eps_mix(weak, strong, x_old)
    rm1, rm_t = r_m(rf), r_m(target)
    print()
    print(f"今: {weak}:{strong} = {ratio}  X_B={x_old:.3f}  ε_mix={e_old:.3f}  Rf={rf}  R_M={rm1:.3f}")
    print(f"目標: Rf={target}  R_M={rm_t:.3f}  (ΔR_M={rm1 - rm_t:+.3f})")
    print(f"移す先: {new_weak}:{new_strong}  (弱 ε°={EPS_SILICA[new_weak]}, 強 ε°={EPS_SILICA[new_strong]})")

    if ratio2:
        e2 = eps_mix(weak, strong, mole_fraction_strong(weak, strong, ratio2))
        print(f"2 枚目: {weak}:{strong} = {ratio2}  ε_mix={e2:.3f}  Rf={rf2}  R_M={r_m(rf2):.3f}")
        if abs(e2 - e_old) < 1e-9 or rf2 == rf:
            print("2 枚目が 1 枚目と同じ比・Rf なので、s は決まりません。")
            return
        s_measured = (rm1 - r_m(rf2)) / (e2 - e_old)
        print(f"  → 実測の s = {s_measured:.1f}（経験則の範囲は {S_RANGE[0]:.1f}〜{S_RANGE[1]:.1f}）")
        if s_measured <= 0:
            print("s が 0 以下です（強溶媒を増やしても Rf が上がっていない）。測定を確認してください。")
            return
        s_values = [s_measured]
    else:
        print("2 枚目なし: 経験則の s の範囲で見積もる")
        s_values = S_RANGE

    for s in s_values:
        e_new = e_old + (rm1 - rm_t) / s
        label = f"s={s:.1f}: 必要な ε_mix = {e_new:.3f}"
        x_new = x_for_eps(new_weak, new_strong, e_new)
        if x_new is None:
            lo, hi = EPS_SILICA[new_weak], EPS_SILICA[new_strong]
            print(f"{label}  → この溶媒系では届かない（ε_mix は {lo}〜{hi} の間）")
            continue
        weak_per_strong, phi = volume_ratio_from_x(new_weak, new_strong, x_new)
        print(f"{label}  → {new_weak}:{new_strong} ≒ {weak_per_strong:.2f}:1"
              f"  (強溶媒の体積分率 {phi * 100:.0f}%)")

def main():
    print("TLC: Rf と溶媒比から、目標の Rf になる比を見積もる（同じ化合物、シリカ）")
    while True:
        print()
        print("1: 同じ溶媒系で比を変える   2: 別の溶媒系に移す")
        mode = ask("モード", "1")
        print()
        if mode == "2":
            other_system()
        else:
            same_system()
        print()
        if input("もう一度? [y/N]: ").strip().lower() != "y":
            break


if __name__ == "__main__":
    main()
