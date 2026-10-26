import iconsUrl from "../assets/icons.svg?url";
import type { Module } from "../types";

// 有独立图标的页面（icons.svg 里的 symbol id = 模块 id）；其余页面用小类图标
const MODULE_ICONS = new Set([
  // Detectability
  "platform_fit", "separation", "ionization", "adducts", "non_ms", "derivatization", "known_methods",
  // MS1
  "ms1_search", "ms1_theory",
  // Inorganic MS
  "icp_isotopes", "icp_interferences", "icp_idms", "icp_massbias", "icp_sp",
]);

export function moduleIcon(m: Pick<Module, "id" | "subgroup">) {
  return MODULE_ICONS.has(m.id) ? m.id : m.subgroup ?? "module";
}

export function Icon({ name, className = "" }: { name: string; className?: string }) {
  return (
    <svg className={className} viewBox="0 0 24 24" aria-hidden="true" focusable="false">
      <use href={`${iconsUrl}#${name}`} />
    </svg>
  );
}
