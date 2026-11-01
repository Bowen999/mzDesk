import iconsUrl from "../assets/icons.svg?url";
import type { Module } from "../types";

// 有独立图标的页面（icons.svg 里的 symbol id = 模块 id）；其余页面用小类图标
const MODULE_ICONS = new Set([
  // Detectability
  "platform_fit", "separation", "ionization", "adducts", "non_ms", "derivatization", "known_methods",
  // MS1
  "ms1_search", "ms1_theory",
  // Spectrum
  "ms2_clean", "ms2_similarity", "ms2_library", "ms2_structure", "ms1_formula", "ms2_denovo", "ms2_ei", "ms2_esi", "ms2_fragments",
  // Inorganic MS
  "icp_isotopes", "icp_interferences", "icp_idms", "icp_massbias", "icp_sp",
  // Molecular
  "identity", "structure_formats", "properties", "lipid_short_name",
  // GNPS2
  "gnps_fastmasst", "gnps_redu", "gnps_massql",
  // Utilities
  "exclusivity", "mrm", "calculator", "report",
]);

export function moduleIcon(m: Pick<Module, "id" | "subgroup">) {
  return MODULE_ICONS.has(m.id) ? m.id : m.subgroup ?? "module";
}

export function pipelineIcon(id: string) {
  return id === "esi_ms2_annotation" ? id : "pipeline";
}

export function Icon({ name, className = "" }: { name: string; className?: string }) {
  return (
    <svg className={className} viewBox="0 0 24 24" aria-hidden="true" focusable="false">
      <use href={`${iconsUrl}#${name}`} />
    </svg>
  );
}
