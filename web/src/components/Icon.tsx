import iconsUrl from "../assets/icons.svg?url";
import type { Module } from "../types";

export function moduleIcon(m: Pick<Module, "id" | "subgroup">) {
  return m.subgroup === "detectability" ? m.id : m.subgroup ?? "module";
}

export function Icon({ name, className = "" }: { name: string; className?: string }) {
  return (
    <svg className={className} viewBox="0 0 24 24" aria-hidden="true" focusable="false">
      <use href={`${iconsUrl}#${name}`} />
    </svg>
  );
}
